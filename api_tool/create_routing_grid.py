"""Create a routing grid from study-area bounds."""

import os

import geopandas as gpd
import numpy as np
from shapely.geometry import box

from api_tool.select_study_area import prompt_float, prompt_int, select_study_area


CRS_LAT_LON = "EPSG:4326"
CRS_METERS = "EPSG:3857"


def request_grid_settings_from_user():
    """Ask the user for routing-grid settings."""
    city_name = input("City or study area name: ").strip()
    while city_name == "":
        city_name = input("City or study area name is required: ").strip()

    bounds = {
        "west": prompt_float("West longitude: ", required=True),
        "south": prompt_float("South latitude: ", required=True),
        "east": prompt_float("East longitude: ", required=True),
        "north": prompt_float("North latitude: ", required=True),
    }
    cell_size_m = prompt_int("Grid cell size in meters: ", required=True)
    output_path = input("Grid output path (optional): ").strip() or None
    return city_name, bounds, cell_size_m, output_path


def validate_bounds(bounds):
    """Check that bounds contain west, south, east, and north values."""
    required_fields = ["west", "south", "east", "north"]
    missing_fields = [field for field in required_fields if field not in bounds]
    if missing_fields:
        missing_text = ", ".join(missing_fields)
        raise ValueError(f"Bounds are missing: {missing_text}")

    if bounds["west"] >= bounds["east"]:
        raise ValueError("west must be less than east.")
    if bounds["south"] >= bounds["north"]:
        raise ValueError("south must be less than north.")


def get_projected_crs(study_area):
    """Return a local projected CRS for meter-based grid creation."""
    projected_crs = study_area.estimate_utm_crs()
    if projected_crs is None:
        return CRS_METERS
    return projected_crs


def create_grid_from_bounds(bounds, cell_size_m, output_crs=CRS_LAT_LON, projected_crs=None):
    """Create a clipped rectangular grid from latitude/longitude bounds."""
    validate_bounds(bounds)
    if cell_size_m <= 0:
        raise ValueError("cell_size_m must be greater than 0.")

    study_area_lat_lon = gpd.GeoDataFrame(
        geometry=[
            box(
                bounds["west"],
                bounds["south"],
                bounds["east"],
                bounds["north"],
            )
        ],
        crs=CRS_LAT_LON,
    )
    projected_crs = projected_crs or get_projected_crs(study_area_lat_lon)
    study_area = study_area_lat_lon.to_crs(projected_crs)

    xmin, ymin, xmax, ymax = study_area.total_bounds
    x_values = np.arange(xmin, xmax + cell_size_m, cell_size_m)
    y_values = np.arange(ymin, ymax + cell_size_m, cell_size_m)

    cells = []
    for x_value in x_values[:-1]:
        for y_value in y_values[:-1]:
            cells.append(
                box(
                    x_value,
                    y_value,
                    x_value + cell_size_m,
                    y_value + cell_size_m,
                )
            )

    grid = gpd.GeoDataFrame(geometry=cells, crs=projected_crs)
    grid = gpd.clip(grid, study_area).reset_index(drop=True)
    grid["cell_id"] = range(len(grid))
    grid["cell_size_m"] = cell_size_m

    if output_crs:
        grid = grid.to_crs(output_crs)

    return grid


def save_grid(grid, output_path):
    """Save a routing grid to disk."""
    output_folder = os.path.dirname(output_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)
    grid.to_file(output_path, driver="GeoJSON")


def create_routing_grid(
    city_name=None,
    bounds=None,
    cell_size_m=None,
    output_path=None,
    output_crs=CRS_LAT_LON,
    projected_crs=None,
):
    """Create and optionally save a routing grid."""
    if city_name is None:
        city_name, bounds, cell_size_m, output_path = request_grid_settings_from_user()

    grid = create_grid_from_bounds(
        bounds,
        cell_size_m,
        output_crs=output_crs,
        projected_crs=projected_crs,
    )

    if output_path:
        save_grid(grid, output_path)

    return {
        "study_area": select_study_area(
            city_name=city_name,
            cell_size_m=cell_size_m,
            bounds=bounds,
            risk_grid_path=output_path,
        ),
        "grid": grid,
        "cell_count": len(grid),
        "bounds": bounds,
        "cell_size_m": cell_size_m,
        "output_path": output_path,
        "output_exists": os.path.exists(output_path) if output_path else False,
    }


if __name__ == "__main__":
    result = create_routing_grid()
    print(
        {
            "city": result["study_area"]["city"],
            "cell_count": result["cell_count"],
            "output_path": result["output_path"],
        }
    )
