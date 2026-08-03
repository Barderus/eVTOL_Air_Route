"""Score routing-grid cells with population, airspace, and traffic costs."""

import os

import geopandas as gpd
import numpy as np
import pandas as pd

from api_tool.build_airspace_cost_layer import score_airspace_sites
from api_tool.select_study_area import select_study_area


CRS_LAT_LON = "EPSG:4326"
CRS_METERS = "EPSG:3857"


def load_grid(grid_path):
    """Load a routing grid from disk."""
    if not os.path.exists(grid_path):
        raise FileNotFoundError(f"Routing grid not found: {grid_path}")
    return gpd.read_file(grid_path).reset_index(drop=True)


def classify_values(values, medium_threshold, high_threshold):
    """Convert numeric values into Low, Medium, and High classes."""
    classes = np.full(len(values), "Low", dtype=object)
    classes[values > medium_threshold] = "Medium"
    classes[values > high_threshold] = "High"
    return classes


def add_population_score(
    grid,
    population_path=None,
    density_field="density_p_km2",
    risk_divisor=10.0,
    max_city_risk=120.0,
):
    """Assign population density and city_risk to each grid-cell centroid."""
    scored_grid = grid.copy()
    scored_grid["population_data_available"] = False
    scored_grid["density_p_km2"] = np.nan
    scored_grid["density_risk"] = None
    scored_grid["city_risk"] = 0.0
    scored_grid["pop_class"] = "No Data"

    if not population_path:
        return scored_grid
    if not os.path.exists(population_path):
        raise FileNotFoundError(f"Population layer not found: {population_path}")

    projected_grid = scored_grid.to_crs(CRS_METERS)
    population = gpd.read_file(population_path).to_crs(CRS_METERS)
    if density_field not in population.columns:
        raise ValueError(f"Population layer is missing density field: {density_field}")

    centroids = gpd.GeoDataFrame(
        {"cell_id": np.arange(len(scored_grid))},
        geometry=projected_grid.geometry.centroid,
        crs=projected_grid.crs,
    )
    joined = gpd.sjoin(
        centroids,
        population[[density_field, "geometry"]],
        how="left",
        predicate="within",
    )
    joined = joined.drop_duplicates("cell_id")
    joined = joined.set_index("cell_id").reindex(range(len(scored_grid)))

    density = pd.to_numeric(joined[density_field], errors="coerce")
    scored_grid["population_data_available"] = density.notna()
    scored_grid["density_p_km2"] = density
    scored_grid["city_risk"] = (density / risk_divisor).clip(0, max_city_risk).fillna(0.0)
    scored_grid["pop_class"] = classify_values(
        scored_grid["city_risk"].to_numpy(dtype=float),
        medium_threshold=40,
        high_threshold=90,
    )
    scored_grid.loc[~scored_grid["population_data_available"], "pop_class"] = "No Data"
    return scored_grid


def add_airspace_score(grid, airspace_sites=None):
    """Assign airport_risk_combined to grid cells."""
    if not airspace_sites:
        scored_grid = grid.copy()
        scored_grid["airport_risk"] = 0.0
        scored_grid["airport_risk_combined"] = 0.0
        scored_grid["air_class"] = "Low"
        scored_grid["airspace_source"] = "None"
        scored_grid["airspace_vertical_range"] = "None"
        return scored_grid

    return score_airspace_sites(grid, airspace_sites)


def add_traffic_score(grid, traffic_csv_path=None):
    """Count traffic observations in each grid cell and add traffic_risk."""
    scored_grid = grid.copy()
    scored_grid["traffic_count"] = 0.0
    scored_grid["traffic_risk"] = 0.0

    if not traffic_csv_path:
        return scored_grid
    if not os.path.exists(traffic_csv_path):
        raise FileNotFoundError(f"Traffic CSV not found: {traffic_csv_path}")

    traffic_data = pd.read_csv(traffic_csv_path)
    required_columns = {"lat", "lon"}
    if not required_columns.issubset(traffic_data.columns):
        raise ValueError(f"Traffic CSV must contain lat and lon columns: {traffic_csv_path}")

    traffic_data["lat"] = pd.to_numeric(traffic_data["lat"], errors="coerce")
    traffic_data["lon"] = pd.to_numeric(traffic_data["lon"], errors="coerce")
    traffic_data = traffic_data.dropna(subset=["lat", "lon"])

    if traffic_data.empty:
        return scored_grid

    traffic_points = gpd.GeoDataFrame(
        traffic_data,
        geometry=gpd.points_from_xy(traffic_data["lon"], traffic_data["lat"]),
        crs=CRS_LAT_LON,
    ).to_crs(scored_grid.crs)
    joined = gpd.sjoin(
        traffic_points[["geometry"]],
        scored_grid[["geometry"]],
        how="inner",
        predicate="within",
    )

    counts = joined.groupby("index_right").size()
    traffic_count = counts.reindex(range(len(scored_grid)), fill_value=0).to_numpy(dtype=float)
    scored_grid["traffic_count"] = traffic_count

    if traffic_count.max() > 0:
        scored_grid["traffic_risk"] = traffic_count / traffic_count.max()

    return scored_grid


def add_combined_risk(
    grid,
    population_weight=0.9,
    airspace_weight=1.4,
):
    """Combine scored cell costs into risk_cost and display classes."""
    scored_grid = grid.copy()
    if "city_risk" not in scored_grid.columns:
        scored_grid["city_risk"] = 0.0
    if "airport_risk_combined" not in scored_grid.columns:
        scored_grid["airport_risk_combined"] = 0.0

    scored_grid["risk_cost"] = (
        population_weight * scored_grid["city_risk"].fillna(0.0)
        + airspace_weight * scored_grid["airport_risk_combined"].fillna(0.0)
    )
    scored_grid["risk_class"] = classify_values(
        scored_grid["risk_cost"].to_numpy(dtype=float),
        medium_threshold=30,
        high_threshold=70,
    )

    scored_grid["density_type"] = "low"
    population_dominant = scored_grid["city_risk"].fillna(0.0) >= scored_grid[
        "airport_risk_combined"
    ].fillna(0.0)
    scored_grid.loc[
        (scored_grid["risk_class"] != "Low") & population_dominant,
        "density_type",
    ] = "population"
    scored_grid.loc[
        (scored_grid["risk_class"] != "Low") & ~population_dominant,
        "density_type",
    ] = "airspace"
    return scored_grid


def save_scored_grid(grid, output_path):
    """Save a scored grid as GeoJSON."""
    output_folder = os.path.dirname(output_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)
    grid.to_crs(CRS_LAT_LON).to_file(output_path, driver="GeoJSON")


def score_grid_cells(
    city_name=None,
    grid=None,
    grid_path=None,
    population_path=None,
    population_density_field="density_p_km2",
    airspace_sites=None,
    traffic_csv_path=None,
    population_weight=0.9,
    airspace_weight=1.4,
    output_path=None,
):
    """Score a routing grid and optionally save the result."""
    if city_name is None:
        city_name = input("City or study area name: ").strip()
        grid_path = input("Input grid path: ").strip()
        population_path = input("Population layer path (optional): ").strip() or None
        traffic_csv_path = input("Traffic CSV path (optional): ").strip() or None
        output_path = input("Scored grid output path (optional): ").strip() or None

    if grid is None:
        grid = load_grid(grid_path)

    grid = add_population_score(
        grid,
        population_path=population_path,
        density_field=population_density_field,
    )
    grid = add_airspace_score(grid, airspace_sites=airspace_sites)
    grid = add_traffic_score(grid, traffic_csv_path=traffic_csv_path)
    grid = add_combined_risk(
        grid,
        population_weight=population_weight,
        airspace_weight=airspace_weight,
    )

    if output_path:
        save_scored_grid(grid, output_path)

    return {
        "study_area": select_study_area(city_name=city_name),
        "grid": grid,
        "cell_count": len(grid),
        "population_path": population_path,
        "traffic_csv_path": traffic_csv_path,
        "airspace_site_count": len(airspace_sites or []),
        "output_path": output_path,
        "output_exists": os.path.exists(output_path) if output_path else False,
    }


if __name__ == "__main__":
    result = score_grid_cells()
    print(
        {
            "city": result["study_area"]["city"],
            "cell_count": result["cell_count"],
            "output_path": result["output_path"],
        }
    )
