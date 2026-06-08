"""Build the St. Louis routing grid and map data."""

import os

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point, box

import settings


CRS_LAT_LON = "EPSG:4326"
CRS_METERS = "EPSG:26915"
NM_TO_M = 1852.0
AIRSPACE_MEDIUM_RISK = 60.0
AIRSPACE_HIGH_RISK = 100.0


def build_base_grid():
    """Create the rectangular St. Louis study grid."""
    settings.validate_grid_settings()

    study_area = gpd.GeoDataFrame(
        geometry=[
            box(
                settings.WEST,
                settings.SOUTH,
                settings.EAST,
                settings.NORTH,
            )
        ],
        crs=CRS_LAT_LON,
    ).to_crs(CRS_METERS)

    xmin, ymin, xmax, ymax = study_area.total_bounds
    x_values = np.arange(xmin, xmax + settings.CELL_SIZE_M, settings.CELL_SIZE_M)
    y_values = np.arange(ymin, ymax + settings.CELL_SIZE_M, settings.CELL_SIZE_M)

    cells = []
    for x_value in x_values[:-1]:
        for y_value in y_values[:-1]:
            cells.append(
                box(
                    x_value,
                    y_value,
                    x_value + settings.CELL_SIZE_M,
                    y_value + settings.CELL_SIZE_M,
                )
            )

    grid = gpd.GeoDataFrame(geometry=cells, crs=CRS_METERS)
    return gpd.clip(grid, study_area).reset_index(drop=True)


def add_population_risk(grid):
    """Assign block-group population density to each grid-cell centroid."""
    population = gpd.read_file(settings.POPULATION_GEOJSON).to_crs(CRS_METERS)
    required_columns = {"density_p_km2", "density_risk"}
    missing_columns = required_columns - set(population.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Population GeoJSON is missing columns: {missing_text}")

    centroids = gpd.GeoDataFrame(
        {"cell_id": np.arange(len(grid))},
        geometry=grid.geometry.centroid,
        crs=grid.crs,
    )
    population_join = gpd.sjoin(
        centroids,
        population[["density_p_km2", "density_risk", "geometry"]],
        how="left",
        predicate="within",
    )
    population_join = population_join.drop_duplicates("cell_id")
    population_join = population_join.set_index("cell_id").reindex(
        range(len(grid))
    )

    population_density = pd.to_numeric(
        population_join["density_p_km2"],
        errors="coerce",
    )
    grid["population_data_available"] = population_density.notna()
    grid["density_p_km2"] = population_density
    grid["density_risk"] = population_join["density_risk"]
    grid["city_risk"] = (population_density / 10.0).clip(0, 120)
    grid["pop_class"] = grid["density_risk"].str.title()
    grid.loc[~grid["population_data_available"], "pop_class"] = "No Data"

    return grid


def radius_to_meters(radius, unit):
    """Convert a supplied airspace radius to meters."""
    if unit == "NM":
        return radius * NM_TO_M
    raise ValueError(f"Unsupported airspace radius unit: {unit}")


def add_airspace_risk(grid):
    """Assign airport airspace risk from the configured radial boundaries."""
    centroids = grid.geometry.centroid
    airspace_risk = np.zeros(len(grid), dtype=float)
    airspace_sources = [[] for _ in range(len(grid))]
    vertical_ranges = [[] for _ in range(len(grid))]

    for site in settings.AIRSPACE_SITES:
        latitude, longitude = site["location"]
        airport_point = gpd.GeoSeries(
            [Point(longitude, latitude)],
            crs=CRS_LAT_LON,
        ).to_crs(CRS_METERS).iloc[0]
        distances = centroids.distance(airport_point).to_numpy()
        inner_radius_m = radius_to_meters(
            site["inner_radius"],
            site["inner_unit"],
        )
        inner_mask = distances <= inner_radius_m
        airspace_risk[inner_mask] = np.maximum(
            airspace_risk[inner_mask],
            AIRSPACE_HIGH_RISK,
        )

        for cell_index in np.flatnonzero(inner_mask):
            airspace_sources[cell_index].append(
                f"{site['code']} {site['airspace_type']} core"
            )
            vertical_ranges[cell_index].append(site["vertical_range"])

        if site["outer_radius"] is None:
            continue

        outer_radius_m = radius_to_meters(
            site["outer_radius"],
            site["outer_unit"],
        )
        shelf_mask = (distances > inner_radius_m) & (
            distances <= outer_radius_m
        )
        airspace_risk[shelf_mask] = np.maximum(
            airspace_risk[shelf_mask],
            AIRSPACE_MEDIUM_RISK,
        )

        for cell_index in np.flatnonzero(shelf_mask):
            airspace_sources[cell_index].append(
                f"{site['code']} {site['airspace_type']} shelf"
            )
            vertical_ranges[cell_index].append(site["vertical_range"])

    grid["airport_risk"] = airspace_risk
    grid["airport_risk_combined"] = airspace_risk
    grid["air_class"] = "Low"
    grid.loc[grid["airport_risk_combined"] > 30, "air_class"] = "Medium"
    grid.loc[grid["airport_risk_combined"] > 70, "air_class"] = "High"
    grid["airspace_source"] = [
        ", ".join(source_names) if source_names else "None"
        for source_names in airspace_sources
    ]
    grid["airspace_vertical_range"] = [
        ", ".join(dict.fromkeys(range_names)) if range_names else "None"
        for range_names in vertical_ranges
    ]
    return grid


def add_combined_risk(grid):
    """Combine population and airspace costs using the route weights."""
    grid["risk_cost"] = (
        settings.POPULATION_WEIGHT * grid["city_risk"].fillna(0.0)
        + settings.AIRSPACE_WEIGHT * grid["airport_risk_combined"]
    )
    grid["risk_class"] = "Low"
    grid.loc[grid["risk_cost"] > 30, "risk_class"] = "Medium"
    grid.loc[grid["risk_cost"] > 70, "risk_class"] = "High"

    grid["density_type"] = "low"
    city_risk = grid["city_risk"].fillna(0.0)
    population_dominant = (
        grid["population_data_available"]
        & ((grid["pop_class"] != "Low") | (grid["air_class"] != "Low"))
        & (city_risk >= grid["airport_risk_combined"])
    )
    airspace_dominant = (
        (grid["air_class"] != "Low")
        & (grid["airport_risk_combined"] > city_risk)
    )
    grid.loc[population_dominant, "density_type"] = "population"
    grid.loc[airspace_dominant, "density_type"] = "airspace"
    grid.loc[
        ~grid["population_data_available"] & (grid["air_class"] == "Low"),
        "density_type",
    ] = "population data unavailable"
    return grid


def main():
    grid = build_base_grid()

    if not os.path.exists(settings.POPULATION_GEOJSON):
        raise FileNotFoundError(
            f"Population file not found: {settings.POPULATION_GEOJSON}"
        )

    grid = add_population_risk(grid)
    grid = add_airspace_risk(grid)
    grid = add_combined_risk(grid)
    grid.to_crs(CRS_LAT_LON).to_file(
        settings.RISK_GRID_GEOJSON,
        driver="GeoJSON",
    )

    print(f"Grid cells: {len(grid)}")
    print(grid["pop_class"].value_counts().to_string())
    missing_population = (~grid["population_data_available"]).sum()
    print(f"Cells without population data: {missing_population}")
    print(grid["air_class"].value_counts().to_string())
    print(grid["risk_class"].value_counts().to_string())
    print(f"Saved St. Louis risk grid: {settings.RISK_GRID_GEOJSON}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error
