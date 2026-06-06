"""Build the St. Louis routing grid and map data."""

import os

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

import settings


CRS_LAT_LON = "EPSG:4326"
CRS_METERS = "EPSG:3857"


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

    grid["density_p_km2"] = pd.to_numeric(
        population_join["density_p_km2"],
        errors="coerce",
    ).fillna(0.0)
    grid["density_risk"] = population_join["density_risk"].fillna("low")
    grid["city_risk"] = (grid["density_p_km2"] / 10.0).clip(0, 120)
    grid["pop_class"] = grid["density_risk"].str.title()

    # Preserve the Chicago grid schema until airport risk is added.
    grid["airport_risk"] = 0.0
    grid["airport_risk_combined"] = 0.0
    grid["air_class"] = "Low"
    grid["density_type"] = "population"
    grid["risk_cost"] = settings.POPULATION_WEIGHT * grid["city_risk"]
    grid["risk_class"] = grid["pop_class"]
    return grid


def main():
    grid = build_base_grid()

    if not os.path.exists(settings.POPULATION_GEOJSON):
        raise FileNotFoundError(
            f"Population file not found: {settings.POPULATION_GEOJSON}"
        )

    grid = add_population_risk(grid)
    grid.to_crs(CRS_LAT_LON).to_file(
        settings.RISK_GRID_GEOJSON,
        driver="GeoJSON",
    )

    print(f"Grid cells: {len(grid)}")
    print(grid["pop_class"].value_counts().to_string())
    print(f"Saved St. Louis risk grid: {settings.RISK_GRID_GEOJSON}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error
