"""Prepare population density data for the St. Louis study area."""

import os

import geopandas as gpd
import pandas as pd

import settings


REQUIRED_POPULATION_COLUMNS = {"GEO_ID", "B01003_001E"}
REQUIRED_BLOCK_GROUP_COLUMNS = {"GEOID", "ALAND"}


def load_population_data():
    """Load and clean ACS block-group population totals."""
    population_data = pd.read_csv(settings.POPULATION_CSV, dtype=str)
    missing_columns = REQUIRED_POPULATION_COLUMNS - set(population_data.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Population CSV is missing columns: {missing_text}")

    population_data = population_data[
        population_data["GEO_ID"].str.startswith("1500000US", na=False)
    ].copy()
    population_data["GEOID"] = population_data["GEO_ID"].str.removeprefix(
        "1500000US"
    )
    population_data["population"] = pd.to_numeric(
        population_data["B01003_001E"],
        errors="coerce",
    )

    if population_data["population"].isna().any():
        invalid_count = population_data["population"].isna().sum()
        raise ValueError(
            f"Population CSV contains {invalid_count} invalid population values."
        )

    return population_data[["GEOID", "population"]]


def load_block_groups():
    """Load Missouri TIGER/Line block-group boundaries."""
    block_groups = gpd.read_file(settings.BLOCK_GROUP_SHAPEFILE)
    missing_columns = REQUIRED_BLOCK_GROUP_COLUMNS - set(block_groups.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Block-group shapefile is missing columns: {missing_text}")

    return block_groups


def calculate_population_density(block_groups, population_data):
    """Join population totals and calculate people per square kilometer."""
    merged = block_groups.merge(
        population_data,
        on="GEOID",
        how="left",
        validate="one_to_one",
    )

    missing_population = merged["population"].isna().sum()
    if missing_population:
        raise ValueError(
            f"{missing_population} block groups have no matching population value."
        )

    land_area_km2 = merged["ALAND"] / 1_000_000
    merged["density_p_km2"] = merged["population"].div(
        land_area_km2.where(land_area_km2 > 0)
    )

    valid_density = merged["density_p_km2"].dropna()
    medium_threshold, high_threshold = valid_density.quantile([0.70, 0.90])

    merged["density_risk"] = "low"
    merged.loc[
        merged["density_p_km2"] >= medium_threshold,
        "density_risk",
    ] = "medium"
    merged.loc[
        merged["density_p_km2"] >= high_threshold,
        "density_risk",
    ] = "high"

    return merged, medium_threshold, high_threshold


def main():
    required_paths = [
        settings.POPULATION_CSV,
        settings.BLOCK_GROUP_SHAPEFILE,
    ]
    for required_path in required_paths:
        if not os.path.exists(required_path):
            raise FileNotFoundError(f"Population source not found: {required_path}")

    population_data = load_population_data()
    block_groups = load_block_groups()
    result, medium_threshold, high_threshold = calculate_population_density(
        block_groups,
        population_data,
    )

    result = result.to_crs("EPSG:4326")
    result.to_file(settings.POPULATION_GEOJSON, driver="GeoJSON")

    print(f"Population rows: {len(population_data)}")
    print(f"Block groups: {len(result)}")
    print(f"Medium-risk threshold: {medium_threshold:.2f} people/km2")
    print(f"High-risk threshold: {high_threshold:.2f} people/km2")
    print(f"Saved: {settings.POPULATION_GEOJSON}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error
