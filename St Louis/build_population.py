"""Prepare population density data for the St. Louis study area."""

import os

import geopandas as gpd
import pandas as pd
from shapely.geometry import box

import settings


REQUIRED_POPULATION_COLUMNS = {"GEO_ID", "B01003_001E"}
REQUIRED_BLOCK_GROUP_COLUMNS = {"GEOID", "ALAND"}


def load_population_data(population_csv):
    """Load and clean ACS block-group population totals."""
    population_data = pd.read_csv(population_csv, dtype=str)
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


def load_block_groups(block_group_shapefile, state_code):
    """Load TIGER/Line block-group boundaries for one state."""
    block_groups = gpd.read_file(block_group_shapefile)
    missing_columns = REQUIRED_BLOCK_GROUP_COLUMNS - set(block_groups.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Block-group shapefile is missing columns: {missing_text}")

    block_groups = block_groups.copy()
    block_groups["source_state"] = state_code
    return block_groups


def load_population_sources():
    """Load Missouri and Illinois block groups with their population totals."""
    source_frames = []
    source_crs = None

    for source in settings.POPULATION_SOURCES:
        population_csv = source["population_csv"]
        block_group_shapefile = source["block_group_shapefile"]
        state_code = source["state"]

        print("Current working directory:", os.getcwd())
        print("Looking for population file:", population_csv)
        print("Absolute path:", os.path.abspath(population_csv))

        if not os.path.exists(population_csv):
            raise FileNotFoundError(f"Population source not found: {population_csv}")
        if not os.path.exists(block_group_shapefile):
            raise FileNotFoundError(
                f"Block-group shapefile not found: {block_group_shapefile}"
            )

        population_data = load_population_data(population_csv)
        block_groups = load_block_groups(block_group_shapefile, state_code)
        merged = block_groups.merge(
            population_data,
            on="GEOID",
            how="left",
            validate="one_to_one",
        )

        missing_population = merged["population"].isna().sum()
        if missing_population:
            raise ValueError(
                f"{missing_population} {state_code} block groups have no matching population value."
            )

        if source_crs is None:
            source_crs = merged.crs
        merged["source_state"] = state_code
        source_frames.append(merged)

    combined = gpd.GeoDataFrame(
        pd.concat(source_frames, ignore_index=True),
        crs=source_crs,
    )
    return combined


def calculate_population_density(block_groups):
    """Calculate density and assign study-area population risk bands."""
    land_area_km2 = block_groups["ALAND"] / 1_000_000
    block_groups["density_p_km2"] = block_groups["population"].div(
        land_area_km2.where(land_area_km2 > 0)
    )

    study_area = gpd.GeoSeries(
        [
            box(
                settings.WEST,
                settings.SOUTH,
                settings.EAST,
                settings.NORTH,
            )
        ],
        crs="EPSG:4326",
    ).to_crs(block_groups.crs).iloc[0]
    study_block_groups = block_groups[
        block_groups.geometry.intersects(study_area)
    ].copy()
    if study_block_groups.empty:
        raise ValueError("No population block groups intersect the study area.")

    valid_density = study_block_groups["density_p_km2"].dropna()
    medium_threshold, high_threshold = valid_density.quantile([0.70, 0.90])

    study_block_groups["density_risk"] = "low"
    study_block_groups.loc[
        study_block_groups["density_p_km2"] >= medium_threshold,
        "density_risk",
    ] = "medium"
    study_block_groups.loc[
        study_block_groups["density_p_km2"] >= high_threshold,
        "density_risk",
    ] = "high"

    return study_block_groups, medium_threshold, high_threshold


def main():
    population_sources = settings.POPULATION_SOURCES
    if not population_sources:
        raise ValueError("Add at least one population source in St Louis/settings.py.")

    combined_block_groups = load_population_sources()
    result, medium_threshold, high_threshold = calculate_population_density(
        combined_block_groups,
    )

    result = result.to_crs("EPSG:4326")
    result.to_file(settings.POPULATION_GEOJSON, driver="GeoJSON")

    print(f"Population rows: {len(result)}")
    print(f"Source states: {', '.join(sorted(result['source_state'].unique()))}")
    print(f"Block groups: {len(result)}")
    print(f"Medium-risk threshold: {medium_threshold:.2f} people/km2")
    print(f"High-risk threshold: {high_threshold:.2f} people/km2")
    print(f"Saved: {settings.POPULATION_GEOJSON}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error
