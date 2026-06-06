"""Build the St. Louis routing grid and map data."""

import os

import geopandas as gpd
import numpy as np
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


def main():
    grid = build_base_grid()

    if not os.path.exists(settings.POPULATION_GEOJSON):
        raise FileNotFoundError(
            f"Population file not found: {settings.POPULATION_GEOJSON}"
        )

    raise NotImplementedError(
        "The St. Louis base grid is configured, but its population join and "
        "airspace risk model still need verified local data before the grid "
        "can be saved."
    )


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, NotImplementedError, ValueError) as error:
        raise SystemExit(str(error)) from error
