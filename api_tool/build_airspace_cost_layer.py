"""Build a generic radial airspace cost layer on a routing grid."""

import os

import geopandas as gpd
import numpy as np
from shapely.geometry import Point

from api_tool.select_study_area import select_study_area


CRS_LAT_LON = "EPSG:4326"
CRS_METERS = "EPSG:3857"
NM_TO_M = 1852.0
MI_TO_M = 1609.344
AIRSPACE_MEDIUM_RISK = 60.0
AIRSPACE_HIGH_RISK = 100.0


def radius_to_meters(radius, unit):
    """Convert an airspace radius to meters."""
    if radius is None:
        return None

    unit_upper = unit.upper()
    if unit_upper in ["M", "METER", "METERS"]:
        return float(radius)
    if unit_upper in ["KM", "KILOMETER", "KILOMETERS"]:
        return float(radius) * 1000.0
    if unit_upper in ["NM", "NAUTICAL_MILE", "NAUTICAL_MILES"]:
        return float(radius) * NM_TO_M
    if unit_upper in ["MI", "MILE", "MILES"]:
        return float(radius) * MI_TO_M

    raise ValueError(f"Unsupported radius unit: {unit}")


def load_grid(grid_path):
    """Load a routing grid from disk."""
    if not os.path.exists(grid_path):
        raise FileNotFoundError(f"Routing grid not found: {grid_path}")
    return gpd.read_file(grid_path).reset_index(drop=True)


def score_airspace_sites(grid, airspace_sites):
    """
    Score grid cells from radial airspace site definitions.

    Each site can include label/code, lat, lon, inner_radius, outer_radius, and
    radius_unit. The inner radius receives high risk and the outer shelf
    receives medium risk.
    """
    grid_m = grid.to_crs(CRS_METERS).copy()
    centroids = grid_m.geometry.centroid
    airspace_risk = np.zeros(len(grid_m), dtype=float)
    airspace_sources = [[] for _ in range(len(grid_m))]
    vertical_ranges = [[] for _ in range(len(grid_m))]

    for site in airspace_sites:
        label = site.get("label") or site.get("airport") or site.get("code") or "airspace site"
        latitude = site.get("lat")
        longitude = site.get("lon")
        if latitude is None or longitude is None:
            raise ValueError(f"Airspace site is missing lat/lon: {label}")

        radius_unit = site.get("radius_unit") or site.get("inner_unit") or "NM"
        inner_radius_m = radius_to_meters(site.get("inner_radius"), radius_unit)
        outer_radius_m = radius_to_meters(
            site.get("outer_radius"),
            site.get("outer_unit") or radius_unit,
        )

        if inner_radius_m is None and outer_radius_m is None:
            raise ValueError(f"Airspace site needs an inner or outer radius: {label}")

        airport_point = gpd.GeoSeries(
            [Point(float(longitude), float(latitude))],
            crs=CRS_LAT_LON,
        ).to_crs(grid_m.crs).iloc[0]
        distances = centroids.distance(airport_point).to_numpy()

        site_affected_mask = np.zeros(len(grid_m), dtype=bool)

        if inner_radius_m is not None:
            inner_mask = distances <= inner_radius_m
            site_affected_mask = site_affected_mask | inner_mask
            airspace_risk[inner_mask] = np.maximum(
                airspace_risk[inner_mask],
                site.get("high_risk", AIRSPACE_HIGH_RISK),
            )
            for cell_index in np.flatnonzero(inner_mask):
                airspace_sources[cell_index].append(f"{label} core")

        if outer_radius_m is not None:
            if inner_radius_m is not None and outer_radius_m <= inner_radius_m:
                raise ValueError(f"{label} outer radius must be larger than inner radius.")

            if inner_radius_m is None:
                shelf_mask = distances <= outer_radius_m
            else:
                shelf_mask = (distances > inner_radius_m) & (distances <= outer_radius_m)

            site_affected_mask = site_affected_mask | shelf_mask
            airspace_risk[shelf_mask] = np.maximum(
                airspace_risk[shelf_mask],
                site.get("medium_risk", AIRSPACE_MEDIUM_RISK),
            )
            for cell_index in np.flatnonzero(shelf_mask):
                airspace_sources[cell_index].append(f"{label} shelf")

        vertical_range = site.get("vertical_range")
        if vertical_range:
            for cell_index in np.flatnonzero(site_affected_mask):
                if airspace_sources[cell_index]:
                    vertical_ranges[cell_index].append(vertical_range)

    scored = grid.copy()
    scored["airport_risk"] = airspace_risk
    scored["airport_risk_combined"] = airspace_risk
    scored["air_class"] = "Low"
    scored.loc[scored["airport_risk_combined"] > 30, "air_class"] = "Medium"
    scored.loc[scored["airport_risk_combined"] > 70, "air_class"] = "High"
    scored["airspace_source"] = [
        ", ".join(source_names) if source_names else "None"
        for source_names in airspace_sources
    ]
    scored["airspace_vertical_range"] = [
        ", ".join(dict.fromkeys(range_names)) if range_names else "None"
        for range_names in vertical_ranges
    ]
    return scored


def save_airspace_layer(grid, output_path):
    """Save the airspace-scored grid to disk."""
    output_folder = os.path.dirname(output_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)
    grid.to_file(output_path, driver="GeoJSON")


def build_airspace_cost_layer(
    city_name=None,
    grid=None,
    grid_path=None,
    airspace_sites=None,
    output_path=None,
):
    """Build and optionally save a generic airspace cost layer."""
    if city_name is None:
        city_name = input("City or study area name: ").strip()
        grid_path = input("Routing grid path: ").strip()
        output_path = input("Airspace layer output path (optional): ").strip() or None

    if grid is None:
        grid = load_grid(grid_path)

    airspace_sites = airspace_sites or []
    scored_grid = score_airspace_sites(grid, airspace_sites)

    if output_path:
        save_airspace_layer(scored_grid, output_path)

    return {
        "study_area": select_study_area(city_name=city_name),
        "grid": scored_grid,
        "airspace_sites": airspace_sites,
        "site_count": len(airspace_sites),
        "output_path": output_path,
        "output_exists": os.path.exists(output_path) if output_path else False,
    }


if __name__ == "__main__":
    result = build_airspace_cost_layer()
    print(
        {
            "city": result["study_area"]["city"],
            "site_count": result["site_count"],
            "cell_count": len(result["grid"]),
        }
    )
