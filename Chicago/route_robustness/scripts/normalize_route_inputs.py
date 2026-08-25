"""Normalize route robustness inputs for any city or route pair.

This script prepares the first route robustness artifacts:

- normalized route-factor summary
- normalized weight combinations

It is configured first for Clow International Airport in Bolingbrook to
Chicago Union Station, but the path and field settings near the top can be
edited for another city or route.
"""

import math
import os

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd


# Input files
GRID_PATH = "Chicago/geojson/risk_grid_v7.geojson"
TRAFFIC_CSV_PATH = "Chicago/opensky/output/ohare_2026-03-07_1s_15nm_bbox.csv"

# Output files
OUTPUT_FOLDER = "Chicago/route_robustness/output"
WEIGHT_CONFIGURATIONS_CSV = os.path.join(
    OUTPUT_FOLDER,
    "weight_configurations.csv",
)
NORMALIZATION_SUMMARY_CSV = os.path.join(
    OUTPUT_FOLDER,
    "normalization_summary.csv",
)

# Route settings for the first test case
ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521

# Grid field settings. Edit these for another city if its column names differ.
POPULATION_FIELD = "city_risk"
AIRSPACE_FIELD = "airport_risk_combined"

# Normalization settings
WEIGHT_STEP = 0.10
AIRSPACE_CLIP_PERCENTILE = 99.0
AIRSPACE_POWER = 0.75
TRAFFIC_CLIP_PERCENTILE = 95.0
TRAFFIC_TRANSFORM = "log1p"
TRAFFIC_POWER = 0.75

# Existing project weight ratio, normalized to sum to 1.
CURRENT_DISTANCE_WEIGHT = 0.6
CURRENT_POPULATION_WEIGHT = 0.9
CURRENT_TRAFFIC_WEIGHT = 1.0
CURRENT_AIRSPACE_WEIGHT = 1.4


def detect_csv_encoding(csv_path):
    """Detect common CSV byte-order marks."""
    with open(csv_path, "rb") as file_handle:
        first_bytes = file_handle.read(4)

    if first_bytes.startswith(b"\xff\xfe") or first_bytes.startswith(b"\xfe\xff"):
        return "utf-16"
    if first_bytes.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8"


def normalize_values(
    values,
    clip_percentile=None,
    transform=None,
    power=1.0,
):
    """Scale non-negative values to the 0 to 1 range."""
    raw_values = np.asarray(values, dtype=float)
    transformed_values = np.clip(raw_values, 0.0, None)

    if transform == "log1p":
        transformed_values = np.log1p(transformed_values)
    elif transform is not None:
        raise ValueError(f"Unsupported transform: {transform}")

    if clip_percentile is None:
        scale = float(np.max(transformed_values))
    else:
        positive_values = transformed_values[transformed_values > 0]
        if positive_values.size == 0:
            normalized_values = np.zeros_like(transformed_values)
            return normalized_values, 1.0
        scale = float(np.percentile(positive_values, clip_percentile))

    if scale <= 0:
        normalized_values = np.zeros_like(transformed_values)
        return normalized_values, 1.0

    normalized_values = np.clip(transformed_values / scale, 0.0, 1.0)
    if power != 1.0:
        normalized_values = np.power(normalized_values, power)

    return normalized_values, scale


def require_grid_columns(grid):
    """Stop early if the configured grid fields are missing."""
    required_columns = {POPULATION_FIELD, AIRSPACE_FIELD}
    missing_columns = required_columns.difference(grid.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Routing grid is missing required columns: {missing_text}")


def build_graph(grid):
    """Build a graph from touching grid cells and return edge distances."""
    grid = grid.reset_index(drop=True)
    projected_grid = grid.to_crs("EPSG:3857")
    centroids = projected_grid.geometry.centroid
    centroid_x = centroids.x.to_numpy()
    centroid_y = centroids.y.to_numpy()

    graph = nx.Graph()
    graph.add_nodes_from(range(len(grid)))

    spatial_index = grid.sindex
    for cell_index, geometry in enumerate(grid.geometry):
        candidate_indexes = spatial_index.query(geometry, predicate="touches")
        for candidate_index in candidate_indexes:
            candidate_index = int(candidate_index)
            if candidate_index > cell_index:
                graph.add_edge(cell_index, candidate_index)

    edge_distances = []
    for start_node, end_node in graph.edges:
        distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        edge_distances.append(distance)

    return graph, np.asarray(edge_distances, dtype=float)


def load_traffic_counts(csv_path, grid):
    """Count traffic observations inside each grid cell."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Traffic CSV not found: {csv_path}")

    csv_encoding = detect_csv_encoding(csv_path)
    traffic_data = pd.read_csv(csv_path, encoding=csv_encoding)
    required_columns = {"lat", "lon"}
    if not required_columns.issubset(traffic_data.columns):
        raise ValueError(f"Traffic CSV must contain lat and lon columns: {csv_path}")

    traffic_data = traffic_data.dropna(subset=["lat", "lon"]).copy()
    traffic_data["lat"] = pd.to_numeric(traffic_data["lat"], errors="coerce")
    traffic_data["lon"] = pd.to_numeric(traffic_data["lon"], errors="coerce")
    traffic_data = traffic_data.dropna(subset=["lat", "lon"])

    traffic_points = gpd.GeoDataFrame(
        traffic_data[["lat", "lon"]],
        geometry=gpd.points_from_xy(traffic_data["lon"], traffic_data["lat"]),
        crs="EPSG:4326",
    ).to_crs(grid.crs)

    joined_points = gpd.sjoin(
        traffic_points,
        grid[["geometry"]],
        how="inner",
        predicate="within",
    )
    counts = joined_points.groupby("index_right").size()
    traffic_counts = counts.reindex(range(len(grid)), fill_value=0)
    return traffic_counts.to_numpy(dtype=float)


def make_summary_row(
    factor_name,
    raw_values,
    normalized_values,
    method,
    scale,
    clip_percentile=None,
    transform=None,
    power=1.0,
):
    """Create one normalization summary row."""
    raw_values = np.asarray(raw_values, dtype=float)
    normalized_values = np.asarray(normalized_values, dtype=float)
    return {
        "factor": factor_name,
        "raw_min": float(np.min(raw_values)),
        "raw_max": float(np.max(raw_values)),
        "method": method,
        "scale": float(scale),
        "clip_percentile": clip_percentile,
        "transform": transform,
        "power": power,
        "normalized_min": float(np.min(normalized_values)),
        "normalized_max": float(np.max(normalized_values)),
    }


def build_normalization_summary(grid):
    """Normalize distance, population, traffic, and airspace factors."""
    graph, edge_distances = build_graph(grid)
    traffic_counts = load_traffic_counts(TRAFFIC_CSV_PATH, grid)

    distance_norm, distance_scale = normalize_values(edge_distances)
    population_values = grid[POPULATION_FIELD].to_numpy(dtype=float)
    population_norm, population_scale = normalize_values(population_values)
    traffic_norm, traffic_scale = normalize_values(
        traffic_counts,
        clip_percentile=TRAFFIC_CLIP_PERCENTILE,
        transform=TRAFFIC_TRANSFORM,
        power=TRAFFIC_POWER,
    )
    airspace_values = grid[AIRSPACE_FIELD].to_numpy(dtype=float)
    airspace_norm, airspace_scale = normalize_values(
        airspace_values,
        clip_percentile=AIRSPACE_CLIP_PERCENTILE,
        power=AIRSPACE_POWER,
    )

    summary_rows = [
        make_summary_row(
            "distance",
            edge_distances,
            distance_norm,
            "max edge distance",
            distance_scale,
        ),
        make_summary_row(
            "population",
            population_values,
            population_norm,
            "max value",
            population_scale,
        ),
        make_summary_row(
            "traffic",
            traffic_counts,
            traffic_norm,
            "log1p percentile clipped value",
            traffic_scale,
            clip_percentile=TRAFFIC_CLIP_PERCENTILE,
            transform=TRAFFIC_TRANSFORM,
            power=TRAFFIC_POWER,
        ),
        make_summary_row(
            "airspace",
            airspace_values,
            airspace_norm,
            "percentile clipped value",
            airspace_scale,
            clip_percentile=AIRSPACE_CLIP_PERCENTILE,
            power=AIRSPACE_POWER,
        ),
    ]

    print("Graph nodes:", graph.number_of_nodes())
    print("Graph edges:", graph.number_of_edges())
    print("Traffic observations matched:", int(np.sum(traffic_counts)))

    return pd.DataFrame(summary_rows)


def current_ratio_baseline():
    """Return the existing project weights normalized to sum to 1."""
    total_weight = (
        CURRENT_DISTANCE_WEIGHT
        + CURRENT_POPULATION_WEIGHT
        + CURRENT_TRAFFIC_WEIGHT
        + CURRENT_AIRSPACE_WEIGHT
    )
    return {
        "distance_weight": CURRENT_DISTANCE_WEIGHT / total_weight,
        "population_weight": CURRENT_POPULATION_WEIGHT / total_weight,
        "traffic_weight": CURRENT_TRAFFIC_WEIGHT / total_weight,
        "airspace_weight": CURRENT_AIRSPACE_WEIGHT / total_weight,
    }


def build_weight_configurations():
    """Generate all valid weight combinations for the configured step."""
    step_units = int(round(1.0 / WEIGHT_STEP))
    baseline = current_ratio_baseline()
    rows = []

    for distance_units in range(step_units + 1):
        for population_units in range(step_units + 1):
            for traffic_units in range(step_units + 1):
                airspace_units = (
                    step_units
                    - distance_units
                    - population_units
                    - traffic_units
                )
                if airspace_units < 0:
                    continue

                distance_weight = distance_units / step_units
                population_weight = population_units / step_units
                traffic_weight = traffic_units / step_units
                airspace_weight = airspace_units / step_units
                weight_sum = (
                    distance_weight
                    + population_weight
                    + traffic_weight
                    + airspace_weight
                )

                rows.append(
                    {
                        "weight_id": f"w_{len(rows) + 1:04d}",
                        "distance_weight": distance_weight,
                        "population_weight": population_weight,
                        "traffic_weight": traffic_weight,
                        "airspace_weight": airspace_weight,
                        "weight_sum": round(weight_sum, 10),
                        "weight_set": "grid_0_10",
                        "weight_step": WEIGHT_STEP,
                        "is_grid_weight": True,
                        "is_distance_only": distance_weight == 1.0,
                        "is_population_only": population_weight == 1.0,
                        "is_traffic_only": traffic_weight == 1.0,
                        "is_airspace_only": airspace_weight == 1.0,
                        "is_equal_weight": False,
                        "is_current_ratio_baseline": (
                            abs(distance_weight - baseline["distance_weight"]) < 1e-9
                            and abs(population_weight - baseline["population_weight"]) < 1e-9
                            and abs(traffic_weight - baseline["traffic_weight"]) < 1e-9
                            and abs(airspace_weight - baseline["airspace_weight"]) < 1e-9
                        ),
                    }
                )

    weights = pd.DataFrame(rows)
    if len(weights) != 286 and WEIGHT_STEP == 0.10:
        raise ValueError(f"Expected 286 weight combinations, found {len(weights)}")

    invalid_sums = weights[~np.isclose(weights["weight_sum"], 1.0)]
    if not invalid_sums.empty:
        raise ValueError("At least one weight combination does not sum to 1.")

    equal_weight_row = {
        "weight_id": f"w_{len(weights) + 1:04d}",
        "distance_weight": 0.25,
        "population_weight": 0.25,
        "traffic_weight": 0.25,
        "airspace_weight": 0.25,
        "weight_sum": 1.0,
        "weight_set": "equal_weight",
        "weight_step": "",
        "is_grid_weight": False,
        "is_distance_only": False,
        "is_population_only": False,
        "is_traffic_only": False,
        "is_airspace_only": False,
        "is_equal_weight": True,
        "is_current_ratio_baseline": False,
    }
    weights = pd.concat(
        [weights, pd.DataFrame([equal_weight_row])],
        ignore_index=True,
    )

    return weights


def main():
    """Create normalized input summaries and weight configurations."""
    if not os.path.exists(GRID_PATH):
        raise FileNotFoundError(f"Routing grid not found: {GRID_PATH}")

    grid = gpd.read_file(GRID_PATH).reset_index(drop=True)
    require_grid_columns(grid)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    normalization_summary = build_normalization_summary(grid)
    weight_configurations = build_weight_configurations()

    normalization_summary.to_csv(NORMALIZATION_SUMMARY_CSV, index=False)
    weight_configurations.to_csv(WEIGHT_CONFIGURATIONS_CSV, index=False)

    print(f"Saved normalization summary: {NORMALIZATION_SUMMARY_CSV}")
    print(f"Saved weight configurations: {WEIGHT_CONFIGURATIONS_CSV}")
    print(f"Weight combinations: {len(weight_configurations)}")


if __name__ == "__main__":
    main()
