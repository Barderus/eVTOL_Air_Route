"""Run route robustness A* routes for a 10-route pilot set.

This script uses the first test case:

- Clow International Airport in Bolingbrook
- Chicago Union Station

It reads `Chicago/route_robustness/output/weight_configurations.csv`, selects 10 pilot
weight configurations, and saves route measurements plus route geometry. The
script is still configurable for another city by editing the settings near the
top.
"""

import json
import math
import os
import time

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point


# Input files
GRID_PATH = "Chicago/geojson/risk_grid_v7.geojson"
TRAFFIC_CSV_PATH = "Chicago/opensky/output/ohare_2026-03-07_1s_15nm_bbox.csv"
WEIGHT_CONFIGURATIONS_CSV = "Chicago/route_robustness/output/weight_configurations.csv"

# Output files
OUTPUT_FOLDER = "Chicago/route_robustness/output"
ROUTE_RUNS_CSV = os.path.join(OUTPUT_FOLDER, "pilot_route_runs.csv")
ROUTES_GEOJSON = os.path.join(OUTPUT_FOLDER, "pilot_routes.geojson")

# Route settings for the first test case
ROUTE_PAIR = "clow_to_chicago_union_station"
ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521
TRAFFIC_DATASET = "ohare_2026-03-07_1s_15nm_bbox"

# Grid field settings. Edit these for another city if its column names differ.
POPULATION_FIELD = "city_risk"
AIRSPACE_FIELD = "airport_risk_combined"

# Normalization settings. Keep these aligned with normalize_route_inputs.py.
AIRSPACE_CLIP_PERCENTILE = 99.0
AIRSPACE_POWER = 0.75
TRAFFIC_CLIP_PERCENTILE = 95.0
TRAFFIC_TRANSFORM = "log1p"
TRAFFIC_POWER = 0.75

PILOT_WEIGHT_CONFIGURATIONS = [
    {
        "pilot_name": "distance_only",
        "distance_weight": 1.0,
        "population_weight": 0.0,
        "traffic_weight": 0.0,
        "airspace_weight": 0.0,
    },
    {
        "pilot_name": "population_only",
        "distance_weight": 0.0,
        "population_weight": 1.0,
        "traffic_weight": 0.0,
        "airspace_weight": 0.0,
    },
    {
        "pilot_name": "traffic_only",
        "distance_weight": 0.0,
        "population_weight": 0.0,
        "traffic_weight": 1.0,
        "airspace_weight": 0.0,
    },
    {
        "pilot_name": "airspace_only",
        "distance_weight": 0.0,
        "population_weight": 0.0,
        "traffic_weight": 0.0,
        "airspace_weight": 1.0,
    },
    {
        "pilot_name": "equal_weight",
        "distance_weight": 0.25,
        "population_weight": 0.25,
        "traffic_weight": 0.25,
        "airspace_weight": 0.25,
    },
    {
        "pilot_name": "distance_heavy",
        "distance_weight": 0.7,
        "population_weight": 0.1,
        "traffic_weight": 0.1,
        "airspace_weight": 0.1,
    },
    {
        "pilot_name": "population_heavy",
        "distance_weight": 0.1,
        "population_weight": 0.7,
        "traffic_weight": 0.1,
        "airspace_weight": 0.1,
    },
    {
        "pilot_name": "traffic_heavy",
        "distance_weight": 0.1,
        "population_weight": 0.1,
        "traffic_weight": 0.7,
        "airspace_weight": 0.1,
    },
    {
        "pilot_name": "airspace_heavy",
        "distance_weight": 0.1,
        "population_weight": 0.1,
        "traffic_weight": 0.1,
        "airspace_weight": 0.7,
    },
    {
        "pilot_name": "mixed_balanced",
        "distance_weight": 0.2,
        "population_weight": 0.2,
        "traffic_weight": 0.3,
        "airspace_weight": 0.3,
    },
]


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
    values = np.asarray(values, dtype=float)
    values = np.clip(values, 0.0, None)

    if transform == "log1p":
        values = np.log1p(values)
    elif transform is not None:
        raise ValueError(f"Unsupported transform: {transform}")

    if clip_percentile is None:
        scale = float(np.max(values))
    else:
        positive_values = values[values > 0]
        if positive_values.size == 0:
            return np.zeros_like(values), 1.0
        scale = float(np.percentile(positive_values, clip_percentile))

    if scale <= 0:
        return np.zeros_like(values), 1.0

    normalized_values = np.clip(values / scale, 0.0, 1.0)
    if power != 1.0:
        normalized_values = np.power(normalized_values, power)

    return normalized_values, scale


def require_input_files():
    """Stop early if the route robustness inputs are missing."""
    required_files = [
        GRID_PATH,
        TRAFFIC_CSV_PATH,
        WEIGHT_CONFIGURATIONS_CSV,
    ]
    for file_path in required_files:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Required input file not found: {file_path}")


def require_grid_columns(grid):
    """Stop early if the configured grid fields are missing."""
    required_columns = {POPULATION_FIELD, AIRSPACE_FIELD}
    missing_columns = required_columns.difference(grid.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Routing grid is missing required columns: {missing_text}")


def load_traffic_counts(csv_path, grid):
    """Count traffic observations inside each grid cell."""
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


def build_graph(grid):
    """Build a graph from touching grid cells and assign normalized edge costs."""
    projected_grid = grid.to_crs("EPSG:3857")
    centroids = projected_grid.geometry.centroid
    centroid_x = centroids.x.to_numpy()
    centroid_y = centroids.y.to_numpy()

    population_values = grid[POPULATION_FIELD].to_numpy(dtype=float)
    population_norm, population_scale = normalize_values(population_values)
    airspace_values = grid[AIRSPACE_FIELD].to_numpy(dtype=float)
    airspace_norm, airspace_scale = normalize_values(
        airspace_values,
        clip_percentile=AIRSPACE_CLIP_PERCENTILE,
        power=AIRSPACE_POWER,
    )
    traffic_counts = load_traffic_counts(TRAFFIC_CSV_PATH, grid)
    traffic_norm, traffic_scale = normalize_values(
        traffic_counts,
        clip_percentile=TRAFFIC_CLIP_PERCENTILE,
        transform=TRAFFIC_TRANSFORM,
        power=TRAFFIC_POWER,
    )

    graph = nx.Graph()
    graph.add_nodes_from(range(len(grid)))

    edge_pairs = []
    edge_distances = []
    spatial_index = grid.sindex
    for cell_index, geometry in enumerate(grid.geometry):
        candidate_indexes = spatial_index.query(geometry, predicate="touches")
        for candidate_index in candidate_indexes:
            candidate_index = int(candidate_index)
            if candidate_index > cell_index:
                distance_m = math.hypot(
                    centroid_x[cell_index] - centroid_x[candidate_index],
                    centroid_y[cell_index] - centroid_y[candidate_index],
                )
                edge_pairs.append((cell_index, candidate_index))
                edge_distances.append(distance_m)

    max_edge_distance = max(edge_distances)

    for (start_node, end_node), distance_m in zip(edge_pairs, edge_distances):
        distance_norm = distance_m / max_edge_distance
        population_edge_norm = (
            population_norm[start_node] + population_norm[end_node]
        ) / 2.0
        traffic_edge_norm = (
            traffic_norm[start_node] + traffic_norm[end_node]
        ) / 2.0
        airspace_edge_norm = (
            airspace_norm[start_node] + airspace_norm[end_node]
        ) / 2.0

        graph.add_edge(
            start_node,
            end_node,
            distance_m=distance_m,
            distance_norm=distance_norm,
            population_norm=population_edge_norm,
            traffic_norm=traffic_edge_norm,
            airspace_norm=airspace_edge_norm,
        )

    normalization_scales = {
        "distance_scale": max_edge_distance,
        "population_scale": population_scale,
        "traffic_scale": traffic_scale,
        "airspace_scale": airspace_scale,
        "traffic_observations_matched": int(np.sum(traffic_counts)),
    }

    return graph, projected_grid, centroid_x, centroid_y, normalization_scales


def find_grid_cell(latitude, longitude, grid, projected_grid):
    """Find the grid cell containing or nearest to a point."""
    point = Point(longitude, latitude)
    containing_cells = list(grid.sindex.query(point, predicate="within"))
    if containing_cells:
        return int(containing_cells[0])

    projected_point = gpd.GeoSeries([point], crs=grid.crs).to_crs(
        projected_grid.crs
    ).iloc[0]
    distances = projected_grid.geometry.centroid.distance(projected_point)
    return int(distances.idxmin())


def make_heuristic(distance_weight, centroid_x, centroid_y, max_edge_distance):
    """Create an admissible A* heuristic from normalized straight-line distance."""
    if distance_weight <= 0:
        return lambda _start_node, _end_node: 0.0

    def heuristic(start_node, end_node):
        distance_m = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        return distance_weight * (distance_m / max_edge_distance)

    return heuristic


def make_weight_function(weights):
    """Create a NetworkX weight function for one weight configuration."""
    distance_weight = float(weights["distance_weight"])
    population_weight = float(weights["population_weight"])
    traffic_weight = float(weights["traffic_weight"])
    airspace_weight = float(weights["airspace_weight"])

    def edge_weight(_start_node, _end_node, edge_data):
        return (
            distance_weight * edge_data["distance_norm"]
            + population_weight * edge_data["population_norm"]
            + traffic_weight * edge_data["traffic_norm"]
            + airspace_weight * edge_data["airspace_norm"]
        )

    return edge_weight


def calculate_path_metrics(graph, path, weights):
    """Calculate route distance, normalized factor sums, and weighted score."""
    distance_m = 0.0
    distance_norm_sum = 0.0
    population_norm_sum = 0.0
    traffic_norm_sum = 0.0
    airspace_norm_sum = 0.0

    for start_node, end_node in zip(path, path[1:]):
        edge_data = graph.edges[start_node, end_node]
        distance_m += edge_data["distance_m"]
        distance_norm_sum += edge_data["distance_norm"]
        population_norm_sum += edge_data["population_norm"]
        traffic_norm_sum += edge_data["traffic_norm"]
        airspace_norm_sum += edge_data["airspace_norm"]

    weighted_score = (
        float(weights["distance_weight"]) * distance_norm_sum
        + float(weights["population_weight"]) * population_norm_sum
        + float(weights["traffic_weight"]) * traffic_norm_sum
        + float(weights["airspace_weight"]) * airspace_norm_sum
    )

    return {
        "route_distance_km": distance_m / 1000.0,
        "distance_norm_sum": distance_norm_sum,
        "population_norm_sum": population_norm_sum,
        "traffic_norm_sum": traffic_norm_sum,
        "airspace_norm_sum": airspace_norm_sum,
        "total_weighted_score": weighted_score,
    }


def make_route_geometry(projected_grid, path):
    """Create one route line from selected grid-cell centroids."""
    route_points = [projected_grid.geometry.iloc[index].centroid for index in path]
    route_line = gpd.GeoSeries(
        [LineString(route_points)],
        crs=projected_grid.crs,
    ).to_crs("EPSG:4326")
    return route_line.iloc[0]


def load_weight_configurations():
    """Load and validate normalized weight configurations."""
    weights = pd.read_csv(WEIGHT_CONFIGURATIONS_CSV)
    required_columns = {
        "weight_id",
        "distance_weight",
        "population_weight",
        "traffic_weight",
        "airspace_weight",
        "weight_sum",
    }
    missing_columns = required_columns.difference(weights.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Weight table is missing columns: {missing_text}")

    invalid_weights = weights[
        (weights["distance_weight"] < 0)
        | (weights["distance_weight"] > 1)
        | (weights["population_weight"] < 0)
        | (weights["population_weight"] > 1)
        | (weights["traffic_weight"] < 0)
        | (weights["traffic_weight"] > 1)
        | (weights["airspace_weight"] < 0)
        | (weights["airspace_weight"] > 1)
        | (~np.isclose(weights["weight_sum"], 1.0))
    ]
    if not invalid_weights.empty:
        raise ValueError("At least one weight row is outside the valid range.")

    return weights


def select_pilot_weights(weights):
    """Select the configured 10-route pilot set from the weight table."""
    selected_rows = []

    for pilot_configuration in PILOT_WEIGHT_CONFIGURATIONS:
        matching_rows = weights[
            np.isclose(
                weights["distance_weight"],
                pilot_configuration["distance_weight"],
            )
            & np.isclose(
                weights["population_weight"],
                pilot_configuration["population_weight"],
            )
            & np.isclose(
                weights["traffic_weight"],
                pilot_configuration["traffic_weight"],
            )
            & np.isclose(
                weights["airspace_weight"],
                pilot_configuration["airspace_weight"],
            )
        ]

        if matching_rows.empty:
            raise ValueError(
                "Pilot weight configuration was not found in the weight table: "
                f"{pilot_configuration['pilot_name']}"
            )

        selected_row = matching_rows.iloc[0].copy()
        selected_row["pilot_name"] = pilot_configuration["pilot_name"]
        selected_rows.append(selected_row)

    return pd.DataFrame(selected_rows).reset_index(drop=True)


def main():
    """Run every weighted route and save route-run artifacts."""
    require_input_files()
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    grid = gpd.read_file(GRID_PATH).reset_index(drop=True)
    require_grid_columns(grid)

    weights = select_pilot_weights(load_weight_configurations())
    graph, projected_grid, centroid_x, centroid_y, normalization_scales = build_graph(grid)
    start_node = find_grid_cell(ORIGIN_LAT, ORIGIN_LON, grid, projected_grid)
    end_node = find_grid_cell(DESTINATION_LAT, DESTINATION_LON, grid, projected_grid)

    print("Graph nodes:", graph.number_of_nodes())
    print("Graph edges:", graph.number_of_edges())
    print("Start node:", start_node)
    print("End node:", end_node)
    print("Pilot weight configurations:", len(weights))

    route_rows = []
    route_geometries = []

    for row_index, weights_row in weights.iterrows():
        route_run_id = f"rr_{row_index + 1:04d}"
        start_time = time.perf_counter()
        heuristic = make_heuristic(
            float(weights_row["distance_weight"]),
            centroid_x,
            centroid_y,
            normalization_scales["distance_scale"],
        )
        weight_function = make_weight_function(weights_row)

        status = "success"
        failure_reason = ""
        path = []
        geometry = None
        metrics = {
            "route_distance_km": np.nan,
            "distance_norm_sum": np.nan,
            "population_norm_sum": np.nan,
            "traffic_norm_sum": np.nan,
            "airspace_norm_sum": np.nan,
            "total_weighted_score": np.nan,
        }

        try:
            path = nx.astar_path(
                graph,
                start_node,
                end_node,
                heuristic=heuristic,
                weight=weight_function,
            )
            metrics = calculate_path_metrics(graph, path, weights_row)
            geometry = make_route_geometry(projected_grid, path)
        except nx.NetworkXNoPath as error:
            status = "failed"
            failure_reason = str(error)

        elapsed_seconds = time.perf_counter() - start_time

        route_row = {
            "route_run_id": route_run_id,
            "weight_id": weights_row["weight_id"],
            "pilot_name": weights_row["pilot_name"],
            "weight_set": weights_row.get("weight_set", ""),
            "is_grid_weight": weights_row.get("is_grid_weight", ""),
            "is_equal_weight": weights_row.get("is_equal_weight", ""),
            "route_pair": ROUTE_PAIR,
            "origin_label": ORIGIN_LABEL,
            "destination_label": DESTINATION_LABEL,
            "origin_lat": ORIGIN_LAT,
            "origin_lon": ORIGIN_LON,
            "destination_lat": DESTINATION_LAT,
            "destination_lon": DESTINATION_LON,
            "traffic_dataset": TRAFFIC_DATASET,
            "distance_weight": float(weights_row["distance_weight"]),
            "population_weight": float(weights_row["population_weight"]),
            "traffic_weight": float(weights_row["traffic_weight"]),
            "airspace_weight": float(weights_row["airspace_weight"]),
            "path_nodes": len(path),
            "path_node_ids": json.dumps(path),
            "elapsed_seconds": elapsed_seconds,
            "status": status,
            "failure_reason": failure_reason,
        }
        route_row.update(metrics)
        route_rows.append(route_row)

        if geometry is not None:
            route_geometries.append(geometry)
        else:
            route_geometries.append(LineString())

        print(
            route_run_id,
            weights_row["weight_id"],
            status,
            f"score={metrics['total_weighted_score']:.4f}",
            f"seconds={elapsed_seconds:.3f}",
        )

    route_runs = pd.DataFrame(route_rows)
    routes = gpd.GeoDataFrame(route_runs, geometry=route_geometries, crs="EPSG:4326")

    route_runs.to_csv(ROUTE_RUNS_CSV, index=False)
    routes.to_file(ROUTES_GEOJSON, driver="GeoJSON")

    print(f"Saved route runs: {ROUTE_RUNS_CSV}")
    print(f"Saved route geometries: {ROUTES_GEOJSON}")


if __name__ == "__main__":
    main()
