"""Run St. Louis route robustness A* routes for normalized weight grids."""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString


REPO_ROOT = Path(__file__).resolve().parents[3]
ST_LOUIS_DIR = REPO_ROOT / "St Louis"
sys.path.insert(0, str(ST_LOUIS_DIR))
sys.path.insert(0, str(REPO_ROOT))

import astar_routes  # noqa: E402
import settings  # noqa: E402
from api_tool.generate_weight_configurations import generate_weight_configurations  # noqa: E402


OUTPUT_FOLDER = ST_LOUIS_DIR / "route_robustness" / "output"
ROUTE_RUNS_CSV = OUTPUT_FOLDER / "st_louis_weighted_route_runs.csv"
ROUTES_GEOJSON = OUTPUT_FOLDER / "st_louis_weighted_routes.geojson"
WEIGHT_CONFIGURATIONS_CSV = OUTPUT_FOLDER / "st_louis_weight_configurations.csv"


def slugify(value):
    """Convert a display label into a stable file-friendly identifier."""
    return (
        value.lower()
        .replace("st. ", "st_")
        .replace(".", "")
        .replace("&", "and")
        .replace("/", "_")
        .replace("-", "_")
        .replace(" ", "_")
    )


def make_weight_function(weights):
    """Create a NetworkX edge-weight function for one normalized weight row."""
    distance_weight = float(weights["distance_weight"])
    population_weight = float(weights["population_weight"])
    traffic_weight = float(weights["traffic_weight"])
    airspace_weight = float(weights["airspace_weight"])

    def edge_weight(_start_node, _end_node, edge_data):
        return (
            distance_weight * edge_data["distance_cost"]
            + population_weight * edge_data["population_cost"]
            + traffic_weight * edge_data["traffic_cost"]
            + airspace_weight * edge_data["airspace_cost"]
        )

    return edge_weight


def make_heuristic(distance_weight, centroid_x, centroid_y, maximum_distance):
    """Create an admissible straight-line-distance heuristic for A*."""
    if distance_weight <= 0:
        return lambda _start_node, _end_node: 0.0

    def heuristic(start_node, end_node):
        distance_m = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        return distance_weight * distance_m / maximum_distance

    return heuristic


def add_normalized_edge_costs(
    graph,
    centroid_x,
    centroid_y,
    population_risk,
    airspace_risk,
    traffic_risk,
):
    """Store individual normalized edge costs used by weighted experiments."""
    edge_distances = []
    for start_node, end_node in graph.edges:
        edge_distances.append(
            math.hypot(
                centroid_x[start_node] - centroid_x[end_node],
                centroid_y[start_node] - centroid_y[end_node],
            )
        )

    maximum_distance = max(edge_distances)

    for start_node, end_node in graph.edges:
        distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        graph.edges[start_node, end_node]["distance_m"] = distance
        graph.edges[start_node, end_node]["distance_cost"] = distance / maximum_distance
        graph.edges[start_node, end_node]["population_cost"] = (
            population_risk[start_node] + population_risk[end_node]
        ) / 2.0
        graph.edges[start_node, end_node]["airspace_cost"] = (
            airspace_risk[start_node] + airspace_risk[end_node]
        ) / 2.0
        graph.edges[start_node, end_node]["traffic_cost"] = (
            traffic_risk[start_node] + traffic_risk[end_node]
        ) / 2.0

    return maximum_distance


def calculate_path_metrics(graph, path, weights):
    """Calculate distance, factor sums, and weighted route score."""
    distance_m = 0.0
    distance_cost_sum = 0.0
    population_cost_sum = 0.0
    traffic_cost_sum = 0.0
    airspace_cost_sum = 0.0

    for start_node, end_node in zip(path, path[1:]):
        edge_data = graph.edges[start_node, end_node]
        distance_m += edge_data["distance_m"]
        distance_cost_sum += edge_data["distance_cost"]
        population_cost_sum += edge_data["population_cost"]
        traffic_cost_sum += edge_data["traffic_cost"]
        airspace_cost_sum += edge_data["airspace_cost"]

    total_weighted_score = (
        float(weights["distance_weight"]) * distance_cost_sum
        + float(weights["population_weight"]) * population_cost_sum
        + float(weights["traffic_weight"]) * traffic_cost_sum
        + float(weights["airspace_weight"]) * airspace_cost_sum
    )

    return {
        "route_distance_km": distance_m / 1000.0,
        "distance_cost_sum": distance_cost_sum,
        "population_cost_sum": population_cost_sum,
        "traffic_cost_sum": traffic_cost_sum,
        "airspace_cost_sum": airspace_cost_sum,
        "total_weighted_score": total_weighted_score,
    }


def load_weights(step):
    """Generate and save St. Louis normalized weight configurations."""
    result = generate_weight_configurations(
        step=step,
        output_path=str(WEIGHT_CONFIGURATIONS_CSV),
    )
    return pd.DataFrame(result["weight_configurations"])


def select_traffic_datasets(traffic_date):
    """Select one configured St. Louis traffic dataset by date."""
    if traffic_date is None:
        return [settings.TRAFFIC_DATASETS[0]]

    datasets = [
        dataset for dataset in settings.TRAFFIC_DATASETS if dataset["date"] == traffic_date
    ]
    if not datasets:
        valid_dates = ", ".join(dataset["date"] for dataset in settings.TRAFFIC_DATASETS)
        raise ValueError(f"Unknown traffic date {traffic_date}. Valid dates: {valid_dates}")
    return datasets


def main():
    """Run every configured route pair over the St. Louis weight grid."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--traffic-date",
        default=None,
        help="Configured St. Louis traffic date to use. Defaults to the first configured date.",
    )
    parser.add_argument(
        "--step",
        type=float,
        default=0.10,
        help="Weight-grid step size. The Chicago baseline uses 0.10.",
    )
    args = parser.parse_args()

    settings.validate_route_settings()
    os.chdir(ST_LOUIS_DIR)
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(settings.RISK_GRID_GEOJSON):
        raise FileNotFoundError(f"St. Louis risk grid not found: {settings.RISK_GRID_GEOJSON}")

    grid = gpd.read_file(settings.RISK_GRID_GEOJSON).reset_index(drop=True)
    required_columns = {"city_risk", "airport_risk_combined"}
    missing_columns = required_columns.difference(grid.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"St. Louis risk grid is missing columns: {missing_text}")

    graph, projected_grid, centroid_x, centroid_y = astar_routes.build_graph(grid)
    population_risk = astar_routes.normalize(grid["city_risk"])
    airspace_risk = astar_routes.normalize(grid["airport_risk_combined"])
    weights = load_weights(args.step)
    traffic_datasets = select_traffic_datasets(args.traffic_date)

    route_rows = []
    route_geometries = []

    print("Graph nodes:", graph.number_of_nodes())
    print("Graph edges:", graph.number_of_edges())
    print("Route pairs:", len(settings.ROUTES))
    print("Weight configurations:", len(weights))

    for dataset in traffic_datasets:
        csv_path = os.path.join(settings.TRAFFIC_FOLDER, dataset["filename"])
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Traffic file not found: {csv_path}")

        traffic_counts = astar_routes.load_traffic_counts(csv_path, grid)
        traffic_risk = astar_routes.normalize(
            traffic_counts,
            clip_percentile=95.0,
            transform="log1p",
            power=0.75,
        )
        maximum_distance = add_normalized_edge_costs(
            graph,
            centroid_x,
            centroid_y,
            population_risk,
            airspace_risk,
            traffic_risk,
        )

        for route in settings.ROUTES:
            start_latitude, start_longitude = route["start"]
            destination_latitude, destination_longitude = route["destination"]
            start_node = astar_routes.find_grid_cell(
                start_latitude,
                start_longitude,
                grid,
                projected_grid,
            )
            end_node = astar_routes.find_grid_cell(
                destination_latitude,
                destination_longitude,
                grid,
                projected_grid,
            )
            route_pair = slugify(route["label"])
            route_pair_label = route["label"]

            for row_index, weights_row in weights.iterrows():
                route_run_id = f"stl_{route_pair}_{dataset['date']}_{row_index + 1:04d}"
                start_time = time.perf_counter()
                path = []
                geometry = LineString()
                status = "success"
                failure_reason = ""
                metrics = {
                    "route_distance_km": np.nan,
                    "distance_cost_sum": np.nan,
                    "population_cost_sum": np.nan,
                    "traffic_cost_sum": np.nan,
                    "airspace_cost_sum": np.nan,
                    "total_weighted_score": np.nan,
                }

                try:
                    path = nx.astar_path(
                        graph,
                        start_node,
                        end_node,
                        heuristic=make_heuristic(
                            float(weights_row["distance_weight"]),
                            centroid_x,
                            centroid_y,
                            maximum_distance,
                        ),
                        weight=make_weight_function(weights_row),
                    )
                    metrics = calculate_path_metrics(graph, path, weights_row)
                    geometry = astar_routes.make_route_geometry(projected_grid, path)
                except nx.NetworkXNoPath as error:
                    status = "failed"
                    failure_reason = str(error)

                elapsed_seconds = time.perf_counter() - start_time
                route_row = {
                    "route_run_id": route_run_id,
                    "route_pair": route_pair,
                    "route_pair_label": route_pair_label,
                    "origin_label": route["start_label"],
                    "destination_label": route["destination_label"],
                    "origin_lat": start_latitude,
                    "origin_lon": start_longitude,
                    "destination_lat": destination_latitude,
                    "destination_lon": destination_longitude,
                    "traffic_date": dataset["date"],
                    "traffic_dataset": dataset["label"],
                    "weight_id": weights_row["weight_id"],
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
                route_geometries.append(geometry)

            print(f"Completed {route_pair_label} for {dataset['date']}")

    route_runs = pd.DataFrame(route_rows)
    routes = gpd.GeoDataFrame(route_runs, geometry=route_geometries, crs="EPSG:4326")
    route_runs.to_csv(ROUTE_RUNS_CSV, index=False)
    routes.to_file(ROUTES_GEOJSON, driver="GeoJSON")

    print(f"Saved weight configurations: {WEIGHT_CONFIGURATIONS_CSV}")
    print(f"Saved weighted route runs: {ROUTE_RUNS_CSV}")
    print(f"Saved weighted route geometries: {ROUTES_GEOJSON}")


if __name__ == "__main__":
    main()
