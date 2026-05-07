from __future__ import annotations

import argparse
import math
import runpy
from pathlib import Path

import geopandas as gpd
import networkx as nx
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
GENERATE_SCRIPT = ROOT / "2D" / "generate_astar_toggle_pages.py"
DEFAULT_OUTPUT = ROOT / "geojson" / "combined_route_costs.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recompute combined route costs for each destination/date and write them to CSV."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output CSV path.",
    )
    return parser.parse_args()


def load_helpers() -> dict[str, object]:
    return runpy.run_path(str(GENERATE_SCRIPT))


def main() -> None:
    helpers = load_helpers()
    grid_path = helpers["GRID_PATH"]
    destinations = helpers["DESTINATIONS"]
    start = helpers["START"]
    datasets = helpers["DATASETS"]
    normalize = helpers["normalize"]
    build_graph = helpers["build_graph"]
    node_for_point = helpers["node_for_point"]
    load_traffic_counts = helpers["load_traffic_counts"]
    assign_route_edge_weights = helpers["assign_route_edge_weights"]
    make_heuristic = helpers["make_heuristic"]
    path_cost = helpers["path_cost"]
    route_specs = helpers["ROUTE_SPECS"]

    combined_spec = next(spec for spec in route_specs if spec["name"] == "combined")
    grid = gpd.read_file(grid_path)
    graph, centroids_m, cent_x, cent_y = build_graph(grid)
    sindex = grid.sindex
    geoms = grid.geometry

    city_risk = grid["city_risk"].to_numpy(dtype=float)
    airspace_risk = grid["airport_risk_combined"].to_numpy(dtype=float)
    city_risk_norm, _ = normalize(city_risk)
    airspace_risk_norm, _ = normalize(airspace_risk)
    start_node = node_for_point(start["lat"], start["lon"], sindex, geoms, centroids_m)

    rows: list[dict[str, object]] = []
    for dataset in datasets:
        csv_path = dataset["csv_path"]
        if not csv_path.exists():
            raise FileNotFoundError(f"Traffic CSV not found: {csv_path}")

        traffic_counts = load_traffic_counts(csv_path, grid)
        traffic_risk_norm, _ = normalize(traffic_counts)
        max_edge_distance = assign_route_edge_weights(
            graph,
            cent_x,
            cent_y,
            city_risk_norm,
            airspace_risk_norm,
            traffic_risk_norm,
        )
        heuristic = make_heuristic(combined_spec["distance_weight"], cent_x, cent_y, max_edge_distance)

        for destination in destinations:
            end_node = node_for_point(destination["lat"], destination["lon"], sindex, geoms, centroids_m)
            path = nx.astar_path(
                graph,
                start_node,
                end_node,
                heuristic=heuristic,
                weight="weight_combined",
            )
            rows.append(
                {
                    "destination": destination["label"],
                    "date": dataset["slug"],
                    "cost": float(path_cost(graph, path, "weight_combined")),
                }
            )

    output_df = pd.DataFrame(rows, columns=["destination", "date", "cost"])
    output_df = output_df.sort_values(["destination", "date"]).reset_index(drop=True)

    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(args.output, index=False)
    print(f"Saved combined route costs: {args.output}")
    print(f"Rows written: {len(output_df)}")


if __name__ == "__main__":
    main()
