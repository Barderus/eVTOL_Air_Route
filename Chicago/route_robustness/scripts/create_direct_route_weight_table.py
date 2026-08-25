"""Create the direct-route weight configuration table."""

import ast
import os

import pandas as pd


ALL_ROUTE_RUNS_CSV = "route_robustness/output/all_route_runs.csv"
ALL_ROUTE_CLUSTERS_CSV = "route_robustness/output/all_route_clusters.csv"
OUTPUT_FOLDER = "route_robustness/output"
OUTPUT_CSV = "route_robustness/output/direct_route_weight_configurations.csv"

DIRECT_SOURCE_CLUSTERS = {
    "cluster_001",
    "cluster_003",
    "cluster_004",
    "cluster_006",
    "cluster_012",
}


def load_direct_routes():
    """Load routes from the selected direct-route source clusters."""
    route_runs = pd.read_csv(ALL_ROUTE_RUNS_CSV)
    clusters = pd.read_csv(ALL_ROUTE_CLUSTERS_CSV)

    direct_clusters = clusters[
        clusters["cluster_id"].isin(DIRECT_SOURCE_CLUSTERS)
    ].copy()
    direct_clusters = direct_clusters.rename(
        columns={
            "cluster_id": "direct_source_cluster_id",
            "cluster_size": "direct_source_cluster_size",
        }
    )

    direct_routes = route_runs.merge(
        direct_clusters[
            [
                "route_run_id",
                "direct_source_cluster_id",
                "direct_source_cluster_size",
            ]
        ],
        on="route_run_id",
        how="inner",
    )
    direct_routes = direct_routes[direct_routes["status"] == "success"].copy()
    direct_routes = direct_routes.sort_values("weight_id").reset_index(drop=True)
    direct_routes["path_sequence"] = direct_routes["path_node_ids"].apply(
        lambda value: ast.literal_eval(str(value))
    )
    return direct_routes


def main():
    """Save the direct-route weight table."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    direct_routes = load_direct_routes()

    output_columns = [
        "route_run_id",
        "weight_id",
        "weight_set",
        "is_equal_weight",
        "direct_source_cluster_id",
        "distance_weight",
        "population_weight",
        "traffic_weight",
        "airspace_weight",
        "route_distance_km",
        "total_weighted_score",
        "path_nodes",
        "path_node_ids",
    ]
    direct_routes[output_columns].to_csv(OUTPUT_CSV, index=False)

    print(f"Saved direct route weight configurations: {OUTPUT_CSV}")
    print("Direct route count:", len(direct_routes))


if __name__ == "__main__":
    main()
