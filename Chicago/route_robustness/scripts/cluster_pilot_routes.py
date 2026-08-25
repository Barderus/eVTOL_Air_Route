"""Cluster pilot routes using shared-node similarity."""

import os

import networkx as nx
import pandas as pd


PILOT_ROUTE_RUNS_CSV = "route_robustness/output/pilot_route_runs.csv"
SIMILARITY_PAIRS_CSV = "route_robustness/output/pilot_route_similarity_pairs.csv"
OUTPUT_FOLDER = "route_robustness/output"
PILOT_ROUTE_CLUSTERS_CSV = os.path.join(
    OUTPUT_FOLDER,
    "pilot_route_clusters.csv",
)

SIMILARITY_THRESHOLD = 0.25


def load_inputs():
    """Load pilot route runs and pairwise similarity scores."""
    if not os.path.exists(PILOT_ROUTE_RUNS_CSV):
        raise FileNotFoundError(f"Pilot route runs not found: {PILOT_ROUTE_RUNS_CSV}")
    if not os.path.exists(SIMILARITY_PAIRS_CSV):
        raise FileNotFoundError(f"Similarity pairs not found: {SIMILARITY_PAIRS_CSV}")

    route_runs = pd.read_csv(PILOT_ROUTE_RUNS_CSV)
    pairs = pd.read_csv(SIMILARITY_PAIRS_CSV)

    required_route_columns = {
        "pilot_name",
        "route_run_id",
        "weight_id",
        "route_distance_km",
        "total_weighted_score",
        "status",
    }
    missing_route_columns = required_route_columns.difference(route_runs.columns)
    if missing_route_columns:
        missing_text = ", ".join(sorted(missing_route_columns))
        raise ValueError(f"Pilot route table is missing columns: {missing_text}")

    required_pair_columns = {"route_a", "route_b", "similarity"}
    missing_pair_columns = required_pair_columns.difference(pairs.columns)
    if missing_pair_columns:
        missing_text = ", ".join(sorted(missing_pair_columns))
        raise ValueError(f"Similarity pair table is missing columns: {missing_text}")

    route_runs = route_runs[route_runs["status"] == "success"].copy()
    return route_runs, pairs


def assign_clusters(route_runs, pairs):
    """Assign cluster IDs using connected components above the threshold."""
    graph = nx.Graph()
    graph.add_nodes_from(route_runs["pilot_name"].tolist())

    for row in pairs.itertuples(index=False):
        if float(row.similarity) >= SIMILARITY_THRESHOLD:
            graph.add_edge(row.route_a, row.route_b, similarity=float(row.similarity))

    cluster_rows = []
    components = sorted(
        nx.connected_components(graph),
        key=lambda component: sorted(component)[0],
    )

    for cluster_index, component in enumerate(components, start=1):
        cluster_id = f"cluster_{cluster_index:02d}"
        for pilot_name in sorted(component):
            cluster_rows.append(
                {
                    "pilot_name": pilot_name,
                    "cluster_id": cluster_id,
                    "similarity_threshold": SIMILARITY_THRESHOLD,
                    "cluster_size": len(component),
                }
            )

    clusters = pd.DataFrame(cluster_rows)
    output = route_runs.merge(clusters, on="pilot_name", how="left")
    output = output[
        [
            "cluster_id",
            "cluster_size",
            "pilot_name",
            "route_run_id",
            "weight_id",
            "distance_weight",
            "population_weight",
            "traffic_weight",
            "airspace_weight",
            "route_distance_km",
            "total_weighted_score",
            "path_nodes",
            "similarity_threshold",
        ]
    ].sort_values(
        by=["cluster_id", "pilot_name"],
    )
    return output


def main():
    """Cluster pilot routes and save the route cluster table."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    route_runs, pairs = load_inputs()
    clusters = assign_clusters(route_runs, pairs)
    clusters.to_csv(PILOT_ROUTE_CLUSTERS_CSV, index=False)

    print(f"Saved pilot route clusters: {PILOT_ROUTE_CLUSTERS_CSV}")
    print(clusters[["cluster_id", "cluster_size", "pilot_name"]].to_string(index=False))


if __name__ == "__main__":
    main()
