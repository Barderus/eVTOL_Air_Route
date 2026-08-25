"""Cluster all weighted routes using shared-node similarity."""

import os

import networkx as nx
import pandas as pd


ALL_ROUTE_RUNS_CSV = "Chicago/route_robustness/output/all_route_runs.csv"
SIMILARITY_PAIRS_CSV = "Chicago/route_robustness/output/all_route_similarity_pairs.csv"
OUTPUT_FOLDER = "Chicago/route_robustness/output"
ALL_ROUTE_CLUSTERS_CSV = os.path.join(
    OUTPUT_FOLDER,
    "all_route_clusters.csv",
)

SIMILARITY_THRESHOLD = 0.25


def load_inputs():
    """Load full route runs and pairwise similarity scores."""
    if not os.path.exists(ALL_ROUTE_RUNS_CSV):
        raise FileNotFoundError(f"All route runs not found: {ALL_ROUTE_RUNS_CSV}")
    if not os.path.exists(SIMILARITY_PAIRS_CSV):
        raise FileNotFoundError(f"Similarity pairs not found: {SIMILARITY_PAIRS_CSV}")

    route_runs = pd.read_csv(ALL_ROUTE_RUNS_CSV)
    pairs = pd.read_csv(SIMILARITY_PAIRS_CSV)

    required_route_columns = {
        "route_run_id",
        "weight_id",
        "weight_set",
        "distance_weight",
        "population_weight",
        "traffic_weight",
        "airspace_weight",
        "route_distance_km",
        "total_weighted_score",
        "path_nodes",
        "status",
    }
    missing_route_columns = required_route_columns.difference(route_runs.columns)
    if missing_route_columns:
        missing_text = ", ".join(sorted(missing_route_columns))
        raise ValueError(f"All-route table is missing columns: {missing_text}")

    required_pair_columns = {"route_run_a", "route_run_b", "similarity"}
    missing_pair_columns = required_pair_columns.difference(pairs.columns)
    if missing_pair_columns:
        missing_text = ", ".join(sorted(missing_pair_columns))
        raise ValueError(f"Similarity pair table is missing columns: {missing_text}")

    route_runs = route_runs[route_runs["status"] == "success"].copy()
    return route_runs, pairs


def dominant_weight_category(row):
    """Classify a route by its largest weight value."""
    if str(row.get("is_equal_weight", "")).lower() == "true":
        return "equal_weight"

    weights = {
        "distance_dominant": float(row["distance_weight"]),
        "population_dominant": float(row["population_weight"]),
        "traffic_dominant": float(row["traffic_weight"]),
        "airspace_dominant": float(row["airspace_weight"]),
    }
    maximum_weight = max(weights.values())
    matching_categories = [
        category
        for category, value in weights.items()
        if value == maximum_weight
    ]
    if len(matching_categories) != 1:
        return "tied_dominant"
    return matching_categories[0]


def assign_clusters(route_runs, pairs):
    """Assign cluster IDs using connected components above the threshold."""
    graph = nx.Graph()
    graph.add_nodes_from(route_runs["route_run_id"].tolist())

    for row in pairs.itertuples(index=False):
        if float(row.similarity) >= SIMILARITY_THRESHOLD:
            graph.add_edge(
                row.route_run_a,
                row.route_run_b,
                similarity=float(row.similarity),
            )

    cluster_rows = []
    components = sorted(
        nx.connected_components(graph),
        key=lambda component: (-len(component), sorted(component)[0]),
    )

    for cluster_index, component in enumerate(components, start=1):
        cluster_id = f"cluster_{cluster_index:03d}"
        for route_run_id in sorted(component):
            cluster_rows.append(
                {
                    "route_run_id": route_run_id,
                    "cluster_id": cluster_id,
                    "similarity_threshold": SIMILARITY_THRESHOLD,
                    "cluster_size": len(component),
                }
            )

    clusters = pd.DataFrame(cluster_rows)
    output = route_runs.merge(clusters, on="route_run_id", how="left")
    output["dominant_weight_category"] = output.apply(
        dominant_weight_category,
        axis=1,
    )
    output = output[
        [
            "cluster_id",
            "cluster_size",
            "route_run_id",
            "weight_id",
            "weight_set",
            "dominant_weight_category",
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
        by=["cluster_size", "cluster_id", "weight_id"],
        ascending=[False, True, True],
    )
    return output


def main():
    """Cluster all weighted routes and save the cluster table."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    route_runs, pairs = load_inputs()
    clusters = assign_clusters(route_runs, pairs)
    clusters.to_csv(ALL_ROUTE_CLUSTERS_CSV, index=False)

    summary = (
        clusters.groupby("cluster_id")
        .agg(
            cluster_size=("route_run_id", "count"),
            min_distance_km=("route_distance_km", "min"),
            max_distance_km=("route_distance_km", "max"),
            mean_score=("total_weighted_score", "mean"),
        )
        .sort_values("cluster_size", ascending=False)
        .reset_index()
    )

    print(f"Saved all-route clusters: {ALL_ROUTE_CLUSTERS_CSV}")
    print("Cluster count:", len(summary))
    print(summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
