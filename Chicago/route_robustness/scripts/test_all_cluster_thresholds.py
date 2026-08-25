"""Test alternate full-route clustering thresholds."""

import os

import networkx as nx
import pandas as pd


ALL_ROUTE_RUNS_CSV = "route_robustness/output/all_route_runs.csv"
SIMILARITY_PAIRS_CSV = "route_robustness/output/all_route_similarity_pairs.csv"
OUTPUT_FOLDER = "route_robustness/output"
SUMMARY_CSV = os.path.join(
    OUTPUT_FOLDER,
    "threshold_comparison_summary.csv",
)

THRESHOLDS = [0.20, 0.15, 0.10]


def threshold_label(threshold):
    """Convert a threshold value into a filename-safe label."""
    return f"{int(round(threshold * 100)):03d}"


def load_inputs():
    """Load route runs and pairwise similarity scores."""
    if not os.path.exists(ALL_ROUTE_RUNS_CSV):
        raise FileNotFoundError(f"All route runs not found: {ALL_ROUTE_RUNS_CSV}")
    if not os.path.exists(SIMILARITY_PAIRS_CSV):
        raise FileNotFoundError(f"Similarity pairs not found: {SIMILARITY_PAIRS_CSV}")

    route_runs = pd.read_csv(ALL_ROUTE_RUNS_CSV)
    pairs = pd.read_csv(SIMILARITY_PAIRS_CSV)
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


def assign_clusters(route_runs, pairs, threshold):
    """Assign route clusters for one threshold."""
    graph = nx.Graph()
    graph.add_nodes_from(route_runs["route_run_id"].tolist())

    for row in pairs.itertuples(index=False):
        if float(row.similarity) >= threshold:
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
                    "similarity_threshold": threshold,
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


def summarize_threshold(clusters, threshold):
    """Create a one-row summary for a threshold."""
    cluster_sizes = clusters.groupby("cluster_id").size().sort_values(ascending=False)
    singleton_count = int((cluster_sizes == 1).sum())
    largest_cluster_size = int(cluster_sizes.iloc[0])
    second_largest_cluster_size = int(cluster_sizes.iloc[1]) if len(cluster_sizes) > 1 else 0
    return {
        "threshold": threshold,
        "cluster_count": int(len(cluster_sizes)),
        "largest_cluster_size": largest_cluster_size,
        "largest_cluster_percent": 100.0 * largest_cluster_size / len(clusters),
        "second_largest_cluster_size": second_largest_cluster_size,
        "singleton_cluster_count": singleton_count,
    }


def main():
    """Create alternate cluster outputs and a threshold comparison summary."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    route_runs, pairs = load_inputs()
    summary_rows = []

    for threshold in THRESHOLDS:
        clusters = assign_clusters(route_runs, pairs, threshold)
        output_path = os.path.join(
            OUTPUT_FOLDER,
            f"all_route_clusters_threshold_{threshold_label(threshold)}.csv",
        )
        clusters.to_csv(output_path, index=False)
        summary_rows.append(summarize_threshold(clusters, threshold))
        print(f"Saved threshold clusters: {output_path}")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(SUMMARY_CSV, index=False)
    print(f"Saved threshold comparison summary: {SUMMARY_CSV}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
