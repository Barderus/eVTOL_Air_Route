"""Cluster direct routes using hierarchical clustering with Jaccard distance."""

import ast
import math
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "route_robustness/output/direct_route_weight_configurations.csv"
OUTPUT_FOLDER = "route_robustness/output"
OUTPUT_CSV = "route_robustness/output/direct_route_clusters_hierarchical_jaccard.csv"

HIERARCHICAL_JACCARD_THRESHOLD = 0.75


def load_routes():
    """Load direct routes and parse path node sequences."""
    routes = pd.read_csv(DIRECT_ROUTE_WEIGHTS_CSV)
    routes["path_sequence"] = routes["path_node_ids"].apply(
        lambda value: ast.literal_eval(str(value))
    )
    return routes


def build_jaccard_distance_matrix(sequences):
    """Build a Jaccard distance matrix from route node sequences."""
    route_sets = [set(sequence) for sequence in sequences]
    route_count = len(route_sets)
    matrix = [[0.0 for _ in range(route_count)] for _ in range(route_count)]

    for i in range(route_count):
        for j in range(i + 1, route_count):
            shared_nodes = len(route_sets[i].intersection(route_sets[j]))
            union_nodes = len(route_sets[i].union(route_sets[j]))
            similarity = shared_nodes / union_nodes if union_nodes else 0.0
            distance = 1.0 - similarity
            matrix[i][j] = distance
            matrix[j][i] = distance
    return matrix


def hierarchical_average(distance_matrix):
    """Run average-linkage hierarchical clustering."""
    clusters = [{index} for index in range(len(distance_matrix))]

    def average_distance(cluster_a, cluster_b):
        distances = [
            distance_matrix[index_a][index_b]
            for index_a in cluster_a
            for index_b in cluster_b
        ]
        return sum(distances) / len(distances)

    while True:
        best_pair = None
        best_distance = math.inf

        for index_a in range(len(clusters)):
            for index_b in range(index_a + 1, len(clusters)):
                distance = average_distance(clusters[index_a], clusters[index_b])
                if distance < best_distance:
                    best_distance = distance
                    best_pair = (index_a, index_b)

        if best_pair is None or best_distance > HIERARCHICAL_JACCARD_THRESHOLD:
            break

        index_a, index_b = best_pair
        clusters[index_a] = clusters[index_a].union(clusters[index_b])
        del clusters[index_b]

    clusters = sorted(
        [sorted(cluster) for cluster in clusters],
        key=lambda values: (-len(values), values[0]),
    )

    labels = [None] * len(distance_matrix)
    for cluster_index, cluster in enumerate(clusters, start=1):
        label = f"hier_jaccard_{cluster_index:03d}"
        for route_index in cluster:
            labels[route_index] = label
    return labels


def main():
    """Save hierarchical Jaccard direct-route clusters."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes = load_routes()
    matrix = build_jaccard_distance_matrix(routes["path_sequence"].tolist())
    routes["hierarchical_jaccard_cluster"] = hierarchical_average(matrix)

    output = routes[
        [
            "route_run_id",
            "weight_id",
            "direct_source_cluster_id",
            "hierarchical_jaccard_cluster",
            "distance_weight",
            "population_weight",
            "traffic_weight",
            "airspace_weight",
            "route_distance_km",
            "total_weighted_score",
            "path_nodes",
        ]
    ]
    output.to_csv(OUTPUT_CSV, index=False)

    summary = output.groupby("hierarchical_jaccard_cluster").size().sort_values(ascending=False)
    print(f"Saved hierarchical Jaccard clusters: {OUTPUT_CSV}")
    print("Cluster count:", len(summary))
    print(summary.head(10).to_string())


if __name__ == "__main__":
    main()
