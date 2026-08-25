"""Cluster direct routes using hierarchical clustering with edit distance."""

import ast
import math
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "Chicago/route_robustness/output/direct_route_weight_configurations.csv"
OUTPUT_FOLDER = "Chicago/route_robustness/output"
OUTPUT_CSV = "Chicago/route_robustness/output/direct_route_clusters_edit_distance.csv"

EDIT_DISTANCE_THRESHOLD = 0.35


def load_routes():
    """Load direct routes and parse path node sequences."""
    routes = pd.read_csv(DIRECT_ROUTE_WEIGHTS_CSV)
    routes["path_sequence"] = routes["path_node_ids"].apply(
        lambda value: ast.literal_eval(str(value))
    )
    return routes


def normalized_edit_distance(sequence_a, sequence_b):
    """Calculate Levenshtein edit distance normalized by longer route length."""
    if sequence_a == sequence_b:
        return 0.0

    if len(sequence_a) < len(sequence_b):
        shorter = sequence_a
        longer = sequence_b
    else:
        shorter = sequence_b
        longer = sequence_a

    previous_row = list(range(len(shorter) + 1))
    for longer_index, longer_value in enumerate(longer, start=1):
        current_row = [longer_index]
        for shorter_index, shorter_value in enumerate(shorter, start=1):
            insert_cost = current_row[shorter_index - 1] + 1
            delete_cost = previous_row[shorter_index] + 1
            replace_cost = previous_row[shorter_index - 1]
            if longer_value != shorter_value:
                replace_cost += 1
            current_row.append(min(insert_cost, delete_cost, replace_cost))
        previous_row = current_row

    return previous_row[-1] / max(len(sequence_a), len(sequence_b))


def build_edit_distance_matrix(sequences):
    """Build a normalized route-node edit distance matrix."""
    route_count = len(sequences)
    matrix = [[0.0 for _ in range(route_count)] for _ in range(route_count)]

    for i in range(route_count):
        for j in range(i + 1, route_count):
            distance = normalized_edit_distance(sequences[i], sequences[j])
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

        if best_pair is None or best_distance > EDIT_DISTANCE_THRESHOLD:
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
        label = f"edit_{cluster_index:03d}"
        for route_index in cluster:
            labels[route_index] = label
    return labels


def main():
    """Save edit-distance direct-route clusters."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes = load_routes()
    matrix = build_edit_distance_matrix(routes["path_sequence"].tolist())
    routes["edit_distance_cluster"] = hierarchical_average(matrix)

    output = routes[
        [
            "route_run_id",
            "weight_id",
            "direct_source_cluster_id",
            "edit_distance_cluster",
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

    summary = output.groupby("edit_distance_cluster").size().sort_values(ascending=False)
    print(f"Saved edit-distance clusters: {OUTPUT_CSV}")
    print("Cluster count:", len(summary))
    print(summary.head(10).to_string())


if __name__ == "__main__":
    main()
