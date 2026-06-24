"""Cluster direct routes using DBSCAN with Jaccard distance."""

import ast
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "route_robustness/output/direct_route_weight_configurations.csv"
OUTPUT_FOLDER = "route_robustness/output"
OUTPUT_CSV = "route_robustness/output/direct_route_clusters_dbscan.csv"

DBSCAN_EPS = 0.75
DBSCAN_MIN_SAMPLES = 4


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


def dbscan(distance_matrix):
    """Run DBSCAN with a precomputed distance matrix."""
    route_count = len(distance_matrix)
    labels = [None] * route_count
    visited = [False] * route_count
    cluster_number = 0
    noise_number = 0

    def neighbors(route_index):
        return [
            index
            for index, distance in enumerate(distance_matrix[route_index])
            if distance <= DBSCAN_EPS
        ]

    for route_index in range(route_count):
        if visited[route_index]:
            continue

        visited[route_index] = True
        route_neighbors = neighbors(route_index)

        if len(route_neighbors) < DBSCAN_MIN_SAMPLES:
            noise_number += 1
            labels[route_index] = f"dbscan_noise_{noise_number:03d}"
            continue

        cluster_number += 1
        cluster_label = f"dbscan_{cluster_number:03d}"
        labels[route_index] = cluster_label

        seeds = list(route_neighbors)
        seed_position = 0
        while seed_position < len(seeds):
            neighbor_index = seeds[seed_position]

            if not visited[neighbor_index]:
                visited[neighbor_index] = True
                neighbor_neighbors = neighbors(neighbor_index)
                if len(neighbor_neighbors) >= DBSCAN_MIN_SAMPLES:
                    for candidate in neighbor_neighbors:
                        if candidate not in seeds:
                            seeds.append(candidate)

            if labels[neighbor_index] is None or "_noise_" in labels[neighbor_index]:
                labels[neighbor_index] = cluster_label

            seed_position += 1

    return labels


def main():
    """Save DBSCAN direct-route clusters."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes = load_routes()
    matrix = build_jaccard_distance_matrix(routes["path_sequence"].tolist())
    routes["dbscan_cluster"] = dbscan(matrix)

    output = routes[
        [
            "route_run_id",
            "weight_id",
            "direct_source_cluster_id",
            "dbscan_cluster",
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

    summary = output.groupby("dbscan_cluster").size().sort_values(ascending=False)
    print(f"Saved DBSCAN clusters: {OUTPUT_CSV}")
    print("Cluster count:", len(summary))
    print(summary.head(10).to_string())


if __name__ == "__main__":
    main()
