"""Cluster direct routes using an OPTICS-style Jaccard workflow."""

import ast
import heapq
import math
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "route_robustness/output/direct_route_weight_configurations.csv"
OUTPUT_FOLDER = "route_robustness/output"
OUTPUT_CSV = "route_robustness/output/direct_route_clusters_optics.csv"

OPTICS_EPS = 0.75
OPTICS_MIN_SAMPLES = 4


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


def optics(distance_matrix):
    """Run OPTICS-style ordering and extract clusters at OPTICS_EPS."""
    route_count = len(distance_matrix)
    processed = [False] * route_count
    reachability = [math.inf] * route_count
    ordering = []

    def sorted_neighbors(route_index):
        distances = [
            (distance, index)
            for index, distance in enumerate(distance_matrix[route_index])
            if index != route_index
        ]
        return sorted(distances)

    def core_distance(route_index):
        nearby = sorted_neighbors(route_index)
        if len(nearby) < OPTICS_MIN_SAMPLES - 1:
            return math.inf
        return nearby[OPTICS_MIN_SAMPLES - 2][0]

    core_distances = [core_distance(index) for index in range(route_count)]

    for start_index in range(route_count):
        if processed[start_index]:
            continue

        seeds = [(0.0, start_index)]
        while seeds:
            _, route_index = heapq.heappop(seeds)
            if processed[route_index]:
                continue

            processed[route_index] = True
            ordering.append(route_index)

            route_core_distance = core_distances[route_index]
            if math.isinf(route_core_distance):
                continue

            for distance, neighbor_index in sorted_neighbors(route_index):
                if processed[neighbor_index]:
                    continue

                new_reachability = max(route_core_distance, distance)
                if new_reachability < reachability[neighbor_index]:
                    reachability[neighbor_index] = new_reachability
                    heapq.heappush(seeds, (new_reachability, neighbor_index))

    labels = [None] * route_count
    cluster_number = 0
    noise_number = 0

    for route_index in ordering:
        if reachability[route_index] > OPTICS_EPS:
            if core_distances[route_index] <= OPTICS_EPS:
                cluster_number += 1
                labels[route_index] = f"optics_{cluster_number:03d}"
            else:
                noise_number += 1
                labels[route_index] = f"optics_noise_{noise_number:03d}"
        else:
            if cluster_number == 0:
                cluster_number = 1
            labels[route_index] = f"optics_{cluster_number:03d}"

    return labels


def main():
    """Save OPTICS direct-route clusters."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes = load_routes()
    matrix = build_jaccard_distance_matrix(routes["path_sequence"].tolist())
    routes["optics_cluster"] = optics(matrix)

    output = routes[
        [
            "route_run_id",
            "weight_id",
            "direct_source_cluster_id",
            "optics_cluster",
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

    summary = output.groupby("optics_cluster").size().sort_values(ascending=False)
    print(f"Saved OPTICS clusters: {OUTPUT_CSV}")
    print("Cluster count:", len(summary))
    print(summary.head(10).to_string())


if __name__ == "__main__":
    main()
