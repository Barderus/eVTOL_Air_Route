"""Cluster direct routes using shared-node Jaccard similarity."""

import ast
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "Chicago/route_robustness/output/direct_route_weight_configurations.csv"
OUTPUT_FOLDER = "Chicago/route_robustness/output"
OUTPUT_CSV = "Chicago/route_robustness/output/direct_route_clusters_jaccard.csv"

JACCARD_DISTANCE_THRESHOLD = 0.75


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


def connected_components(distance_matrix):
    """Assign connected-component clusters under the distance threshold."""
    route_count = len(distance_matrix)
    visited = [False] * route_count
    labels = [None] * route_count
    components = []

    for start_index in range(route_count):
        if visited[start_index]:
            continue

        stack = [start_index]
        visited[start_index] = True
        component = []

        while stack:
            current_index = stack.pop()
            component.append(current_index)

            for next_index, distance in enumerate(distance_matrix[current_index]):
                if not visited[next_index] and distance <= JACCARD_DISTANCE_THRESHOLD:
                    visited[next_index] = True
                    stack.append(next_index)

        components.append(sorted(component))

    components = sorted(components, key=lambda values: (-len(values), values[0]))
    for cluster_index, component in enumerate(components, start=1):
        label = f"jaccard_{cluster_index:03d}"
        for route_index in component:
            labels[route_index] = label
    return labels


def main():
    """Save Jaccard direct-route clusters."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes = load_routes()
    matrix = build_jaccard_distance_matrix(routes["path_sequence"].tolist())
    routes["jaccard_cluster"] = connected_components(matrix)

    output = routes[
        [
            "route_run_id",
            "weight_id",
            "direct_source_cluster_id",
            "jaccard_cluster",
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

    summary = output.groupby("jaccard_cluster").size().sort_values(ascending=False)
    print(f"Saved Jaccard clusters: {OUTPUT_CSV}")
    print("Cluster count:", len(summary))
    print(summary.head(10).to_string())


if __name__ == "__main__":
    main()
