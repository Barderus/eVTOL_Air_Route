"""Cluster Chicago weighted routes with multiple route-distance methods."""

import ast
import json
import math
from pathlib import Path

import pandas as pd


OUTPUT_FOLDER = Path("Chicago") / "route_robustness" / "output"
ROUTE_RUNS_CSV = OUTPUT_FOLDER / "additional_route_runs.csv"
ROUTES_GEOJSON = OUTPUT_FOLDER / "additional_routes.geojson"

DBSCAN_OUTPUT_CSV = OUTPUT_FOLDER / "additional_route_clusters_dbscan.csv"
EDIT_OUTPUT_CSV = OUTPUT_FOLDER / "additional_route_clusters_edit_distance.csv"
FRECHET_OUTPUT_CSV = OUTPUT_FOLDER / "additional_route_clusters_frechet.csv"
JACCARD_OUTPUT_CSV = OUTPUT_FOLDER / "additional_route_clusters_hierarchical_jaccard.csv"

DBSCAN_EPS_KM = 2.5
DBSCAN_MIN_SAMPLES = 4
EDIT_DISTANCE_THRESHOLD = 0.35
FRECHET_THRESHOLD_KM = 2.5
HIERARCHICAL_JACCARD_THRESHOLD = 0.75
MAX_FRECHET_POINTS = 80

BASE_OUTPUT_COLUMNS = [
    "route_run_id",
    "route_pair",
    "route_pair_label",
    "weight_id",
    "distance_weight",
    "population_weight",
    "traffic_weight",
    "airspace_weight",
    "route_distance_km",
    "total_weighted_score",
    "path_nodes",
]


def load_routes():
    """Load successful weighted route rows and parsed path sequences."""
    routes = pd.read_csv(ROUTE_RUNS_CSV)
    routes = routes[routes["status"] == "success"].copy()
    if routes.empty:
        raise ValueError("No successful Chicago weighted routes found.")
    routes["path_sequence"] = routes["path_node_ids"].apply(
        lambda value: ast.literal_eval(str(value))
    )
    routes["path_signature"] = routes["path_node_ids"].astype(str)
    return routes


def load_route_features():
    """Load route GeoJSON features by route_run_id."""
    with ROUTES_GEOJSON.open("r", encoding="utf-8") as file_handle:
        geojson = json.load(file_handle)

    return {
        feature["properties"]["route_run_id"]: feature
        for feature in geojson["features"]
    }


def unique_routes_for_pair(routes, route_pair):
    """Return one representative row for each exact path in one route pair."""
    subset = routes[routes["route_pair"] == route_pair].copy()
    return subset.drop_duplicates("path_signature").reset_index(drop=True)


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


def build_jaccard_distance_matrix(sequences):
    """Build a Jaccard distance matrix from route node sets."""
    route_sets = [set(sequence) for sequence in sequences]
    route_count = len(route_sets)
    matrix = [[0.0 for _ in range(route_count)] for _ in range(route_count)]

    for i in range(route_count):
        for j in range(i + 1, route_count):
            shared_nodes = len(route_sets[i].intersection(route_sets[j]))
            union_nodes = len(route_sets[i].union(route_sets[j]))
            similarity = shared_nodes / union_nodes if union_nodes else 0.0
            matrix[i][j] = 1.0 - similarity
            matrix[j][i] = matrix[i][j]
    return matrix


def lonlat_to_km(coordinates):
    """Convert lon/lat coordinates to approximate local kilometers."""
    latitudes = [lat for _, lat in coordinates]
    center_lat = math.radians(sum(latitudes) / len(latitudes))
    return [
        (lon * 111.320 * math.cos(center_lat), lat * 110.574)
        for lon, lat in coordinates
    ]


def resample_points(points):
    """Resample a route geometry to a smaller fixed number of points."""
    if len(points) <= MAX_FRECHET_POINTS:
        return points

    distances = [0.0]
    for index in range(1, len(points)):
        previous_x, previous_y = points[index - 1]
        current_x, current_y = points[index]
        distances.append(distances[-1] + math.hypot(current_x - previous_x, current_y - previous_y))

    total_distance = distances[-1]
    if total_distance == 0:
        return [points[0]]

    resampled = []
    segment_index = 1
    for target_index in range(MAX_FRECHET_POINTS):
        target_distance = total_distance * target_index / (MAX_FRECHET_POINTS - 1)
        while segment_index < len(distances) - 1 and distances[segment_index] < target_distance:
            segment_index += 1

        previous_distance = distances[segment_index - 1]
        current_distance = distances[segment_index]
        previous_point = points[segment_index - 1]
        current_point = points[segment_index]

        if current_distance == previous_distance:
            resampled.append(current_point)
            continue

        ratio = (target_distance - previous_distance) / (current_distance - previous_distance)
        resampled.append(
            (
                previous_point[0] + ratio * (current_point[0] - previous_point[0]),
                previous_point[1] + ratio * (current_point[1] - previous_point[1]),
            )
        )

    return resampled


def discrete_frechet_distance(points_a, points_b):
    """Calculate discrete Frechet distance between two point sequences."""
    previous_row = [math.inf] * len(points_b)

    for index_a, point_a in enumerate(points_a):
        current_row = [math.inf] * len(points_b)
        for index_b, point_b in enumerate(points_b):
            distance = math.hypot(point_a[0] - point_b[0], point_a[1] - point_b[1])
            if index_a == 0 and index_b == 0:
                current_row[index_b] = distance
            elif index_a == 0:
                current_row[index_b] = max(current_row[index_b - 1], distance)
            elif index_b == 0:
                current_row[index_b] = max(previous_row[index_b], distance)
            else:
                current_row[index_b] = max(
                    min(
                        previous_row[index_b],
                        previous_row[index_b - 1],
                        current_row[index_b - 1],
                    ),
                    distance,
                )
        previous_row = current_row

    return previous_row[-1]


def build_frechet_distance_matrix(unique_routes, route_features):
    """Build a discrete Frechet distance matrix from route geometries."""
    route_points = []
    for route_id in unique_routes["route_run_id"].tolist():
        coordinates = route_features[route_id]["geometry"]["coordinates"]
        route_points.append(resample_points(lonlat_to_km(coordinates)))

    route_count = len(route_points)
    matrix = [[0.0 for _ in range(route_count)] for _ in range(route_count)]

    for i in range(route_count):
        for j in range(i + 1, route_count):
            distance = discrete_frechet_distance(route_points[i], route_points[j])
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
            if distance <= DBSCAN_EPS_KM
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


def hierarchical_average(distance_matrix, threshold, label_prefix):
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

        if best_pair is None or best_distance > threshold:
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
        label = f"{label_prefix}_{cluster_index:03d}"
        for route_index in cluster:
            labels[route_index] = label
    return labels


def assign_unique_labels_to_all_routes(routes, unique_routes, unique_labels, cluster_column):
    """Map cluster labels from unique path signatures back to all weights."""
    signature_to_label = dict(zip(unique_routes["path_signature"], unique_labels))
    routes = routes.copy()
    routes[cluster_column] = routes["path_signature"].map(signature_to_label)
    return routes


def relabel_clusters_by_size(routes, cluster_column, label_prefix):
    """Renumber cluster labels by displayed size within one route pair."""
    routes = routes.copy()
    cluster_counts = routes.groupby(cluster_column).size().sort_values(
        ascending=False,
        kind="stable",
    )
    label_map = {
        cluster_id: f"{label_prefix}_{index:03d}"
        for index, cluster_id in enumerate(cluster_counts.index, start=1)
    }
    routes[f"{cluster_column}_source"] = routes[cluster_column]
    routes[cluster_column] = routes[cluster_column].map(label_map)
    return routes


def build_cluster_output(routes, cluster_column):
    """Select common output columns plus the method cluster label."""
    source_column = f"{cluster_column}_source"
    return routes[BASE_OUTPUT_COLUMNS + [cluster_column, source_column]]


def cluster_method_for_all_pairs(
    routes,
    method_name,
    cluster_column,
    label_prefix,
    matrix_builder,
    clusterer,
):
    """Run one clustering method independently for every route pair."""
    outputs = []
    for route_pair in sorted(routes["route_pair"].unique()):
        unique_routes = unique_routes_for_pair(routes, route_pair)
        matrix = matrix_builder(unique_routes)
        labels = clusterer(matrix)
        paired_routes = routes[routes["route_pair"] == route_pair].copy()
        paired_routes = assign_unique_labels_to_all_routes(
            paired_routes,
            unique_routes,
            labels,
            cluster_column,
        )
        paired_routes = relabel_clusters_by_size(
            paired_routes,
            cluster_column,
            label_prefix,
        )
        outputs.append(paired_routes)
        summary = paired_routes.groupby(cluster_column).size().sort_values(ascending=False)
        print(
            f"{method_name} {route_pair}: "
            f"{len(summary)} clusters across {len(paired_routes)} weights"
        )

    return build_cluster_output(pd.concat(outputs, ignore_index=True), cluster_column)


def main():
    """Save Chicago cluster CSVs for all requested route-distance methods."""
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
    routes = load_routes()
    route_features = load_route_features()

    dbscan_output = cluster_method_for_all_pairs(
        routes,
        "DBSCAN",
        "dbscan_cluster",
        "dbscan",
        lambda unique_routes: build_frechet_distance_matrix(unique_routes, route_features),
        dbscan,
    )
    dbscan_output.to_csv(DBSCAN_OUTPUT_CSV, index=False)
    print(f"Saved DBSCAN clusters: {DBSCAN_OUTPUT_CSV}")

    edit_output = cluster_method_for_all_pairs(
        routes,
        "Edit distance",
        "edit_distance_cluster",
        "edit",
        lambda unique_routes: build_edit_distance_matrix(unique_routes["path_sequence"].tolist()),
        lambda matrix: hierarchical_average(matrix, EDIT_DISTANCE_THRESHOLD, "edit"),
    )
    edit_output.to_csv(EDIT_OUTPUT_CSV, index=False)
    print(f"Saved edit-distance clusters: {EDIT_OUTPUT_CSV}")

    frechet_output = cluster_method_for_all_pairs(
        routes,
        "Frechet",
        "frechet_cluster",
        "frechet",
        lambda unique_routes: build_frechet_distance_matrix(unique_routes, route_features),
        lambda matrix: hierarchical_average(matrix, FRECHET_THRESHOLD_KM, "frechet"),
    )
    frechet_output.to_csv(FRECHET_OUTPUT_CSV, index=False)
    print(f"Saved Frechet clusters: {FRECHET_OUTPUT_CSV}")

    jaccard_output = cluster_method_for_all_pairs(
        routes,
        "Hierarchical Jaccard",
        "hierarchical_jaccard_cluster",
        "hier_jaccard",
        lambda unique_routes: build_jaccard_distance_matrix(unique_routes["path_sequence"].tolist()),
        lambda matrix: hierarchical_average(
            matrix,
            HIERARCHICAL_JACCARD_THRESHOLD,
            "hier_jaccard",
        ),
    )
    jaccard_output.to_csv(JACCARD_OUTPUT_CSV, index=False)
    print(f"Saved hierarchical Jaccard clusters: {JACCARD_OUTPUT_CSV}")


if __name__ == "__main__":
    main()
