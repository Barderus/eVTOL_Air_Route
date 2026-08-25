"""Cluster direct routes using hierarchical clustering with Frechet distance."""

import json
import math
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "Chicago/route_robustness/output/direct_route_weight_configurations.csv"
ALL_ROUTES_GEOJSON = "Chicago/route_robustness/output/all_routes.geojson"
OUTPUT_FOLDER = "Chicago/route_robustness/output"
OUTPUT_CSV = "Chicago/route_robustness/output/direct_route_clusters_frechet.csv"

FRECHET_THRESHOLD_KM = 2.5
MAX_FRECHET_POINTS = 80


def load_routes():
    """Load direct route records."""
    return pd.read_csv(DIRECT_ROUTE_WEIGHTS_CSV)


def load_route_features():
    """Load all route GeoJSON features by route_run_id."""
    with open(ALL_ROUTES_GEOJSON, "r", encoding="utf-8") as file:
        geojson = json.load(file)

    features = {}
    for feature in geojson["features"]:
        route_id = feature["properties"]["route_run_id"]
        features[route_id] = feature
    return features


def lonlat_to_km(coordinates):
    """Convert lon/lat coordinates to approximate local kilometers."""
    latitudes = [lat for _, lat in coordinates]
    center_lat = math.radians(sum(latitudes) / len(latitudes))
    points = []

    for lon, lat in coordinates:
        x = lon * 111.320 * math.cos(center_lat)
        y = lat * 110.574
        points.append((x, y))
    return points


def resample_points(points):
    """Resample a route geometry to a smaller fixed number of points."""
    if len(points) <= MAX_FRECHET_POINTS:
        return points

    distances = [0.0]
    for index in range(1, len(points)):
        previous_x, previous_y = points[index - 1]
        current_x, current_y = points[index]
        step = math.hypot(current_x - previous_x, current_y - previous_y)
        distances.append(distances[-1] + step)

    total_distance = distances[-1]
    if total_distance == 0:
        return [points[0]]

    resampled = []
    segment_index = 1
    for target_index in range(MAX_FRECHET_POINTS):
        target_distance = total_distance * target_index / (MAX_FRECHET_POINTS - 1)
        while (
            segment_index < len(distances) - 1
            and distances[segment_index] < target_distance
        ):
            segment_index += 1

        previous_distance = distances[segment_index - 1]
        current_distance = distances[segment_index]
        previous_point = points[segment_index - 1]
        current_point = points[segment_index]

        if current_distance == previous_distance:
            resampled.append(current_point)
            continue

        ratio = (target_distance - previous_distance) / (
            current_distance - previous_distance
        )
        x = previous_point[0] + ratio * (current_point[0] - previous_point[0])
        y = previous_point[1] + ratio * (current_point[1] - previous_point[1])
        resampled.append((x, y))

    return resampled


def point_distance(point_a, point_b):
    """Calculate Euclidean distance between two points."""
    return math.hypot(point_a[0] - point_b[0], point_a[1] - point_b[1])


def discrete_frechet_distance(points_a, points_b):
    """Calculate discrete Frechet distance between two point sequences."""
    previous_row = [math.inf] * len(points_b)

    for index_a, point_a in enumerate(points_a):
        current_row = [math.inf] * len(points_b)
        for index_b, point_b in enumerate(points_b):
            distance = point_distance(point_a, point_b)
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


def build_frechet_distance_matrix(routes, route_features):
    """Build a discrete Frechet distance matrix from route geometries."""
    route_ids = routes["route_run_id"].tolist()
    route_points = []

    for route_id in route_ids:
        coordinates = route_features[route_id]["geometry"]["coordinates"]
        points_km = lonlat_to_km(coordinates)
        route_points.append(resample_points(points_km))

    route_count = len(route_points)
    matrix = [[0.0 for _ in range(route_count)] for _ in range(route_count)]

    for i in range(route_count):
        for j in range(i + 1, route_count):
            distance = discrete_frechet_distance(route_points[i], route_points[j])
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

        if best_pair is None or best_distance > FRECHET_THRESHOLD_KM:
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
        label = f"frechet_{cluster_index:03d}"
        for route_index in cluster:
            labels[route_index] = label
    return labels


def main():
    """Save Frechet direct-route clusters."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes = load_routes()
    route_features = load_route_features()
    matrix = build_frechet_distance_matrix(routes, route_features)
    routes["frechet_cluster"] = hierarchical_average(matrix)

    output = routes[
        [
            "route_run_id",
            "weight_id",
            "direct_source_cluster_id",
            "frechet_cluster",
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

    summary = output.groupby("frechet_cluster").size().sort_values(ascending=False)
    print(f"Saved Frechet clusters: {OUTPUT_CSV}")
    print("Cluster count:", len(summary))
    print(summary.head(10).to_string())


if __name__ == "__main__":
    main()
