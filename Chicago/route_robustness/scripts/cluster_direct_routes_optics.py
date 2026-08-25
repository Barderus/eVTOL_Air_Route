"""Cluster direct routes using an OPTICS-style Frechet workflow."""

import heapq
import json
import math
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "Chicago/route_robustness/output/direct_route_weight_configurations.csv"
ALL_ROUTES_GEOJSON = "Chicago/route_robustness/output/all_routes.geojson"
OUTPUT_FOLDER = "Chicago/route_robustness/output"
OUTPUT_CSV = "Chicago/route_robustness/output/direct_route_clusters_optics.csv"

OPTICS_EPS_KM = 2.5
OPTICS_MIN_SAMPLES = 4
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


def optics(distance_matrix):
    """Run OPTICS-style ordering and extract clusters at OPTICS_EPS_KM."""
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
        if reachability[route_index] > OPTICS_EPS_KM:
            if core_distances[route_index] <= OPTICS_EPS_KM:
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
    route_features = load_route_features()
    matrix = build_frechet_distance_matrix(routes, route_features)
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
