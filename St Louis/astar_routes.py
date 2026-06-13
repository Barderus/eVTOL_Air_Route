"""Generate St. Louis routes with A* search."""

import math
import os

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point

import settings


def normalize(
    values,
    clip_percentile=None,
    transform=None,
    power=1.0,
):
    """Scale non-negative values to the 0 to 1 range."""
    values = np.asarray(values, dtype=float)
    values = np.clip(values, 0.0, None)

    if transform == "log1p":
        values = np.log1p(values)
    elif transform is not None:
        raise ValueError(f"Unsupported transform: {transform}")

    if clip_percentile is None:
        scale = float(np.max(values))
    else:
        positive_values = values[values > 0]
        if positive_values.size == 0:
            return np.zeros_like(values)
        scale = float(np.percentile(positive_values, clip_percentile))

    if scale <= 0:
        return np.zeros_like(values)

    normalized_values = np.clip(values / scale, 0.0, 1.0)
    if power != 1.0:
        normalized_values = np.power(normalized_values, power)
    return normalized_values


def build_graph(grid):
    """Build a graph from touching grid cells."""
    grid = grid.reset_index(drop=True)
    projected_grid = grid.to_crs("EPSG:3857")
    centroids = projected_grid.geometry.centroid
    centroid_x = centroids.x.to_numpy()
    centroid_y = centroids.y.to_numpy()

    graph = nx.Graph()
    graph.add_nodes_from(range(len(grid)))

    spatial_index = grid.sindex
    for cell_index, geometry in enumerate(grid.geometry):
        candidate_indexes = spatial_index.query(geometry, predicate="touches")
        for candidate_index in candidate_indexes:
            candidate_index = int(candidate_index)
            if candidate_index > cell_index:
                graph.add_edge(cell_index, candidate_index)

    return graph, projected_grid, centroid_x, centroid_y


def find_grid_cell(latitude, longitude, grid, projected_grid):
    """Find the grid cell containing or nearest to a point."""
    point = Point(longitude, latitude)
    containing_cells = list(grid.sindex.query(point, predicate="within"))
    if containing_cells:
        return int(containing_cells[0])

    projected_point = gpd.GeoSeries([point], crs=grid.crs).to_crs(
        projected_grid.crs
    ).iloc[0]
    distances = projected_grid.geometry.centroid.distance(projected_point)
    return int(distances.idxmin())


def load_traffic_counts(csv_path, grid):
    """Count traffic observations inside each grid cell."""
    traffic_data = pd.read_csv(csv_path)
    required_columns = {"lat", "lon"}
    if not required_columns.issubset(traffic_data.columns):
        raise ValueError(f"Traffic file must contain lat and lon columns: {csv_path}")

    traffic_points = gpd.GeoDataFrame(
        traffic_data,
        geometry=gpd.points_from_xy(traffic_data["lon"], traffic_data["lat"]),
        crs="EPSG:4326",
    ).to_crs(grid.crs)

    joined_points = gpd.sjoin(
        traffic_points,
        grid[["geometry"]],
        how="inner",
        predicate="within",
    )
    counts = joined_points.groupby("index_right").size()
    return counts.reindex(range(len(grid)), fill_value=0).to_numpy(dtype=float)


def add_edge_weights(
    graph,
    centroid_x,
    centroid_y,
    population_risk,
    airspace_risk,
    traffic_risk,
):
    """Add the configured combined cost to every graph edge."""
    edge_distances = []
    for start_node, end_node in graph.edges:
        distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        edge_distances.append(distance)

    maximum_distance = max(edge_distances)

    for start_node, end_node in graph.edges:
        distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        distance_cost = distance / maximum_distance
        population_cost = (
            population_risk[start_node] + population_risk[end_node]
        ) / 2.0
        airspace_cost = (
            airspace_risk[start_node] + airspace_risk[end_node]
        ) / 2.0
        traffic_cost = (
            traffic_risk[start_node] + traffic_risk[end_node]
        ) / 2.0

        graph.edges[start_node, end_node]["weight"] = (
            settings.DISTANCE_WEIGHT * distance_cost
            + settings.POPULATION_WEIGHT * population_cost
            + settings.AIRSPACE_WEIGHT * airspace_cost
            + settings.TRAFFIC_WEIGHT * traffic_cost
        )

    return maximum_distance


def calculate_path_costs(
    path,
    centroid_x,
    centroid_y,
    population_risk,
    airspace_risk,
    traffic_risk,
    maximum_distance,
):
    """Calculate weighted cost contributions for one route path."""
    distance_cost = 0.0
    population_cost = 0.0
    airspace_cost = 0.0
    traffic_cost = 0.0
    route_distance_m = 0.0

    for start_node, end_node in zip(path, path[1:]):
        edge_distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        route_distance_m += edge_distance

        distance_cost += (
            settings.DISTANCE_WEIGHT * edge_distance / maximum_distance
        )
        population_cost += settings.POPULATION_WEIGHT * (
            population_risk[start_node] + population_risk[end_node]
        ) / 2.0
        airspace_cost += settings.AIRSPACE_WEIGHT * (
            airspace_risk[start_node] + airspace_risk[end_node]
        ) / 2.0
        traffic_cost += settings.TRAFFIC_WEIGHT * (
            traffic_risk[start_node] + traffic_risk[end_node]
        ) / 2.0

    total_cost = (
        distance_cost
        + population_cost
        + airspace_cost
        + traffic_cost
    )

    if total_cost > 0:
        distance_percent = 100.0 * distance_cost / total_cost
        population_percent = 100.0 * population_cost / total_cost
        airspace_percent = 100.0 * airspace_cost / total_cost
        traffic_percent = 100.0 * traffic_cost / total_cost
    else:
        distance_percent = 0.0
        population_percent = 0.0
        airspace_percent = 0.0
        traffic_percent = 0.0

    return {
        "total_cost": total_cost,
        "route_distance_km": route_distance_m / 1000.0,
        "distance_cost": distance_cost,
        "population_cost": population_cost,
        "airspace_cost": airspace_cost,
        "traffic_cost": traffic_cost,
        "distance_percent": distance_percent,
        "population_percent": population_percent,
        "airspace_percent": airspace_percent,
        "traffic_percent": traffic_percent,
    }


def make_route_geometry(projected_grid, path):
    """Create one route line from the selected grid cells."""
    route_points = [projected_grid.geometry.iloc[index].centroid for index in path]
    route_line = gpd.GeoSeries(
        [LineString(route_points)],
        crs=projected_grid.crs,
    ).to_crs("EPSG:4326")
    return route_line.iloc[0]


def main():
    settings.validate_route_settings()

    if not os.path.exists(settings.RISK_GRID_GEOJSON):
        raise FileNotFoundError(
            f"St. Louis risk grid not found: {settings.RISK_GRID_GEOJSON}"
        )

    grid = gpd.read_file(settings.RISK_GRID_GEOJSON).reset_index(drop=True)
    required_columns = {"city_risk", "airport_risk_combined"}
    missing_columns = required_columns.difference(grid.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"St. Louis risk grid is missing columns: {missing_text}")

    graph, projected_grid, centroid_x, centroid_y = build_graph(grid)
    population_risk = normalize(grid["city_risk"])
    airspace_risk = normalize(grid["airport_risk_combined"])
    route_rows = []
    analysis_rows = []
    route_geometries = []

    for dataset in settings.TRAFFIC_DATASETS:
        csv_path = os.path.join(settings.TRAFFIC_FOLDER, dataset["filename"])
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Traffic file not found: {csv_path}")

        traffic_counts = load_traffic_counts(csv_path, grid)
        traffic_risk = normalize(
            traffic_counts,
            clip_percentile=95.0,
            transform="log1p",
            power=0.75,
        )
        maximum_distance = add_edge_weights(
            graph,
            centroid_x,
            centroid_y,
            population_risk,
            airspace_risk,
            traffic_risk,
        )

        for route in settings.ROUTES:
            start_latitude, start_longitude = route["start"]
            destination_latitude, destination_longitude = route["destination"]
            start_node = find_grid_cell(
                start_latitude,
                start_longitude,
                grid,
                projected_grid,
            )
            end_node = find_grid_cell(
                destination_latitude,
                destination_longitude,
                grid,
                projected_grid,
            )
            path = nx.astar_path(graph, start_node, end_node, weight="weight")
            path_costs = calculate_path_costs(
                path,
                centroid_x,
                centroid_y,
                population_risk,
                airspace_risk,
                traffic_risk,
                maximum_distance,
            )

            route_rows.append(
                {
                    "route": route["label"],
                    "dataset": dataset["label"],
                    "date": dataset["date"],
                    "origin": route["start_label"],
                    "destination": route["destination_label"],
                    "path_nodes": len(path),
                }
            )
            analysis_row = {
                "route": route["label"],
                "dataset": dataset["label"],
                "date": dataset["date"],
                "origin": route["start_label"],
                "destination": route["destination_label"],
                "path_nodes": len(path),
            }
            analysis_row.update(path_costs)
            analysis_rows.append(analysis_row)
            route_geometries.append(make_route_geometry(projected_grid, path))

    routes = gpd.GeoDataFrame(
        route_rows,
        geometry=route_geometries,
        crs="EPSG:4326",
    )

    os.makedirs(settings.ROUTE_GEOJSON_FOLDER, exist_ok=True)
    output_path = os.path.join(
        settings.ROUTE_GEOJSON_FOLDER,
        "st_louis_astar_routes.geojson",
    )
    routes.to_file(output_path, driver="GeoJSON")
    print(f"Saved St. Louis routes: {output_path}")

    analysis_columns = [
        "route",
        "dataset",
        "date",
        "origin",
        "destination",
        "path_nodes",
        "route_distance_km",
        "total_cost",
        "distance_cost",
        "distance_percent",
        "population_cost",
        "population_percent",
        "airspace_cost",
        "airspace_percent",
        "traffic_cost",
        "traffic_percent",
    ]
    analysis = pd.DataFrame(analysis_rows)[analysis_columns]
    analysis.to_csv(settings.ROUTE_ANALYSIS_CSV, index=False)
    print(f"Saved route cost analysis: {settings.ROUTE_ANALYSIS_CSV}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error
