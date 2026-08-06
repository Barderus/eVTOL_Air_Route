"""Run A* routing against an API routing graph package."""

import math
import os
import pickle

import geopandas as gpd
import networkx as nx
from shapely.geometry import LineString, Point

from api_tool.build_routing_graph import build_routing_graph
from api_tool.validate_route_weights import (
    request_route_weights_from_user,
    validate_route_weights,
)


def load_graph_package(graph_path):
    """Load a graph package saved by build_routing_graph."""
    if not os.path.exists(graph_path):
        raise FileNotFoundError(f"Routing graph package not found: {graph_path}")

    with open(graph_path, "rb") as input_file:
        return pickle.load(input_file)


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


def add_weighted_cost(graph, weights, weight_field="weighted_cost"):
    """Add one weighted edge cost field to a graph."""
    for start_node, end_node in graph.edges:
        edge = graph.edges[start_node, end_node]
        edge[weight_field] = (
            weights["distance_weight"] * edge.get("distance_norm", 0.0)
            + weights["population_weight"] * edge.get("population_norm", 0.0)
            + weights["traffic_weight"] * edge.get("traffic_norm", 0.0)
            + weights["airspace_weight"] * edge.get("airspace_norm", 0.0)
        )
    return weight_field


def make_heuristic(distance_weight, centroid_x, centroid_y, maximum_distance):
    """Build an admissible A* heuristic using only the distance term."""
    if distance_weight <= 0:
        return lambda _start_node, _end_node: 0.0

    def heuristic(start_node, end_node):
        distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        return distance_weight * (distance / maximum_distance)

    return heuristic


def calculate_path_costs(graph, path, weights):
    """Calculate total and per-factor route costs for one path."""
    totals = {
        "distance_norm_sum": 0.0,
        "population_norm_sum": 0.0,
        "traffic_norm_sum": 0.0,
        "airspace_norm_sum": 0.0,
        "route_distance_m": 0.0,
    }

    for start_node, end_node in zip(path, path[1:]):
        edge = graph.edges[start_node, end_node]
        totals["distance_norm_sum"] += edge.get("distance_norm", 0.0)
        totals["population_norm_sum"] += edge.get("population_norm", 0.0)
        totals["traffic_norm_sum"] += edge.get("traffic_norm", 0.0)
        totals["airspace_norm_sum"] += edge.get("airspace_norm", 0.0)
        totals["route_distance_m"] += edge.get("distance_m", 0.0)

    distance_cost = weights["distance_weight"] * totals["distance_norm_sum"]
    population_cost = weights["population_weight"] * totals["population_norm_sum"]
    traffic_cost = weights["traffic_weight"] * totals["traffic_norm_sum"]
    airspace_cost = weights["airspace_weight"] * totals["airspace_norm_sum"]
    total_cost = distance_cost + population_cost + traffic_cost + airspace_cost

    return {
        "total_cost": total_cost,
        "route_distance_km": totals["route_distance_m"] / 1000.0,
        "distance_cost": distance_cost,
        "population_cost": population_cost,
        "traffic_cost": traffic_cost,
        "airspace_cost": airspace_cost,
        "distance_norm_sum": totals["distance_norm_sum"],
        "population_norm_sum": totals["population_norm_sum"],
        "traffic_norm_sum": totals["traffic_norm_sum"],
        "airspace_norm_sum": totals["airspace_norm_sum"],
    }


def make_route_geometry(projected_grid, path):
    """Create one route line from selected grid-cell centroids."""
    if len(path) < 2:
        raise ValueError(
            "Route path must contain at least two grid cells. "
            "Check that the origin and destination are different locations."
        )

    route_points = [projected_grid.geometry.iloc[index].centroid for index in path]
    route_line = gpd.GeoSeries(
        [LineString(route_points)],
        crs=projected_grid.crs,
    ).to_crs("EPSG:4326")
    return route_line.iloc[0]


def route_geometry_to_path(geometry):
    """Convert a route LineString to JSON-safe lon/lat path coordinates."""
    return [
        {"lon": float(longitude), "lat": float(latitude)}
        for longitude, latitude in geometry.coords
    ]


def get_lat_lon(point):
    """Return latitude and longitude from a location dictionary or tuple."""
    if isinstance(point, dict):
        return point["lat"], point["lon"]
    return point[0], point[1]


def run_astar_route(
    city_name=None,
    graph_package=None,
    graph_path=None,
    grid_path=None,
    origin=None,
    destination=None,
    weights=None,
    route_id="route_001",
    interactive=False,
):
    """
    Run A* for one origin, destination, and weight configuration.

    Pass graph_package directly, load graph_path, or pass grid_path to build a
    graph before routing.
    """
    if city_name is None and interactive:
        city_name = input("City or study area name: ").strip()
        grid_path = input("Scored grid path: ").strip()

    if graph_package is None:
        if graph_path:
            graph_package = load_graph_package(graph_path)
        elif grid_path:
            if city_name is None:
                raise ValueError("city_name is required when building from grid_path.")
            graph_package = build_routing_graph(city_name=city_name, grid_path=grid_path)
        else:
            raise ValueError("Pass graph_package, graph_path, or grid_path.")

    if origin is None or destination is None:
        raise ValueError("origin and destination are required.")

    if weights is None and interactive:
        weights = request_route_weights_from_user(require_sum_to_one=True)
    elif weights is None:
        raise ValueError("weights are required for non-interactive A* routing.")

    weight_validation = validate_route_weights(
        weights["distance_weight"],
        weights["population_weight"],
        weights["traffic_weight"],
        weights["airspace_weight"],
    )
    if not weight_validation["is_valid"]:
        error_text = "; ".join(weight_validation["errors"])
        raise ValueError(error_text)

    graph = graph_package["graph"]
    grid = graph_package["grid"]
    projected_grid = graph_package["projected_grid"]
    centroid_x = graph_package["centroid_x"]
    centroid_y = graph_package["centroid_y"]
    maximum_distance = graph_package["maximum_distance"]

    origin_lat, origin_lon = get_lat_lon(origin)
    destination_lat, destination_lon = get_lat_lon(destination)
    origin_node = find_grid_cell(origin_lat, origin_lon, grid, projected_grid)
    destination_node = find_grid_cell(destination_lat, destination_lon, grid, projected_grid)
    if origin_node == destination_node:
        raise ValueError(
            "Origin and destination resolve to the same routing grid cell. "
            "Use different coordinates or a smaller grid cell size."
        )

    weight_field = add_weighted_cost(graph, weights)
    heuristic = make_heuristic(
        weights["distance_weight"],
        centroid_x,
        centroid_y,
        maximum_distance,
    )
    path = nx.astar_path(
        graph,
        origin_node,
        destination_node,
        heuristic=heuristic,
        weight=weight_field,
    )
    path_costs = calculate_path_costs(graph, path, weights)
    geometry = make_route_geometry(projected_grid, path)

    return {
        "route_id": route_id,
        "origin": origin,
        "destination": destination,
        "origin_node": origin_node,
        "destination_node": destination_node,
        "path_node_ids": path,
        "path_nodes": len(path),
        "weights": weights,
        "costs": path_costs,
        "route_path": route_geometry_to_path(geometry),
        "route_wkt": geometry.wkt,
        "geometry": geometry,
        "status": "success",
    }


if __name__ == "__main__":
    result = run_astar_route(interactive=True)
    printable = result.copy()
    printable["geometry"] = result["geometry"].wkt
    print(printable)
