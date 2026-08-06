"""Agent-safe tool wrapper for generating one A* route."""

from numbers import Integral, Real

from shapely.geometry import mapping

from api_tool.run_astar_route import run_astar_route
from api_tool.save_route_outputs import save_route_outputs


ASTAR_ROUTE_TOOL = {
    "name": "generate_astar_route",
    "description": (
        "Generate one eVTOL A* route from a saved routing graph package, "
        "requested origin/destination coordinates, and route cost weights. "
        "The output separates requested coordinates from snapped routing-grid "
        "centroids so callers do not confuse user input with grid-cell geometry."
    ),
    "inputs": {
        "graph_path": {
            "type": "string",
            "description": "Path to a saved routing graph .pkl package.",
            "required": True,
        },
        "origin": {
            "type": "object",
            "description": "Origin coordinate as {'lat': number, 'lon': number}.",
            "required": True,
        },
        "destination": {
            "type": "object",
            "description": "Destination coordinate as {'lat': number, 'lon': number}.",
            "required": True,
        },
        "weights": {
            "type": "object",
            "description": (
                "Route cost weights with distance_weight, population_weight, "
                "traffic_weight, and airspace_weight. Values must sum to 1.0."
            ),
            "required": True,
        },
        "route_id": {
            "type": "string",
            "description": "Caller-supplied route identifier.",
            "required": True,
        },
        "output_geojson_path": {
            "type": "string",
            "description": "Optional path where the generated route GeoJSON should be saved.",
            "required": False,
        },
    },
    "output": {
        "type": "object",
        "description": (
            "Dictionary with status, requested coordinates, snapped grid-cell "
            "coordinates, route start/end geometry, distance, cost components, "
            "optional output paths, and failure details when routing fails."
        ),
    },
}


def to_json_safe(value):
    """Convert NumPy-like scalar values inside route results to plain Python."""
    if isinstance(value, dict):
        return {key: to_json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [to_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [to_json_safe(item) for item in value]
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        return float(value)
    return value


def validate_tool_inputs(graph_path, origin, destination, weights, route_id):
    """Validate the wrapper's strict tool-facing argument contract."""
    if not graph_path:
        raise ValueError("graph_path is required.")
    if not route_id:
        raise ValueError("route_id is required.")

    for label, point in (("origin", origin), ("destination", destination)):
        if not isinstance(point, dict):
            raise ValueError(f"{label} must be a dictionary with lat and lon.")
        if "lat" not in point or "lon" not in point:
            raise ValueError(f"{label} must include lat and lon.")

    required_weights = {
        "distance_weight",
        "population_weight",
        "traffic_weight",
        "airspace_weight",
    }
    if not isinstance(weights, dict):
        raise ValueError("weights must be a dictionary.")
    missing_weights = sorted(required_weights - set(weights))
    if missing_weights:
        raise ValueError(f"weights missing required keys: {missing_weights}.")


def serialize_route_result(route_result, output_geojson_path=None):
    """Convert run_astar_route output into a simple JSON-friendly dictionary."""
    geometry = route_result["geometry"]
    route_coordinates = route_result["route_path"]
    return to_json_safe(
        {
            "route_id": route_result["route_id"],
            "status": route_result["status"],
            "requested": {
                "origin": route_result["origin"],
                "destination": route_result["destination"],
            },
            "snapped": {
                "origin": {
                    "node_id": route_result["origin_node"],
                    "coordinate": route_coordinates[0],
                },
                "destination": {
                    "node_id": route_result["destination_node"],
                    "coordinate": route_coordinates[-1],
                },
            },
            "path_node_ids": route_result["path_node_ids"],
            "path_nodes": route_result["path_nodes"],
            "route": {
                "start_coordinate": route_coordinates[0],
                "end_coordinate": route_coordinates[-1],
                "path_coordinates": route_coordinates,
                "geometry_geojson": mapping(geometry),
                "geometry_wkt": route_result["route_wkt"],
                "distance_km": route_result["costs"]["route_distance_km"],
            },
            "weights": route_result["weights"],
            "costs": {
                "total_cost": route_result["costs"]["total_cost"],
                "distance_cost": route_result["costs"]["distance_cost"],
                "population_cost": route_result["costs"]["population_cost"],
                "traffic_cost": route_result["costs"]["traffic_cost"],
                "airspace_cost": route_result["costs"]["airspace_cost"],
                "normalized_sums": {
                    "distance": route_result["costs"]["distance_norm_sum"],
                    "population": route_result["costs"]["population_norm_sum"],
                    "traffic": route_result["costs"]["traffic_norm_sum"],
                    "airspace": route_result["costs"]["airspace_norm_sum"],
                },
            },
            "outputs": {
                "route_geojson_path": output_geojson_path,
                "route_map_path": None,
            },
        }
    )


def generate_astar_route(
    graph_path,
    origin,
    destination,
    weights,
    route_id,
    output_geojson_path=None,
):
    """
    Tool wrapper for one A* route generation call.

    This wrapper accepts only explicit machine-provided inputs and never asks
    for terminal input.
    """
    try:
        validate_tool_inputs(graph_path, origin, destination, weights, route_id)
        route_result = run_astar_route(
            graph_path=graph_path,
            origin=origin,
            destination=destination,
            weights=weights,
            route_id=route_id,
            interactive=False,
        )
        save_route_outputs([route_result], geojson_path=output_geojson_path)
        return {
            "status": "success",
            "tool": ASTAR_ROUTE_TOOL["name"],
            "route": serialize_route_result(
                route_result,
                output_geojson_path=output_geojson_path,
            ),
            "error": None,
        }
    except Exception as error:
        return {
            "status": "failure",
            "tool": ASTAR_ROUTE_TOOL["name"],
            "route": None,
            "error": {
                "type": type(error).__name__,
                "message": str(error),
            },
        }
