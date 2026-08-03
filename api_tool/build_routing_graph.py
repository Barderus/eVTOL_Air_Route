"""Build a routing graph from a scored geospatial grid."""

import math
import os
import pickle

import geopandas as gpd
import networkx as nx
import numpy as np

from api_tool.select_study_area import select_study_area


DEFAULT_COST_FIELDS = {
    "population": "city_risk",
    "airspace": "airport_risk_combined",
    "traffic": "traffic_risk",
}


def normalize_values(values, clip_percentile=None, transform=None, power=1.0):
    """Scale non-negative values to the 0 to 1 range."""
    values = np.asarray(values, dtype=float)
    values = np.clip(values, 0.0, None)

    if transform == "log1p":
        values = np.log1p(values)
    elif transform is not None:
        raise ValueError(f"Unsupported transform: {transform}")

    if clip_percentile is None:
        scale = float(np.max(values)) if values.size else 0.0
    else:
        positive_values = values[values > 0]
        if positive_values.size == 0:
            return np.zeros_like(values), 1.0
        scale = float(np.percentile(positive_values, clip_percentile))

    if scale <= 0:
        return np.zeros_like(values), 1.0

    normalized_values = np.clip(values / scale, 0.0, 1.0)
    if power != 1.0:
        normalized_values = np.power(normalized_values, power)

    return normalized_values, scale


def load_grid(grid_path):
    """Load a GeoJSON or other geospatial grid file."""
    if not os.path.exists(grid_path):
        raise FileNotFoundError(f"Routing grid not found: {grid_path}")
    return gpd.read_file(grid_path).reset_index(drop=True)


def get_field_values(grid, field_name):
    """Return one grid field as a numeric array, or zeros if not available."""
    if field_name is None or field_name not in grid.columns:
        return np.zeros(len(grid), dtype=float)
    return grid[field_name].fillna(0).to_numpy(dtype=float)


def normalize_cost_layers(grid, cost_fields=None, normalization_settings=None):
    """Normalize population, airspace, and traffic cost fields."""
    cost_fields = cost_fields or DEFAULT_COST_FIELDS
    normalization_settings = normalization_settings or {}
    cost_layers = {}
    normalization_summary = {}

    for cost_name, field_name in cost_fields.items():
        settings = normalization_settings.get(cost_name, {})
        values = get_field_values(grid, field_name)
        normalized_values, scale = normalize_values(
            values,
            clip_percentile=settings.get("clip_percentile"),
            transform=settings.get("transform"),
            power=settings.get("power", 1.0),
        )
        cost_layers[cost_name] = normalized_values
        normalization_summary[cost_name] = {
            "field": field_name,
            "scale": scale,
            "clip_percentile": settings.get("clip_percentile"),
            "transform": settings.get("transform"),
            "power": settings.get("power", 1.0),
        }

    return cost_layers, normalization_summary


def build_touching_cell_graph(grid, allow_diagonal=True):
    """Build graph nodes from grid rows and edges from touching cells."""
    graph = nx.Graph()
    projected_grid = grid.to_crs("EPSG:3857")
    centroids = projected_grid.geometry.centroid
    centroid_x = centroids.x.to_numpy()
    centroid_y = centroids.y.to_numpy()

    for index, row in grid.iterrows():
        graph.add_node(index, geometry=row.geometry)

    spatial_index = grid.sindex
    for cell_index, geometry in enumerate(grid.geometry):
        candidate_indexes = spatial_index.query(geometry, predicate="touches")
        for candidate_index in candidate_indexes:
            candidate_index = int(candidate_index)
            if candidate_index <= cell_index:
                continue

            if not allow_diagonal:
                intersection = geometry.intersection(grid.geometry.iloc[candidate_index])
                if intersection.length <= 0:
                    continue

            graph.add_edge(cell_index, candidate_index)

    return graph, projected_grid, centroid_x, centroid_y


def add_edge_costs(graph, centroid_x, centroid_y, cost_layers):
    """Add normalized edge cost terms to every graph edge."""
    edge_distances = []
    for start_node, end_node in graph.edges:
        distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        edge_distances.append(distance)

    maximum_distance = max(edge_distances) if edge_distances else 1.0

    for start_node, end_node in graph.edges:
        edge_distance = math.hypot(
            centroid_x[start_node] - centroid_x[end_node],
            centroid_y[start_node] - centroid_y[end_node],
        )
        graph.edges[start_node, end_node]["distance_norm"] = edge_distance / maximum_distance
        graph.edges[start_node, end_node]["distance_m"] = edge_distance

        for cost_name, values in cost_layers.items():
            graph.edges[start_node, end_node][f"{cost_name}_norm"] = (
                values[start_node] + values[end_node]
            ) / 2.0

    return maximum_distance


def save_graph_package(graph_package, output_path):
    """Save the graph package as a pickle file."""
    output_folder = os.path.dirname(output_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    with open(output_path, "wb") as output_file:
        pickle.dump(graph_package, output_file)


def build_routing_graph(
    city_name=None,
    grid_path=None,
    grid=None,
    cost_fields=None,
    normalization_settings=None,
    allow_diagonal=True,
    output_path=None,
):
    """
    Build a reusable routing graph from a scored grid.

    The returned graph package can be passed directly to run_astar_route.
    """
    if city_name is None:
        city_name = input("City or study area name: ").strip()
        grid_path = input("Scored grid path: ").strip()
        output_path = input("Graph pickle output path (optional): ").strip() or None

    if grid is None:
        grid = load_grid(grid_path)
    else:
        grid = grid.reset_index(drop=True)

    cost_layers, normalization_summary = normalize_cost_layers(
        grid,
        cost_fields=cost_fields,
        normalization_settings=normalization_settings,
    )
    graph, projected_grid, centroid_x, centroid_y = build_touching_cell_graph(
        grid,
        allow_diagonal=allow_diagonal,
    )
    maximum_distance = add_edge_costs(graph, centroid_x, centroid_y, cost_layers)

    graph_package = {
        "study_area": select_study_area(city_name=city_name),
        "graph": graph,
        "grid": grid,
        "projected_grid": projected_grid,
        "centroid_x": centroid_x,
        "centroid_y": centroid_y,
        "maximum_distance": maximum_distance,
        "cost_fields": cost_fields or DEFAULT_COST_FIELDS,
        "normalization_summary": normalization_summary,
        "allow_diagonal": allow_diagonal,
    }

    if output_path:
        save_graph_package(graph_package, output_path)

    return graph_package


if __name__ == "__main__":
    result = build_routing_graph()
    print(
        {
            "nodes": result["graph"].number_of_nodes(),
            "edges": result["graph"].number_of_edges(),
            "maximum_distance": result["maximum_distance"],
        }
    )
