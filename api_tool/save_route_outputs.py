"""Save route outputs for an API-style workflow."""

import csv
import json
import os

import geopandas as gpd


def ensure_output_folder(output_path):
    """Create the output folder when the path includes one."""
    output_folder = os.path.dirname(output_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)


def route_result_to_properties(route_result):
    """Convert one route result into GeoJSON-safe properties."""
    costs = route_result.get("costs", {})
    weights = route_result.get("weights", {})
    properties = {
        "route_id": route_result.get("route_id"),
        "origin_node": route_result.get("origin_node"),
        "destination_node": route_result.get("destination_node"),
        "path_nodes": route_result.get("path_nodes"),
        "path_node_ids": json.dumps(route_result.get("path_node_ids", [])),
        "status": route_result.get("status"),
    }
    properties.update(weights)
    properties.update(costs)
    return properties


def save_route_geojson(route_results, output_path):
    """Save route results with geometry as GeoJSON."""
    ensure_output_folder(output_path)
    rows = []
    geometries = []

    for route_result in route_results:
        rows.append(route_result_to_properties(route_result))
        geometries.append(route_result["geometry"])

    routes_gdf = gpd.GeoDataFrame(rows, geometry=geometries, crs="EPSG:4326")
    routes_gdf.to_file(output_path, driver="GeoJSON")


def save_route_csv(route_results, output_path):
    """Save route result properties as CSV."""
    ensure_output_folder(output_path)
    rows = [route_result_to_properties(route_result) for route_result in route_results]

    if not rows:
        with open(output_path, "w", encoding="utf-8") as output_file:
            output_file.write("")
        return

    fieldnames = list(rows[0].keys())
    with open(output_path, "w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_json(data, output_path):
    """Save data as JSON."""
    ensure_output_folder(output_path)
    with open(output_path, "w", encoding="utf-8") as output_file:
        json.dump(data, output_file, indent=2)


def save_route_outputs(route_results=None, geojson_path=None, csv_path=None, metadata_path=None):
    """Save API route results as GeoJSON, CSV, or metadata JSON."""
    route_results = route_results or []
    saved_outputs = []

    if geojson_path:
        save_route_geojson(route_results, geojson_path)
        saved_outputs.append({"type": "geojson", "path": geojson_path})

    if csv_path:
        save_route_csv(route_results, csv_path)
        saved_outputs.append({"type": "csv", "path": csv_path})

    if metadata_path:
        metadata = {
            "route_count": len(route_results),
            "route_ids": [route_result.get("route_id") for route_result in route_results],
        }
        save_json(metadata, metadata_path)
        saved_outputs.append({"type": "metadata", "path": metadata_path})

    return {"saved_outputs": saved_outputs}


if __name__ == "__main__":
    print("Pass route_results from Python to save API route outputs.")
