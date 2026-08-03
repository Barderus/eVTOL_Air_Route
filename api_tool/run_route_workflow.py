"""Run the route workflow from grid creation through saved route outputs."""

import os

from api_tool.build_routing_graph import build_routing_graph
from api_tool.create_routing_grid import create_routing_grid
from api_tool.make_density_map import make_density_map
from api_tool.make_route_map import make_route_map
from api_tool.run_astar_route import run_astar_route
from api_tool.save_route_outputs import save_route_outputs
from api_tool.score_grid_cells import score_grid_cells
from api_tool.select_study_area import select_study_area
from api_tool.set_route_locations import set_route_locations
from api_tool.validate_route_weights import (
    request_route_weights_from_user,
    validate_route_weights,
)


def request_workflow_settings_from_user():
    """Ask for the minimum settings needed to run a route workflow."""
    city_name = input("City or study area name: ").strip()
    while city_name == "":
        city_name = input("City or study area name is required: ").strip()

    print("Study area bounds")
    west = float(input("West longitude: ").strip())
    south = float(input("South latitude: ").strip())
    east = float(input("East longitude: ").strip())
    north = float(input("North latitude: ").strip())
    cell_size_m = int(input("Grid cell size in meters: ").strip())

    print("Pre-route density inputs")
    population_path = clean_optional_path(
        input("Population data path, GeoJSON or geospatial file (optional): ")
    )
    traffic_csv_path = clean_optional_path(input("Air traffic CSV path (optional): "))
    airspace_sites = request_airspace_sites_from_user()

    print("Drone takeoff location")
    origin_label = input("Takeoff label: ").strip()
    origin_lon = float(input("Takeoff longitude: ").strip())
    origin_lat = float(input("Takeoff latitude: ").strip())

    print("Destination location")
    destination_label = input("Destination label: ").strip()
    destination_lon = float(input("Destination longitude: ").strip())
    destination_lat = float(input("Destination latitude: ").strip())

    print("Route weights")
    weights = request_route_weights_from_user(require_sum_to_one=True)

    output_folder = clean_optional_path(
        input("Folder where workflow outputs should be saved: ")
    )
    while not output_folder:
        output_folder = clean_optional_path(input("Output folder is required: "))

    locations = {
        "origin": {
            "label": origin_label,
            "lat": origin_lat,
            "lon": origin_lon,
            "type": "origin",
        },
        "destination": {
            "label": destination_label,
            "lat": destination_lat,
            "lon": destination_lon,
            "type": "destination",
        },
    }
    routes = [
        {
            "route_id": "route_001",
            "label": f"{origin_label} to {destination_label}",
            "origin_key": "origin",
            "destination_key": "destination",
        }
    ]

    return {
        "city_name": city_name,
        "bounds": {
            "west": west,
            "south": south,
            "east": east,
            "north": north,
        },
        "cell_size_m": cell_size_m,
        "locations": locations,
        "routes": routes,
        "weights": weights,
        "population_path": population_path,
        "traffic_csv_path": traffic_csv_path,
        "airspace_sites": airspace_sites,
        "output_folder": output_folder,
    }


def parse_optional_float(value):
    """Return a float for non-empty text, otherwise None."""
    value = value.strip()
    if value == "":
        return None
    return float(value)


def clean_optional_path(value):
    """Normalize an optional path typed into the terminal."""
    value = value.strip().strip('"').strip("'").strip()
    if value == "":
        return None
    return value


def request_airspace_sites_from_user():
    """Ask for radial airspace sites used by the pre-route airspace map."""
    site_count_text = input("How many airports or radial airspace sites? [0]: ").strip()
    site_count = int(site_count_text) if site_count_text else 0
    airspace_sites = []

    for index in range(site_count):
        print(f"Airport or airspace site {index + 1}")
        label = input("  Airport/site name: ").strip()
        while label == "":
            label = input("  Airport/site name is required: ").strip()

        airspace_sites.append(
            {
                "label": label,
                "lon": float(input("  Longitude: ").strip()),
                "lat": float(input("  Latitude: ").strip()),
                "inner_radius": float(input("  Inner radius: ").strip()),
                "inner_unit": input("  Inner radius unit [NM]: ").strip() or "NM",
                "outer_radius": parse_optional_float(input("  Outer radius (optional): ")),
                "outer_unit": input("  Outer radius unit [same as inner]: ").strip() or None,
                "vertical_range": input("  Vertical range (optional): ").strip() or None,
            }
        )

    return airspace_sites


def get_default_output_paths(output_folder):
    """Return standard workflow output paths inside one folder."""
    if not output_folder:
        return {}

    return {
        "output_grid_path": os.path.join(output_folder, "grid.geojson"),
        "output_scored_grid_path": os.path.join(output_folder, "scored_grid.geojson"),
        "output_graph_path": os.path.join(output_folder, "routing_graph.pkl"),
        "output_geojson_path": os.path.join(output_folder, "routes.geojson"),
        "output_csv_path": os.path.join(output_folder, "routes.csv"),
        "output_route_map_path": os.path.join(output_folder, "route_map.html"),
        "output_population_map_path": os.path.join(output_folder, "population_density_map.html"),
        "output_airspace_map_path": os.path.join(output_folder, "airspace_density_map.html"),
        "output_traffic_map_path": os.path.join(output_folder, "traffic_density_map.html"),
    }


def run_route_workflow(
    city_name=None,
    bounds=None,
    cell_size_m=None,
    locations=None,
    routes=None,
    weights=None,
    population_path=None,
    population_density_field="density_p_km2",
    airspace_sites=None,
    traffic_csv_path=None,
    grid=None,
    scored_grid=None,
    output_folder=None,
    output_grid_path=None,
    output_scored_grid_path=None,
    output_graph_path=None,
    output_geojson_path=None,
    output_csv_path=None,
    output_route_map_path=None,
    output_population_map_path=None,
    output_airspace_map_path=None,
    output_traffic_map_path=None,
    make_map=True,
    make_pre_route_maps=True,
):
    """
    Run the main route workflow using API-tool functions.

    A caller can pass a raw grid, a pre-scored grid, or bounds and cell size.
    """
    if city_name is None:
        settings = request_workflow_settings_from_user()
        city_name = settings["city_name"]
        bounds = settings["bounds"]
        cell_size_m = settings["cell_size_m"]
        locations = settings["locations"]
        routes = settings["routes"]
        weights = settings["weights"]
        population_path = settings["population_path"]
        traffic_csv_path = settings["traffic_csv_path"]
        airspace_sites = settings["airspace_sites"]
        output_folder = settings["output_folder"]

    default_output_paths = get_default_output_paths(output_folder)
    output_grid_path = output_grid_path or default_output_paths.get("output_grid_path")
    output_scored_grid_path = output_scored_grid_path or default_output_paths.get(
        "output_scored_grid_path"
    )
    output_graph_path = output_graph_path or default_output_paths.get("output_graph_path")
    output_geojson_path = output_geojson_path or default_output_paths.get(
        "output_geojson_path"
    )
    output_csv_path = output_csv_path or default_output_paths.get("output_csv_path")
    output_route_map_path = output_route_map_path or default_output_paths.get(
        "output_route_map_path"
    )
    output_population_map_path = output_population_map_path or default_output_paths.get(
        "output_population_map_path"
    )
    output_airspace_map_path = output_airspace_map_path or default_output_paths.get(
        "output_airspace_map_path"
    )
    output_traffic_map_path = output_traffic_map_path or default_output_paths.get(
        "output_traffic_map_path"
    )

    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    study_area = select_study_area(
        city_name=city_name,
        bounds=bounds,
        cell_size_m=cell_size_m,
        risk_grid_path=output_scored_grid_path or output_grid_path,
        population_output_path=population_path,
        traffic_output_folder=None,
        route_output_folder=None,
    )

    if weights is None:
        weights = request_route_weights_from_user(require_sum_to_one=True)

    weight_validation = validate_route_weights(
        weights["distance_weight"],
        weights["population_weight"],
        weights["traffic_weight"],
        weights["airspace_weight"],
    )
    if not weight_validation["is_valid"]:
        error_text = "; ".join(weight_validation["errors"])
        raise ValueError(error_text)

    grid_result = None
    if scored_grid is None:
        if grid is None:
            grid_result = create_routing_grid(
                city_name=city_name,
                bounds=bounds,
                cell_size_m=cell_size_m,
                output_path=output_grid_path,
            )
            grid = grid_result["grid"]

        score_result = score_grid_cells(
            city_name=city_name,
            grid=grid,
            population_path=population_path,
            population_density_field=population_density_field,
            airspace_sites=airspace_sites,
            traffic_csv_path=traffic_csv_path,
            output_path=output_scored_grid_path,
        )
        scored_grid = score_result["grid"]
    else:
        score_result = {
            "grid": scored_grid,
            "cell_count": len(scored_grid),
            "output_path": output_scored_grid_path,
        }

    pre_route_map_results = []
    if make_pre_route_maps:
        map_specs = [
            ("Population", "density_p_km2", output_population_map_path),
            ("Airspace", "airport_risk_combined", output_airspace_map_path),
            ("Air Traffic", "traffic_risk", output_traffic_map_path),
        ]
        for map_label, value_field, output_path in map_specs:
            if output_path:
                pre_route_map_results.append(
                    make_density_map(
                        city_name=city_name,
                        grid=scored_grid,
                        grid_path=output_scored_grid_path,
                        value_field=value_field,
                        output_html_path=output_path,
                        title=f"{city_name} {map_label} Density",
                    )
                )

    graph_package = build_routing_graph(
        city_name=city_name,
        grid=scored_grid,
        output_path=output_graph_path,
    )

    if locations is None or routes is None:
        route_location_result = set_route_locations(city_name=city_name)
        locations = route_location_result["locations"]
        routes = route_location_result["routes"]
    else:
        route_location_result = set_route_locations(
            city_name=city_name,
            locations=locations,
            routes=routes,
        )
        routes = route_location_result["routes"]

    route_results = []
    for route in routes:
        route_result = run_astar_route(
            city_name=city_name,
            graph_package=graph_package,
            origin=route["origin"],
            destination=route["destination"],
            weights=weights,
            route_id=route["route_id"],
        )
        route_result["route_label"] = route["label"]
        route_results.append(route_result)

    save_result = save_route_outputs(
        route_results,
        geojson_path=output_geojson_path,
        csv_path=output_csv_path,
    )

    map_result = None
    if make_map and output_geojson_path and output_route_map_path:
        map_result = make_route_map(
            route_geojson_path=output_geojson_path,
            scored_grid_path=output_scored_grid_path,
            locations=locations,
            output_html_path=output_route_map_path,
            title=f"{city_name} eVTOL Route Map",
        )

    return {
        "study_area": study_area,
        "grid_result": grid_result,
        "score_result": score_result,
        "graph_summary": {
            "nodes": graph_package["graph"].number_of_nodes(),
            "edges": graph_package["graph"].number_of_edges(),
            "maximum_distance": graph_package["maximum_distance"],
        },
        "route_location_result": route_location_result,
        "weight_validation": weight_validation,
        "route_results": route_results,
        "save_result": save_result,
        "map_result": map_result,
        "pre_route_map_results": pre_route_map_results,
        "output_folder": output_folder,
        "output_paths": {
            "grid": output_grid_path,
            "scored_grid": output_scored_grid_path,
            "graph": output_graph_path,
            "routes_geojson": output_geojson_path,
            "routes_csv": output_csv_path,
            "route_map": output_route_map_path,
            "population_density_map": output_population_map_path,
            "airspace_density_map": output_airspace_map_path,
            "traffic_density_map": output_traffic_map_path,
        },
    }


if __name__ == "__main__":
    result = run_route_workflow()
    print(
        {
            "city": result["study_area"]["city"],
            "routes": len(result["route_results"]),
            "graph": result["graph_summary"],
            "saved_outputs": result["save_result"]["saved_outputs"],
            "route_map": result["output_paths"]["route_map"],
        }
    )
