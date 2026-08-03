"""Terminal menu for running API-tool workflows."""

import os
import sys


if __package__ in [None, ""]:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api_tool.generate_weight_configurations import generate_weight_configurations
from api_tool.make_route_map import make_route_map
from api_tool.run_route_workflow import run_route_workflow
from api_tool.validate_route_weights import (
    request_route_weights_from_user,
    validate_route_weights,
)


def prompt_choice(prompt_text, valid_choices):
    """Ask for a menu choice until the user enters a valid option."""
    while True:
        choice = input(prompt_text).strip()
        if choice in valid_choices:
            return choice
        print("Choose one of:", ", ".join(valid_choices))


def prompt_float(prompt_text):
    """Ask for a required floating-point value."""
    while True:
        value = input(prompt_text).strip()
        try:
            return float(value)
        except ValueError:
            print("Enter a valid number.")


def prompt_optional_path(prompt_text):
    """Ask for an optional path."""
    value = input(prompt_text).strip()
    if value == "":
        return None
    return value


def print_header(title):
    """Print a simple terminal section header."""
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


def print_main_menu():
    """Print the main terminal UI menu."""
    print_header("eVTOL API Tool Terminal UI")
    print("1. Build density maps and routes")
    print("2. Generate weight configurations")
    print("3. Validate route weights")
    print("4. Show workflow output filenames")
    print("5. Make route map")
    print("0. Exit")


def run_full_route_workflow_menu():
    """Run the interactive full route workflow."""
    print_header("Build Density Maps And Routes")
    print("This builds population, airspace, and traffic maps before routing.")
    try:
        result = run_route_workflow()
    except (FileNotFoundError, OSError, ValueError) as error:
        print()
        print("Workflow failed:")
        print(error)
        return

    print_header("Workflow Complete")
    print("City:", result["study_area"]["city"])
    print("Routes:", len(result["route_results"]))
    print("Graph nodes:", result["graph_summary"]["nodes"])
    print("Graph edges:", result["graph_summary"]["edges"])
    print("Pre-route maps:")
    for map_result in result["pre_route_map_results"]:
        print("-", map_result["value_field"], map_result["output_html_path"])
    print("Saved outputs:")
    for saved_output in result["save_result"]["saved_outputs"]:
        print("-", saved_output["type"], saved_output["path"])


def generate_weight_configurations_menu():
    """Prompt for a weight step and optionally save weight configurations."""
    print_header("Generate Weight Configurations")
    step = prompt_float("Weight step, for example 0.10: ")
    output_path = prompt_optional_path("CSV output path, blank to only print summary: ")

    result = generate_weight_configurations(step=step, output_path=output_path)

    print("Generated configurations:", result["count"])
    print("Resolved step:", result["step"])
    if result["output_path"]:
        print("Saved:", result["output_path"])
    else:
        print("First five configurations:")
        for row in result["weight_configurations"][:5]:
            print(row)


def validate_route_weights_menu():
    """Prompt for route weights and print validation result."""
    print_header("Validate Route Weights")
    weights = request_route_weights_from_user(require_sum_to_one=True)
    result = validate_route_weights(**weights)

    print("Weight sum:", result["weight_sum"])
    print("Valid:", result["is_valid"])
    if result["errors"]:
        print("Errors:")
        for error in result["errors"]:
            print("-", error)


def show_workflow_output_filenames_menu():
    """Show the standard filenames created by the workflow output folder."""
    print_header("Workflow Output Filenames")
    output_folder = input("Output folder: ").strip()
    while output_folder == "":
        output_folder = input("Output folder is required: ").strip()

    print("The workflow will save pre-route maps first:")
    print(f"- {output_folder}/grid.geojson")
    print(f"- {output_folder}/scored_grid.geojson")
    print(f"- {output_folder}/population_density_map.html")
    print(f"- {output_folder}/airspace_density_map.html")
    print(f"- {output_folder}/traffic_density_map.html")
    print("Then it will save route outputs:")
    print(f"- {output_folder}/routing_graph.pkl")
    print(f"- {output_folder}/routes.geojson")
    print(f"- {output_folder}/routes.csv")
    print(f"- {output_folder}/route_map.html")


def make_route_map_menu():
    """Prompt for route map inputs and create a Leaflet HTML map."""
    print_header("Make Route Map")
    route_geojson_path = input("Route GeoJSON path: ").strip()
    scored_grid_path = input("Scored grid GeoJSON path (optional): ").strip() or None
    output_html_path = input("Route map HTML output path: ").strip()
    title = input("Map title, blank for default: ").strip() or "eVTOL Route Map"

    result = make_route_map(
        route_geojson_path=route_geojson_path,
        scored_grid_path=scored_grid_path,
        output_html_path=output_html_path,
        title=title,
    )
    print("Saved:", result["output_html_path"])
    print("Routes:", result["route_count"])
    print("Grid cells:", result["grid_cell_count"])


def run_terminal_ui():
    """Run the API-tool terminal menu."""
    while True:
        print_main_menu()
        choice = prompt_choice("Choose an option: ", ["0", "1", "2", "3", "4", "5"])

        if choice == "0":
            print("Exiting.")
            return
        if choice == "1":
            run_full_route_workflow_menu()
        elif choice == "2":
            generate_weight_configurations_menu()
        elif choice == "3":
            validate_route_weights_menu()
        elif choice == "4":
            show_workflow_output_filenames_menu()
        elif choice == "5":
            make_route_map_menu()

        input("\nPress Enter to return to the main menu.")


if __name__ == "__main__":
    run_terminal_ui()
