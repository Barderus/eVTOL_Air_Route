"""Collect route-robustness map settings for an API-style workflow."""

import os
import subprocess

from api_tool.select_study_area import select_study_area


def prompt_yes_no(prompt_text, default=False):
    """Ask for a yes/no value."""
    default_text = "y" if default else "n"
    value = input(f"{prompt_text} [y/n, default {default_text}]: ").strip().lower()
    if value == "":
        return default
    return value in ["y", "yes"]


def make_route_robustness_maps(
    city_name=None,
    routes_path=None,
    clusters_path=None,
    output_html_path=None,
    script_path=None,
    working_directory=".",
    run_script=False,
):
    """Return route-robustness map metadata supplied by a user or API caller."""
    if city_name is None:
        city_name = input("City or study area name: ").strip()
        routes_path = input("Routes GeoJSON path (optional): ").strip() or None
        clusters_path = input("Clusters CSV path (optional): ").strip() or None
        output_html_path = input("Map HTML output path (optional): ").strip() or None
        script_path = input("Map script path (optional): ").strip() or None
        working_directory = input("Working directory, blank for project root: ").strip() or "."
        run_script = prompt_yes_no("Run the map script now?", default=False)

    run_result = None
    if run_script:
        if not script_path:
            raise ValueError("script_path is required when run_script is True.")
        command = ["uv", "run", "python", script_path]
        result = subprocess.run(
            command,
            cwd=working_directory,
            check=True,
            capture_output=True,
            text=True,
        )
        run_result = {"command": command, "stdout": result.stdout, "stderr": result.stderr}

    return {
        "study_area": select_study_area(city_name=city_name),
        "routes_path": routes_path,
        "routes_exists": os.path.exists(routes_path) if routes_path else False,
        "clusters_path": clusters_path,
        "clusters_exists": os.path.exists(clusters_path) if clusters_path else False,
        "output_html_path": output_html_path,
        "output_exists": os.path.exists(output_html_path) if output_html_path else False,
        "script_path": script_path,
        "working_directory": working_directory,
        "run_script": run_script,
        "run_result": run_result,
    }


if __name__ == "__main__":
    result = make_route_robustness_maps()
    print(result)

