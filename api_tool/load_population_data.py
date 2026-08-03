"""Load or describe population data for an API-style workflow."""

import os
import subprocess

import geopandas as gpd

from api_tool.select_study_area import select_study_area


def prompt_yes_no(prompt_text, default=False):
    """Ask for a yes/no value."""
    default_text = "y" if default else "n"
    value = input(f"{prompt_text} [y/n, default {default_text}]: ").strip().lower()
    if value == "":
        return default
    return value in ["y", "yes"]


def prompt_path_list():
    """Ask for population source paths until the user leaves a blank line."""
    paths = []
    print("Enter population source paths. Leave blank when done.")
    while True:
        path = input("Source path: ").strip()
        if path == "":
            break
        paths.append(path)
    return paths


def get_path_status(paths):
    """Return exists flags for a list of paths."""
    path_status = []
    for path in paths:
        path_status.append({"path": path, "exists": os.path.exists(path)})
    return path_status


def summarize_population_output(output_path):
    """Read a population spatial output and return a small metadata summary."""
    if not output_path:
        return {
            "output_path": None,
            "exists": False,
            "feature_count": 0,
            "columns": [],
            "bounds": None,
        }

    if not os.path.exists(output_path):
        return {
            "output_path": output_path,
            "exists": False,
            "feature_count": 0,
            "columns": [],
            "bounds": None,
        }

    population_gdf = gpd.read_file(output_path)
    bounds = population_gdf.total_bounds.tolist()

    return {
        "output_path": output_path,
        "exists": True,
        "feature_count": int(len(population_gdf)),
        "columns": list(population_gdf.columns),
        "bounds": {
            "west": bounds[0],
            "south": bounds[1],
            "east": bounds[2],
            "north": bounds[3],
        },
    }


def run_population_script(script_path, working_directory):
    """Run a user-supplied population script."""
    command = ["uv", "run", "python", script_path]
    result = subprocess.run(
        command,
        cwd=working_directory or ".",
        check=True,
        capture_output=True,
        text=True,
    )
    return {
        "command": command,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def request_population_settings_from_user():
    """Ask the user for population workflow settings."""
    city_name = input("City or study area name: ").strip()
    while city_name == "":
        city_name = input("City or study area name is required: ").strip()

    source_paths = prompt_path_list()
    output_path = input("Population output path to summarize (optional): ").strip()
    if output_path == "":
        output_path = None

    run_script = prompt_yes_no("Run a population preprocessing script now?", default=False)
    script_path = None
    working_directory = "."
    if run_script:
        script_path = input("Population script path: ").strip()
        working_directory = input("Working directory, blank for project root: ").strip()
        if working_directory == "":
            working_directory = "."

    return city_name, source_paths, output_path, run_script, script_path, working_directory


def load_population_data(
    city_name=None,
    source_paths=None,
    output_path=None,
    run_script=False,
    script_path=None,
    working_directory=".",
):
    """
    Return population data metadata supplied by a user or API caller.

    If run_script is True, script_path is executed before output_path is summarized.
    """
    if city_name is None:
        (
            city_name,
            source_paths,
            output_path,
            run_script,
            script_path,
            working_directory,
        ) = request_population_settings_from_user()

    source_paths = source_paths or []

    run_result = None
    if run_script:
        if not script_path:
            raise ValueError("script_path is required when run_script is True.")
        run_result = run_population_script(script_path, working_directory)

    return {
        "study_area": select_study_area(city_name=city_name),
        "source_status": get_path_status(source_paths),
        "population_output": summarize_population_output(output_path),
        "run_script": run_script,
        "script_path": script_path,
        "working_directory": working_directory,
        "run_result": run_result,
    }


if __name__ == "__main__":
    population_result = load_population_data()
    print(population_result)
