"""Collect route-smoothing settings for an API-style workflow."""

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


def smooth_route(
    city_name=None,
    route_path=None,
    method="turn_penalty",
    output_path=None,
    settings=None,
    script_path=None,
    working_directory=".",
    run_script=False,
):
    """Return route-smoothing metadata supplied by a user or API caller."""
    if city_name is None:
        city_name = input("City or study area name: ").strip()
        route_path = input("Route path: ").strip()
        method = input("Smoothing method, blank for turn_penalty: ").strip() or "turn_penalty"
        output_path = input("Smoothed route output path (optional): ").strip() or None
        script_path = input("Smoothing script path (optional): ").strip() or None
        working_directory = input("Working directory, blank for project root: ").strip() or "."
        run_script = prompt_yes_no("Run the smoothing script now?", default=False)

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
        "route_path": route_path,
        "route_exists": os.path.exists(route_path) if route_path else False,
        "method": method,
        "settings": settings or {},
        "output_path": output_path,
        "output_exists": os.path.exists(output_path) if output_path else False,
        "script_path": script_path,
        "working_directory": working_directory,
        "run_script": run_script,
        "run_result": run_result,
    }


if __name__ == "__main__":
    result = smooth_route()
    print(result)

