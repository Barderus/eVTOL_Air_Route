"""Prepare or run OpenSky traffic exports for an API-style workflow."""

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


def prompt_optional(prompt_text):
    """Ask for optional text and return None when left blank."""
    value = input(prompt_text).strip()
    if value == "":
        return None
    return value


def list_existing_exports(output_folder):
    """Return CSV files that already exist in the traffic output folder."""
    if not output_folder or not os.path.isdir(output_folder):
        return []

    exports = []
    for file_name in sorted(os.listdir(output_folder)):
        if file_name.lower().endswith(".csv"):
            file_path = os.path.join(output_folder, file_name)
            exports.append(
                {
                    "file_name": file_name,
                    "path": file_path,
                    "size_bytes": os.path.getsize(file_path),
                }
            )
    return exports


def build_export_command(script_path, export_date, trino_user=None, trino_path=None):
    """Build a generic legacy OpenSky export command."""
    command = ["uv", "run", "python", script_path]

    if export_date:
        command.extend(["--date", export_date])
    if trino_user:
        command.extend(["--user", trino_user])
    if trino_path:
        command.extend(["--trino-path", trino_path])

    return command


def run_export_command(command, working_directory, trino_user):
    """Run a user-supplied OpenSky export command."""
    run_environment = os.environ.copy()
    if trino_user:
        run_environment["TRINO_USER"] = trino_user

    result = subprocess.run(
        command,
        cwd=working_directory or ".",
        env=run_environment,
        check=True,
        capture_output=True,
        text=True,
    )
    return {
        "command": command,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def request_traffic_settings_from_user():
    """Ask the user for OpenSky export settings."""
    city_name = input("City or study area name: ").strip()
    while city_name == "":
        city_name = input("City or study area name is required: ").strip()

    export_date = prompt_optional("Export date, YYYY-MM-DD (optional): ")
    query_path = prompt_optional("SQL query path (optional): ")
    output_folder = prompt_optional("Traffic output folder (optional): ")
    script_path = prompt_optional("OpenSky export script path (optional): ")
    working_directory = prompt_optional("Working directory, blank for project root: ") or "."
    trino_user = prompt_optional("Trino user, blank to use environment/default: ")
    trino_path = prompt_optional("Trino CLI JAR path (optional): ")
    run_export = prompt_yes_no("Run the OpenSky export now?", default=False)

    return {
        "city_name": city_name,
        "export_date": export_date,
        "query_path": query_path,
        "output_folder": output_folder,
        "script_path": script_path,
        "working_directory": working_directory,
        "trino_user": trino_user,
        "trino_path": trino_path,
        "run_export": run_export,
    }


def export_opensky_traffic(
    city_name=None,
    export_date=None,
    query_path=None,
    output_folder=None,
    script_path=None,
    working_directory=".",
    trino_user=None,
    trino_path=None,
    run_export=False,
):
    """
    Return OpenSky export metadata, and optionally run a user-supplied export.

    Network access and Trino authentication are only used when run_export is True.
    """
    if city_name is None:
        settings = request_traffic_settings_from_user()
        city_name = settings["city_name"]
        export_date = settings["export_date"]
        query_path = settings["query_path"]
        output_folder = settings["output_folder"]
        script_path = settings["script_path"]
        working_directory = settings["working_directory"]
        trino_user = settings["trino_user"]
        trino_path = settings["trino_path"]
        run_export = settings["run_export"]

    command = None
    if script_path:
        command = build_export_command(script_path, export_date, trino_user, trino_path)

    run_result = None
    if run_export:
        if not command:
            raise ValueError("script_path is required when run_export is True.")
        run_result = run_export_command(command, working_directory, trino_user)

    return {
        "study_area": select_study_area(city_name=city_name),
        "query_path": query_path,
        "query_exists": os.path.exists(query_path) if query_path else False,
        "output_folder": output_folder,
        "output_folder_exists": os.path.isdir(output_folder) if output_folder else False,
        "existing_exports": list_existing_exports(output_folder),
        "export_date": export_date,
        "script_path": script_path,
        "working_directory": working_directory,
        "command": command,
        "run_export": run_export,
        "run_result": run_result,
    }


if __name__ == "__main__":
    traffic_result = export_opensky_traffic()
    print(traffic_result)
