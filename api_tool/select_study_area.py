"""Collect study-area settings for an API-style workflow."""

import os


def normalize_name(value):
    """Convert a user-facing name into a stable key."""
    return value.strip().lower().replace("-", "_").replace(" ", "_")


def prompt_optional(prompt_text):
    """Ask for optional text and return None when left blank."""
    value = input(prompt_text).strip()
    if value == "":
        return None
    return value


def prompt_float(prompt_text, required=False):
    """Ask for a floating-point value."""
    while True:
        value = input(prompt_text).strip()
        if value == "" and not required:
            return None
        try:
            return float(value)
        except ValueError:
            print("Enter a number, or leave it blank if optional.")


def prompt_int(prompt_text, required=False):
    """Ask for an integer value."""
    while True:
        value = input(prompt_text).strip()
        if value == "" and not required:
            return None
        try:
            return int(value)
        except ValueError:
            print("Enter a whole number, or leave it blank if optional.")


def add_path_status(study_area):
    """Add exists flags for important local paths."""
    checked_area = study_area.copy()
    path_status = {}

    for field_name in [
        "folder",
        "risk_grid_path",
        "population_output_path",
        "traffic_output_folder",
        "route_output_folder",
        "map_output_folder",
    ]:
        path_value = checked_area.get(field_name)
        if path_value:
            path_status[field_name] = {
                "path": path_value,
                "exists": os.path.exists(path_value),
            }

    checked_area["path_status"] = path_status
    return checked_area


def request_study_area_from_user():
    """Ask the user for study-area settings."""
    city_name = input("City or study area name: ").strip()
    while city_name == "":
        city_name = input("City or study area name is required: ").strip()

    bounds = {
        "west": prompt_float("West longitude (optional): "),
        "south": prompt_float("South latitude (optional): "),
        "east": prompt_float("East longitude (optional): "),
        "north": prompt_float("North latitude (optional): "),
    }
    bounds = {key: value for key, value in bounds.items() if value is not None}

    return {
        "city": city_name,
        "city_key": normalize_name(city_name),
        "folder": prompt_optional("Project folder for this area (optional): "),
        "timezone": prompt_optional("Timezone, for example America/Chicago (optional): "),
        "cell_size_m": prompt_int("Grid cell size in meters (optional): "),
        "bounds": bounds or None,
        "risk_grid_path": prompt_optional("Risk grid path (optional): "),
        "population_output_path": prompt_optional("Population output path (optional): "),
        "traffic_output_folder": prompt_optional("Traffic output folder (optional): "),
        "route_output_folder": prompt_optional("Route output folder (optional): "),
        "map_output_folder": prompt_optional("Map output folder (optional): "),
        "notes": prompt_optional("Notes (optional): "),
    }


def select_study_area(
    city_name=None,
    folder=None,
    timezone=None,
    cell_size_m=None,
    bounds=None,
    risk_grid_path=None,
    population_output_path=None,
    traffic_output_folder=None,
    route_output_folder=None,
    map_output_folder=None,
    notes=None,
):
    """
    Return study-area metadata supplied by a user or API caller.

    If city_name is not supplied, this function asks for values interactively.
    """
    if city_name is None:
        return add_path_status(request_study_area_from_user())

    study_area = {
        "city": city_name,
        "city_key": normalize_name(city_name),
        "folder": folder,
        "timezone": timezone,
        "cell_size_m": cell_size_m,
        "bounds": bounds,
        "risk_grid_path": risk_grid_path,
        "population_output_path": population_output_path,
        "traffic_output_folder": traffic_output_folder,
        "route_output_folder": route_output_folder,
        "map_output_folder": map_output_folder,
        "notes": notes,
    }
    return add_path_status(study_area)


if __name__ == "__main__":
    selected_area = select_study_area()
    print(selected_area)
