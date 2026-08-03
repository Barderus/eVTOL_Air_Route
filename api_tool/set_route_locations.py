"""Collect route locations for an API-style workflow."""

from api_tool.select_study_area import normalize_name, select_study_area


def prompt_float(prompt_text):
    """Ask for a required floating-point value."""
    while True:
        value = input(prompt_text).strip()
        try:
            return float(value)
        except ValueError:
            print("Enter a valid number.")


def prompt_int(prompt_text):
    """Ask for a required integer value."""
    while True:
        value = input(prompt_text).strip()
        try:
            return int(value)
        except ValueError:
            print("Enter a whole number.")


def make_location(label, lat, lon, location_type=None):
    """Create one normalized location record."""
    return {
        "label": label,
        "lat": float(lat),
        "lon": float(lon),
        "type": location_type,
    }


def make_route(route_id, label, origin_key, destination_key, locations):
    """Create one route record with attached origin and destination data."""
    if origin_key not in locations:
        raise ValueError(f"Unknown origin location key: {origin_key}")
    if destination_key not in locations:
        raise ValueError(f"Unknown destination location key: {destination_key}")

    return {
        "route_id": route_id,
        "label": label,
        "origin_key": origin_key,
        "destination_key": destination_key,
        "origin": locations[origin_key],
        "destination": locations[destination_key],
    }


def request_locations_from_user():
    """Ask the user for location and route-pair settings."""
    city_name = input("City or study area name: ").strip()
    while city_name == "":
        city_name = input("City or study area name is required: ").strip()

    location_count = prompt_int("How many locations do you want to enter? ")
    locations = {}

    for index in range(location_count):
        print(f"Location {index + 1}")
        label = input("  Label: ").strip()
        while label == "":
            label = input("  Label is required: ").strip()

        location_key = normalize_name(label)
        lat = prompt_float("  Latitude: ")
        lon = prompt_float("  Longitude: ")
        location_type = input("  Type, for example airport/origin/destination (optional): ").strip()
        if location_type == "":
            location_type = None

        locations[location_key] = make_location(label, lat, lon, location_type)

    route_count = prompt_int("How many route pairs do you want to enter? ")
    routes = []

    for index in range(route_count):
        print(f"Route {index + 1}")
        label = input("  Route label: ").strip()
        while label == "":
            label = input("  Route label is required: ").strip()

        print("  Available location keys:", ", ".join(locations.keys()))
        origin_key = input("  Origin location key: ").strip()
        destination_key = input("  Destination location key: ").strip()
        route_id = normalize_name(label)
        routes.append(make_route(route_id, label, origin_key, destination_key, locations))

    return city_name, locations, routes


def set_route_locations(city_name=None, locations=None, routes=None):
    """
    Return route locations and route pairs supplied by a user or API caller.

    locations should be a dictionary keyed by location ID.
    routes should be a list with origin_key and destination_key values.
    """
    if city_name is None or locations is None or routes is None:
        city_name, locations, prepared_routes = request_locations_from_user()
    else:
        prepared_routes = []
        for route in routes:
            route_id = route.get("route_id") or normalize_name(route["label"])
            prepared_routes.append(
                make_route(
                    route_id,
                    route["label"],
                    route["origin_key"],
                    route["destination_key"],
                    locations,
                )
            )

    return {
        "study_area": select_study_area(city_name=city_name),
        "locations": locations,
        "routes": prepared_routes,
    }


if __name__ == "__main__":
    route_locations = set_route_locations()
    print(route_locations)
