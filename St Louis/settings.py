"""Editable settings for the St. Louis workflow."""


# Study area
WEST = None
SOUTH = None
EAST = None
NORTH = None
MAP_CENTER_LAT = None
MAP_CENTER_LON = None
MAP_ZOOM = 10
CELL_SIZE_M = 500

# Source data
POPULATION_SOURCE = "St Louis/population/source_population_data"
POPULATION_GEOJSON = "St Louis/population/st_louis_population_density.geojson"

# Generated outputs
RISK_GRID_GEOJSON = "St Louis/maps/st_louis_risk_grid.geojson"
TRAFFIC_FOLDER = "St Louis/traffic/output"
ROUTE_GEOJSON_FOLDER = "St Louis/routes/output"
ROUTE_HTML_FOLDER = "St Louis/maps"

# Route endpoints
START = None
DESTINATIONS = []

# Traffic datasets used during routing
TRAFFIC_DATASETS = []

# OpenSky export settings
TRAFFIC_CENTER_NAME = None
TRAFFIC_CENTER_LAT = None
TRAFFIC_CENTER_LON = None
TRAFFIC_RADIUS_NM = 15.0
TRAFFIC_SOUTH = None
TRAFFIC_NORTH = None
TRAFFIC_WEST = None
TRAFFIC_EAST = None
LOCAL_TIMEZONE = "America/Chicago"

# Route weights
DISTANCE_WEIGHT = 0.6
POPULATION_WEIGHT = 0.9
AIRSPACE_WEIGHT = 1.4
TRAFFIC_WEIGHT = 1.0


def validate_grid_settings():
    """Check that the study area is ready for grid development."""
    required_values = [WEST, SOUTH, EAST, NORTH]
    if any(value is None for value in required_values):
        raise ValueError(
            "Set WEST, SOUTH, EAST, and NORTH in St Louis/settings.py "
            "before building the St. Louis grid."
        )


def validate_traffic_settings():
    """Check that the OpenSky query area is configured."""
    if TRAFFIC_CENTER_NAME is None:
        raise ValueError("Set TRAFFIC_CENTER_NAME in St Louis/settings.py.")
    if TRAFFIC_CENTER_LAT is None or TRAFFIC_CENTER_LON is None:
        raise ValueError(
            "Set TRAFFIC_CENTER_LAT and TRAFFIC_CENTER_LON in "
            "St Louis/settings.py before exporting traffic data."
        )


def validate_route_settings():
    """Check that routing inputs and endpoints are configured."""
    if START is None:
        raise ValueError("Set START in St Louis/settings.py before generating routes.")
    if not DESTINATIONS:
        raise ValueError("Add at least one destination in St Louis/settings.py.")
    if not TRAFFIC_DATASETS:
        raise ValueError("Add at least one traffic dataset in St Louis/settings.py.")
