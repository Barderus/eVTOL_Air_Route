"""Editable settings for the St. Louis workflow."""


# Study area
WEST = -90.839767
SOUTH = 38.432831
EAST = -89.678330
NORTH = 38.902497
MAP_CENTER_LAT = 38.667664
MAP_CENTER_LON = -90.259049
MAP_ZOOM = 10
CELL_SIZE_M = 500

# Required map landmarks
LOCATIONS = {
    "st_louis_downtown": (38.627003, -90.199402),
    "st_louis_downtown_airport": (38.5703611, -90.1550833),
    "st_louis_lambert_airport": (38.7486982, -90.3700257),
    "midamerica_st_louis_airport": (38.5451731, -89.8351856),
}

# Simplified low-altitude airspace assumptions
AIRSPACE_SITES = [
    {
        "code": "STL",
        "airport": "St. Louis Lambert International Airport",
        "location": LOCATIONS["st_louis_lambert_airport"],
        "airspace_type": "Class B",
        "inner_radius": 6.0,
        "inner_unit": "NM",
        "outer_radius": 15.0,
        "outer_unit": "NM",
        "vertical_range": "surface to 4000 ft MSL",
    },
    {
        "code": "CPS",
        "airport": "St. Louis Downtown-Parks Airport",
        "location": LOCATIONS["st_louis_downtown_airport"],
        "airspace_type": "Class D",
        "inner_radius": 3.475905,
        "inner_unit": "NM",
        "outer_radius": None,
        "outer_unit": None,
        "vertical_range": "surface to 2900 ft MSL",
    },
    {
        "code": "BLV",
        "airport": "Scott AFB / MidAmerica St. Louis Airport",
        "location": LOCATIONS["midamerica_st_louis_airport"],
        "airspace_type": "Class D",
        "inner_radius": 4.257984,
        "inner_unit": "NM",
        "outer_radius": None,
        "outer_unit": None,
        "vertical_range": "surface to 3000 ft MSL",
    },
]

# Source data
POPULATION_SOURCES = [
    {
        "state": "MO",
        "population_csv": "population/population_data.csv",
        "block_group_shapefile": "population/tl_2025_29_bg/tl_2025_29_bg.shp",
    },
    {
        "state": "IL",
        "population_csv": "../Chicago/population/population_data.csv",
        "block_group_shapefile": "../Chicago/population/tl_2025_17_bg/tl_2025_17_bg.shp",
    },
]
POPULATION_GEOJSON = "population/st_louis_population_density.geojson"

# Generated outputs
RISK_GRID_GEOJSON = "maps/st_louis_risk_grid.geojson"
MAP_HTML = "maps/st_louis_map.html"
TRAFFIC_FOLDER = "traffic/output"
ROUTE_GEOJSON_FOLDER = "routes/output"
ROUTE_HTML_FOLDER = "maps"

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
