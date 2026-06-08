"""Editable settings for the St. Louis workflow."""


# Study area
WEST = -90.501400
SOUTH = 38.432831
EAST = -89.678330
NORTH = 38.863100
MAP_CENTER_LAT = 38.719886
MAP_CENTER_LON = -90.259049
TRAFFIC_MAP_CENTER_LAT = 38.719886
TRAFFIC_MAP_CENTER_LON = -90.340000
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
TRAFFIC_MAP_HTML = "maps/st_louis_traffic_leaflet.html"
ROUTE_MAP_HTML = "maps/st_louis_routes_map.html"
ROUTE_ANALYSIS_CSV = "routes/output/st_louis_route_cost_analysis.csv"
ROUTE_ANALYSIS_FOLDER = "routes/analysis"
TRAFFIC_FOLDER = "traffic/output"
ROUTE_GEOJSON_FOLDER = "routes/output"
ROUTE_HTML_FOLDER = "maps"

# Route endpoints
ROUTES = [
    {
        "label": "MidAmerica to Downtown St. Louis",
        "start_label": "MidAmerica St. Louis Airport",
        "start": LOCATIONS["midamerica_st_louis_airport"],
        "destination_label": "Downtown St. Louis",
        "destination": LOCATIONS["st_louis_downtown"],
    },
    {
        "label": "MidAmerica to St. Louis Lambert",
        "start_label": "MidAmerica St. Louis Airport",
        "start": LOCATIONS["midamerica_st_louis_airport"],
        "destination_label": "St. Louis Lambert International Airport",
        "destination": LOCATIONS["st_louis_lambert_airport"],
    },
    {
        "label": "St. Louis Downtown Airport to St. Louis Lambert",
        "start_label": "St. Louis Downtown Airport",
        "start": LOCATIONS["st_louis_downtown_airport"],
        "destination_label": "St. Louis Lambert International Airport",
        "destination": LOCATIONS["st_louis_lambert_airport"],
    },
]

# Traffic datasets used during routing
TRAFFIC_DATASETS = [
    {
        "date": "2026-03-07",
        "label": "March 7th - Saturday",
        "filename": "st_louis_2026-03-07_1s_15nm_bbox.csv",
    },
    {
        "date": "2026-03-09",
        "label": "March 9th - Monday",
        "filename": "st_louis_2026-03-09_1s_15nm_bbox.csv",
    },
    {
        "date": "2026-01-10",
        "label": "January 10th - Saturday",
        "filename": "st_louis_2026-01-10_1s_15nm_bbox.csv",
    },
    {
        "date": "2026-01-12",
        "label": "January 12th - Monday",
        "filename": "st_louis_2026-01-12_1s_15nm_bbox.csv",
    },
    {
        "date": "2025-07-12",
        "label": "July 12th - Saturday",
        "filename": "st_louis_2025-07-12_1s_15nm_bbox.csv",
    },
    {
        "date": "2025-07-14",
        "label": "July 14th - Monday",
        "filename": "st_louis_2025-07-14_1s_15nm_bbox.csv",
    },
]

# OpenSky export settings
TRAFFIC_CENTER_NAME = "st_louis"
TRAFFIC_CENTER_LAT = MAP_CENTER_LAT
TRAFFIC_CENTER_LON = MAP_CENTER_LON
TRAFFIC_RADIUS_NM = 15.0
TRAFFIC_SOUTH = SOUTH
TRAFFIC_NORTH = NORTH
TRAFFIC_WEST = WEST
TRAFFIC_EAST = EAST
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
    if not ROUTES:
        raise ValueError("Add at least one route in St Louis/settings.py.")
    if not TRAFFIC_DATASETS:
        raise ValueError("Add at least one traffic dataset in St Louis/settings.py.")
