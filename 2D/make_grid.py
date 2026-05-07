import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from shapely.geometry import box, Point, Polygon

# Initial set up
CELL_SIZE_M = 500
BASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = BASE_DIR.parent

CRS_LL = "EPSG:4326"
CRS_M  = "EPSG:3857"

# Conversion constants
NM_TO_M = 1852.0
MI_TO_M = 1609.344

# Study bounds anchored to Aurora, Frankfort, and Lake Zurich.
# West bound from 41°45′50″N 88°17′24″W
# North bound from 42°07′25″N 87°55′15″W
# Aurora west, Frankfort south, Lake Zurich north; east preserves Chicago Loop.
WEST, SOUTH, EAST, NORTH = -88.34, 41.48, -87.52, 42.21

# Airports
# name, lat, lon, heading_deg, A (peak risk), L_m (Length), W_m (Width)
AIRPORT_SITES = [
    ("O'Hare", 41.97807408541273, -87.90902412382081, 90, 120.0, 7000, 2000),
    ("Midway", 41.7856116663475,  -87.75331135429448, 130, 100.0, 6000, 2000),
    # ("Clow International Airport", 41.695923717435235, -88.12876224517822, 90, 80.0, 4500, 1500),
    ("Lewis University Airport", 41.60690586957971, -88.09526487573747, 90, 90, 3500, 1500),
    ("Brookeridge Airpark", 41.73268507397754, -87.9972105460431, 90, 90, 4000, 2000),
]

# FAA-based airspace radius (2D widest extent) — for the airspace "stay-away" layer
# Airspace radii (2D) — using only the two innermost layers (from SkyVector rings)
# name -> (inner_radius_m, outer_radius_m)
AIRSPACE_RADIAL_SHELVES = {
}

# NO-FLY ZONES
# name, lat, lon, radius_m
NO_FLY_SITES = []


def bearing_degrees(center, point):
    lat1 = np.deg2rad(center[0])
    lat2 = np.deg2rad(point[0])
    dlon = np.deg2rad(point[1] - center[1])
    y = np.sin(dlon) * np.cos(lat2)
    x = (np.cos(lat1) * np.sin(lat2)) - (np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
    return (np.rad2deg(np.arctan2(y, x)) + 360.0) % 360.0


def destination_point(center, radius_m, bearing_deg):
    earth_radius_m = 6371000.0
    bearing = np.deg2rad(bearing_deg)
    lat1 = np.deg2rad(center[0])
    lon1 = np.deg2rad(center[1])
    angular_distance = radius_m / earth_radius_m

    lat2 = np.arcsin(
        (np.sin(lat1) * np.cos(angular_distance)) +
        (np.cos(lat1) * np.sin(angular_distance) * np.cos(bearing))
    )
    lon2 = lon1 + np.arctan2(
        np.sin(bearing) * np.sin(angular_distance) * np.cos(lat1),
        np.cos(angular_distance) - (np.sin(lat1) * np.sin(lat2))
    )
    return (float(np.rad2deg(lat2)), float(np.rad2deg(lon2)))


def circle_arc(center, radius_m, start_point, end_point, direction, steps):
    start_bearing = bearing_degrees(center, start_point)
    end_bearing = bearing_degrees(center, end_point)

    if direction == "cw":
        sweep = end_bearing - start_bearing if end_bearing >= start_bearing else (360.0 - start_bearing) + end_bearing
    else:
        sweep = start_bearing - end_bearing if start_bearing >= end_bearing else start_bearing + (360.0 - end_bearing)

    points = []
    for i in range(steps + 1):
        fraction = i / steps
        if direction == "cw":
            bearing = (start_bearing + (sweep * fraction)) % 360.0
        else:
            bearing = (start_bearing - (sweep * fraction) + 360.0) % 360.0
        points.append(destination_point(center, radius_m, bearing))
    return points


def circle_polygon(center, radius_m, steps=96):
    return [destination_point(center, radius_m, bearing) for bearing in np.linspace(0.0, 360.0, steps, endpoint=False)]


def sector_polygon(center, radius_m, start_bearing_deg, end_bearing_deg, steps=32):
    if end_bearing_deg < start_bearing_deg:
        end_bearing_deg += 360.0
    arc_points = [
        destination_point(center, radius_m, bearing % 360.0)
        for bearing in np.linspace(start_bearing_deg, end_bearing_deg, steps + 1)
    ]
    return [center, *arc_points, center]


def ll_polygon_to_m(points_ll):
    polygon_ll = Polygon([(lon, lat) for lat, lon in points_ll]).buffer(0)
    return gpd.GeoSeries([polygon_ll], crs=CRS_LL).to_crs(CRS_M).iloc[0].buffer(0)


def to_local_xy(point, lat0_deg):
    cos_lat = np.cos(np.deg2rad(lat0_deg))
    return (
        point[1] * 111320.0 * cos_lat,
        point[0] * 111320.0,
    )


def from_local_xy(point_xy, lat0_deg):
    cos_lat = np.cos(np.deg2rad(lat0_deg))
    return (
        point_xy[1] / 111320.0,
        point_xy[0] / (111320.0 * cos_lat),
    )


def east_line_circle_intersection(start_point, center, radius_m):
    lat0 = (start_point[0] + center[0]) / 2.0
    sx, sy = to_local_xy(start_point, lat0)
    cx_local, cy_local = to_local_xy(center, lat0)
    dy = sy - cy_local
    inside = (radius_m * radius_m) - (dy * dy)
    if inside < 0:
        raise ValueError("No east-line circle intersection for provided geometry.")
    dx = np.sqrt(inside)
    candidate1 = from_local_xy((cx_local + dx, sy), lat0)
    candidate2 = from_local_xy((cx_local - dx, sy), lat0)
    return candidate1 if candidate1[1] > candidate2[1] else candidate2


def build_ord_class_b_shelves():
    ord_poo = (41.9877777778, -87.9047222222)
    mdw_legal = (41.7861111111, -87.7525000000)

    poo5 = 5.0 * NM_TO_M
    poo6 = 6.0 * NM_TO_M
    poo10 = 10.0 * NM_TO_M
    poo10_5 = 10.5 * NM_TO_M
    mdw5 = 5.0 * NM_TO_M
    mdw10 = 10.0 * NM_TO_M

    area_a_start = (42.0694444444, -87.9252777778)
    area_a_east_arc = (41.9875000000, -87.7930555556)
    area_a_east_line = (41.9875000000, -87.7708333333)
    area_a_i290_south = (41.9533333333, -88.0322222222)
    area_a_i290_north = (42.0222222222, -88.0308333333)
    area_a_us12 = (42.0841666667, -87.9405555556)

    area_a_points = [
        area_a_start,
        *circle_arc(ord_poo, poo5, area_a_start, area_a_east_arc, "cw", 28)[1:],
        area_a_east_line,
        *circle_arc(ord_poo, poo6, area_a_east_line, area_a_i290_south, "cw", 40)[1:],
        area_a_i290_north,
        *circle_arc(ord_poo, poo6, area_a_i290_north, area_a_us12, "cw", 22)[1:],
    ]

    area_b_start = (42.0658333333, -87.8652777778)
    area_b_willow = (42.1055555556, -87.8272222222)
    area_b_willow_east = (42.1011111111, -87.7405555556)
    area_b_mdw5_on10 = (41.8261111111, -87.8500000000)
    area_b_mdw5_on10_5 = (41.8163888889, -87.8561111111)
    area_b_mdw10_on10_5 = (41.8197222222, -87.9705555556)
    area_b_mdw10_on10 = (41.8277777778, -87.9680555556)
    area_b_us12 = (42.1338888889, -88.0122222222)
    area_b_back_to5 = (42.0694444444, -87.9252777778)

    area_b_points = [
        area_b_start,
        area_b_willow,
        area_b_willow_east,
        *circle_arc(ord_poo, poo10, area_b_willow_east, area_b_mdw5_on10, "cw", 44)[1:],
        *circle_arc(mdw_legal, mdw5, area_b_mdw5_on10, area_b_mdw5_on10_5, "ccw", 12)[1:],
        *circle_arc(ord_poo, poo10_5, area_b_mdw5_on10_5, area_b_mdw10_on10_5, "cw", 22)[1:],
        *circle_arc(mdw_legal, mdw10, area_b_mdw10_on10_5, area_b_mdw10_on10, "cw", 14)[1:],
        *circle_arc(ord_poo, poo10, area_b_mdw10_on10, area_b_us12, "cw", 50)[1:],
        area_b_back_to5,
        *circle_arc(ord_poo, poo5, area_b_back_to5, area_b_start, "cw", 10)[1:],
    ]

    return ll_polygon_to_m(area_a_points), ll_polygon_to_m(area_b_points)


ORD_CLASS_B_AREA_A_M, ORD_CLASS_B_AREA_B_M = build_ord_class_b_shelves()
ORD_CLASS_B_FOOTPRINT_M = ORD_CLASS_B_AREA_A_M.union(ORD_CLASS_B_AREA_B_M).buffer(0)


def build_midway_class_c_shelves():
    ord_poo = (41.9877777778, -87.9047222222)
    mdw = (41.7861111111, -87.7525000000)

    ord10 = 10.0 * NM_TO_M
    ord10_5 = 10.5 * NM_TO_M
    mdw5 = 5.0 * NM_TO_M
    mdw10 = 10.0 * NM_TO_M

    land_start = (41.8261111111, -87.8500000000)
    shoreline_north = (41.8691666667, -87.6163888889)
    shoreline_south = (41.7497222222, -87.5350000000)
    mdw10_on_ord10_5 = (41.8197222222, -87.9705555556)
    mdw5_on_ord10_5 = (41.8163888889, -87.8561111111)

    lake_east = east_line_circle_intersection(shoreline_north, mdw, mdw10)

    land_shelf_points = [
        land_start,
        shoreline_north,
        shoreline_south,
        *circle_arc(mdw, mdw10, shoreline_south, mdw10_on_ord10_5, "cw", 48)[1:],
        *circle_arc(ord_poo, ord10_5, mdw10_on_ord10_5, mdw5_on_ord10_5, "ccw", 22)[1:],
        *circle_arc(mdw, mdw5, mdw5_on_ord10_5, land_start, "ccw", 16)[1:],
    ]

    lake_shelf_points = [
        shoreline_north,
        lake_east,
        *circle_arc(mdw, mdw10, lake_east, shoreline_south, "cw", 24)[1:],
        shoreline_south,
    ]

    core_polygon = ll_polygon_to_m(circle_polygon(mdw, mdw5))
    land_shelf_polygon = ll_polygon_to_m(land_shelf_points)
    lake_shelf_polygon = ll_polygon_to_m(lake_shelf_points)
    shelf_polygon = land_shelf_polygon.union(lake_shelf_polygon).buffer(0)

    # The legal Class C excludes overlapping Chicago Class B airspace.
    core_polygon = core_polygon.difference(ORD_CLASS_B_FOOTPRINT_M).buffer(0)
    shelf_polygon = shelf_polygon.difference(ORD_CLASS_B_FOOTPRINT_M).buffer(0)
    return core_polygon, shelf_polygon


MIDWAY_CLASS_C_CORE_M, MIDWAY_CLASS_C_SHELF_M = build_midway_class_c_shelves()


def build_lewis_class_d_polygon():
    lewis = (41.60690586957971, -88.09526487573747)
    radius = 4.1 * NM_TO_M

    # FAA LTA-C90-17 describes Lewis as a 4.1 NM Class D with a cutout along
    # I-55 for Clow traffic. Without published polygon vertices in the repo,
    # model that cutout as a northwest-facing wedge toward Clow.
    full_circle = ll_polygon_to_m(circle_polygon(lewis, radius, steps=128))
    cutout_wedge = ll_polygon_to_m(sector_polygon(lewis, radius, 322.0, 18.0, steps=28))
    return full_circle.difference(cutout_wedge).buffer(0)


LEWIS_CLASS_D_M = build_lewis_class_d_polygon()

# Build the study grid once in projected meters so every downstream cost uses
# consistent cell sizes and distance measurements.
study_ll = gpd.GeoDataFrame(geometry=[box(WEST, SOUTH, EAST, NORTH)], crs=CRS_LL)
study_m = study_ll.to_crs(CRS_M)

xmin, ymin, xmax, ymax = study_m.total_bounds

xs = np.arange(xmin, xmax + CELL_SIZE_M, CELL_SIZE_M)
ys = np.arange(ymin, ymax + CELL_SIZE_M, CELL_SIZE_M)

cells = [box(x, y, x + CELL_SIZE_M, y + CELL_SIZE_M)
         for x in xs[:-1] for y in ys[:-1]]

grid_m = gpd.GeoDataFrame({"geometry": cells}, crs=CRS_M)
grid_m = gpd.clip(grid_m, study_m).reset_index(drop=True)

grid_m["centroid"] = grid_m.geometry.centroid
cx = grid_m["centroid"].x.to_numpy()
cy = grid_m["centroid"].y.to_numpy()

# Join population density by cell centroid so each cell gets a single dominant
# census block group instead of blending neighboring polygons.
BG_GEOJSON_PATH = REPO_ROOT / "geojson" / "il_blockgroups_population_density.geojson"
centroids = gpd.GeoDataFrame(
    {"cell_id": np.arange(len(grid_m))},
    geometry=grid_m["centroid"],
    crs=CRS_M
)

if BG_GEOJSON_PATH.exists():
    bg = gpd.read_file(BG_GEOJSON_PATH).to_crs(CRS_M)
    bg = gpd.clip(bg, study_m)
    join = gpd.sjoin(
        centroids,
        bg[["density_p_km2", "density_risk", "geometry"]],
        how="left",
        predicate="within"
    )
    join = join.set_index("cell_id").reindex(np.arange(len(grid_m)))
    grid_m["density_p_km2"] = join["density_p_km2"].astype(float).fillna(0.0)
    grid_m["density_risk"] = join["density_risk"].fillna("low")
else:
    prior_grid = gpd.read_file(REPO_ROOT / "geojson" / "risk_grid_v6.geojson").to_crs(CRS_M)
    prior_join = gpd.sjoin(
        centroids,
        prior_grid[["density_p_km2", "density_risk", "city_risk", "pop_class", "geometry"]],
        how="left",
        predicate="within"
    )
    prior_join = prior_join.set_index("cell_id").reindex(np.arange(len(grid_m)))
    grid_m["density_p_km2"] = prior_join["density_p_km2"].astype(float).fillna(0.0)
    grid_m["density_risk"] = prior_join["density_risk"].fillna("low")

# CITY RISK (from density)
grid_m["city_risk"] = 0.0
DENSITY_SCALE = 10.0
grid_m["city_risk"] = (grid_m["density_p_km2"] / DENSITY_SCALE).clip(0, 120)

# Model runway-aligned approach/departure corridors separately from the radial
# airspace rings so both effects can be tuned independently.
def runway_uv(dx, dy, heading_deg):
    theta = np.deg2rad(heading_deg - 90.0)
    u = np.cos(theta) * dx + np.sin(theta) * dy
    v = -np.sin(theta) * dx + np.cos(theta) * dy
    return u, v

def corridor_risk(u, v, A, L, W):
    # risk = peak * (along-runway decay) * (sideways decay)
    return A * np.exp(-np.abs(u)/L) * np.exp(-(v/W)**2)

grid_m["airport_risk"] = 0.0

for name, lat, lon, heading, A, L, W in AIRPORT_SITES:
    pt = gpd.GeoSeries([Point(lon, lat)], crs=CRS_LL).to_crs(CRS_M).iloc[0]
    dx = cx - pt.x
    dy = cy - pt.y
    u, v = runway_uv(dx, dy, heading)
    grid_m["airport_risk"] += corridor_risk(u, v, A, L, W)

# Thresholds + classes (existing)
grid_m["pop_class"] = "Low"
grid_m["air_class"] = "Low"
grid_m["density_type"] = "low"

POP_MED_TH = 40
POP_HIGH_TH = 90

AIR_MED_TH = 30
AIR_HIGH_TH = 70

grid_m.loc[grid_m["city_risk"] > POP_MED_TH, "pop_class"] = "Medium"
grid_m.loc[grid_m["city_risk"] > POP_HIGH_TH, "pop_class"] = "High"

# Add airport airspace structure on top of the corridor risk:
# ORD uses FAA-like Class B shelf polygons; Midway and Lewis keep radial
# approximations but follow the same high / medium / outside-decay scoring.
def radial_decay(d_m, A, S_m):
    return A * np.exp(-d_m / S_m)

AIRSPACE_HIGH_VAL = 1.0
AIRSPACE_MED_PEAK = 0.6
AIRSPACE_DECAY_S_M = 8000.0
AIRSPACE_DECAY_MAX_DISTANCE_M = 3 * CELL_SIZE_M

grid_m["airport_airspace_high"] = 0.0
grid_m["airport_airspace_med"]  = 0.0

centroid_series = gpd.GeoSeries(grid_m["centroid"], crs=CRS_M)

# ORD Class B: Area A is the high-cost core, Area B is the medium-cost shelf,
# and cost decays outward from the Area B outer boundary.
in_ord_area_a = centroid_series.within(ORD_CLASS_B_AREA_A_M) | centroid_series.touches(ORD_CLASS_B_AREA_A_M)
in_ord_area_b = centroid_series.within(ORD_CLASS_B_AREA_B_M) | centroid_series.touches(ORD_CLASS_B_AREA_B_M)
in_ord_medium = in_ord_area_b & ~in_ord_area_a
outside_ord_area_b = ~in_ord_area_b
ord_boundary_distance = centroid_series.distance(ORD_CLASS_B_AREA_B_M.boundary).to_numpy()

grid_m.loc[in_ord_area_a, "airport_airspace_high"] = AIRSPACE_HIGH_VAL
grid_m.loc[in_ord_medium, "airport_airspace_med"] = np.maximum(
    grid_m.loc[in_ord_medium, "airport_airspace_med"].to_numpy(),
    AIRSPACE_MED_PEAK
)

ord_decay = np.zeros(len(grid_m), dtype=float)
ord_decay_mask = outside_ord_area_b.to_numpy() & (ord_boundary_distance <= AIRSPACE_DECAY_MAX_DISTANCE_M)
ord_decay[ord_decay_mask] = radial_decay(
    ord_boundary_distance[ord_decay_mask],
    A=AIRSPACE_MED_PEAK,
    S_m=AIRSPACE_DECAY_S_M
)
grid_m["airport_airspace_med"] = np.maximum(grid_m["airport_airspace_med"].to_numpy(), ord_decay)

# Midway Class C: 5 NM core is high cost, the shelf is medium cost, and cost
# decays outward from the outer shelf boundary.
in_mdw_core = centroid_series.within(MIDWAY_CLASS_C_CORE_M) | centroid_series.touches(MIDWAY_CLASS_C_CORE_M)
in_mdw_shelf = centroid_series.within(MIDWAY_CLASS_C_SHELF_M) | centroid_series.touches(MIDWAY_CLASS_C_SHELF_M)
in_mdw_medium = in_mdw_shelf & ~in_mdw_core
outside_mdw_shelf = ~in_mdw_shelf
mdw_boundary_distance = centroid_series.distance(MIDWAY_CLASS_C_SHELF_M.boundary).to_numpy()

grid_m.loc[in_mdw_core, "airport_airspace_high"] = AIRSPACE_HIGH_VAL
grid_m.loc[in_mdw_medium, "airport_airspace_med"] = np.maximum(
    grid_m.loc[in_mdw_medium, "airport_airspace_med"].to_numpy(),
    AIRSPACE_MED_PEAK
)

mdw_decay = np.zeros(len(grid_m), dtype=float)
mdw_decay_mask = outside_mdw_shelf.to_numpy() & (mdw_boundary_distance <= AIRSPACE_DECAY_MAX_DISTANCE_M)
mdw_decay[mdw_decay_mask] = radial_decay(
    mdw_boundary_distance[mdw_decay_mask],
    A=AIRSPACE_MED_PEAK,
    S_m=AIRSPACE_DECAY_S_M
)
grid_m["airport_airspace_med"] = np.maximum(grid_m["airport_airspace_med"].to_numpy(), mdw_decay)

# Lewis Class D: single controlled polygon with a Clow cutout. Treat the
# interior as high cost and decay outward from the Class D boundary.
in_lewis_class_d = centroid_series.within(LEWIS_CLASS_D_M) | centroid_series.touches(LEWIS_CLASS_D_M)
outside_lewis_class_d = ~in_lewis_class_d
lewis_boundary_distance = centroid_series.distance(LEWIS_CLASS_D_M.boundary).to_numpy()

grid_m.loc[in_lewis_class_d, "airport_airspace_high"] = AIRSPACE_HIGH_VAL

lewis_decay = np.zeros(len(grid_m), dtype=float)
lewis_decay_mask = outside_lewis_class_d.to_numpy() & (lewis_boundary_distance <= AIRSPACE_DECAY_MAX_DISTANCE_M)
lewis_decay[lewis_decay_mask] = radial_decay(
    lewis_boundary_distance[lewis_decay_mask],
    A=AIRSPACE_MED_PEAK,
    S_m=AIRSPACE_DECAY_S_M
)
grid_m["airport_airspace_med"] = np.maximum(grid_m["airport_airspace_med"].to_numpy(), lewis_decay)

for name, lat, lon, heading, A, L, W in AIRPORT_SITES:
    if name not in AIRSPACE_RADIAL_SHELVES:
        continue

    r_inner_m, r_outer_m = AIRSPACE_RADIAL_SHELVES[name]
    r_inner_m = float(r_inner_m)
    r_outer_m = float(r_outer_m)

    pt = gpd.GeoSeries([Point(lon, lat)], crs=CRS_LL).to_crs(CRS_M).iloc[0]

    dx = cx - pt.x
    dy = cy - pt.y
    d  = np.sqrt(dx*dx + dy*dy)

    # HIGH: inside inner radius
    in_high = d <= r_inner_m
    grid_m.loc[in_high, "airport_airspace_high"] = np.maximum(
        grid_m.loc[in_high, "airport_airspace_high"].to_numpy(),
        AIRSPACE_HIGH_VAL
    )

    # MEDIUM: between inner and outer radius
    in_med_band = (d > r_inner_m) & (d <= r_outer_m)
    med_val = np.zeros_like(d, dtype=float)
    med_val[in_med_band] = AIRSPACE_MED_PEAK

    # Outside the medium shelf, decay away from the outer boundary.
    outside_outer = d > r_outer_m
    decay_distance = d - r_outer_m
    decay_mask = outside_outer & (decay_distance <= AIRSPACE_DECAY_MAX_DISTANCE_M)
    med_val[decay_mask] = radial_decay(decay_distance[decay_mask], A=AIRSPACE_MED_PEAK, S_m=AIRSPACE_DECAY_S_M)

    grid_m["airport_airspace_med"] = np.maximum(
        grid_m["airport_airspace_med"].to_numpy(),
        med_val
    )

grid_m["airport_airspace_risk"] = np.maximum(
    grid_m["airport_airspace_high"].to_numpy(),
    grid_m["airport_airspace_med"].to_numpy()
)


# Combine runway corridor risk + airspace radius risk
# Scale airspace layer to your current airport thresholds:
# inside FAA radius should be strong enough to register as HIGH airspace.
AIRSPACE_WEIGHT = AIR_HIGH_TH  # 70

grid_m["airport_risk_combined"] = grid_m["airport_risk"] + AIRSPACE_WEIGHT * grid_m["airport_airspace_risk"]

# Now classify airspace based on the combined airport score
grid_m.loc[grid_m["airport_risk_combined"] > AIR_MED_TH, "air_class"] = "Medium"
grid_m.loc[grid_m["airport_risk_combined"] > AIR_HIGH_TH, "air_class"] = "High"

# Keep the dominant source so the map can explain whether a cell is driven by
# population exposure or airport constraints.
dominant_mask = (grid_m["pop_class"] != "Low") | (grid_m["air_class"] != "Low")
grid_m.loc[dominant_mask & (grid_m["city_risk"] >= grid_m["airport_risk_combined"]), "density_type"] = "population"
grid_m.loc[dominant_mask & (grid_m["airport_risk_combined"] > grid_m["city_risk"]), "density_type"] = "airspace"

# Total risk cost
grid_m["risk_cost"] = grid_m["airport_risk_combined"] + grid_m["city_risk"]

# NO-FLY OVERRIDE
grid_m["risk_class"] = "Low"

for name, lat, lon, radius in NO_FLY_SITES:
    pt = gpd.GeoSeries([Point(lon, lat)], crs=CRS_LL).to_crs(CRS_M).iloc[0]
    mask = grid_m["centroid"].distance(pt) <= radius
    grid_m.loc[mask, "risk_class"] = "No-Fly"
    grid_m.loc[mask, "density_type"] = "no-fly"
    grid_m.loc[mask, "pop_class"] = "Low"
    grid_m.loc[mask, "air_class"] = "Low"

# FINAL CLASSIFICATION
grid_m.loc[grid_m["risk_cost"] > 30, "risk_class"] = "Medium"
grid_m.loc[grid_m["risk_cost"] > 70, "risk_class"] = "High"

# Freeze no-fly cells at an overwhelming cost so route search will never prefer
# them over ordinary risk tradeoffs.
grid_m.loc[grid_m["risk_class"] == "No-Fly", "risk_cost"] = 9999

# EXPORT
OUTPUT_GEOJSON_PATH = REPO_ROOT / "geojson" / "risk_grid_v6.geojson"
grid_ll = grid_m.drop(columns=["centroid"]).to_crs(CRS_LL)
grid_ll.to_file(OUTPUT_GEOJSON_PATH, driver="GeoJSON")
print("Saved:", OUTPUT_GEOJSON_PATH)
