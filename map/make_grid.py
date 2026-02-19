import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box, Point

# Initial set up
CELL_SIZE_M = 500

CRS_LL = "EPSG:4326"
CRS_M  = "EPSG:3857"

# Conversion constants
NM_TO_M = 1852.0
MI_TO_M = 1609.344

# Bounds coords
WEST, SOUTH, EAST, NORTH = -88.89589344859755, 41.104190944576466, -87.52534901500378, 42.21224516288584

# Airports
# name, lat, lon, heading_deg, A (peak risk), L_m (Length), W_m (Width)
AIRPORT_SITES = [
    ("O'Hare", 41.97807408541273, -87.90902412382081, 90, 120.0, 7000, 2000),
    ("Midway", 41.7856116663475,  -87.75331135429448, 130, 100.0, 6000, 2000),
    # ("Clow International Airport", 41.695923717435235, -88.12876224517822, 90, 80.0, 4500, 1500),
    ("Lewis University Airport", 41.60690586957971, -88.09526487573747, 90, 90, 3500, 1500),
    ("Brookeridge Airpark", 41.73268507397754, -87.9972105460431, 90, 90, 2500, 1500),
]

# FAA-based airspace radius (2D widest extent) — for the airspace "stay-away" layer
# Airspace radii (2D) — using only the two innermost layers (from SkyVector rings)
# name -> (inner_radius_m, outer_radius_m)
AIRSPACE_RADII_M = {
    "O'Hare": (10.0 * NM_TO_M, 15.0 * NM_TO_M),
    "Midway": (5.0 * NM_TO_M, 10.0 * NM_TO_M),
    # Lewis only has one controlled radius in our simplified model
    "Lewis University Airport": (4.1 * MI_TO_M, 4.1 * MI_TO_M),
}

# NO-FLY ZONES
# name, lat, lon, radius_m
NO_FLY_SITES = [
    ("Argonne National Lab", 41.7150, -87.9830, 2500),
    ("Fermilab", 41.8407, -88.2620, 2500),
]

# Build the grid using geopandas
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

# Load + join population density to grid (centroid-in-polygon)
BG_GEOJSON_PATH = "../population/il_blockgroups_population_density.geojson"
bg = gpd.read_file(BG_GEOJSON_PATH).to_crs(CRS_M)

# Clip block groups to the study area
bg = gpd.clip(bg, study_m)

centroids = gpd.GeoDataFrame(
    {"cell_id": np.arange(len(grid_m))},
    geometry=grid_m["centroid"],
    crs=CRS_M
)

join = gpd.sjoin(
    centroids,
    bg[["density_p_km2", "density_risk", "geometry"]],
    how="left",
    predicate="within"
)

join = join.set_index("cell_id").reindex(np.arange(len(grid_m)))

grid_m["density_p_km2"] = join["density_p_km2"].astype(float).fillna(0.0)
grid_m["density_risk"]  = join["density_risk"].fillna("low")

# CITY RISK (from density)
grid_m["city_risk"] = 0.0
DENSITY_SCALE = 10.0
grid_m["city_risk"] = (grid_m["density_p_km2"] / DENSITY_SCALE).clip(0, 120)

# Runway corridor model (existing)
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

# 2D Airspace radius risk (two innermost rings)
# High inside inner radius, medium decays between inner and outer
def radial_decay(d_m, A, S_m):
    return A * np.exp(-d_m / S_m)

AIRSPACE_HIGH_VAL = 1.0
AIRSPACE_MED_PEAK = 0.6
AIRSPACE_DECAY_S_M = 8000.0

grid_m["airport_airspace_high"] = 0.0
grid_m["airport_airspace_med"]  = 0.0

for name, lat, lon, heading, A, L, W in AIRPORT_SITES:
    if name not in AIRSPACE_RADII_M:
        continue

    r_inner_m, r_outer_m = AIRSPACE_RADII_M[name]
    r_inner_m = float(r_inner_m)
    r_outer_m = float(r_outer_m)

    pt = gpd.GeoSeries([Point(lon, lat)], crs=CRS_LL).to_crs(CRS_M).iloc[0]

    dx = cx - pt.x
    dy = cy - pt.y
    d  = np.sqrt(dx*dx + dy*dy)

    # HIGH: inside inner radius
    in_high = d <= r_inner_m
    grid_m.loc[in_high, "airport_airspace_high"] = AIRSPACE_HIGH_VAL

    # MEDIUM: only between inner and outer radius
    in_med_band = (d > r_inner_m) & (d <= r_outer_m)
    d_from_inner = d - r_inner_m

    med_val = np.zeros_like(d, dtype=float)
    med_val[in_med_band] = radial_decay(d_from_inner[in_med_band], A=AIRSPACE_MED_PEAK, S_m=AIRSPACE_DECAY_S_M)

    # Combine airports via MAX
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

# Dominant type selection (population vs airspace)
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

# No-fly overrides everything
grid_m.loc[grid_m["risk_class"] == "No-Fly", "risk_cost"] = 9999

# Debug prints
print("Density stats (grid):")
print(grid_m["density_risk"].value_counts())

print("Airspace high cells:", int((grid_m["airport_airspace_high"] > 0).sum()))
print("Airspace med  cells:", int((grid_m["airport_airspace_med"] > 0).sum()))
print("Combined airport risk stats:")
print(grid_m["airport_risk_combined"].describe())

# EXPORT
grid_ll = grid_m.drop(columns=["centroid"]).to_crs(CRS_LL)
grid_ll.to_file("risk_grid_v5.geojson", driver="GeoJSON")
print("Saved: risk_grid_v5.geojson")
