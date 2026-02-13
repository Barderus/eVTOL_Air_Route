import numpy as np
import geopandas as gpd
from shapely.geometry import box, Point

# Initial set up
CELL_SIZE_M = 500

CRS_LL = "EPSG:4326"
CRS_M  = "EPSG:3857"

# Bounds coords
WEST, SOUTH, EAST, NORTH = -88.89589344859755, 41.104190944576466, -87.52534901500378, 42.21224516288584

# AIRPORTS (Directional Corridor)
# name, lat, lon, heading_deg, A (peak risk), L_m (Length), W_m (Width)
AIRPORT_SITES = [
    ("O'Hare", 41.97807408541273, -87.90902412382081, 90, 120.0, 7000, 2000),
    ("Midway", 41.7856116663475,  -87.75331135429448, 130, 100.0, 6000, 2000),
    ("Clow International Airport", 41.695923717435235, -88.12876224517822, 90, 80.0, 4500, 1500),
    ("Lewis University Airport", 41.60690586957971, -88.09526487573747,90, 90, 3500, 1500),
    ("Brookeridge Airpark", 41.73268507397754, -87.9972105460431, 90, 90, 2500, 1500),
]

# HIGH DENSITY CITIES (High Risk Core)
# name, lat, lon, B, S_m
HIGH_CITY_SITES = [
    ("Chicago Loop", 41.878772143057425, -87.6336991787736, 120.0, 3000),
    ("Cicero", 41.85192566701417, -87.74822044539016, 90.0, 2500),
    ("Bolingbrook", 41.705768800035415, -88.08616587761664, 90.0, 2500),
    ("Schaumburg", 42.02725118627193, -88.10469329746029, 90.0, 2500),
    ("Joliet", 41.527304324777546, -88.0896799216755, 85.0, 3000),
    ("Naperville", 41.774956808956375, -88.14718471433731, 95.0, 3200),
    ("Aurora", 41.75775499242817, -88.31395741142339, 85.0, 2800),
    ("Evanston", 42.04841129067057, -87.68331671106188, 85.0, 2000),
]

# MEDIUM DENSITY AREAS
# name, lat, lon, B, S_m
MEDIUM_CITY_SITES = [
    ("Medium 1", 41.853488482488636, -87.66546887595683, 60.0, 2000),
    ("Medium 2", 41.90126747541342, -87.67576522277758, 60.0, 2000),
    ("Medium 3", 41.880065165414, -87.72209878347078, 60.0, 2000),
    ("Medium 4", 41.83201463509129, -87.6448761823154, 60.0, 2000),
    ("Medium 5", 41.93599339613172, -87.66821456844237, 60.0, 2000),
]

# NO-FLY ZONES
# name, lat, lon, radius_m
# TODO: Look for other no-fly zones in Chicago region
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

# Airport runnaway corridor model
# Rotate the runnaway so it's more realistic
def runway_uv(dx, dy, heading_deg):
    theta = np.deg2rad(heading_deg - 90.0)
    u = np.cos(theta) * dx + np.sin(theta) * dy
    v = -np.sin(theta) * dx + np.cos(theta) * dy
    return u, v

# Corridor risk computes the risk value of the corridor
def corridor_risk(u, v, A, L, W):
        # risk = peak * (along-runway decay) * (sideways decay)
    return A * np.exp(-np.abs(u)/L) * np.exp(-(v/W)**2)

# Force the airport to be risk = 0 since drones are class B
grid_m["airport_risk"] = 0.0

# Add each airport/runnaway corridor risk into the grid cell
for name, lat, lon, heading, A, L, W in AIRPORT_SITES:
    pt = gpd.GeoSeries([Point(lon, lat)], crs=CRS_LL).to_crs(CRS_M).iloc[0]
    dx = cx - pt.x
    dy = cy - pt.y
    u, v = runway_uv(dx, dy, heading) # Rotate the corridor
    grid_m["airport_risk"] += corridor_risk(u, v, A, L, W) # Evaluates the corridor risk function

# CITY RADIAL DECAY MODEL

# Let's make sure the defaulty city risk is 0
grid_m["city_risk"] = 0.0

# horizontal offset, vertical offset, peak risk, decay scale
def radial_decay(dx, dy, B, S):
    d = np.sqrt(dx**2 + dy**2) # Calculate the euclidean distance
    return B * np.exp(-d/S) # Exponential decay

# Same thing we did previous with the airport, we are doing here
# Add the risk to each high risk city and add the decay function to it
for name, lat, lon, B, S in HIGH_CITY_SITES:
    pt = gpd.GeoSeries([Point(lon, lat)], crs=CRS_LL).to_crs(CRS_M).iloc[0]
    dx = cx - pt.x
    dy = cy - pt.y
    grid_m["city_risk"] += radial_decay(dx, dy, B, S)

# Totalk risk
grid_m["risk_cost"] = grid_m["airport_risk"] + grid_m["city_risk"]

# NO-FLY OVERRIDE
grid_m["risk_class"] = "Low"

# Same thing we did for high risk and medium risk for the airport and cities
for name, lat, lon, radius in NO_FLY_SITES:
    pt = gpd.GeoSeries([Point(lon, lat)], crs=CRS_LL).to_crs(CRS_M).iloc[0]
    mask = grid_m["centroid"].distance(pt) <= radius
    grid_m.loc[mask, "risk_class"] = "No-Fly"

# FINAL CLASSIFICATION
grid_m.loc[grid_m["risk_cost"] > 30, "risk_class"] = "Medium"
grid_m.loc[grid_m["risk_cost"] > 70, "risk_class"] = "High"

# No-fly overrides everything
grid_m.loc[grid_m["risk_class"] == "No-Fly", "risk_cost"] = 9999

# EXPORT
grid_ll = grid_m.drop(columns=["centroid"]).to_crs(CRS_LL)
grid_ll.to_file("risk_grid_v2.geojson", driver="GeoJSON")
print("Saved: risk_grid_v2.geojson")
