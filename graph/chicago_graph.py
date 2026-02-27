import geopandas as gpd
import networkx as nx
from shapely.geometry import Point, LineString
import time


# Load grid
grid = gpd.read_file("../map/risk_grid_v5.geojson")
grid = grid.reset_index(drop=True)  # ensure 0..n-1 index

# Build graph nodes
G = nx.Graph()

for i, row in grid.iterrows():
    G.add_node(
        i,
        risk_cost=float(row["risk_cost"]),
        risk_class=row["risk_class"],
        geometry=row.geometry
    )
print("Nodes:", G.number_of_nodes())

# Add edges as Dr. Xu mentioned, they must be diagonal too
sindex = grid.sindex
geoms = grid.geometry

for i, geom in enumerate(geoms):
    cand = list(sindex.intersection(geom.bounds))
    for j in cand:
        if j <= i:
            continue
        if geom.touches(geoms.iloc[j]):
            G.add_edge(i, j)

print("Edges:", G.number_of_edges())

# Set weights = distance (meters) * (1 + risk)
grid_m = grid.to_crs("EPSG:3857")
centroids_m = grid_m.geometry.centroid

for u, v in G.edges():
    d = centroids_m.iloc[u].distance(centroids_m.iloc[v])  # meters
    r = 0.5 * (grid.loc[u, "risk_cost"] + grid.loc[v, "risk_cost"])
    G[u][v]["weight"] = float(d * (1.0 + r))

# Nearest-node lookup
def node_for_point(lat, lon):
    pt = Point(lon, lat)  # lon, lat

    cand = list(sindex.intersection(pt.bounds))

    for i in cand:
        if geoms.iloc[i].contains(pt):
            return int(i)

    # fallback: nearest centroid in meters
    pt_m = gpd.GeoSeries([pt], crs="EPSG:4326").to_crs("EPSG:3857").iloc[0]
    return int(centroids_m.distance(pt_m).idxmin())

start_lat, start_lon = 41.695923717435235, -88.12876224517822  # Clow
end_lat, end_lon     = 41.87838051825937,  -87.63905525207521  # Union Station
start_node = node_for_point(start_lat, start_lon)
end_node   = node_for_point(end_lat, end_lon)

print("Start node:", start_node, "End node:", end_node)
dij_path = nx.dijkstra_path(G, start_node, end_node, weight="weight")
print("Path nodes:", len(dij_path))

# A* algorithm
def heuristic(u, v):
    return centroids_m.iloc[u].distance(centroids_m.iloc[v])  # meters

a_path = nx.astar_path(G, start_node, end_node, heuristic=heuristic, weight="weight")
print("A* path nodes:", len(a_path))

# Export route GeoJSON
dij_line = LineString([grid.loc[n].geometry.centroid for n in dij_path])
a_line   = LineString([grid.loc[n].geometry.centroid for n in a_path])

routes = gpd.GeoDataFrame(
    {"name": ["dijkstra", "astar"]},
    geometry=[dij_line, a_line],
    crs=grid.crs
)

routes.to_file("routes.geojson", driver="GeoJSON")
print("Saved routes.geojson")

### PERFORMANCE TEST
def path_cost(G, path, weight="weight"):
    return sum(G[path[i]][path[i+1]][weight] for i in range(len(path)-1))


t0 = time.perf_counter()
p1 = nx.dijkstra_path(G, start_node, end_node, weight="weight")
t1 = time.perf_counter()

t2 = time.perf_counter()
p2 = nx.astar_path(G, start_node, end_node, heuristic=heuristic, weight="weight")
t3 = time.perf_counter()

print("Dijkstra seconds:", t1 - t0)
print("A* seconds:", t3 - t2)
print("Same path?", p1 == p2)
print("Dijkstra cost:", path_cost(G, p1))
print("A* cost:", path_cost(G, p2))

"""
Output:
Dijkstra seconds: 0.26059980000718497
A* seconds: 0.8576807000063127
Same path? False
Dijkstra cost: 3212055.1838673027
A* cost: 3212055.1838673023
"""