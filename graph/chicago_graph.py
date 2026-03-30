import math
import time
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
from shapely.geometry import LineString, Point


GRID_PATH = Path(__file__).resolve().parent.parent / "map" / "risk_grid_v5.geojson"
OUTPUT_PATH = Path(__file__).resolve().parent / "routes.geojson"

# Keep all normalized terms in a comparable range, but bias the combined route
# toward avoiding constrained airspace more than distance or population.
DISTANCE_WEIGHT = 0.8
POPULATION_WEIGHT = 0.9
AIRSPACE_WEIGHT = 1.4

ROUTE_SPECS = [
    {
        "name": "combined",
        "label": "Combined Cost",
        "distance_weight": DISTANCE_WEIGHT,
        "population_weight": POPULATION_WEIGHT,
        "airspace_weight": AIRSPACE_WEIGHT,
    },
    {
        "name": "population_only",
        "label": "Population Only",
        "distance_weight": 0.0,
        "population_weight": 1.0,
        "airspace_weight": 0.0,
    },
    {
        "name": "distance_only",
        "label": "Distance Only",
        "distance_weight": 1.0,
        "population_weight": 0.0,
        "airspace_weight": 0.0,
    },
    {
        "name": "airspace_only",
        "label": "Airspace Only",
        "distance_weight": 0.0,
        "population_weight": 0.0,
        "airspace_weight": 1.0,
    },
]


def normalize(values: np.ndarray) -> tuple[np.ndarray, float]:
    maximum = float(np.max(values))
    if maximum <= 0:
        return np.zeros_like(values, dtype=float), 1.0
    return values / maximum, maximum


def path_cost(graph: nx.Graph, path: list[int], weight: str) -> float:
    return sum(graph[path[i]][path[i + 1]][weight] for i in range(len(path) - 1))


def build_graph(grid: gpd.GeoDataFrame) -> tuple[nx.Graph, gpd.GeoSeries, np.ndarray, np.ndarray]:
    graph = nx.Graph()
    for i, row in grid.iterrows():
        graph.add_node(
            i,
            risk_cost=float(row["risk_cost"]),
            city_risk=float(row["city_risk"]),
            airport_risk_combined=float(row["airport_risk_combined"]),
            risk_class=row["risk_class"],
            geometry=row.geometry,
        )

    sindex = grid.sindex
    geoms = grid.geometry
    for i, geom in enumerate(geoms):
        for j in sindex.intersection(geom.bounds):
            if j <= i:
                continue
            if geom.touches(geoms.iloc[j]):
                graph.add_edge(i, j)

    grid_m = grid.to_crs("EPSG:3857")
    centroids_m = grid_m.geometry.centroid
    cent_x = centroids_m.x.to_numpy()
    cent_y = centroids_m.y.to_numpy()
    return graph, centroids_m, cent_x, cent_y


def assign_edge_weights(
    graph: nx.Graph,
    cent_x: np.ndarray,
    cent_y: np.ndarray,
    city_risk_norm: np.ndarray,
    airspace_risk_norm: np.ndarray,
) -> float:
    neighbor_distances = []
    for u, v in graph.edges():
        dx = cent_x[u] - cent_x[v]
        dy = cent_y[u] - cent_y[v]
        neighbor_distances.append(math.hypot(dx, dy))

    max_edge_distance = max(neighbor_distances)

    for u, v in graph.edges():
        dx = cent_x[u] - cent_x[v]
        dy = cent_y[u] - cent_y[v]
        distance_norm = math.hypot(dx, dy) / max_edge_distance
        population_norm = 0.5 * (city_risk_norm[u] + city_risk_norm[v])
        airspace_norm = 0.5 * (airspace_risk_norm[u] + airspace_risk_norm[v])

        for spec in ROUTE_SPECS:
            graph[u][v][f"weight_{spec['name']}"] = (
                spec["distance_weight"] * distance_norm
                + spec["population_weight"] * population_norm
                + spec["airspace_weight"] * airspace_norm
            )

    return max_edge_distance


def node_for_point(
    lat: float,
    lon: float,
    sindex,
    geoms: gpd.GeoSeries,
    centroids_m: gpd.GeoSeries,
) -> int:
    pt = Point(lon, lat)
    for i in sindex.intersection(pt.bounds):
        if geoms.iloc[i].contains(pt):
            return int(i)

    pt_m = gpd.GeoSeries([pt], crs="EPSG:4326").to_crs("EPSG:3857").iloc[0]
    return int(centroids_m.distance(pt_m).idxmin())


def make_heuristic(distance_weight: float, cent_x: np.ndarray, cent_y: np.ndarray, max_edge_distance: float):
    if distance_weight <= 0:
        return lambda _u, _v: 0.0

    def heuristic(u: int, v: int) -> float:
        dx = cent_x[u] - cent_x[v]
        dy = cent_y[u] - cent_y[v]
        return distance_weight * (math.hypot(dx, dy) / max_edge_distance)

    return heuristic


def route_line(grid: gpd.GeoDataFrame, path: list[int]) -> LineString:
    return LineString([grid.loc[n].geometry.centroid for n in path])


def main() -> None:
    grid = gpd.read_file(GRID_PATH).reset_index(drop=True)
    graph, centroids_m, cent_x, cent_y = build_graph(grid)

    city_risk = grid["city_risk"].to_numpy(dtype=float)
    airspace_risk = grid["airport_risk_combined"].to_numpy(dtype=float)
    city_risk_norm, city_risk_max = normalize(city_risk)
    airspace_risk_norm, airspace_risk_max = normalize(airspace_risk)
    max_edge_distance = assign_edge_weights(graph, cent_x, cent_y, city_risk_norm, airspace_risk_norm)

    print("Nodes:", graph.number_of_nodes())
    print("Edges:", graph.number_of_edges())
    print("Normalization maxima:")
    print("  distance max edge (m):", max_edge_distance)
    print("  city_risk max:", city_risk_max)
    print("  airport_risk_combined max:", airspace_risk_max)
    print("Route weights:")
    print("  combined:", DISTANCE_WEIGHT, POPULATION_WEIGHT, AIRSPACE_WEIGHT)

    sindex = grid.sindex
    geoms = grid.geometry

    start_lat, start_lon = 41.695923717435235, -88.12876224517822  # Clow
    end_lat, end_lon = 41.87838051825937, -87.63905525207521  # Union Station
    start_node = node_for_point(start_lat, start_lon, sindex, geoms, centroids_m)
    end_node = node_for_point(end_lat, end_lon, sindex, geoms, centroids_m)
    print("Start node:", start_node, "End node:", end_node)

    route_rows = []
    route_geometries = []

    for spec in ROUTE_SPECS:
        weight_key = f"weight_{spec['name']}"
        heuristic = make_heuristic(spec["distance_weight"], cent_x, cent_y, max_edge_distance)

        t0 = time.perf_counter()
        path = nx.astar_path(graph, start_node, end_node, heuristic=heuristic, weight=weight_key)
        elapsed = time.perf_counter() - t0
        total_cost = path_cost(graph, path, weight_key)

        route_rows.append(
            {
                "name": spec["name"],
                "label": spec["label"],
                "distance_weight": spec["distance_weight"],
                "population_weight": spec["population_weight"],
                "airspace_weight": spec["airspace_weight"],
                "path_nodes": len(path),
                "total_cost": total_cost,
                "algorithm": "astar",
            }
        )
        route_geometries.append(route_line(grid, path))

        print(
            f"{spec['name']}:",
            f"nodes={len(path)}",
            f"cost={total_cost:.4f}",
            f"seconds={elapsed:.4f}",
        )

    routes = gpd.GeoDataFrame(route_rows, geometry=route_geometries, crs=grid.crs)
    routes.to_file(OUTPUT_PATH, driver="GeoJSON")
    print(f"Saved {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
