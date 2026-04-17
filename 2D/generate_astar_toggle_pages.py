from __future__ import annotations

import json
import math
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point


ROOT = Path(__file__).resolve().parent.parent
GRID_PATH = ROOT / "geojson" / "risk_grid_v5.geojson"
OPENSKY_OUTPUT = ROOT / "opensky" / "output"
GEOJSON_OUTPUT = ROOT / "geojson"
HTML_OUTPUT = ROOT / "html"

DISTANCE_WEIGHT = 0.8
POPULATION_WEIGHT = 0.9
AIRSPACE_WEIGHT = 1.4
TRAFFIC_WEIGHT = 1.0

# These specs drive both the A* edge weighting and the layer metadata emitted
# into the self-contained HTML pages.
ROUTE_SPECS = [
    {
        "name": "combined",
        "label": "Combined",
        "distance_weight": DISTANCE_WEIGHT,
        "population_weight": POPULATION_WEIGHT,
        "airspace_weight": AIRSPACE_WEIGHT,
        "traffic_weight": TRAFFIC_WEIGHT,
        "color": "#dc2626",
    },
    {
        "name": "airspace_only",
        "label": "Airspace Only",
        "distance_weight": 0.0,
        "population_weight": 0.0,
        "airspace_weight": 1.0,
        "traffic_weight": 0.0,
        "color": "#f59e0b",
    },
    {
        "name": "flight_density_only",
        "label": "Flight Density Only",
        "distance_weight": 0.0,
        "population_weight": 0.0,
        "airspace_weight": 0.0,
        "traffic_weight": 1.0,
        "color": "#7c3aed",
    },
    {
        "name": "population_only",
        "label": "Population Only",
        "distance_weight": 0.0,
        "population_weight": 1.0,
        "airspace_weight": 0.0,
        "traffic_weight": 0.0,
        "color": "#2563eb",
    },
    {
        "name": "distance_only",
        "label": "Distance Only",
        "distance_weight": 1.0,
        "population_weight": 0.0,
        "airspace_weight": 0.0,
        "traffic_weight": 0.0,
        "color": "#16a34a",
    },
]

START = {
    "label": "Clow International Airport",
    "lat": 41.695923717435235,
    "lon": -88.12876224517822,
}

DESTINATIONS = [
    {
        "slug": "union_station",
        "label": "Union Station",
        "lat": 41.87838051825937,
        "lon": -87.63905525207521,
    },
    {
        "slug": "ohare",
        "label": "O'Hare International Airport",
        "lat": 41.97807408541273,
        "lon": -87.90902412382081,
    },
    {
        "slug": "midway",
        "label": "Midway International Airport",
        "lat": 41.7856116663475,
        "lon": -87.75331135429448,
    },
]

DATASETS = [
    {
        "slug": "2026-03-07",
        "label": "2026-03-07",
        "csv_path": OPENSKY_OUTPUT / "ohare_2026-03-07_1s_15nm_bbox.csv",
    },
    {
        "slug": "2026-03-09",
        "label": "2026-03-09",
        "csv_path": OPENSKY_OUTPUT / "ohare_2026-03-09_1s_15nm_bbox.csv",
    },
    {
        "slug": "2026-01-10",
        "label": "2026-01-10",
        "csv_path": OPENSKY_OUTPUT / "ohare_2026-01-10_1s_15nm_bbox.csv",
    },
    {
        "slug": "2026-01-12",
        "label": "2026-01-12",
        "csv_path": OPENSKY_OUTPUT / "ohare_2026-01-12_1s_15nm_bbox.csv",
    },
    {
        "slug": "2025-07-14",
        "label": "2025-07-14",
        "csv_path": OPENSKY_OUTPUT / "ohare_2025-07-14_1s_15nm_bbox.csv",
    },
    {
        "slug": "2025-07-12",
        "label": "2025-07-12",
        "csv_path": OPENSKY_OUTPUT / "ohare_2025-07-12_1s_15nm_bbox.csv",
    },
]

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{page_title}</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
    integrity="sha256-p4NxAoJBhIIN+hmNHrzRCf9tD/miZyoHS5obTRR9BMY="
    crossorigin=""
  />
  <style>
    html, body {{ height: 100%; margin: 0; }}
    #map {{ height: 100%; width: 100%; }}

    .date-panel {{
      position: absolute;
      top: 12px;
      left: 12px;
      z-index: 1000;
      background: white;
      padding: 10px 12px;
      border-radius: 8px;
      box-shadow: 0 1px 8px rgba(0,0,0,0.2);
      font-family: Arial, sans-serif;
      font-size: 13px;
      line-height: 1.4;
      max-width: 360px;
    }}

    .date-toggle {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 8px;
    }}

    .date-toggle button {{
      border: 0;
      border-radius: 999px;
      padding: 6px 10px;
      background: #cbd5e1;
      color: #0f172a;
      font-size: 12px;
      cursor: pointer;
    }}

    .date-toggle button.active {{
      background: #0f172a;
      color: #f8fafc;
    }}

    .legend-box {{
      background: white;
      padding: 10px 12px;
      border-radius: 8px;
      box-shadow: 0 1px 8px rgba(0,0,0,0.2);
      font-family: Arial, sans-serif;
      font-size: 13px;
      line-height: 1.4;
    }}

    .legend-row {{
      display: flex;
      align-items: center;
      gap: 8px;
      margin: 4px 0;
    }}

    .swatch-line {{
      width: 34px;
      height: 0;
      border-top: 4px solid #000;
      display: inline-block;
    }}
  </style>
</head>
<body>
  <div id="map"></div>
  <div class="date-panel">
    <div style="font-weight:700;">{page_title}</div>
    <div>Select the traffic dataset below. Route layers stay as Leaflet checkboxes on the map.</div>
    <div class="date-toggle" id="dateToggle"></div>
  </div>

  <script
    src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"
    integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo="
    crossorigin=""
  ></script>
  <script>
    const map = L.map("map").setView([41.88, -87.63], 10);
    // Embed every precomputed route so the page can switch datasets instantly
    // without fetching separate GeoJSON files.
    const routeData = {route_data};
    const datasetOrder = {dataset_order};
    const startPoint = {start_point};
    const destinationPoint = {destination_point};

    L.tileLayer("https://{{s}}.tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors"
    }}).addTo(map);

    const startMarker = L.marker([startPoint.lat, startPoint.lon]).bindPopup(`<b>Start</b><br>${{startPoint.label}}`);
    const destinationMarker = L.marker([destinationPoint.lat, destinationPoint.lon]).bindPopup(`<b>Destination</b><br>${{destinationPoint.label}}`);
    startMarker.addTo(map);
    destinationMarker.addTo(map);

    const allBounds = L.latLngBounds([startPoint.lat, startPoint.lon], [destinationPoint.lat, destinationPoint.lon]);
    const routeDefinitions = {route_definitions};
    const activeRouteLayers = {{}};
    // Index features by dataset slug and route name so a date switch only swaps
    // the visible route geometry, not the whole map.
    const datasetLayers = new Map();
    const dateToggle = document.getElementById("dateToggle");
    const weekdayLabels = ["SUN", "MON", "TUE", "WED", "THU", "FRI", "SAT"];

    function formatDatasetLabel(dateLabel) {{
      const parsedDate = new Date(`${{dateLabel}}T00:00:00`);
      if (Number.isNaN(parsedDate.getTime())) {{
        return dateLabel;
      }}
      return `${{dateLabel}}-${{weekdayLabels[parsedDate.getDay()]}}`;
    }}

    routeData.features.forEach((feature) => {{
      const properties = feature.properties || {{}};
      if (!datasetLayers.has(properties.dataset_slug)) {{
        datasetLayers.set(properties.dataset_slug, {{}});
      }}
      datasetLayers.get(properties.dataset_slug)[properties.route_name] = feature;
    }});

    routeDefinitions.forEach((routeDefinition) => {{
      const layer = L.geoJSON(null, {{
        style: (feature) => {{
        const p = feature.properties || {{}};
        return {{
          color: p.color || "#dc2626",
          weight: 6,
          opacity: 0.92
        }};
        }},
        onEachFeature: (feature, routeLayer) => {{
          const p = feature.properties || {{}};
          routeLayer.bindPopup(
            `<b>${{p.destination_label}}</b><br>` +
            `<b>Route:</b> ${{p.route_label}}<br>` +
            `<b>Dataset:</b> ${{formatDatasetLabel(p.dataset_label)}}<br>` +
            `<b>Path nodes:</b> ${{p.path_nodes}}<br>` +
            `<b>Total cost:</b> ${{Number(p.total_cost).toFixed(4)}}<br>` +
            `<b>Weights:</b> D=${{p.distance_weight}}, P=${{p.population_weight}}, A=${{p.airspace_weight}}, T=${{p.traffic_weight}}`
          );
        }}
      }});
      activeRouteLayers[routeDefinition.name] = layer;
    }});

    const layerControlEntries = Object.fromEntries(
      routeDefinitions.map((routeDefinition) => [routeDefinition.label, activeRouteLayers[routeDefinition.name]])
    );

    const layerControl = L.control.layers(null, layerControlEntries, {{
      collapsed: false
    }}).addTo(map);

    function setActiveDateButton(activeSlug) {{
      Array.from(dateToggle.querySelectorAll("button")).forEach((button) => {{
        button.classList.toggle("active", button.dataset.slug === activeSlug);
      }});
    }}

    function showDataset(datasetSlug) {{
      const routeLayerMap = datasetLayers.get(datasetSlug) || {{}};
      Object.entries(activeRouteLayers).forEach(([routeName, layer]) => {{
        const feature = routeLayerMap[routeName];
        const wasVisible = map.hasLayer(layer);
        layer.clearLayers();
        if (feature) {{
          layer.addData(feature);
          allBounds.extend(layer.getBounds());
        }}
        if (wasVisible && feature) {{
          layer.addTo(map);
        }}
      }});

      if (!map.hasLayer(activeRouteLayers["combined"])) {{
        activeRouteLayers["combined"].addTo(map);
      }}
      setActiveDateButton(datasetSlug);
    }}

    datasetOrder.forEach((dataset) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.slug = dataset.slug;
      button.textContent = formatDatasetLabel(dataset.label);
      button.addEventListener("click", () => showDataset(dataset.slug));
      dateToggle.appendChild(button);
    }});

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "legend-box");
      const legendRows = routeDefinitions.map((routeDefinition) => `
        <div class="legend-row">
          <span class="swatch-line" style="border-top-color:${{routeDefinition.color}};"></span>
          <span>${{routeDefinition.label}}</span>
        </div>
      `).join("");
      div.innerHTML = `
        <div style="font-weight:700; margin-bottom:6px;">{page_title}</div>
        ${{legendRows}}
      `;
      return div;
    }};
    legend.addTo(map);

    map.fitBounds(allBounds, {{ padding: [20, 20] }});
    showDataset(datasetOrder[0].slug);
  </script>
</body>
</html>
"""


def normalize(values: np.ndarray) -> tuple[np.ndarray, float]:
    maximum = float(np.max(values))
    if maximum <= 0:
        return np.zeros_like(values, dtype=float), 1.0
    return values / maximum, maximum


def path_cost(graph: nx.Graph, path: list[int], weight: str) -> float:
    return sum(graph[path[i]][path[i + 1]][weight] for i in range(len(path) - 1))


def detect_csv_encoding(csv_path: Path) -> str:
    with csv_path.open("rb") as file_handle:
        first_bytes = file_handle.read(4)

    if first_bytes.startswith(b"\xff\xfe") or first_bytes.startswith(b"\xfe\xff"):
        return "utf-16"
    if first_bytes.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8"


def build_graph(grid: gpd.GeoDataFrame) -> tuple[nx.Graph, gpd.GeoSeries, np.ndarray, np.ndarray]:
    graph = nx.Graph()
    for index, row in grid.iterrows():
        graph.add_node(
            index,
            risk_cost=float(row["risk_cost"]),
            city_risk=float(row["city_risk"]),
            airport_risk_combined=float(row["airport_risk_combined"]),
            geometry=row.geometry,
        )

    sindex = grid.sindex
    geoms = grid.geometry
    for i, geom in enumerate(geoms):
        for j in sindex.intersection(geom.bounds):
            if j <= i:
                continue
            # `touches` keeps both edge-sharing and corner-sharing neighbors, so
            # diagonal steps are legal moves in the routing graph.
            if geom.touches(geoms.iloc[j]):
                graph.add_edge(i, j)

    grid_m = grid.to_crs("EPSG:3857")
    centroids_m = grid_m.geometry.centroid
    cent_x = centroids_m.x.to_numpy()
    cent_y = centroids_m.y.to_numpy()
    return graph, centroids_m, cent_x, cent_y


def node_for_point(
    lat: float,
    lon: float,
    sindex,
    geoms: gpd.GeoSeries,
    centroids_m: gpd.GeoSeries,
) -> int:
    point = Point(lon, lat)
    for index in sindex.intersection(point.bounds):
        if geoms.iloc[index].contains(point):
            return int(index)

    # Origins and destinations can land on a cell boundary, so fall back to the
    # nearest centroid instead of failing the route build.
    point_m = gpd.GeoSeries([point], crs="EPSG:4326").to_crs("EPSG:3857").iloc[0]
    return int(centroids_m.distance(point_m).idxmin())


def make_heuristic(distance_weight: float, cent_x: np.ndarray, cent_y: np.ndarray, max_edge_distance: float):
    if distance_weight <= 0:
        return lambda _u, _v: 0.0

    def heuristic(u: int, v: int) -> float:
        dx = cent_x[u] - cent_x[v]
        dy = cent_y[u] - cent_y[v]
        return distance_weight * (math.hypot(dx, dy) / max_edge_distance)

    return heuristic


def route_line(grid: gpd.GeoDataFrame, path: list[int], simplify_tolerance_m: float = 0.0) -> LineString:
    line = LineString([grid.loc[node].geometry.centroid for node in path])
    if simplify_tolerance_m <= 0:
        return line

    simplified = (
        gpd.GeoSeries([line], crs=grid.crs)
        .to_crs("EPSG:3857")
        .simplify(simplify_tolerance_m, preserve_topology=False)
        .to_crs(grid.crs)
        .iloc[0]
    )
    if isinstance(simplified, LineString):
        return simplified
    return line


def direct_route_line(start: dict[str, object], destination: dict[str, object]) -> LineString:
    return LineString(
        [
            Point(float(start["lon"]), float(start["lat"])),
            Point(float(destination["lon"]), float(destination["lat"])),
        ]
    )


def load_traffic_counts(csv_path: Path, grid: gpd.GeoDataFrame) -> np.ndarray:
    data = pd.read_csv(csv_path, encoding=detect_csv_encoding(csv_path))
    data = data.dropna(subset=["lat", "lon"]).copy()
    data["lat"] = pd.to_numeric(data["lat"], errors="coerce")
    data["lon"] = pd.to_numeric(data["lon"], errors="coerce")
    data = data.dropna(subset=["lat", "lon"])

    # Convert raw OpenSky observations into per-cell traffic counts for the
    # dataset date represented by this route page.
    observations = gpd.GeoDataFrame(
        data[["lat", "lon"]].copy(),
        geometry=gpd.points_from_xy(data["lon"], data["lat"]),
        crs="EPSG:4326",
    )
    joined = gpd.sjoin(
        observations,
        grid[["geometry"]],
        how="left",
        predicate="within",
    )
    counts = joined["index_right"].value_counts().sort_index()
    traffic_counts = np.zeros(len(grid), dtype=float)
    valid_counts = counts[counts.index.notna()]
    if not valid_counts.empty:
        traffic_counts[valid_counts.index.astype(int)] = valid_counts.to_numpy(dtype=float)
    return traffic_counts


def assign_route_edge_weights(
    graph: nx.Graph,
    cent_x: np.ndarray,
    cent_y: np.ndarray,
    city_risk_norm: np.ndarray,
    airspace_risk_norm: np.ndarray,
    traffic_risk_norm: np.ndarray,
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
        # Treat the edge as inheriting the average risk of the two cells it
        # connects, which keeps node-based cost layers compatible with A*.
        population_norm = 0.5 * (city_risk_norm[u] + city_risk_norm[v])
        airspace_norm = 0.5 * (airspace_risk_norm[u] + airspace_risk_norm[v])
        traffic_norm = 0.5 * (traffic_risk_norm[u] + traffic_risk_norm[v])
        for spec in ROUTE_SPECS:
            graph[u][v][f"weight_{spec['name']}"] = (
                spec["distance_weight"] * distance_norm
                + spec["population_weight"] * population_norm
                + spec["airspace_weight"] * airspace_norm
                + spec["traffic_weight"] * traffic_norm
            )

    return max_edge_distance


def build_route_features() -> dict[str, gpd.GeoDataFrame]:
    grid = gpd.read_file(GRID_PATH).reset_index(drop=True)
    graph, centroids_m, cent_x, cent_y = build_graph(grid)
    sindex = grid.sindex
    geoms = grid.geometry

    city_risk = grid["city_risk"].to_numpy(dtype=float)
    airspace_risk = grid["airport_risk_combined"].to_numpy(dtype=float)
    city_risk_norm, _ = normalize(city_risk)
    airspace_risk_norm, _ = normalize(airspace_risk)

    start_node = node_for_point(START["lat"], START["lon"], sindex, geoms, centroids_m)
    destination_nodes = {
        destination["slug"]: node_for_point(destination["lat"], destination["lon"], sindex, geoms, centroids_m)
        for destination in DESTINATIONS
    }

    route_rows_by_destination: dict[str, list[dict[str, object]]] = {destination["slug"]: [] for destination in DESTINATIONS}
    route_geometries_by_destination: dict[str, list[LineString]] = {destination["slug"]: [] for destination in DESTINATIONS}

    # Precompute every destination/date/preset combination once so the HTML
    # pages only need to toggle prebuilt route features.
    for dataset in DATASETS:
        if not dataset["csv_path"].exists():
            raise SystemExit(f"Traffic CSV not found: {dataset['csv_path']}")

        traffic_counts = load_traffic_counts(dataset["csv_path"], grid)
        traffic_risk_norm, traffic_max = normalize(traffic_counts)
        max_edge_distance = assign_route_edge_weights(
            graph,
            cent_x,
            cent_y,
            city_risk_norm,
            airspace_risk_norm,
            traffic_risk_norm,
        )
        for destination in DESTINATIONS:
            end_node = destination_nodes[destination["slug"]]
            for spec in ROUTE_SPECS:
                weight_key = f"weight_{spec['name']}"
                heuristic = make_heuristic(spec["distance_weight"], cent_x, cent_y, max_edge_distance)
                path = nx.astar_path(graph, start_node, end_node, heuristic=heuristic, weight=weight_key)
                total_cost = path_cost(graph, path, weight_key)
                route_rows_by_destination[destination["slug"]].append(
                    {
                        "dataset_slug": dataset["slug"],
                        "dataset_label": dataset["label"],
                        "destination_slug": destination["slug"],
                        "destination_label": destination["label"],
                        "route_name": spec["name"],
                        "route_label": spec["label"],
                        "path_nodes": len(path),
                        "total_cost": total_cost,
                        "distance_weight": spec["distance_weight"],
                        "population_weight": spec["population_weight"],
                        "airspace_weight": spec["airspace_weight"],
                        "traffic_weight": spec["traffic_weight"],
                        "traffic_samples_max": float(traffic_max),
                        "algorithm": "astar",
                        "color": spec["color"],
                    }
                )
                geometry = route_line(grid, path)
                route_geometries_by_destination[destination["slug"]].append(geometry)
                print(
                    f"{destination['slug']} | {dataset['slug']} | {spec['name']}:",
                    f"nodes={len(path)}",
                    f"cost={total_cost:.4f}",
                    f"traffic_max={traffic_max:.0f}",
                )

    return {
        destination["slug"]: gpd.GeoDataFrame(
            route_rows_by_destination[destination["slug"]],
            geometry=route_geometries_by_destination[destination["slug"]],
            crs=grid.crs,
        )
        for destination in DESTINATIONS
    }


def write_html(destination: dict[str, object], route_geojson: dict[str, object]) -> None:
    page_title = f"Clow To {destination['label']} A* Route"
    output_path = HTML_OUTPUT / f"clow_to_{destination['slug']}_astar.html"
    # Emit a single self-contained HTML artifact per destination.
    html = HTML_TEMPLATE.format(
        page_title=page_title,
        route_data=json.dumps(route_geojson),
        dataset_order=json.dumps([{"slug": item["slug"], "label": item["label"]} for item in DATASETS]),
        route_definitions=json.dumps(
            [{"name": item["name"], "label": item["label"], "color": item["color"]} for item in ROUTE_SPECS]
        ),
        start_point=json.dumps(START),
        destination_point=json.dumps(destination),
        destination_label=destination["label"],
        distance_weight=DISTANCE_WEIGHT,
        population_weight=POPULATION_WEIGHT,
        airspace_weight=AIRSPACE_WEIGHT,
        traffic_weight=TRAFFIC_WEIGHT,
    )
    output_path.write_text(html, encoding="utf-8")
    print(f"Saved HTML: {output_path}")


def main() -> None:
    HTML_OUTPUT.mkdir(parents=True, exist_ok=True)
    GEOJSON_OUTPUT.mkdir(parents=True, exist_ok=True)
    routes_by_destination = build_route_features()

    for destination in DESTINATIONS:
        route_gdf = routes_by_destination[destination["slug"]]
        geojson_path = GEOJSON_OUTPUT / f"clow_to_{destination['slug']}_astar_routes.geojson"
        route_gdf.to_file(geojson_path, driver="GeoJSON")
        print(f"Saved GeoJSON: {geojson_path}")
        write_html(destination, json.loads(route_gdf.to_json()))


if __name__ == "__main__":
    main()
