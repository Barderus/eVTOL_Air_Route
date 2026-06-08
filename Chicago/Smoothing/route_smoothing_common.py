from __future__ import annotations

"""Shared helpers for route-smoothing experiments.

The smoothing scripts reuse the base A* route outputs, then either modify the
search itself or post-process the displayed line geometry. This module holds the
shared data loading, route evaluation, HTML generation, and cost re-scoring
logic for those experiments.
"""

import csv
import io
import json
import math
import sys
from contextlib import redirect_stdout
from dataclasses import dataclass
from heapq import heappop, heappush
from pathlib import Path
from typing import Any

import geopandas as gpd
import networkx as nx
import numpy as np
from shapely.geometry import LineString, Point, mapping


PROJECT_ROOT = Path(__file__).resolve().parent.parent
HTML_OUTPUT = PROJECT_ROOT / "html"

ROUTES_2D_DIR = PROJECT_ROOT / "2D"
if str(ROUTES_2D_DIR) not in sys.path:
    sys.path.insert(0, str(ROUTES_2D_DIR))

import generate_astar_toggle_pages as base_routes


TURN_PENALTY_WEIGHT = 1.5
DISPLAY_SAMPLE_SPACING_M = 75.0
SPLINE_SAMPLES_PER_SEGMENT = 6

CSV_FIELD_NAMES = [
    "destination_slug",
    "destination",
    "dataset_slug",
    "dataset_label",
    "total_cost",
    "path_nodes",
    "final_combined_cost",
    "final_path_nodes",
]


@dataclass
class DatasetContext:
    """All per-dataset objects needed to generate and score smoothed routes."""
    dataset_slug: str
    dataset_label: str
    traffic_scale: float
    graph: nx.Graph
    centroids_m: gpd.GeoSeries
    cent_x: np.ndarray
    cent_y: np.ndarray
    city_risk_norm: np.ndarray
    airspace_risk_norm: np.ndarray
    traffic_risk_norm: np.ndarray
    max_edge_distance: float
    start_node: int
    destination_nodes: dict[str, int]
    grid: gpd.GeoDataFrame
    geoms: gpd.GeoSeries
    sindex: Any


@dataclass
class RouteResult:
    """One rendered route plus the metadata needed for HTML and CSV outputs."""
    destination_slug: str
    destination_label: str
    dataset_slug: str
    dataset_label: str
    total_cost: float
    path_nodes: int
    final_combined_cost: float
    final_path_nodes: int
    display_line: LineString
    raw_path: list[int]
    metadata: dict[str, Any]


def destination_display_name(destination_slug: str, original_label: str) -> str:
    """Keep display labels consistent with the existing handoff artifacts."""
    if destination_slug == "union_station":
        return "Chicago"
    return original_label


def sort_rows(rows: list[dict[str, Any]]) -> None:
    """Sort CSV rows into a stable destination/date order."""
    destination_positions = {
        destination["slug"]: index
        for index, destination in enumerate(base_routes.DESTINATIONS)
    }
    dataset_positions = {
        dataset["slug"]: index
        for index, dataset in enumerate(base_routes.DATASETS)
    }
    rows.sort(
        key=lambda row: (
            destination_positions[row["destination_slug"]],
            dataset_positions[row["dataset_slug"]],
        )
    )


def write_cost_csv(output_path: Path, rows: list[dict[str, Any]]) -> None:
    """Write one summary CSV for a smoothing experiment."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=CSV_FIELD_NAMES)
        writer.writeheader()
        writer.writerows(rows)


def base_combined_routes_by_destination() -> dict[str, gpd.GeoDataFrame]:
    """Load the baseline combined A* routes without printing their build log."""
    with redirect_stdout(io.StringIO()):
        routes_by_destination = base_routes.build_route_features()

    return {
        destination_slug: route_gdf[route_gdf["route_name"] == "combined"].copy()
        for destination_slug, route_gdf in routes_by_destination.items()
    }


def build_dataset_contexts() -> dict[str, DatasetContext]:
    """Build reusable per-dataset routing contexts for all experiments."""
    grid = gpd.read_file(base_routes.GRID_PATH).reset_index(drop=True)
    graph_template, centroids_m, cent_x, cent_y = base_routes.build_graph(grid)
    sindex = grid.sindex
    geoms = grid.geometry

    city_risk = grid["city_risk"].to_numpy(dtype=float)
    airspace_risk = grid["airport_risk_combined"].to_numpy(dtype=float)
    city_risk_norm, _ = base_routes.normalize(city_risk)
    airspace_risk_norm, _ = base_routes.normalize(airspace_risk, clip_percentile=99.0, power=0.75)

    start_node = base_routes.node_for_point(
        base_routes.START["lat"],
        base_routes.START["lon"],
        sindex,
        geoms,
        centroids_m,
    )
    destination_nodes = {
        destination["slug"]: base_routes.node_for_point(
            destination["lat"],
            destination["lon"],
            sindex,
            geoms,
            centroids_m,
        )
        for destination in base_routes.DESTINATIONS
    }

    contexts: dict[str, DatasetContext] = {}
    for dataset in base_routes.DATASETS:
        graph = nx.Graph(graph_template)
        traffic_counts = base_routes.load_traffic_counts(dataset["csv_path"], grid)
        traffic_risk_norm, traffic_scale = base_routes.normalize(
            traffic_counts,
            clip_percentile=95.0,
            transform="log1p",
            power=0.75,
        )
        max_edge_distance = base_routes.assign_route_edge_weights(
            graph,
            cent_x,
            cent_y,
            city_risk_norm,
            airspace_risk_norm,
            traffic_risk_norm,
        )
        contexts[dataset["slug"]] = DatasetContext(
            dataset_slug=dataset["slug"],
            dataset_label=dataset["label"],
            traffic_scale=float(traffic_scale),
            graph=graph,
            centroids_m=centroids_m,
            cent_x=cent_x,
            cent_y=cent_y,
            city_risk_norm=city_risk_norm,
            airspace_risk_norm=airspace_risk_norm,
            traffic_risk_norm=traffic_risk_norm,
            max_edge_distance=max_edge_distance,
            start_node=start_node,
            destination_nodes=destination_nodes,
            grid=grid,
            geoms=geoms,
            sindex=sindex,
        )

    return contexts


def _segment_angle_degrees(
    prev_node: int | None,
    current_node: int,
    next_node: int,
    cent_x: np.ndarray,
    cent_y: np.ndarray,
) -> float:
    """Measure the heading change between two consecutive route segments."""
    if prev_node is None:
        return 0.0

    incoming = np.array(
        [
            cent_x[current_node] - cent_x[prev_node],
            cent_y[current_node] - cent_y[prev_node],
        ],
        dtype=float,
    )
    outgoing = np.array(
        [
            cent_x[next_node] - cent_x[current_node],
            cent_y[next_node] - cent_y[current_node],
        ],
        dtype=float,
    )
    incoming_norm = np.linalg.norm(incoming)
    outgoing_norm = np.linalg.norm(outgoing)
    if incoming_norm <= 0 or outgoing_norm <= 0:
        return 0.0

    cosine = float(np.dot(incoming, outgoing) / (incoming_norm * outgoing_norm))
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def turn_penalty_cost(
    prev_node: int | None,
    current_node: int,
    next_node: int,
    cent_x: np.ndarray,
    cent_y: np.ndarray,
    weight: float = TURN_PENALTY_WEIGHT,
) -> float:
    """Convert a heading change into the quadratic penalty used by approach 1."""
    angle_degrees = _segment_angle_degrees(prev_node, current_node, next_node, cent_x, cent_y)
    if angle_degrees <= 0:
        return 0.0
    return weight * ((angle_degrees / 180.0) ** 2)


def astar_with_turn_penalty(
    context: DatasetContext,
    destination_slug: str,
    weight_key: str = "weight_combined",
) -> tuple[list[int], float]:
    """Run A* over an expanded state space that remembers the previous node."""
    goal_node = context.destination_nodes[destination_slug]
    heuristic = base_routes.make_heuristic(
        base_routes.DISTANCE_WEIGHT,
        context.cent_x,
        context.cent_y,
        context.max_edge_distance,
    )

    start_state = (context.start_node, None)
    open_heap: list[tuple[float, float, tuple[int, int | None]]] = [
        (heuristic(context.start_node, goal_node), 0.0, start_state)
    ]
    g_score: dict[tuple[int, int | None], float] = {start_state: 0.0}
    came_from: dict[tuple[int, int | None], tuple[int, int | None]] = {}

    best_goal_state: tuple[int, int | None] | None = None
    best_goal_cost = math.inf

    while open_heap:
        _, current_cost, state = heappop(open_heap)
        current_node, prev_node = state

        if current_cost > g_score.get(state, math.inf):
            continue

        if current_node == goal_node:
            best_goal_state = state
            best_goal_cost = current_cost
            break

        for neighbor in context.graph.neighbors(current_node):
            edge_cost = float(context.graph[current_node][neighbor][weight_key])
            penalty = turn_penalty_cost(
                prev_node,
                current_node,
                neighbor,
                context.cent_x,
                context.cent_y,
            )
            tentative_g = current_cost + edge_cost + penalty
            neighbor_state = (neighbor, current_node)
            if tentative_g >= g_score.get(neighbor_state, math.inf):
                continue

            came_from[neighbor_state] = state
            g_score[neighbor_state] = tentative_g
            estimate = tentative_g + heuristic(neighbor, goal_node)
            heappush(open_heap, (estimate, tentative_g, neighbor_state))

    if best_goal_state is None:
        raise RuntimeError(f"Failed to find turn-penalty route for {destination_slug}")

    ordered_path = [best_goal_state[0]]
    state = best_goal_state
    while state in came_from:
        state = came_from[state]
        ordered_path.append(state[0])
    ordered_path.reverse()
    return ordered_path, best_goal_cost


def catmull_rom_spline(line: LineString, samples_per_segment: int = SPLINE_SAMPLES_PER_SEGMENT) -> LineString:
    """Smooth a polyline with a Catmull-Rom spline for display purposes."""
    coordinates = list(line.coords)
    if len(coordinates) < 4:
        return line

    points = [np.asarray(coord, dtype=float) for coord in coordinates]
    padded = [points[0], *points, points[-1]]
    smoothed: list[tuple[float, float]] = [tuple(points[0])]

    for index in range(1, len(padded) - 2):
        p0, p1, p2, p3 = padded[index - 1 : index + 3]
        for step in range(1, samples_per_segment + 1):
            t = step / samples_per_segment
            t2 = t * t
            t3 = t2 * t
            point = 0.5 * (
                (2.0 * p1)
                + (-p0 + p2) * t
                + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
                + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3
            )
            smoothed.append((float(point[0]), float(point[1])))

    if smoothed[-1] != tuple(points[-1]):
        smoothed.append(tuple(points[-1]))

    return LineString(smoothed)


def sample_line_to_cell_path(
    line: LineString,
    context: DatasetContext,
    sample_spacing_m: float = DISPLAY_SAMPLE_SPACING_M,
) -> list[int]:
    """Map a displayed line back onto a neighboring sequence of grid cells."""
    if line.is_empty:
        return []

    line_m = gpd.GeoSeries([line], crs=context.grid.crs).to_crs("EPSG:3857").iloc[0]
    distances = list(np.arange(0.0, line_m.length, sample_spacing_m))
    if not distances or distances[-1] != line_m.length:
        distances.append(line_m.length)

    sample_points_m = [line_m.interpolate(distance) for distance in distances]
    sample_points = gpd.GeoSeries(sample_points_m, crs="EPSG:3857").to_crs(context.grid.crs)

    cell_path: list[int] = []
    for point in sample_points:
        node = base_routes.node_for_point(point.y, point.x, context.sindex, context.geoms, context.centroids_m)
        if not cell_path or cell_path[-1] != node:
            cell_path.append(node)

    if not cell_path:
        return []

    expanded_path = [cell_path[0]]
    for node in cell_path[1:]:
        last_node = expanded_path[-1]
        if node == last_node:
            continue
        if context.graph.has_edge(last_node, node):
            expanded_path.append(node)
            continue

        bridge = nx.shortest_path(context.graph, last_node, node, weight="weight_combined")
        expanded_path.extend(bridge[1:])

    return expanded_path


def evaluate_combined_cost_for_line(line: LineString, context: DatasetContext) -> tuple[float, int, list[int]]:
    """Re-score a rendered line against the base combined route weights."""
    cell_path = sample_line_to_cell_path(line, context)
    if len(cell_path) < 2:
        return 0.0, len(cell_path), cell_path

    total_cost = base_routes.path_cost(context.graph, cell_path, "weight_combined")
    return float(total_cost), len(cell_path), cell_path


def route_result_to_row(result: RouteResult) -> dict[str, Any]:
    """Flatten a route result into one CSV row."""
    return {
        "destination_slug": result.destination_slug,
        "destination": destination_display_name(result.destination_slug, result.destination_label),
        "dataset_slug": result.dataset_slug,
        "dataset_label": result.dataset_label,
        "total_cost": float(result.total_cost),
        "path_nodes": int(result.path_nodes),
        "final_combined_cost": float(result.final_combined_cost),
        "final_path_nodes": int(result.final_path_nodes),
    }


def route_results_to_feature_collection(results: list[RouteResult]) -> dict[str, Any]:
    """Convert experiment results into GeoJSON-like data for HTML embedding."""
    features = []
    for index, result in enumerate(results):
        properties = {
            "dataset_slug": result.dataset_slug,
            "dataset_label": result.dataset_label,
            "destination_slug": result.destination_slug,
            "destination_label": result.destination_label,
            "display_destination": destination_display_name(result.destination_slug, result.destination_label),
            "total_cost": result.total_cost,
            "path_nodes": result.path_nodes,
            "final_combined_cost": result.final_combined_cost,
            "final_path_nodes": result.final_path_nodes,
            **result.metadata,
        }
        features.append(
            {
                "id": str(index),
                "type": "Feature",
                "properties": properties,
                "geometry": mapping(result.display_line),
            }
        )
    return {"type": "FeatureCollection", "features": features}


def write_html(output_filename: str, page_title: str, technique_name: str, results: list[RouteResult]) -> Path:
    """Write one self-contained HTML page for a smoothing experiment."""
    route_data = route_results_to_feature_collection(results)
    output_path = HTML_OUTPUT / output_filename
    output_path.parent.mkdir(parents=True, exist_ok=True)

    html = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{page_title}</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin=""
  />
  <style>
    :root {{
      --surface: #ffffff;
      --surface-soft: #f5f7f8;
      --line: #cfd7dd;
      --text: #172026;
      --muted: #5f6f7a;
      --active: #165a72;
      --route: #dc2626;
    }}

    * {{ box-sizing: border-box; }}

    html,
    body {{
      height: 100%;
      margin: 0;
      color: var(--text);
      font-family: Arial, Helvetica, sans-serif;
      overflow: hidden;
    }}

    body {{
      display: flex;
      flex-direction: column;
      background: var(--surface-soft);
    }}

    .toolbar {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      padding: 10px 14px;
      border-bottom: 1px solid var(--line);
      background: var(--surface);
      z-index: 1000;
    }}

    .title-block {{ min-width: 240px; }}

    h1 {{
      margin: 0;
      font-size: 17px;
      font-weight: 700;
      line-height: 1.2;
    }}

    .control-strip {{
      display: flex;
      align-items: center;
      justify-content: flex-end;
      gap: 10px;
      flex-wrap: wrap;
    }}

    .segmented {{
      display: inline-grid;
      min-height: 38px;
      border: 1px solid var(--line);
      border-radius: 7px;
      overflow: hidden;
      background: var(--surface);
    }}

    .destination-toggle {{ grid-template-columns: repeat(3, minmax(128px, 1fr)); }}
    .date-toggle {{ grid-template-columns: repeat(6, minmax(86px, 1fr)); }}

    .segmented button {{
      min-width: 0;
      padding: 8px 12px;
      border: 0;
      border-right: 1px solid var(--line);
      color: var(--text);
      background: transparent;
      font: inherit;
      font-size: 13px;
      cursor: pointer;
    }}

    .segmented button:last-child {{ border-right: 0; }}
    .segmented button:hover {{ background: #eef3f5; }}
    .segmented button.active {{ color: #ffffff; background: var(--active); }}

    #map {{
      flex: 1 1 auto;
      width: 100%;
      height: calc(100vh - 59px);
      min-height: 360px;
    }}

    .legend-box {{
      background: white;
      padding: 10px 12px;
      border-radius: 7px;
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.14);
      font-size: 13px;
      line-height: 1.45;
    }}

    .legend-title {{
      font-weight: 700;
      margin-bottom: 6px;
    }}

    .legend-row {{
      display: flex;
      align-items: center;
      gap: 8px;
      margin: 5px 0;
    }}

    .swatch-line {{
      width: 34px;
      height: 0;
      border-top: 4px solid var(--route);
      display: inline-block;
    }}

    @media (max-width: 980px) {{
      .toolbar {{ align-items: stretch; flex-direction: column; gap: 9px; }}
      .control-strip {{ justify-content: stretch; }}
      .segmented {{ width: 100%; }}
      .destination-toggle {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
      .date-toggle {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
      .segmented button {{ padding-inline: 6px; font-size: 12px; }}
      #map {{ height: calc(100vh - 151px); }}
    }}
  </style>
</head>
<body>
  <header class="toolbar">
    <div class="title-block">
      <h1>{page_title}</h1>
    </div>
    <div class="control-strip">
      <div class="segmented destination-toggle" id="destinationToggle" role="group" aria-label="Destination"></div>
      <div class="segmented date-toggle" id="dateToggle" role="group" aria-label="Traffic dataset"></div>
    </div>
  </header>
  <main id="map" aria-label="{page_title} map"></main>

  <script
    src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""
  ></script>
  <script>
    const STUDY_BOUNDS = L.latLngBounds(
      [41.48, -88.34],
      [42.21, -87.52]
    );
    const routeData = {route_data};
    const datasetOrder = {dataset_order};
    const destinations = {destinations};
    const startPoint = {start_point};
    const destinationToggle = document.getElementById("destinationToggle");
    const dateToggle = document.getElementById("dateToggle");

    let activeDestinationSlug = destinations[0].slug;
    let activeDatasetSlug = datasetOrder[0].slug;
    let activeLayer = null;

    const routeByKey = new Map();
    routeData.features.forEach((feature) => {{
      const p = feature.properties || {{}};
      routeByKey.set(`${{p.destination_slug}}|${{p.dataset_slug}}`, feature);
    }});

    const map = L.map("map", {{
      center: [41.84, -87.93],
      zoom: 10
    }});

    L.tileLayer("https://{{s}}.tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors"
    }}).addTo(map);

    const startMarker = L.marker([startPoint.lat, startPoint.lon]).bindPopup(`<b>Start</b><br>${{startPoint.label}}`);
    startMarker.addTo(map);
    let destinationMarker = null;

    function formatDatasetLabel(dateLabel) {{
      if (/^(Sun|Mon|Tue|Wed|Thu|Fri|Sat)\\s/.test(dateLabel)) {{
        return dateLabel;
      }}
      const parsedDate = new Date(`${{dateLabel}}T00:00:00`);
      if (Number.isNaN(parsedDate.getTime())) {{
        return dateLabel;
      }}
      return parsedDate.toLocaleDateString("en-US", {{
        weekday: "short",
        year: "numeric",
        month: "2-digit",
        day: "2-digit"
      }}).replace(/,\\s*/, " ");
    }}

    function setActiveButtons(toggleEl, field, value) {{
      Array.from(toggleEl.querySelectorAll("button")).forEach((button) => {{
        const isActive = button.dataset[field] === value;
        button.classList.toggle("active", isActive);
        button.setAttribute("aria-pressed", isActive ? "true" : "false");
      }});
    }}

    function syncRoute() {{
      const key = `${{activeDestinationSlug}}|${{activeDatasetSlug}}`;
      const feature = routeByKey.get(key);
      if (!feature) {{
        if (activeLayer) {{
          map.removeLayer(activeLayer);
          activeLayer = null;
        }}
        if (destinationMarker) {{
          map.removeLayer(destinationMarker);
          destinationMarker = null;
        }}
        map.fitBounds(STUDY_BOUNDS, {{ padding: [20, 20] }});
        setActiveButtons(destinationToggle, "destination", activeDestinationSlug);
        setActiveButtons(dateToggle, "slug", activeDatasetSlug);
        return;
      }}

      if (activeLayer) {{
        map.removeLayer(activeLayer);
      }}
      if (destinationMarker) {{
        map.removeLayer(destinationMarker);
      }}

      activeLayer = L.geoJSON(feature, {{
        style: () => ({{
          color: "#dc2626",
          weight: 6,
          opacity: 0.94
        }})
      }}).addTo(map);

      const p = feature.properties || {{}};
      destinationMarker = L.marker([
        p.destination_lat,
        p.destination_lon
      ]).bindPopup(`<b>Destination</b><br>${{p.display_destination}}`).addTo(map);

      const bounds = activeLayer.getBounds();
      if (bounds.isValid()) {{
        map.fitBounds(bounds.pad(0.08), {{ animate: false }});
      }} else {{
        map.fitBounds(STUDY_BOUNDS, {{ padding: [20, 20] }});
      }}

      setActiveButtons(destinationToggle, "destination", activeDestinationSlug);
      setActiveButtons(dateToggle, "slug", activeDatasetSlug);
    }}

    destinations.forEach((destination) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.destination = destination.slug;
      button.textContent = destination.label;
      button.addEventListener("click", () => {{
        activeDestinationSlug = destination.slug;
        syncRoute();
      }});
      destinationToggle.appendChild(button);
    }});

    datasetOrder.forEach((dataset) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.slug = dataset.slug;
      button.textContent = formatDatasetLabel(dataset.label);
      button.addEventListener("click", () => {{
        activeDatasetSlug = dataset.slug;
        syncRoute();
      }});
      dateToggle.appendChild(button);
    }});

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "legend-box");
      div.innerHTML = `
        <div class="legend-title">{page_title}</div>
        <div class="legend-row">
          <span class="swatch-line"></span>
          <span>Displayed route geometry</span>
        </div>
      `;
      return div;
    }};
    legend.addTo(map);

    syncRoute();
  </script>
</body>
</html>
""".format(
        page_title=page_title,
        technique_name=technique_name,
        route_data=json.dumps(route_data),
        dataset_order=json.dumps(
            [{"slug": item["slug"], "label": item["label"]} for item in base_routes.DATASETS]
        ),
        destinations=json.dumps(
            [
                {
                    "slug": item["slug"],
                    "label": destination_display_name(item["slug"], item["label"]),
                }
                for item in base_routes.DESTINATIONS
            ]
        ),
        start_point=json.dumps(base_routes.START),
    )

    output_path.write_text(html, encoding="utf-8")
    return output_path
