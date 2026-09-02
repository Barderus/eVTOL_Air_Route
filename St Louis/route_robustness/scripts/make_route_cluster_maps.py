"""Build Leaflet maps for St. Louis route clusters by route pair."""

import colorsys
import json
from pathlib import Path

import pandas as pd


ST_LOUIS_ROBUSTNESS = Path("St Louis") / "route_robustness"
OUTPUT_FOLDER = ST_LOUIS_ROBUSTNESS / "output"
MAPS_FOLDER = ST_LOUIS_ROBUSTNESS / "route_clusters"
ROUTES_GEOJSON = OUTPUT_FOLDER / "st_louis_weighted_routes.geojson"

METHODS = [
    {
        "key": "dbscan",
        "label": "DBSCAN",
        "file": OUTPUT_FOLDER / "st_louis_route_clusters_dbscan.csv",
        "cluster_column": "dbscan_cluster",
    },
    {
        "key": "edit_distance",
        "label": "Edit Distance",
        "file": OUTPUT_FOLDER / "st_louis_route_clusters_edit_distance.csv",
        "cluster_column": "edit_distance_cluster",
    },
    {
        "key": "frechet",
        "label": "Frechet",
        "file": OUTPUT_FOLDER / "st_louis_route_clusters_frechet.csv",
        "cluster_column": "frechet_cluster",
    },
    {
        "key": "hierarchical_jaccard",
        "label": "Hierarchical Jaccard",
        "file": OUTPUT_FOLDER / "st_louis_route_clusters_hierarchical_jaccard.csv",
        "cluster_column": "hierarchical_jaccard_cluster",
    },
]

ROUTE_PAIR_ORDER = [
    "midamerica_to_st_louis_lambert",
    "midamerica_to_st_louis_union_station",
    "st_louis_downtown_airport_to_st_louis_lambert",
]

BASE_COLORS = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
    "#1f78b4",
    "#b15928",
    "#6a3d9a",
    "#fb9a99",
    "#33a02c",
]


def color_for_index(index, total_count):
    """Return a deterministic route-cluster color."""
    if total_count <= len(BASE_COLORS):
        return BASE_COLORS[index]

    hue = (index * 137.508) % 360.0
    lightness = 38.0 + (index % 4) * 8.0
    red, green, blue = colorsys.hls_to_rgb(hue / 360.0, lightness / 100.0, 0.70)
    return f"#{round(red * 255):02x}{round(green * 255):02x}{round(blue * 255):02x}"


def require_columns(dataframe, required_columns, table_name):
    """Validate that a table has the expected columns."""
    missing_columns = set(required_columns).difference(dataframe.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"{table_name} is missing columns: {missing_text}")


def load_route_features():
    """Load successful weighted route GeoJSON features by route_run_id."""
    with ROUTES_GEOJSON.open("r", encoding="utf-8") as file_handle:
        geojson = json.load(file_handle)

    features_by_id = {}
    for feature in geojson["features"]:
        properties = feature.get("properties", {})
        if properties.get("status") != "success":
            continue
        features_by_id[properties["route_run_id"]] = feature
    return features_by_id


def load_cluster_tables():
    """Load all requested cluster assignments."""
    cluster_tables = {}
    for method in METHODS:
        clusters = pd.read_csv(method["file"])
        require_columns(
            clusters,
            ["route_run_id", "route_pair", method["cluster_column"]],
            str(method["file"]),
        )
        cluster_tables[method["key"]] = clusters
    return cluster_tables


def coordinates_to_latlon(coordinates):
    """Convert GeoJSON lon/lat line coordinates to Leaflet lat/lon points."""
    return [[lat, lon] for lon, lat in coordinates]


def simplify_line(points, max_points=90):
    """Downsample long lines for lighter embedded HTML."""
    if len(points) <= max_points:
        return points

    keep_indexes = {
        round(index * (len(points) - 1) / (max_points - 1))
        for index in range(max_points)
    }
    return [point for index, point in enumerate(points) if index in keep_indexes]


def build_method_payload(method, clusters, features_by_id, route_pair):
    """Build map-ready routes and cluster summaries for one method."""
    cluster_column = method["cluster_column"]
    route_rows = clusters[clusters["route_pair"] == route_pair].copy()
    cluster_counts = route_rows[cluster_column].value_counts()
    cluster_colors = {
        cluster_id: color_for_index(index, len(cluster_counts))
        for index, cluster_id in enumerate(cluster_counts.index)
    }

    routes = []
    bounds_points = []
    for row in route_rows.itertuples(index=False):
        feature = features_by_id[row.route_run_id]
        properties = feature["properties"]
        cluster_id = getattr(row, cluster_column)
        points = simplify_line(
            coordinates_to_latlon(feature["geometry"]["coordinates"])
        )
        bounds_points.extend(points)
        routes.append(
            {
                "route_run_id": row.route_run_id,
                "weight_id": row.weight_id,
                "cluster": cluster_id,
                "color": cluster_colors[cluster_id],
                "points": points,
                "distance_km": round(float(row.route_distance_km), 3),
                "score": round(float(row.total_weighted_score), 3),
                "path_nodes": int(row.path_nodes),
                "weights": {
                    "distance": float(row.distance_weight),
                    "population": float(row.population_weight),
                    "traffic": float(row.traffic_weight),
                    "airspace": float(row.airspace_weight),
                },
            }
        )

    cluster_summary = [
        {
            "cluster": cluster_id,
            "count": int(cluster_counts.loc[cluster_id]),
            "color": cluster_colors[cluster_id],
        }
        for cluster_id in cluster_counts.index
    ]

    return {
        "label": method["label"],
        "routes": routes,
        "clusters": cluster_summary,
        "bounds": bounds_points,
    }


def build_route_pair_payload(route_pair, cluster_tables, features_by_id):
    """Build all method payloads for one route pair."""
    first_table = next(iter(cluster_tables.values()))
    first_row = first_table[first_table["route_pair"] == route_pair].iloc[0]

    payload = {
        "routePair": route_pair,
        "routeLabel": first_row["route_pair_label"],
        "originLabel": "",
        "destinationLabel": "",
        "origin": None,
        "destination": None,
        "methods": {},
    }

    sample_feature = features_by_id[first_row["route_run_id"]]
    sample_properties = sample_feature["properties"]
    payload["originLabel"] = sample_properties["origin_label"]
    payload["destinationLabel"] = sample_properties["destination_label"]
    payload["origin"] = [
        float(sample_properties["origin_lat"]),
        float(sample_properties["origin_lon"]),
    ]
    payload["destination"] = [
        float(sample_properties["destination_lat"]),
        float(sample_properties["destination_lon"]),
    ]

    for method in METHODS:
        payload["methods"][method["key"]] = build_method_payload(
            method,
            cluster_tables[method["key"]],
            features_by_id,
            route_pair,
        )

    return payload


def html_template(payload, title_prefix="St. Louis Route Clusters"):
    """Build one route-pair Leaflet cluster map."""
    payload_json = json.dumps(payload)
    method_order_json = json.dumps(
        [{"key": method["key"], "label": method["label"]} for method in METHODS]
    )
    title = f"{title_prefix} - {payload['routeLabel']}"
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin=""
  />
  <style>
    :root {{
      --surface: #ffffff;
      --line: #d1d5db;
      --text: #172026;
      --muted: #5f6f7a;
      --active: #165a72;
    }}
    * {{ box-sizing: border-box; }}
    html, body {{
      height: 100%;
      margin: 0;
      color: var(--text);
      font-family: Arial, Helvetica, sans-serif;
      overflow: hidden;
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      min-height: 58px;
      padding: 10px 14px;
      border-bottom: 1px solid var(--line);
      background: var(--surface);
      z-index: 1000;
    }}
    h1 {{
      margin: 0;
      font-size: 17px;
      line-height: 1.2;
    }}
    .status {{
      margin-top: 3px;
      color: var(--muted);
      font-size: 12px;
    }}
    .segmented {{
      display: inline-grid;
      grid-template-columns: repeat(4, minmax(120px, 1fr));
      min-height: 38px;
      border: 1px solid var(--line);
      border-radius: 7px;
      overflow: hidden;
      background: var(--surface);
    }}
    .segmented button {{
      padding: 8px 10px;
      border: 0;
      border-right: 1px solid var(--line);
      color: var(--text);
      background: transparent;
      font: inherit;
      font-size: 12px;
      cursor: pointer;
    }}
    .segmented button:last-child {{ border-right: 0; }}
    .segmented button:hover {{ background: #eef3f5; }}
    .segmented button.active {{
      color: #ffffff;
      background: var(--active);
    }}
    #map {{
      width: 100%;
      height: calc(100vh - 58px);
      min-height: 420px;
    }}
    .cluster-panel {{
      width: 258px;
      max-height: calc(100vh - 90px);
      overflow: auto;
      padding: 11px 12px;
      border-radius: 7px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.14);
      font-size: 12px;
      line-height: 1.35;
    }}
    .cluster-panel strong {{
      display: block;
      margin-bottom: 7px;
      font-size: 13px;
    }}
    .cluster-row {{
      appearance: none;
      width: 100%;
      border: 0;
      border-bottom: 1px solid rgba(209, 213, 219, 0.65);
      background: transparent;
      display: grid;
      grid-template-columns: 12px 1fr auto;
      gap: 7px;
      align-items: center;
      padding: 4px 0;
      color: inherit;
      cursor: pointer;
      font: inherit;
      text-align: left;
    }}
    .cluster-row:hover {{ background: rgba(229, 231, 235, 0.55); }}
    .cluster-row.is-hidden {{ opacity: 0.42; }}
    .cluster-swatch {{
      width: 10px;
      height: 10px;
      border-radius: 999px;
      border: 1px solid rgba(17, 24, 39, 0.25);
    }}
    .leaflet-popup-content {{
      min-width: 210px;
      font-size: 12px;
      line-height: 1.45;
    }}
    @media (max-width: 900px) {{
      .toolbar {{ align-items: stretch; flex-direction: column; gap: 9px; }}
      .segmented {{ width: 100%; grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      #map {{ height: calc(100vh - 118px); }}
    }}
  </style>
</head>
<body>
  <header class="toolbar">
    <div>
      <h1>{payload['routeLabel']}</h1>
      <div id="status" class="status">Select a clustering method.</div>
    </div>
    <div class="segmented" id="methodToggle" role="group" aria-label="Cluster method"></div>
  </header>

  <main id="map" aria-label="St. Louis route cluster map"></main>

  <script
    src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""
  ></script>
  <script>
    const routePayload = {payload_json};
    const methodOrder = {method_order_json};
    const map = L.map("map", {{
      center: routePayload.origin,
      zoom: 10,
      preferCanvas: true
    }});

    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);

    const methodToggle = document.getElementById("methodToggle");
    const statusEl = document.getElementById("status");
    const routeLayer = L.layerGroup().addTo(map);
    const endpointLayer = L.layerGroup().addTo(map);
    const clusterControl = L.control({{ position: "topright" }});
    let activeMethod = methodOrder[0].key;
    let hiddenClusters = new Set();
    let clusterLayers = new Map();
    let controlDiv = null;

    clusterControl.onAdd = function() {{
      controlDiv = L.DomUtil.create("div", "cluster-panel");
      L.DomEvent.disableClickPropagation(controlDiv);
      L.DomEvent.disableScrollPropagation(controlDiv);
      return controlDiv;
    }};
    clusterControl.addTo(map);

    function drawEndpoints() {{
      endpointLayer.clearLayers();
      [
        {{ label: "Origin", name: routePayload.originLabel, point: routePayload.origin }},
        {{ label: "Destination", name: routePayload.destinationLabel, point: routePayload.destination }}
      ].forEach((endpoint) => {{
        L.circleMarker(endpoint.point, {{
          radius: 7,
          color: "#172026",
          weight: 2,
          fillColor: endpoint.label === "Origin" ? "#ffffff" : "#facc15",
          fillOpacity: 1
        }}).bindPopup(`<b>${{endpoint.label}}</b><br>${{endpoint.name}}`).addTo(endpointLayer);
      }});
    }}

    function popupHtml(route) {{
      const weights = route.weights;
      return `
        <b>${{route.cluster}}</b><br>
        <b>Weight:</b> ${{route.weight_id}}<br>
        <b>D/P/T/A:</b>
        ${{weights.distance.toFixed(2)}} /
        ${{weights.population.toFixed(2)}} /
        ${{weights.traffic.toFixed(2)}} /
        ${{weights.airspace.toFixed(2)}}<br>
        <b>Distance:</b> ${{route.distance_km.toFixed(2)}} km<br>
        <b>Score:</b> ${{route.score.toFixed(3)}}<br>
        <b>Nodes:</b> ${{route.path_nodes}}
      `;
    }}

    function drawRoutes() {{
      routeLayer.clearLayers();
      clusterLayers = new Map();
      const methodData = routePayload.methods[activeMethod];
      const largestCluster = Math.max(...methodData.clusters.map((cluster) => cluster.count), 1);
      const clusterSizes = new Map(
        methodData.clusters.map((cluster) => [cluster.cluster, cluster.count])
      );

      methodData.routes.forEach((route) => {{
        const clusterSize = clusterSizes.get(route.cluster) || 1;
        const lineWeight = 2 + 5 * Math.sqrt(clusterSize / largestCluster);
        const line = L.polyline(route.points, {{
          color: route.color,
          weight: lineWeight,
          opacity: hiddenClusters.has(route.cluster) ? 0 : 0.42,
          interactive: true
        }}).bindPopup(popupHtml(route));

        if (!clusterLayers.has(route.cluster)) {{
          clusterLayers.set(route.cluster, []);
        }}
        clusterLayers.get(route.cluster).push(line);
        if (!hiddenClusters.has(route.cluster)) {{
          line.addTo(routeLayer);
        }}
      }});

      renderClusterControl();
      statusEl.textContent =
        `${{methodData.label}}: ${{methodData.clusters.length}} clusters, ` +
        `${{methodData.routes.length}} route weights`;
    }}

    function renderClusterControl() {{
      const methodData = routePayload.methods[activeMethod];
      controlDiv.innerHTML = `
        <strong>${{methodData.label}} Clusters (${{methodData.clusters.length}})</strong>
        ${{methodData.clusters.map((cluster) => `
          <button
            class="cluster-row ${{hiddenClusters.has(cluster.cluster) ? "is-hidden" : ""}}"
            type="button"
            data-cluster="${{cluster.cluster}}"
          >
            <span class="cluster-swatch" style="background:${{cluster.color}}"></span>
            <span>${{cluster.cluster}}</span>
            <span>${{cluster.count}}</span>
          </button>
        `).join("")}}
      `;

      controlDiv.querySelectorAll(".cluster-row").forEach((row) => {{
        row.addEventListener("click", () => {{
          const cluster = row.dataset.cluster;
          if (hiddenClusters.has(cluster)) {{
            hiddenClusters.delete(cluster);
          }} else {{
            hiddenClusters.add(cluster);
          }}
          const isHidden = hiddenClusters.has(cluster);
          row.classList.toggle("is-hidden", isHidden);
          (clusterLayers.get(cluster) || []).forEach((layer) => {{
            if (isHidden) {{
              routeLayer.removeLayer(layer);
            }} else {{
              layer.addTo(routeLayer);
            }}
          }});
        }});
      }});
    }}

    function setActiveMethod(methodKey) {{
      activeMethod = methodKey;
      hiddenClusters = new Set();
      Array.from(methodToggle.querySelectorAll("button")).forEach((button) => {{
        const isActive = button.dataset.method === activeMethod;
        button.classList.toggle("active", isActive);
        button.setAttribute("aria-pressed", isActive ? "true" : "false");
      }});
      drawRoutes();
    }}

    methodOrder.forEach((method) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.method = method.key;
      button.textContent = method.label;
      button.addEventListener("click", () => setActiveMethod(method.key));
      methodToggle.appendChild(button);
    }});

    drawEndpoints();
    setActiveMethod(activeMethod);

    const firstBounds = routePayload.methods[activeMethod].bounds;
    if (firstBounds.length > 0) {{
      map.fitBounds(firstBounds, {{ padding: [28, 28] }});
    }}
  </script>
</body>
</html>
"""


def main():
    """Create one Leaflet route-cluster map for each St. Louis route pair."""
    if not ROUTES_GEOJSON.exists():
        raise FileNotFoundError(f"Weighted routes not found: {ROUTES_GEOJSON}")

    MAPS_FOLDER.mkdir(parents=True, exist_ok=True)
    features_by_id = load_route_features()
    cluster_tables = load_cluster_tables()

    route_pairs = [
        route_pair
        for route_pair in ROUTE_PAIR_ORDER
        if route_pair in set(next(iter(cluster_tables.values()))["route_pair"])
    ]

    for route_pair in route_pairs:
        payload = build_route_pair_payload(route_pair, cluster_tables, features_by_id)
        output_path = MAPS_FOLDER / f"{route_pair}_route_clusters.html"
        output_path.write_text(html_template(payload), encoding="utf-8")
        print(f"Saved route-cluster map: {output_path}")


if __name__ == "__main__":
    main()
