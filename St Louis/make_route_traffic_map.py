"""Overlay St. Louis A* routes on the 24-hour traffic-density heatmap."""

import json
import os

import geopandas as gpd
import pandas as pd

import make_traffic_map
import settings
from generate_astar_toggle_pages import ROUTES, ROUTE_SPECS


OUTPUT_PATH = os.path.join(
    settings.ROUTE_HTML_FOLDER,
    "st_louis_routes_24h_traffic.html",
)


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>St. Louis Routes Over 24-Hour Traffic</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin=""
  />
  <style>
    html, body {{
      height: 100%;
      margin: 0;
      font-family: Arial, sans-serif;
      background: #0d1117;
    }}

    #map {{
      width: 100vw;
      height: 100vh;
    }}

    .panel {{
      position: absolute;
      top: 12px;
      left: 12px;
      z-index: 1000;
      width: min(430px, calc(100vw - 24px));
      max-height: calc(100vh - 24px);
      overflow-y: auto;
      padding: 14px 16px;
      border-radius: 10px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 10px 24px rgba(0, 0, 0, 0.18);
      line-height: 1.35;
    }}

    .panel h1 {{
      margin: 0 0 5px;
      font-size: 18px;
    }}

    .panel p {{
      margin: 0 0 8px;
      color: #475569;
      font-size: 12px;
    }}

    .metric {{
      display: inline-block;
      margin: 2px 6px 0 0;
      padding: 4px 8px;
      border-radius: 999px;
      background: #e2e8f0;
      color: #0f172a;
      font-size: 12px;
    }}

    .toggle-label {{
      margin-top: 10px;
      color: #334155;
      font-size: 12px;
      font-weight: 700;
    }}

    .toggle-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 6px;
    }}

    .toggle-row button {{
      border: 0;
      border-radius: 999px;
      padding: 6px 10px;
      background: #cbd5e1;
      color: #0f172a;
      font-size: 12px;
      cursor: pointer;
    }}

    .toggle-row button.active {{
      background: #0f172a;
      color: #f8fafc;
    }}

    .legend-box {{
      padding: 9px 11px;
      border-radius: 7px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.14);
      font-size: 12px;
      line-height: 1.5;
    }}

    .legend-line {{
      display: inline-block;
      width: 30px;
      margin-right: 7px;
      border-top: 5px solid #dc2626;
      vertical-align: middle;
    }}
  </style>
</head>
<body>
  <div id="map"></div>
  <div class="panel">
    <h1>St. Louis Routes Over 24-Hour Traffic</h1>
    <p>The heatmap always represents the full 24-hour traffic dataset.</p>
    <div class="metric" id="datasetMetric"></div>
    <div class="metric" id="densityMetric"></div>
    <div class="toggle-label">Traffic date</div>
    <div class="toggle-row" id="dateToggle"></div>
    <div class="toggle-label">Route pair</div>
    <div class="toggle-row" id="pairToggle"></div>
    <div class="toggle-label">Route weighting</div>
    <div class="toggle-row" id="modeToggle"></div>
  </div>

  <script
    src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""
  ></script>
  <script src="https://unpkg.com/leaflet.heat/dist/leaflet-heat.js"></script>
  <script>
    const trafficDatasets = {traffic_datasets};
    const routeData = {route_data};
    const routePairs = {route_pairs};
    const routeModes = {route_modes};

    const map = L.map("map", {{
      center: [{center_lat}, {center_lon}],
      zoom: {zoom},
      preferCanvas: true,
      zoomSnap: 0.25
    }});

    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO",
      maxZoom: 19
    }}).addTo(map);

    const heatLayer = L.heatLayer([], {{
      radius: 18,
      blur: 14,
      maxZoom: 12,
      minOpacity: 0.35,
      gradient: {{
        0.2: "#1d4ed8",
        0.4: "#06b6d4",
        0.6: "#facc15",
        0.8: "#f97316",
        1.0: "#dc2626"
      }}
    }}).addTo(map);

    const routeLayer = L.geoJSON(null, {{
      style: (feature) => ({{
        color: feature.properties.color || "#dc2626",
        weight: 7,
        opacity: 0.95
      }}),
      onEachFeature: (feature, layer) => {{
        const p = feature.properties || {{}};
        layer.bindPopup(
          `<b>${{p.route_label_full}}</b><br>` +
          `<b>Weighting:</b> ${{p.route_label}}<br>` +
          `<b>Traffic:</b> ${{p.dataset_label}} (24h)<br>` +
          `<b>Path nodes:</b> ${{p.path_nodes}}<br>` +
          `<b>Total cost:</b> ${{Number(p.total_cost).toFixed(4)}}`
        );
      }}
    }}).addTo(map);

    const endpointLayer = L.layerGroup().addTo(map);
    const dateToggle = document.getElementById("dateToggle");
    const pairToggle = document.getElementById("pairToggle");
    const modeToggle = document.getElementById("modeToggle");
    const datasetMetric = document.getElementById("datasetMetric");
    const densityMetric = document.getElementById("densityMetric");

    let activeDate = Object.keys(trafficDatasets)[0];
    let activePair = routePairs[0].slug;
    let activeMode = "combined";

    function drawSelection() {{
      const traffic = trafficDatasets[activeDate];
      const heatMax = traffic.heatmap.reduce(
        (maximum, row) => Math.max(maximum, row.count),
        0
      ) || 1;
      heatLayer.setLatLngs(
        traffic.heatmap.map((row) => [row.lat, row.lon, row.count / heatMax])
      );

      routeLayer.clearLayers();
      const feature = routeData.features.find((candidate) => {{
        const p = candidate.properties || {{}};
        return p.dataset_slug === activeDate
          && p.route_slug === activePair
          && p.route_name === activeMode;
      }});
      if (feature) routeLayer.addData(feature);

      endpointLayer.clearLayers();
      const pair = routePairs.find((candidate) => candidate.slug === activePair);
      [
        {{ label: pair.start.label, point: pair.start, fill: "#ffffff" }},
        {{ label: pair.destination.label, point: pair.destination, fill: "#facc15" }}
      ].forEach((endpoint) => {{
        L.circleMarker([endpoint.point.lat, endpoint.point.lon], {{
          radius: 7,
          color: "#111827",
          weight: 2,
          fillColor: endpoint.fill,
          fillOpacity: 1
        }}).bindPopup(`<b>${{endpoint.label}}</b>`).addTo(endpointLayer);
      }});

      datasetMetric.textContent = `Dataset: ${{traffic.label}}`;
      densityMetric.textContent = `24h density cells: ${{traffic.heatmap.length.toLocaleString()}}`;
      const mode = routeModes.find((candidate) => candidate.name === activeMode);
      const legendLabel = document.getElementById("legendLabel");
      const legendLine = document.querySelector(".legend-line");
      if (legendLabel && mode) legendLabel.textContent = mode.label;
      if (legendLine && mode) legendLine.style.borderTopColor = mode.color;
      updateButtons();
    }}

    function updateButtons() {{
      document.querySelectorAll("#dateToggle button").forEach((button) => {{
        button.classList.toggle("active", button.dataset.value === activeDate);
      }});
      document.querySelectorAll("#pairToggle button").forEach((button) => {{
        button.classList.toggle("active", button.dataset.value === activePair);
      }});
      document.querySelectorAll("#modeToggle button").forEach((button) => {{
        button.classList.toggle("active", button.dataset.value === activeMode);
      }});
    }}

    function addButtons(container, rows, valueKey, labelKey, onSelect) {{
      rows.forEach((row) => {{
        const button = document.createElement("button");
        button.type = "button";
        button.dataset.value = row[valueKey];
        button.textContent = row[labelKey];
        button.addEventListener("click", () => {{
          onSelect(row[valueKey]);
          drawSelection();
        }});
        container.appendChild(button);
      }});
    }}

    addButtons(
      dateToggle,
      Object.entries(trafficDatasets).map(([date, data]) => ({{ date, label: data.label }})),
      "date",
      "label",
      (value) => activeDate = value
    );
    addButtons(pairToggle, routePairs, "slug", "short_label", (value) => activePair = value);
    addButtons(modeToggle, routeModes, "name", "label", (value) => activeMode = value);

    L.control.layers(null, {{
      "24-Hour Traffic Density": heatLayer,
      "Selected A* Route": routeLayer,
      "Route Endpoints": endpointLayer
    }}, {{ collapsed: false }}).addTo(map);

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function() {{
      const div = L.DomUtil.create("div", "legend-box");
      div.innerHTML =
        `<div><b>Selected route</b></div>` +
        `<div><span class="legend-line"></span><span id="legendLabel"></span></div>`;
      return div;
    }};
    legend.addTo(map);

    drawSelection();
  </script>
</body>
</html>
"""


def build_24_hour_traffic():
    """Build only the 24-hour density payload needed by the overlay page."""
    output = {}
    for dataset in settings.TRAFFIC_DATASETS:
        csv_path = os.path.join(
            settings.TRAFFIC_FOLDER,
            dataset["filename"],
        )
        if not os.path.isfile(csv_path):
            raise FileNotFoundError(f"Traffic file not found: {csv_path}")

        data = make_traffic_map.load_flight_data(csv_path)
        heatmap = make_traffic_map.build_windowed_heatmap(data, [24])["24"]
        output[dataset["date"]] = {
            "label": dataset["label"],
            "heatmap": heatmap,
        }
        print(f"Loaded 24h density for {dataset['label']}: {len(heatmap)} cells")
    return output


def load_routes():
    """Load all route-pair comparison files into one feature collection."""
    frames = []
    for route in ROUTES:
        route_path = os.path.join(
            settings.ROUTE_GEOJSON_FOLDER,
            f"{route['slug']}_astar_routes.geojson",
        )
        if not os.path.isfile(route_path):
            raise FileNotFoundError(
                f"Route file not found: {route_path}. "
                "Run generate_astar_toggle_pages.py first."
            )
        frames.append(gpd.read_file(route_path))

    combined = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True),
        crs=frames[0].crs,
    )
    if pd.api.types.is_datetime64_any_dtype(combined["dataset_slug"]):
        combined["dataset_slug"] = combined["dataset_slug"].dt.strftime(
            "%Y-%m-%d"
        )
    return json.loads(combined.to_json())


def main():
    traffic_datasets = build_24_hour_traffic()
    route_data = load_routes()
    route_pairs = [
        {
            **route,
            "short_label": route["label"]
            .replace("St. Louis ", "")
            .replace("MidAmerica To ", "MidAmerica to "),
        }
        for route in ROUTES
    ]
    route_modes = [
        {
            "name": spec["name"],
            "label": spec["label"],
            "color": spec["color"],
        }
        for spec in ROUTE_SPECS
    ]

    html = HTML_TEMPLATE.format(
        traffic_datasets=json.dumps(traffic_datasets),
        route_data=json.dumps(route_data),
        route_pairs=json.dumps(route_pairs),
        route_modes=json.dumps(route_modes),
        center_lat=settings.MAP_CENTER_LAT,
        center_lon=settings.MAP_CENTER_LON,
        zoom=settings.MAP_ZOOM,
    )
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, "w", encoding="utf-8") as output_file:
        output_file.write(html)
    print(f"Saved route and 24h traffic overlay: {OUTPUT_PATH}")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as error:
        raise SystemExit(str(error)) from error
