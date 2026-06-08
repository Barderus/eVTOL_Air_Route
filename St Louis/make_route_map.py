"""Create the interactive St. Louis route demonstration map."""

import json
import os

import geopandas as gpd

import settings


ROUTE_GEOJSON_PATH = os.path.join(
    settings.ROUTE_GEOJSON_FOLDER,
    "st_louis_astar_routes.geojson",
)
OUTPUT_PATH = settings.ROUTE_MAP_HTML

ROUTE_STYLES = {
    "MidAmerica to Downtown St. Louis": {
        "short_label": "MidAmerica to Downtown",
        "color": "#dc2626",
    },
    "MidAmerica to St. Louis Lambert": {
        "short_label": "MidAmerica to Lambert",
        "color": "#2563eb",
    },
    "St. Louis Downtown Airport to St. Louis Lambert": {
        "short_label": "Downtown Airport to Lambert",
        "color": "#16a34a",
    },
}


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>St. Louis eVTOL Routes</title>
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

    h1 {{
      margin: 0;
      font-size: 17px;
      line-height: 1.2;
    }}

    .status {{
      margin-top: 2px;
      color: var(--muted);
      font-size: 12px;
    }}

    .control-strip {{
      display: flex;
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

    .route-toggle {{ grid-template-columns: repeat(3, minmax(150px, 1fr)); }}
    .date-toggle {{ grid-template-columns: repeat(6, minmax(125px, 1fr)); }}

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
      flex: 1 1 auto;
      width: 100%;
      height: calc(100vh - 59px);
      min-height: 360px;
    }}

    .legend-box {{
      padding: 10px 12px;
      border-radius: 7px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.14);
      font-size: 12px;
      line-height: 1.4;
    }}

    .legend-row {{
      display: flex;
      align-items: center;
      gap: 8px;
      margin: 5px 0;
    }}

    .swatch-line {{
      width: 34px;
      border-top: 5px solid #000000;
    }}

    @media (max-width: 1050px) {{
      .toolbar {{ align-items: stretch; flex-direction: column; gap: 9px; }}
      .control-strip {{ justify-content: stretch; }}
      .segmented {{ width: 100%; }}
      .route-toggle {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
      .date-toggle {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
      #map {{ height: calc(100vh - 145px); }}
    }}
  </style>
</head>
<body>
  <header class="toolbar">
    <div>
      <h1>St. Louis eVTOL Routes</h1>
      <div id="status" class="status">Select routes and a traffic date.</div>
    </div>
    <div class="control-strip">
      <div class="segmented route-toggle" id="routeToggle"></div>
      <div class="segmented date-toggle" id="dateToggle"></div>
    </div>
  </header>

  <main id="map" aria-label="St. Louis eVTOL route map"></main>

  <script
    src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""
  ></script>
  <script>
    const routeData = {route_data};
    const routeStyles = {route_styles};
    const routeOrder = {route_order};
    const datasetOrder = {dataset_order};
    const endpointData = {endpoint_data};

    const map = L.map("map", {{
      center: [{center_lat}, {center_lon}],
      zoom: {zoom},
      preferCanvas: true
    }});

    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);

    const routeToggle = document.getElementById("routeToggle");
    const dateToggle = document.getElementById("dateToggle");
    const statusEl = document.getElementById("status");
    const routeLayers = {{}};
    const endpointLayer = L.layerGroup().addTo(map);
    let activeDate = datasetOrder[0].date;
    let activeRoutes = new Set(routeOrder);

    routeOrder.forEach((routeName) => {{
      routeLayers[routeName] = L.geoJSON(null, {{
        style: {{
          color: routeStyles[routeName].color,
          weight: 6,
          opacity: 0.92
        }},
        onEachFeature: (feature, layer) => {{
          const properties = feature.properties || {{}};
          layer.bindPopup(
            `<b>${{properties.route}}</b><br>` +
            `<b>Traffic:</b> ${{properties.dataset}}<br>` +
            `<b>Origin:</b> ${{properties.origin}}<br>` +
            `<b>Destination:</b> ${{properties.destination}}<br>` +
            `<b>Path nodes:</b> ${{properties.path_nodes}}`
          );
        }}
      }});
    }});

    function updateEndpoints() {{
      endpointLayer.clearLayers();
      const endpointKeys = new Set();

      routeOrder.forEach((routeName) => {{
        if (!activeRoutes.has(routeName)) return;
        const endpoints = endpointData[routeName];

        [
          {{ type: "Origin", data: endpoints.start }},
          {{ type: "Destination", data: endpoints.destination }}
        ].forEach((endpoint) => {{
          const key = `${{endpoint.data[0]}},${{endpoint.data[1]}}`;
          if (endpointKeys.has(key)) return;
          endpointKeys.add(key);

          L.circleMarker(endpoint.data, {{
            radius: 7,
            color: "#172026",
            weight: 2,
            fillColor: endpoint.type === "Origin" ? "#ffffff" : "#facc15",
            fillOpacity: 1
          }}).bindPopup(
            `<b>${{endpoint.type}}</b><br>${{endpoints[endpoint.type === "Origin" ? "start_label" : "destination_label"]}}`
          ).addTo(endpointLayer);
        }});
      }});
    }}

    function drawSelection() {{
      Object.values(routeLayers).forEach((layer) => {{
        layer.clearLayers();
        map.removeLayer(layer);
      }});

      routeData.features.forEach((feature) => {{
        const properties = feature.properties || {{}};
        if (
          properties.date === activeDate &&
          activeRoutes.has(properties.route)
        ) {{
          routeLayers[properties.route].addData(feature);
        }}
      }});

      routeOrder.forEach((routeName) => {{
        const layer = routeLayers[routeName];
        if (activeRoutes.has(routeName) && layer.getLayers().length > 0) {{
          layer.addTo(map);
        }}
      }});

      updateEndpoints();
      updateButtons();

      const selectedDataset = datasetOrder.find(
        (dataset) => dataset.date === activeDate
      );
      statusEl.textContent =
        `${{activeRoutes.size}} route(s), ${{selectedDataset.label}} traffic`;
    }}

    function updateButtons() {{
      Array.from(routeToggle.querySelectorAll("button")).forEach((button) => {{
        const isActive = activeRoutes.has(button.dataset.route);
        button.classList.toggle("active", isActive);
        button.setAttribute("aria-pressed", isActive ? "true" : "false");
      }});

      Array.from(dateToggle.querySelectorAll("button")).forEach((button) => {{
        const isActive = button.dataset.date === activeDate;
        button.classList.toggle("active", isActive);
        button.setAttribute("aria-pressed", isActive ? "true" : "false");
      }});
    }}

    routeOrder.forEach((routeName) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.route = routeName;
      button.textContent = routeStyles[routeName].short_label;
      button.addEventListener("click", () => {{
        if (activeRoutes.has(routeName)) {{
          if (activeRoutes.size === 1) return;
          activeRoutes.delete(routeName);
        }} else {{
          activeRoutes.add(routeName);
        }}
        drawSelection();
      }});
      routeToggle.appendChild(button);
    }});

    datasetOrder.forEach((dataset) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.date = dataset.date;
      button.textContent = dataset.label;
      button.addEventListener("click", () => {{
        activeDate = dataset.date;
        drawSelection();
      }});
      dateToggle.appendChild(button);
    }});

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function() {{
      const div = L.DomUtil.create("div", "legend-box");
      const rows = routeOrder.map((routeName) => `
        <div class="legend-row">
          <span class="swatch-line" style="border-top-color:${{routeStyles[routeName].color}}"></span>
          <span>${{routeStyles[routeName].short_label}}</span>
        </div>
      `).join("");
      div.innerHTML = `<b>Routes</b>${{rows}}`;
      return div;
    }};
    legend.addTo(map);

    drawSelection();
  </script>
</body>
</html>
"""


def build_endpoint_data():
    endpoint_data = {}

    for route in settings.ROUTES:
        endpoint_data[route["label"]] = {
            "start_label": route["start_label"],
            "start": list(route["start"]),
            "destination_label": route["destination_label"],
            "destination": list(route["destination"]),
        }

    return endpoint_data


def main():
    if not os.path.isfile(ROUTE_GEOJSON_PATH):
        raise FileNotFoundError(
            f"Route GeoJSON not found: {ROUTE_GEOJSON_PATH}. "
            "Run astar_routes.py first."
        )

    routes = gpd.read_file(ROUTE_GEOJSON_PATH)
    routes["date"] = routes["date"].dt.strftime("%Y-%m-%d")
    route_data = json.loads(routes.to_json())

    dataset_order = []
    for dataset in settings.TRAFFIC_DATASETS:
        dataset_order.append(
            {
                "date": dataset["date"],
                "label": dataset["label"],
            }
        )

    route_order = [route["label"] for route in settings.ROUTES]
    endpoint_data = build_endpoint_data()

    html = HTML_TEMPLATE.format(
        route_data=json.dumps(route_data),
        route_styles=json.dumps(ROUTE_STYLES),
        route_order=json.dumps(route_order),
        dataset_order=json.dumps(dataset_order),
        endpoint_data=json.dumps(endpoint_data),
        center_lat=settings.MAP_CENTER_LAT,
        center_lon=settings.MAP_CENTER_LON,
        zoom=settings.MAP_ZOOM,
    )

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, "w", encoding="utf-8") as output_file:
        output_file.write(html)

    print(f"Saved St. Louis route map: {OUTPUT_PATH}")
    print(f"Routes embedded: {len(routes)}")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as error:
        raise SystemExit(str(error)) from error
