"""Create a standalone Leaflet route map from API route outputs."""

import json
import os

import geopandas as gpd


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{title}</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
    integrity="sha256-p4NxAoJBhIINfQPDQo4O4FwBdxuV6pS2j3u2Bf8e6a4="
    crossorigin=""
  />
  <style>
    :root {{
      color-scheme: light;
      --surface: #ffffff;
      --surface-soft: #f5f7f8;
      --line: #cfd7dd;
      --text: #172026;
      --muted: #5f6f7a;
      --low: #2ca25f;
      --medium-combined: #f4b183;
      --high-combined: #c84c5a;
    }}

    * {{
      box-sizing: border-box;
    }}

    html, body {{
      height: 100%;
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--text);
      overflow: hidden;
      background: var(--surface-soft);
    }}

    body {{
      display: flex;
      flex-direction: column;
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
      font-weight: 700;
      line-height: 1.2;
    }}

    .status {{
      margin-top: 2px;
      color: var(--muted);
      font-size: 12px;
      line-height: 1.25;
      white-space: nowrap;
    }}

    #map {{
      flex: 1 1 auto;
      width: 100%;
      height: calc(100vh - 59px);
      min-height: 360px;
    }}

    .leaflet-container {{
      background: #edf1f2;
    }}

    .legend-box {{
      padding: 10px 12px;
      border: 1px solid rgba(23, 32, 38, 0.15);
      border-radius: 7px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.14);
      color: var(--text);
      font-size: 12px;
      line-height: 1.35;
    }}

    .legend-title {{
      margin-bottom: 7px;
      font-weight: 700;
    }}

    .legend-row {{
      display: grid;
      grid-template-columns: 32px auto;
      align-items: center;
      gap: 7px;
      margin: 5px 0;
      white-space: nowrap;
    }}

    .line {{
      width: 32px;
      height: 0;
      border-top: 4px solid #2563eb;
    }}
  </style>
</head>
<body>
  <header class="toolbar">
    <div>
      <h1>{title}</h1>
      <div id="status" class="status">Loading route map...</div>
    </div>
    <div id="counts" class="status"></div>
  </header>
  <div id="map"></div>

  <script
    src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"
    integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo="
    crossorigin=""
  ></script>
  <script>
    const routeData = {route_data};
    const gridData = {grid_data};
    const locations = {locations};
    const statusEl = document.getElementById("status");
    const countsEl = document.getElementById("counts");

    const map = L.map("map", {{
      preferCanvas: true,
      zoomControl: true
    }});
    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);

    const routeColors = [
      "#2563eb", "#dc2626", "#16a34a", "#9333ea", "#f59e0b",
      "#0891b2", "#be123c", "#4f46e5"
    ];

    function riskColor(properties) {{
      const riskClass = properties.risk_class || "Low";
      if (riskClass === "High") return "#dc2626";
      if (riskClass === "Medium") return "#f59e0b";
      return "#16a34a";
    }}

    let gridLayer = null;
    if (gridData && gridData.features && gridData.features.length > 0) {{
      gridLayer = L.geoJSON(gridData, {{
        renderer: L.canvas({{ padding: 0.35 }}),
        style: (feature) => ({{
          color: riskColor(feature.properties || {{}}),
          weight: 0.35,
          opacity: 0.45,
          fillColor: riskColor(feature.properties || {{}}),
          fillOpacity: 0.36
        }}),
        onEachFeature: (feature, layer) => {{
          const p = feature.properties || {{}};
          layer.bindPopup(
            `<b>Grid Cell</b><br>` +
            `Risk class: ${{p.risk_class || "N/A"}}<br>` +
            `Risk cost: ${{Number(p.risk_cost || 0).toFixed(2)}}<br>` +
            `Population risk: ${{Number(p.city_risk || 0).toFixed(2)}}<br>` +
            `Airspace risk: ${{Number(p.airport_risk_combined || 0).toFixed(2)}}<br>` +
            `Traffic risk: ${{Number(p.traffic_risk || 0).toFixed(2)}}`
          );
        }}
      }}).addTo(map);
    }}

    const routeLayer = L.geoJSON(routeData, {{
      style: (feature) => {{
        const index = routeData.features.indexOf(feature);
        return {{
          color: routeColors[index % routeColors.length],
          weight: 5,
          opacity: 0.92
        }};
      }},
      onEachFeature: (feature, layer) => {{
        const p = feature.properties || {{}};
        layer.bindPopup(
          `<b>${{p.route_id || "Route"}}</b><br>` +
          `Total cost: ${{Number(p.total_cost || 0).toFixed(4)}}<br>` +
          `Distance: ${{Number(p.route_distance_km || 0).toFixed(2)}} km<br>` +
          `Path nodes: ${{p.path_nodes || "N/A"}}`
        );
      }}
    }}).addTo(map);

    const markerLayer = L.layerGroup().addTo(map);
    Object.values(locations || {{}}).forEach((location) => {{
      if (location.lat !== undefined && location.lon !== undefined) {{
        L.marker([location.lat, location.lon])
          .bindPopup(`<b>${{location.label || "Location"}}</b><br>${{location.type || ""}}`)
          .addTo(markerLayer);
      }}
    }});

    const bounds = L.latLngBounds([]);
    if (gridLayer) {{
      bounds.extend(gridLayer.getBounds());
    }}
    bounds.extend(routeLayer.getBounds());
    markerLayer.eachLayer((layer) => bounds.extend(layer.getLatLng()));

    if (bounds.isValid()) {{
      map.fitBounds(bounds.pad(0.08));
      if (gridLayer) {{
        map.setMaxBounds(gridLayer.getBounds().pad(0.04));
      }}
      statusEl.textContent = "Route map loaded";
    }} else {{
      map.setView([39.5, -98.35], 4);
      statusEl.textContent = "Route map has no valid bounds";
    }}

    countsEl.textContent =
      `Routes: ${{routeData.features ? routeData.features.length : 0}} | ` +
      `Grid cells: ${{gridData && gridData.features ? gridData.features.length.toLocaleString() : 0}}`;

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "legend-box");
      div.innerHTML = `
        <div class="legend-title">Layers</div>
        <div class="legend-row"><span class="line"></span><span>Route</span></div>
        <div class="legend-row"><span style="width:16px;height:16px;background:#16a34a;opacity:0.4;display:inline-block;"></span><span>Low risk</span></div>
        <div class="legend-row"><span style="width:16px;height:16px;background:#f59e0b;opacity:0.4;display:inline-block;"></span><span>Medium risk</span></div>
        <div class="legend-row"><span style="width:16px;height:16px;background:#dc2626;opacity:0.4;display:inline-block;"></span><span>High risk</span></div>
      `;
      return div;
    }};
    legend.addTo(map);

    if (gridLayer) {{
      L.control.layers(null, {{
        "Risk grid": gridLayer,
        "Routes": routeLayer,
        "Locations": markerLayer
      }}).addTo(map);
    }}
  </script>
</body>
</html>
"""


def load_geojson(path):
    """Load a geospatial file and return GeoJSON data as a dictionary."""
    if not path or not os.path.exists(path):
        return {"type": "FeatureCollection", "features": []}

    geodata = gpd.read_file(path).to_crs("EPSG:4326")
    return json.loads(geodata.to_json(default=str))


def save_route_map_html(html_text, output_html_path):
    """Save route map HTML."""
    output_folder = os.path.dirname(output_html_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    with open(output_html_path, "w", encoding="utf-8") as output_file:
        output_file.write(html_text)


def make_route_map(
    route_geojson_path=None,
    scored_grid_path=None,
    locations=None,
    output_html_path=None,
    title="eVTOL Route Map",
):
    """Create a standalone Leaflet HTML map for route results."""
    if route_geojson_path is None:
        route_geojson_path = input("Route GeoJSON path: ").strip()
        scored_grid_path = input("Scored grid GeoJSON path (optional): ").strip() or None
        output_html_path = input("Route map HTML output path: ").strip()

    if not output_html_path:
        raise ValueError("output_html_path is required.")
    if not os.path.exists(route_geojson_path):
        raise FileNotFoundError(f"Route GeoJSON not found: {route_geojson_path}")

    route_data = load_geojson(route_geojson_path)
    grid_data = load_geojson(scored_grid_path)
    html_text = HTML_TEMPLATE.format(
        title=title,
        route_data=json.dumps(route_data),
        grid_data=json.dumps(grid_data),
        locations=json.dumps(locations or {}),
    )
    save_route_map_html(html_text, output_html_path)

    return {
        "route_geojson_path": route_geojson_path,
        "scored_grid_path": scored_grid_path,
        "output_html_path": output_html_path,
        "output_exists": os.path.exists(output_html_path),
        "route_count": len(route_data.get("features", [])),
        "grid_cell_count": len(grid_data.get("features", [])),
    }


if __name__ == "__main__":
    result = make_route_map()
    print(result)
