"""Create an HTML map for the 10-route robustness pilot."""

import json
import os

import geopandas as gpd


PILOT_ROUTES_GEOJSON = "route_robustness/output/pilot_routes.geojson"
OUTPUT_FOLDER = "route_robustness/maps"
OUTPUT_HTML = os.path.join(OUTPUT_FOLDER, "pilot_routes.html")

ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521

MAP_CENTER_LAT = 41.79
MAP_CENTER_LON = -87.88
MAP_ZOOM = 10

ROUTE_COLORS = {
    "distance_only": "#1b9e77",
    "population_only": "#d95f02",
    "traffic_only": "#7570b3",
    "airspace_only": "#e7298a",
    "equal_weight": "#000000",
    "distance_heavy": "#66a61e",
    "population_heavy": "#e6ab02",
    "traffic_heavy": "#a6761d",
    "airspace_heavy": "#1f78b4",
    "mixed_balanced": "#b15928",
}


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Route Robustness Pilot Routes</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin=""
  />
  <style>
    html,
    body {{
      height: 100%;
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: #172026;
    }}

    #map {{
      width: 100%;
      height: 100%;
    }}

    .panel {{
      background: #ffffff;
      padding: 10px 12px;
      border-radius: 7px;
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.16);
      font-size: 13px;
      line-height: 1.35;
    }}

    .panel h1 {{
      margin: 0 0 6px;
      font-size: 16px;
      line-height: 1.2;
    }}

    .legend-row {{
      display: flex;
      align-items: center;
      gap: 8px;
      margin: 4px 0;
      white-space: nowrap;
    }}

    .legend-line {{
      display: inline-block;
      width: 30px;
      border-top: 4px solid #000000;
    }}
  </style>
</head>
<body>
  <div id="map"></div>
  <script
    src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""
  ></script>
  <script>
    const map = L.map("map", {{
      center: [{map_center_lat}, {map_center_lon}],
      zoom: {map_zoom}
    }});

    L.tileLayer("https://{{s}}.tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors"
    }}).addTo(map);

    const routeData = {route_data};
    const routeColors = {route_colors};
    const routeLayers = {{}};
    const layerControlItems = {{}};

    L.marker([{origin_lat}, {origin_lon}])
      .bindPopup("<b>Origin</b><br>{origin_label}")
      .addTo(map);

    L.marker([{destination_lat}, {destination_lon}])
      .bindPopup("<b>Destination</b><br>{destination_label}")
      .addTo(map);

    function formatNumber(value, digits) {{
      const numberValue = Number(value);
      if (Number.isNaN(numberValue)) {{
        return "n/a";
      }}
      return numberValue.toFixed(digits);
    }}

    routeData.features.forEach((feature) => {{
      const properties = feature.properties || {{}};
      const pilotName = properties.pilot_name || "route";
      const routeColor = routeColors[pilotName] || "#dc2626";
      const layer = L.geoJSON(feature, {{
        style: {{
          color: routeColor,
          weight: pilotName === "equal_weight" ? 7 : 5,
          opacity: pilotName === "equal_weight" ? 1.0 : 0.86
        }},
        onEachFeature: (routeFeature, routeLayer) => {{
          const p = routeFeature.properties || {{}};
          routeLayer.bindPopup(
            `<b>${{p.pilot_name}}</b><br>` +
            `<b>Weight ID:</b> ${{p.weight_id}}<br>` +
            `<b>Distance:</b> ${{formatNumber(p.route_distance_km, 2)}} km<br>` +
            `<b>Score:</b> ${{formatNumber(p.total_weighted_score, 4)}}<br>` +
            `<b>Path nodes:</b> ${{p.path_nodes}}<br>` +
            `<b>Weights:</b> D=${{p.distance_weight}}, P=${{p.population_weight}}, T=${{p.traffic_weight}}, A=${{p.airspace_weight}}`
          );
        }}
      }});
      routeLayers[pilotName] = layer;
      layerControlItems[pilotName.replaceAll("_", " ")] = layer;
      layer.addTo(map);
    }});

    L.control.layers(null, layerControlItems, {{ collapsed: false }}).addTo(map);

    const bounds = L.geoJSON(routeData).getBounds();
    if (bounds.isValid()) {{
      map.fitBounds(bounds.pad(0.12));
    }}

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "panel");
      const rows = routeData.features.map((feature) => {{
        const p = feature.properties || {{}};
        const pilotName = p.pilot_name || "route";
        const color = routeColors[pilotName] || "#dc2626";
        return `
          <div class="legend-row">
            <span class="legend-line" style="border-top-color:${{color}};"></span>
            <span>${{pilotName.replaceAll("_", " ")}}</span>
          </div>
        `;
      }}).join("");
      div.innerHTML = `
        <h1>Route Robustness Pilot</h1>
        <div>${{rows}}</div>
      `;
      return div;
    }};
    legend.addTo(map);
  </script>
</body>
</html>
"""


def main():
    """Write the pilot route map HTML."""
    if not os.path.exists(PILOT_ROUTES_GEOJSON):
        raise FileNotFoundError(
            f"Pilot route GeoJSON not found: {PILOT_ROUTES_GEOJSON}"
        )

    routes = gpd.read_file(PILOT_ROUTES_GEOJSON)
    if "path_node_ids" in routes.columns:
        routes = routes.drop(columns=["path_node_ids"])
    route_data = json.loads(routes.to_json())

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    html = HTML_TEMPLATE.format(
        route_data=json.dumps(route_data),
        route_colors=json.dumps(ROUTE_COLORS),
        origin_label=ORIGIN_LABEL,
        origin_lat=ORIGIN_LAT,
        origin_lon=ORIGIN_LON,
        destination_label=DESTINATION_LABEL,
        destination_lat=DESTINATION_LAT,
        destination_lon=DESTINATION_LON,
        map_center_lat=MAP_CENTER_LAT,
        map_center_lon=MAP_CENTER_LON,
        map_zoom=MAP_ZOOM,
    )

    with open(OUTPUT_HTML, "w", encoding="utf-8") as output_file:
        output_file.write(html)

    print(f"Saved pilot route map: {OUTPUT_HTML}")


if __name__ == "__main__":
    main()
