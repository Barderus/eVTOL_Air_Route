"""Create an HTML map for 0.15-threshold representative cluster routes."""

import json
import os

import geopandas as gpd


REPRESENTATIVES_GEOJSON = (
    "route_robustness/output/all_cluster_representative_routes_threshold_015.geojson"
)
OUTPUT_FOLDER = "route_robustness/maps"
OUTPUT_HTML = os.path.join(
    OUTPUT_FOLDER,
    "all_cluster_representative_routes_threshold_015.html",
)

ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521

MAP_CENTER_LAT = 41.79
MAP_CENTER_LON = -87.88
MAP_ZOOM = 10

CLUSTER_COLORS = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
]


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Route Robustness Representative Routes Threshold 0.15</title>
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
      font-family: Arial, Helvetica, sans-serif;
      color: #172026;
    }}
    #map {{ width: 100%; height: 100%; }}
    .panel {{
      background: #ffffff;
      padding: 10px 12px;
      border-radius: 7px;
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.16);
      font-size: 13px;
      line-height: 1.35;
      max-width: 390px;
      max-height: 48vh;
      overflow: auto;
    }}
    .panel h1 {{ margin: 0 0 6px; font-size: 16px; line-height: 1.2; }}
    .legend-row {{ display: flex; align-items: flex-start; gap: 8px; margin: 5px 0; }}
    .legend-line {{
      display: inline-block;
      flex: 0 0 auto;
      width: 30px;
      margin-top: 7px;
      border-top: 4px solid #000000;
    }}
    .legend-label {{ min-width: 0; }}
    .legend-label strong {{ display: block; }}
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
    const clusterColors = {cluster_colors};
    const layerControlItems = {{}};

    L.marker([{origin_lat}, {origin_lon}]).bindPopup("<b>Origin</b><br>{origin_label}").addTo(map);
    L.marker([{destination_lat}, {destination_lon}]).bindPopup("<b>Destination</b><br>{destination_label}").addTo(map);

    function formatNumber(value, digits) {{
      const numberValue = Number(value);
      if (Number.isNaN(numberValue)) {{ return "n/a"; }}
      return numberValue.toFixed(digits);
    }}

    routeData.features.forEach((feature) => {{
      const properties = feature.properties || {{}};
      const clusterId = properties.cluster_id || "unclustered";
      const routeColor = clusterColors[clusterId] || "#666666";
      const layer = L.geoJSON(feature, {{
        style: {{
          color: routeColor,
          weight: properties.cluster_size > 20 ? 7 : 5,
          opacity: 0.92
        }},
        onEachFeature: (routeFeature, routeLayer) => {{
          const p = routeFeature.properties || {{}};
          routeLayer.bindPopup(
            `<b>${{p.cluster_id}}</b><br>` +
            `<b>Representative:</b> ${{p.weight_id}}<br>` +
            `<b>Cluster size:</b> ${{p.cluster_size}} routes<br>` +
            `<b>Weight space:</b> ${{formatNumber(p.cluster_weight_space_percent, 1)}}%<br>` +
            `<b>Distance:</b> ${{formatNumber(p.route_distance_km, 2)}} km<br>` +
            `<b>Score:</b> ${{formatNumber(p.total_weighted_score, 4)}}<br>` +
            `<b>Mean similarity:</b> ${{formatNumber(p.representative_mean_similarity, 3)}}<br>` +
            `<b>Weights:</b> D=${{p.distance_weight}}, P=${{p.population_weight}}, T=${{p.traffic_weight}}, A=${{p.airspace_weight}}`
          );
        }}
      }}).addTo(map);
      layerControlItems[`${{clusterId}} (${{properties.cluster_size}} routes)`] = layer;
    }});

    L.control.layers(null, layerControlItems, {{ collapsed: false }}).addTo(map);

    const bounds = L.geoJSON(routeData).getBounds();
    if (bounds.isValid()) {{ map.fitBounds(bounds.pad(0.12)); }}

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "panel");
      const rows = routeData.features.map((feature) => {{
        const p = feature.properties || {{}};
        const color = clusterColors[p.cluster_id] || "#666666";
        return `
          <div class="legend-row">
            <span class="legend-line" style="border-top-color:${{color}};"></span>
            <span class="legend-label">
              <strong>${{p.cluster_id}} (${{p.cluster_size}} routes)</strong>
              representative ${{p.weight_id}},
              ${{formatNumber(p.cluster_weight_space_percent, 1)}}% of tested weights
            </span>
          </div>
        `;
      }}).join("");
      div.innerHTML = `<h1>Representative Routes Threshold 0.15</h1><div>${{rows}}</div>`;
      return div;
    }};
    legend.addTo(map);
  </script>
</body>
</html>
"""


def build_cluster_colors(cluster_ids):
    """Assign stable display colors to clusters."""
    return {
        cluster_id: CLUSTER_COLORS[index % len(CLUSTER_COLORS)]
        for index, cluster_id in enumerate(sorted(cluster_ids))
    }


def main():
    """Write the 0.15-threshold representative route map HTML."""
    if not os.path.exists(REPRESENTATIVES_GEOJSON):
        raise FileNotFoundError(
            f"Representative route GeoJSON not found: {REPRESENTATIVES_GEOJSON}"
        )

    representatives = gpd.read_file(REPRESENTATIVES_GEOJSON)
    cluster_colors = build_cluster_colors(representatives["cluster_id"].dropna().unique())
    route_data = json.loads(representatives.to_json())

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    html = HTML_TEMPLATE.format(
        route_data=json.dumps(route_data),
        cluster_colors=json.dumps(cluster_colors),
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

    print(f"Saved 0.15 representative route map: {OUTPUT_HTML}")


if __name__ == "__main__":
    main()
