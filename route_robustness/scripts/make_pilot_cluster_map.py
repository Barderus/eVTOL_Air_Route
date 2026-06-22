"""Create an HTML map for clustered pilot routes."""

import json
import os

import geopandas as gpd
import pandas as pd


PILOT_ROUTES_GEOJSON = "route_robustness/output/pilot_routes.geojson"
PILOT_ROUTE_CLUSTERS_CSV = "route_robustness/output/pilot_route_clusters.csv"
OUTPUT_FOLDER = "route_robustness/maps"
OUTPUT_HTML = os.path.join(OUTPUT_FOLDER, "pilot_route_clusters.html")

ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521

MAP_CENTER_LAT = 41.79
MAP_CENTER_LON = -87.88
MAP_ZOOM = 10

CLUSTER_COLORS = {
    "cluster_01": "#1b9e77",
    "cluster_02": "#d95f02",
    "cluster_03": "#7570b3",
    "cluster_04": "#e7298a",
    "cluster_05": "#000000",
    "cluster_06": "#66a61e",
    "cluster_07": "#e6ab02",
    "cluster_08": "#a6761d",
    "cluster_09": "#1f78b4",
}


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Route Robustness Pilot Clusters</title>
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
      max-width: 360px;
    }}

    .panel h1 {{
      margin: 0 0 6px;
      font-size: 16px;
      line-height: 1.2;
    }}

    .legend-row {{
      display: flex;
      align-items: flex-start;
      gap: 8px;
      margin: 5px 0;
    }}

    .legend-line {{
      display: inline-block;
      flex: 0 0 auto;
      width: 30px;
      margin-top: 7px;
      border-top: 4px solid #000000;
    }}

    .legend-label {{
      min-width: 0;
    }}

    .legend-label strong {{
      display: block;
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
    const clusterColors = {cluster_colors};
    const clusterRoutes = {cluster_routes};
    const clusterLayers = {{}};
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
      const clusterId = properties.cluster_id || "unclustered";
      const routeColor = clusterColors[clusterId] || "#dc2626";
      if (!clusterLayers[clusterId]) {{
        clusterLayers[clusterId] = L.layerGroup().addTo(map);
        const routeNames = clusterRoutes[clusterId] || [clusterId];
        const label = `${{clusterId}}: ${{routeNames.map((name) => name.replaceAll("_", " ")).join(", ")}}`;
        layerControlItems[label] = clusterLayers[clusterId];
      }}
      const layer = L.geoJSON(feature, {{
        style: {{
          color: routeColor,
          weight: properties.cluster_size > 1 ? 7 : 5,
          opacity: properties.cluster_size > 1 ? 1.0 : 0.86
        }},
        onEachFeature: (routeFeature, routeLayer) => {{
          const p = routeFeature.properties || {{}};
          routeLayer.bindPopup(
            `<b>${{p.pilot_name}}</b><br>` +
            `<b>Cluster:</b> ${{p.cluster_id}}<br>` +
            `<b>Cluster size:</b> ${{p.cluster_size}}<br>` +
            `<b>Distance:</b> ${{formatNumber(p.route_distance_km, 2)}} km<br>` +
            `<b>Score:</b> ${{formatNumber(p.total_weighted_score, 4)}}<br>` +
            `<b>Weights:</b> D=${{p.distance_weight}}, P=${{p.population_weight}}, T=${{p.traffic_weight}}, A=${{p.airspace_weight}}`
          );
        }}
      }});
      layer.addTo(clusterLayers[clusterId]);
    }});

    L.control.layers(null, layerControlItems, {{ collapsed: false }}).addTo(map);

    const bounds = L.geoJSON(routeData).getBounds();
    if (bounds.isValid()) {{
      map.fitBounds(bounds.pad(0.12));
    }}

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "panel");
      const rows = Object.keys(clusterRoutes).sort().map((clusterId) => {{
        const color = clusterColors[clusterId] || "#dc2626";
        const routes = clusterRoutes[clusterId].map((name) => name.replaceAll("_", " ")).join(", ");
        return `
          <div class="legend-row">
            <span class="legend-line" style="border-top-color:${{color}};"></span>
            <span class="legend-label">
              <strong>${{clusterId}}</strong>
              ${{routes}}
            </span>
          </div>
        `;
      }}).join("");
      div.innerHTML = `
        <h1>Pilot Route Clusters</h1>
        <div>${{rows}}</div>
      `;
      return div;
    }};
    legend.addTo(map);
  </script>
</body>
</html>
"""


def load_clustered_routes():
    """Join pilot route geometries with cluster labels."""
    if not os.path.exists(PILOT_ROUTES_GEOJSON):
        raise FileNotFoundError(
            f"Pilot route GeoJSON not found: {PILOT_ROUTES_GEOJSON}"
        )
    if not os.path.exists(PILOT_ROUTE_CLUSTERS_CSV):
        raise FileNotFoundError(
            f"Pilot route clusters not found: {PILOT_ROUTE_CLUSTERS_CSV}"
        )

    routes = gpd.read_file(PILOT_ROUTES_GEOJSON)
    clusters = pd.read_csv(PILOT_ROUTE_CLUSTERS_CSV)

    if "path_node_ids" in routes.columns:
        routes = routes.drop(columns=["path_node_ids"])

    cluster_columns = [
        "pilot_name",
        "cluster_id",
        "cluster_size",
        "similarity_threshold",
    ]
    routes = routes.drop(
        columns=[
            column
            for column in ["cluster_id", "cluster_size", "similarity_threshold"]
            if column in routes.columns
        ]
    )
    return routes.merge(
        clusters[cluster_columns],
        on="pilot_name",
        how="left",
    )


def build_cluster_route_lookup(routes):
    """Create a cluster-to-route-name dictionary for the legend."""
    lookup = {}
    for row in routes[["cluster_id", "pilot_name"]].itertuples(index=False):
        lookup.setdefault(row.cluster_id, []).append(row.pilot_name)
    return lookup


def main():
    """Write the clustered pilot route map HTML."""
    routes = load_clustered_routes()
    route_data = json.loads(routes.to_json())
    cluster_routes = build_cluster_route_lookup(routes)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    html = HTML_TEMPLATE.format(
        route_data=json.dumps(route_data),
        cluster_colors=json.dumps(CLUSTER_COLORS),
        cluster_routes=json.dumps(cluster_routes),
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

    print(f"Saved pilot cluster map: {OUTPUT_HTML}")


if __name__ == "__main__":
    main()
