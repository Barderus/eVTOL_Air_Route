"""Create one HTML map per direct-route clustering method."""

import json
import os

import pandas as pd


DIRECT_ROUTE_WEIGHTS_CSV = "route_robustness/output/direct_route_weight_configurations.csv"
JACCARD_CSV = "route_robustness/output/direct_route_clusters_jaccard.csv"
DBSCAN_CSV = "route_robustness/output/direct_route_clusters_dbscan.csv"
OPTICS_CSV = "route_robustness/output/direct_route_clusters_optics.csv"
HIERARCHICAL_CSV = "route_robustness/output/direct_route_clusters_hierarchical_jaccard.csv"
FRECHET_CSV = "route_robustness/output/direct_route_clusters_frechet.csv"
EDIT_DISTANCE_CSV = "route_robustness/output/direct_route_clusters_edit_distance.csv"
ALL_ROUTES_GEOJSON = "route_robustness/output/all_routes.geojson"
MAP_FOLDER = "route_robustness/maps"

MAP_CENTER_LAT = 41.79
MAP_CENTER_LON = -87.88
MAP_ZOOM = 10

ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521

METHODS = [
    {
        "label": "Jaccard Threshold",
        "file": JACCARD_CSV,
        "cluster_column": "jaccard_cluster",
        "output_html": "direct_route_clusters_jaccard.html",
        "description": "Connected components using shared path-node Jaccard distance <= 0.75.",
    },
    {
        "label": "DBSCAN On Jaccard",
        "file": DBSCAN_CSV,
        "cluster_column": "dbscan_cluster",
        "output_html": "direct_route_clusters_dbscan.html",
        "description": "DBSCAN using Jaccard distance with eps 0.75 and min_samples 4.",
    },
    {
        "label": "OPTICS On Jaccard",
        "file": OPTICS_CSV,
        "cluster_column": "optics_cluster",
        "output_html": "direct_route_clusters_optics.html",
        "description": "OPTICS-style ordering using Jaccard distance with extraction eps 0.75.",
    },
    {
        "label": "Hierarchical On Jaccard",
        "file": HIERARCHICAL_CSV,
        "cluster_column": "hierarchical_jaccard_cluster",
        "output_html": "direct_route_clusters_hierarchical_jaccard.html",
        "description": "Average-linkage hierarchical clustering using Jaccard distance <= 0.75.",
    },
    {
        "label": "Hierarchical On Frechet",
        "file": FRECHET_CSV,
        "cluster_column": "frechet_cluster",
        "output_html": "direct_route_clusters_frechet.html",
        "description": "Average-linkage hierarchical clustering using discrete Frechet distance <= 2.5 km.",
    },
    {
        "label": "Hierarchical On Edit Distance",
        "file": EDIT_DISTANCE_CSV,
        "cluster_column": "edit_distance_cluster",
        "output_html": "direct_route_clusters_edit_distance.html",
        "description": "Average-linkage hierarchical clustering using normalized route-node edit distance <= 0.35.",
    },
]

CLUSTER_COLORS = [
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
    "#e31a1c",
    "#ff7f00",
    "#6d6d6d",
]


def load_route_features():
    """Load full route GeoJSON features by route_run_id."""
    with open(ALL_ROUTES_GEOJSON, "r", encoding="utf-8") as file:
        geojson = json.load(file)

    features = {}
    for feature in geojson["features"]:
        route_id = feature["properties"]["route_run_id"]
        features[route_id] = feature
    return features


def load_method_assignments(method):
    """Load direct routes and one clustering method output."""
    routes = pd.read_csv(DIRECT_ROUTE_WEIGHTS_CSV)
    clusters = pd.read_csv(method["file"])
    assignments = routes.merge(
        clusters[["route_run_id", method["cluster_column"]]],
        on="route_run_id",
        how="left",
    )
    return assignments


def build_cluster_summary(assignments, cluster_column):
    """Summarize clusters for one method."""
    summary = (
        assignments.groupby(cluster_column)
        .agg(
            route_count=("route_run_id", "count"),
            min_distance_km=("route_distance_km", "min"),
            max_distance_km=("route_distance_km", "max"),
            mean_distance_km=("route_distance_km", "mean"),
            mean_score=("total_weighted_score", "mean"),
        )
        .sort_values("route_count", ascending=False)
        .reset_index()
    )

    cluster_colors = {}
    cluster_summary = {}
    for color_index, row in enumerate(summary.itertuples(index=False)):
        cluster_id = getattr(row, cluster_column)
        cluster_colors[cluster_id] = CLUSTER_COLORS[color_index % len(CLUSTER_COLORS)]
        cluster_summary[cluster_id] = {
            "route_count": int(row.route_count),
            "route_share": float(row.route_count / len(assignments)),
            "min_distance_km": float(row.min_distance_km),
            "max_distance_km": float(row.max_distance_km),
            "mean_distance_km": float(row.mean_distance_km),
            "mean_score": float(row.mean_score),
        }

    return cluster_colors, cluster_summary


def build_map_geojson(assignments, route_features):
    """Build GeoJSON for the direct routes."""
    assignment_lookup = assignments.set_index("route_run_id").to_dict("index")
    features = []

    for route_id in assignments["route_run_id"]:
        feature = json.loads(json.dumps(route_features[route_id]))
        properties = feature["properties"]
        for key, value in assignment_lookup[route_id].items():
            if key != "route_run_id":
                properties[key] = value
        features.append(feature)

    return {
        "type": "FeatureCollection",
        "features": features,
    }


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{map_title}</title>
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
      max-width: 410px;
      max-height: 56vh;
      overflow: auto;
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

    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);

    const routeData = {route_data};
    const clusterColumn = "{cluster_column}";
    const clusterColors = {cluster_colors};
    const clusterSummary = {cluster_summary};
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
      const clusterId = properties[clusterColumn] || "unclustered";
      const routeColor = clusterColors[clusterId] || "#666666";

      if (!clusterLayers[clusterId]) {{
        clusterLayers[clusterId] = L.layerGroup().addTo(map);
        const summary = clusterSummary[clusterId] || {{}};
        const count = summary.route_count || 0;
        layerControlItems[`${{clusterId}} (${{count}} routes)`] = clusterLayers[clusterId];
      }}

      const layer = L.geoJSON(feature, {{
        style: {{
          color: routeColor,
          weight: 3,
          opacity: 0.58
        }},
        onEachFeature: (routeFeature, routeLayer) => {{
          const p = routeFeature.properties || {{}};
          routeLayer.bindPopup(
            `<b>${{p.weight_id}}</b><br>` +
            `<b>Cluster:</b> ${{p[clusterColumn]}}<br>` +
            `<b>Source 0.25 cluster:</b> ${{p.direct_source_cluster_id}}<br>` +
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
      const rows = Object.keys(clusterSummary).map((clusterId) => {{
        const summary = clusterSummary[clusterId];
        const color = clusterColors[clusterId] || "#666666";
        return `
          <div class="legend-row">
            <span class="legend-line" style="border-top-color:${{color}};"></span>
            <span class="legend-label">
              <strong>${{clusterId}} (${{summary.route_count}} routes)</strong>
              ${{formatNumber(summary.route_share * 100, 1)}}%,
              ${{formatNumber(summary.min_distance_km, 1)}}-${{formatNumber(summary.max_distance_km, 1)}} km
            </span>
          </div>
        `;
      }}).join("");
      div.innerHTML = `
        <h1>Clusters</h1>
        <div>${{rows}}</div>
      `;
      return div;
    }};
    legend.addTo(map);
  </script>
</body>
</html>
"""


def write_method_map(method, route_features):
    """Write one method-specific cluster map."""
    assignments = load_method_assignments(method)
    cluster_colors, cluster_summary = build_cluster_summary(
        assignments,
        method["cluster_column"],
    )
    route_data = build_map_geojson(assignments, route_features)

    html = HTML_TEMPLATE.format(
        map_title=f"Direct Routes - {method['label']}",
        description=method["description"],
        map_center_lat=MAP_CENTER_LAT,
        map_center_lon=MAP_CENTER_LON,
        map_zoom=MAP_ZOOM,
        origin_lat=ORIGIN_LAT,
        origin_lon=ORIGIN_LON,
        origin_label=ORIGIN_LABEL,
        destination_lat=DESTINATION_LAT,
        destination_lon=DESTINATION_LON,
        destination_label=DESTINATION_LABEL,
        route_data=json.dumps(route_data),
        cluster_column=method["cluster_column"],
        cluster_colors=json.dumps(cluster_colors),
        cluster_summary=json.dumps(cluster_summary),
    )

    output_path = os.path.join(MAP_FOLDER, method["output_html"])
    with open(output_path, "w", encoding="utf-8") as file:
        file.write(html)

    return output_path, len(cluster_summary)


def main():
    """Create one HTML map for each direct-route clustering method."""
    os.makedirs(MAP_FOLDER, exist_ok=True)
    route_features = load_route_features()

    for method in METHODS:
        output_path, cluster_count = write_method_map(method, route_features)
        print(f"Saved {method['label']} map: {output_path}")
        print("Cluster count:", cluster_count)


if __name__ == "__main__":
    main()
