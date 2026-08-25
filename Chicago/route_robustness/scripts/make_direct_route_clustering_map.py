"""Create an HTML map for comparing direct-route clustering methods."""

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
OUTPUT_HTML = "route_robustness/maps/direct_route_clustering.html"

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
        "key": "jaccard_cluster",
        "label": "Jaccard threshold",
        "file": JACCARD_CSV,
        "cluster_column": "jaccard_cluster",
        "description": "Connected components using shared path-node Jaccard distance <= 0.75.",
    },
    {
        "key": "dbscan_cluster",
        "label": "DBSCAN on Frechet",
        "file": DBSCAN_CSV,
        "cluster_column": "dbscan_cluster",
        "description": "DBSCAN using discrete Frechet distance with eps 2.5 km and min_samples 4.",
    },
    {
        "key": "optics_cluster",
        "label": "OPTICS on Frechet",
        "file": OPTICS_CSV,
        "cluster_column": "optics_cluster",
        "description": "OPTICS-style ordering using discrete Frechet distance with extraction eps 2.5 km.",
    },
    {
        "key": "hierarchical_jaccard_cluster",
        "label": "Hierarchical on Jaccard",
        "file": HIERARCHICAL_CSV,
        "cluster_column": "hierarchical_jaccard_cluster",
        "description": "Average-linkage hierarchical clustering using Jaccard distance <= 0.75.",
    },
    {
        "key": "frechet_cluster",
        "label": "Hierarchical on Frechet",
        "file": FRECHET_CSV,
        "cluster_column": "frechet_cluster",
        "description": "Average-linkage hierarchical clustering using discrete Frechet distance <= 2.5 km.",
    },
    {
        "key": "edit_distance_cluster",
        "label": "Hierarchical on edit distance",
        "file": EDIT_DISTANCE_CSV,
        "cluster_column": "edit_distance_cluster",
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


def load_assignments():
    """Load direct routes and merge each clustering output."""
    assignments = pd.read_csv(DIRECT_ROUTE_WEIGHTS_CSV)

    for method in METHODS:
        method_table = pd.read_csv(method["file"])
        assignments = assignments.merge(
            method_table[["route_run_id", method["cluster_column"]]],
            on="route_run_id",
            how="left",
        )
        assignments = assignments.rename(
            columns={method["cluster_column"]: method["key"]}
        )

    return assignments


def load_route_features():
    """Load full route GeoJSON features by route_run_id."""
    with open(ALL_ROUTES_GEOJSON, "r", encoding="utf-8") as file:
        geojson = json.load(file)

    features = {}
    for feature in geojson["features"]:
        route_id = feature["properties"]["route_run_id"]
        features[route_id] = feature
    return features


def build_map_geojson(assignments, route_features):
    """Build GeoJSON with direct-route cluster assignments attached."""
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


def build_method_summary(assignments):
    """Build cluster summary data used by the HTML legend."""
    method_summary = {}

    for method in METHODS:
        cluster_column = method["key"]
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
        clusters = {}
        for color_index, row in enumerate(summary.itertuples(index=False)):
            cluster_id = getattr(row, cluster_column)
            cluster_colors[cluster_id] = CLUSTER_COLORS[color_index % len(CLUSTER_COLORS)]
            clusters[cluster_id] = {
                "route_count": int(row.route_count),
                "route_share": float(row.route_count / len(assignments)),
                "min_distance_km": float(row.min_distance_km),
                "max_distance_km": float(row.max_distance_km),
                "mean_distance_km": float(row.mean_distance_km),
                "mean_score": float(row.mean_score),
            }

        method_summary[cluster_column] = {
            "label": method["label"],
            "description": method["description"],
            "cluster_colors": cluster_colors,
            "clusters": clusters,
        }

    return method_summary


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Direct Route Clustering</title>
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

    .method-option {{
      display: block;
      margin: 6px 0;
      cursor: pointer;
    }}

    .method-option input {{
      margin-right: 6px;
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
    const methodSummary = {method_summary};
    const methods = {methods};
    const routeLayers = [];
    let selectedMethod = methods[0].key;
    let legendControl = null;

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

    function routeStyle(feature) {{
      const properties = feature.properties || {{}};
      const clusterId = properties[selectedMethod] || "unclustered";
      const colorMap = methodSummary[selectedMethod].cluster_colors || {{}};
      const color = colorMap[clusterId] || "#666666";
      return {{
        color: color,
        weight: 3,
        opacity: 0.58
      }};
    }}

    function popupText(properties) {{
      const clusterId = properties[selectedMethod] || "unclustered";
      return (
        `<b>${{properties.weight_id}}</b><br>` +
        `<b>Method:</b> ${{methodSummary[selectedMethod].label}}<br>` +
        `<b>Cluster:</b> ${{clusterId}}<br>` +
        `<b>Source 0.25 cluster:</b> ${{properties.direct_source_cluster_id}}<br>` +
        `<b>Distance:</b> ${{formatNumber(properties.route_distance_km, 2)}} km<br>` +
        `<b>Score:</b> ${{formatNumber(properties.total_weighted_score, 4)}}<br>` +
        `<b>Weights:</b> D=${{properties.distance_weight}}, P=${{properties.population_weight}}, T=${{properties.traffic_weight}}, A=${{properties.airspace_weight}}`
      );
    }}

    routeData.features.forEach((feature) => {{
      const layer = L.geoJSON(feature, {{
        style: routeStyle,
        onEachFeature: (routeFeature, routeLayer) => {{
          routeLayer.bindPopup(popupText(routeFeature.properties || {{}}));
        }}
      }}).addTo(map);
      routeLayers.push(layer);
    }});

    const bounds = L.geoJSON(routeData).getBounds();
    if (bounds.isValid()) {{
      map.fitBounds(bounds.pad(0.12));
    }}

    function refreshRouteStyles() {{
      routeLayers.forEach((layer) => {{
        layer.eachLayer((subLayer) => {{
          subLayer.setStyle(routeStyle(subLayer.feature));
          subLayer.bindPopup(popupText(subLayer.feature.properties || {{}}));
        }});
      }});
    }}

    const methodControl = L.control({{ position: "topright" }});
    methodControl.onAdd = function () {{
      const div = L.DomUtil.create("div", "panel");
      const rows = methods.map((method, index) => {{
        const checked = index === 0 ? "checked" : "";
        return `
          <label class="method-option">
            <input type="radio" name="method" value="${{method.key}}" ${{checked}}>
            ${{method.label}}
          </label>
        `;
      }}).join("");
      div.innerHTML = `
        <h1>Direct Route Clustering</h1>
        ${{rows}}
      `;
      L.DomEvent.disableClickPropagation(div);
      return div;
    }};
    methodControl.addTo(map);

    function drawLegend() {{
      if (legendControl) {{
        map.removeControl(legendControl);
      }}

      legendControl = L.control({{ position: "bottomright" }});
      legendControl.onAdd = function () {{
        const div = L.DomUtil.create("div", "panel");
        const selectedSummary = methodSummary[selectedMethod];
        const clusters = selectedSummary.clusters || {{}};
        const rows = Object.keys(clusters).map((clusterId) => {{
          const summary = clusters[clusterId];
          const color = selectedSummary.cluster_colors[clusterId] || "#666666";
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
          <h1>${{selectedSummary.label}}</h1>
          <div>${{rows}}</div>
        `;
        return div;
      }};
      legendControl.addTo(map);
    }}

    document.addEventListener("change", (event) => {{
      if (event.target && event.target.name === "method") {{
        selectedMethod = event.target.value;
        refreshRouteStyles();
        drawLegend();
      }}
    }});

    drawLegend();
  </script>
</body>
</html>
"""


def write_html(route_data, method_summary):
    """Write the direct route clustering HTML file."""
    os.makedirs(MAP_FOLDER, exist_ok=True)
    methods_for_js = [
        {
            "key": method["key"],
            "label": method["label"],
        }
        for method in METHODS
    ]

    html = HTML_TEMPLATE.format(
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
        method_summary=json.dumps(method_summary),
        methods=json.dumps(methods_for_js),
    )

    with open(OUTPUT_HTML, "w", encoding="utf-8") as file:
        file.write(html)


def main():
    """Create the direct route clustering map."""
    assignments = load_assignments()
    route_features = load_route_features()
    route_data = build_map_geojson(assignments, route_features)
    method_summary = build_method_summary(assignments)
    write_html(route_data, method_summary)

    print(f"Saved direct route clustering map: {OUTPUT_HTML}")
    print("Direct route count:", len(assignments))
    for method in METHODS:
        cluster_count = assignments[method["key"]].nunique()
        print(f"{method['label']}: {cluster_count} clusters")


if __name__ == "__main__":
    main()
