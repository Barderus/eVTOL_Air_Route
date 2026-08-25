"""Select representative routes for DBSCAN Frechet clusters."""

import json
import os

import pandas as pd

from cluster_direct_routes_dbscan import (
    ALL_ROUTES_GEOJSON,
    build_frechet_distance_matrix,
    load_route_features,
)
from make_weight_barycentric_plot import CLUSTER_COLORS


DBSCAN_CLUSTERS_CSV = "route_robustness/output/direct_route_clusters_dbscan.csv"
OUTPUT_FOLDER = "route_robustness/output"
OUTPUT_CSV = os.path.join(OUTPUT_FOLDER, "representative_route.csv")
OUTPUT_GEOJSON = os.path.join(OUTPUT_FOLDER, "representative_route.geojson")
OUTPUT_HTML = os.path.join(OUTPUT_FOLDER, "representative_route.html")

MAP_CENTER_LAT = 41.79
MAP_CENTER_LON = -87.88
MAP_ZOOM = 10

ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>DBSCAN Frechet Representative Routes</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin=""
  />
  <style>
    html,
    body {
      height: 100%;
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: #172026;
    }

    #map {
      width: 100%;
      height: 100%;
    }

    .panel {
      background: #ffffff;
      padding: 10px 12px;
      border-radius: 7px;
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.16);
      font-size: 13px;
      line-height: 1.35;
      max-width: 410px;
      max-height: 54vh;
      overflow: auto;
    }

    .panel h1 {
      margin: 0 0 6px;
      font-size: 16px;
      line-height: 1.2;
    }

    .legend-row {
      display: flex;
      align-items: flex-start;
      gap: 8px;
      margin: 6px 0;
    }

    .legend-line {
      display: inline-block;
      flex: 0 0 auto;
      width: 30px;
      margin-top: 7px;
      border-top: 5px solid #000000;
    }

    .legend-label {
      min-width: 0;
    }

    .legend-label strong {
      display: block;
    }
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
    const map = L.map("map", {
      center: [__MAP_CENTER_LAT__, __MAP_CENTER_LON__],
      zoom: __MAP_ZOOM__
    });

    L.tileLayer("https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png", {
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }).addTo(map);

    const routeData = __ROUTE_DATA__;
    const clusterColors = __CLUSTER_COLORS__;

    L.marker([__ORIGIN_LAT__, __ORIGIN_LON__])
      .bindPopup("<b>Origin</b><br>__ORIGIN_LABEL__")
      .addTo(map);

    L.marker([__DESTINATION_LAT__, __DESTINATION_LON__])
      .bindPopup("<b>Destination</b><br>__DESTINATION_LABEL__")
      .addTo(map);

    function formatNumber(value, digits) {
      const numberValue = Number(value);
      if (Number.isNaN(numberValue)) {
        return "n/a";
      }
      return numberValue.toFixed(digits);
    }

    const layerControlItems = {};

    routeData.features.forEach((feature) => {
      const properties = feature.properties || {};
      const clusterId = properties.dbscan_cluster || "unclustered";
      const routeColor = clusterColors[clusterId] || "#666666";
      const layer = L.geoJSON(feature, {
        style: {
          color: routeColor,
          weight: 7,
          opacity: 0.94
        },
        onEachFeature: (routeFeature, routeLayer) => {
          const p = routeFeature.properties || {};
          routeLayer.bindPopup(
            `<b>${p.dbscan_cluster}</b><br>` +
            `<b>Representative route:</b> ${p.representative_route_run_id}<br>` +
            `<b>Weight ID:</b> ${p.weight_id}<br>` +
            `<b>Cluster size:</b> ${p.cluster_size} routes<br>` +
            `<b>Weight space:</b> ${formatNumber(p.cluster_weight_space_percent, 1)}%<br>` +
            `<b>Total Frechet distance:</b> ${formatNumber(p.medoid_total_frechet_km, 3)} km<br>` +
            `<b>Mean Frechet distance:</b> ${formatNumber(p.medoid_mean_frechet_km, 3)} km<br>` +
            `<b>Route distance:</b> ${formatNumber(p.route_distance_km, 2)} km<br>` +
            `<b>Score:</b> ${formatNumber(p.total_weighted_score, 4)}<br>` +
            `<b>Weights:</b> D=${p.distance_weight}, P=${p.population_weight}, T=${p.traffic_weight}, A=${p.airspace_weight}`
          );
        }
      }).addTo(map);
      layerControlItems[`${clusterId} (${properties.cluster_size} routes)`] = layer;
    });

    L.control.layers(null, layerControlItems, { collapsed: false }).addTo(map);

    const bounds = L.geoJSON(routeData).getBounds();
    if (bounds.isValid()) {
      map.fitBounds(bounds.pad(0.12));
    }

    const legend = L.control({ position: "bottomright" });
    legend.onAdd = function () {
      const div = L.DomUtil.create("div", "panel");
      const rows = routeData.features.map((feature) => {
        const p = feature.properties || {};
        const color = clusterColors[p.dbscan_cluster] || "#666666";
        return `
          <div class="legend-row">
            <span class="legend-line" style="border-top-color:${color};"></span>
            <span class="legend-label">
              <strong>${p.dbscan_cluster} (${p.cluster_size} routes)</strong>
              representative ${p.weight_id},
              mean Frechet ${formatNumber(p.medoid_mean_frechet_km, 3)} km
            </span>
          </div>
        `;
      }).join("");
      div.innerHTML = `
        <h1>DBSCAN Frechet Representatives</h1>
        <div>${rows}</div>
      `;
      return div;
    };
    legend.addTo(map);
  </script>
</body>
</html>
"""


def load_clustered_routes():
    """Load DBSCAN route assignments."""
    if not os.path.exists(DBSCAN_CLUSTERS_CSV):
        raise FileNotFoundError(f"Required input file not found: {DBSCAN_CLUSTERS_CSV}")

    return pd.read_csv(DBSCAN_CLUSTERS_CSV)


def is_noise_cluster(cluster_id):
    """Return true when a DBSCAN label is a noise label."""
    return "noise" in str(cluster_id).lower()


def select_representatives(routes, distance_matrix):
    """Select the Frechet medoid route for each DBSCAN cluster."""
    representatives = []
    total_routes = len(routes)

    for cluster_id, cluster_routes in routes.groupby("dbscan_cluster", sort=True):
        if is_noise_cluster(cluster_id):
            continue

        cluster_indices = cluster_routes.index.tolist()
        best_route_index = None
        best_total_distance = None

        for route_index in cluster_indices:
            total_distance = sum(
                distance_matrix[route_index][other_index]
                for other_index in cluster_indices
            )
            if best_total_distance is None or total_distance < best_total_distance:
                best_total_distance = total_distance
                best_route_index = route_index

        representative = routes.loc[best_route_index].copy()
        cluster_size = len(cluster_indices)
        representative["representative_route_run_id"] = representative["route_run_id"]
        representative["cluster_size"] = cluster_size
        representative["cluster_weight_space_percent"] = 100.0 * cluster_size / total_routes
        representative["medoid_total_frechet_km"] = float(best_total_distance)
        if cluster_size > 1:
            representative["medoid_mean_frechet_km"] = float(
                best_total_distance / (cluster_size - 1)
            )
        else:
            representative["medoid_mean_frechet_km"] = 0.0
        representatives.append(representative)

    output = pd.DataFrame(representatives)
    output = output.sort_values(
        by=["cluster_size", "dbscan_cluster"],
        ascending=[False, True],
    )
    return output


def build_geojson(representatives, route_features):
    """Build a GeoJSON feature collection for representative routes."""
    features = []

    for row in representatives.to_dict("records"):
        route_id = row["route_run_id"]
        feature = json.loads(json.dumps(route_features[route_id]))
        properties = feature["properties"]
        properties.update(row)
        features.append(feature)

    return {
        "type": "FeatureCollection",
        "features": features,
    }


def build_cluster_colors(representatives):
    """Assign colors using the same size-ordered palette as the tetrahedron."""
    cluster_ids = representatives["dbscan_cluster"].tolist()
    return {
        cluster_id: CLUSTER_COLORS[index % len(CLUSTER_COLORS)]
        for index, cluster_id in enumerate(cluster_ids)
    }


def write_html(route_data, cluster_colors):
    """Write the representative route HTML map."""
    html = (
        HTML_TEMPLATE.replace("__ROUTE_DATA__", json.dumps(route_data))
        .replace("__CLUSTER_COLORS__", json.dumps(cluster_colors))
        .replace("__MAP_CENTER_LAT__", str(MAP_CENTER_LAT))
        .replace("__MAP_CENTER_LON__", str(MAP_CENTER_LON))
        .replace("__MAP_ZOOM__", str(MAP_ZOOM))
        .replace("__ORIGIN_LAT__", str(ORIGIN_LAT))
        .replace("__ORIGIN_LON__", str(ORIGIN_LON))
        .replace("__ORIGIN_LABEL__", ORIGIN_LABEL)
        .replace("__DESTINATION_LAT__", str(DESTINATION_LAT))
        .replace("__DESTINATION_LON__", str(DESTINATION_LON))
        .replace("__DESTINATION_LABEL__", DESTINATION_LABEL)
    )

    with open(OUTPUT_HTML, "w", encoding="utf-8") as file:
        file.write(html)


def main():
    """Save DBSCAN Frechet representative routes."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    routes = load_clustered_routes()
    route_features = load_route_features()
    distance_matrix = build_frechet_distance_matrix(routes, route_features)
    representatives = select_representatives(routes, distance_matrix)
    route_data = build_geojson(representatives, route_features)
    cluster_colors = build_cluster_colors(representatives)

    representatives.to_csv(OUTPUT_CSV, index=False)
    with open(OUTPUT_GEOJSON, "w", encoding="utf-8") as file:
        json.dump(route_data, file)
    write_html(route_data, cluster_colors)

    print(f"Saved representative route CSV: {OUTPUT_CSV}")
    print(f"Saved representative route GeoJSON: {OUTPUT_GEOJSON}")
    print(f"Saved representative route HTML: {OUTPUT_HTML}")
    print(
        representatives[
            [
                "dbscan_cluster",
                "cluster_size",
                "representative_route_run_id",
                "weight_id",
                "medoid_mean_frechet_km",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
