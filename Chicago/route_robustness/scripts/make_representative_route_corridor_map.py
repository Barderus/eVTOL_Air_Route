"""Create a corridor-width map for DBSCAN Frechet representative routes."""

import json
import os

import geopandas as gpd

from make_weight_barycentric_plot import CLUSTER_COLORS


REPRESENTATIVE_ROUTES_GEOJSON = "Chicago/route_robustness/output/representative_route.geojson"
OUTPUT_FOLDER = "Chicago/route_robustness/output"
OUTPUT_GEOJSON = os.path.join(OUTPUT_FOLDER, "representative_route_corridor.geojson")
OUTPUT_HTML = os.path.join(OUTPUT_FOLDER, "representative_route_corridor.html")

CORRIDOR_WIDTH_M = 4828.0
CORRIDOR_BUFFER_M = CORRIDOR_WIDTH_M / 2.0
PROJECTED_CRS = "EPSG:32616"

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
  <title>DBSCAN Frechet Representative Route Corridors</title>
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
      max-width: 420px;
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
      border-top: 7px solid #000000;
      opacity: 0.72;
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

    const corridorData = __CORRIDOR_DATA__;
    const routeData = __ROUTE_DATA__;
    const clusterColors = __CLUSTER_COLORS__;
    const corridorWidthM = __CORRIDOR_WIDTH_M__;

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

    corridorData.features.forEach((feature) => {
      const properties = feature.properties || {};
      const clusterId = properties.dbscan_cluster || "unclustered";
      const routeColor = clusterColors[clusterId] || "#666666";
      const layer = L.geoJSON(feature, {
        style: {
          color: routeColor,
          weight: 2,
          opacity: 0.72,
          fillColor: routeColor,
          fillOpacity: 0.24
        },
        onEachFeature: (routeFeature, routeLayer) => {
          const p = routeFeature.properties || {};
          routeLayer.bindPopup(
            `<b>${p.dbscan_cluster}</b><br>` +
            `<b>Representative route:</b> ${p.representative_route_run_id}<br>` +
            `<b>Corridor width:</b> ${formatNumber(p.corridor_width_m, 0)} m<br>` +
            `<b>Buffer each side:</b> ${formatNumber(p.corridor_buffer_m, 0)} m<br>` +
            `<b>Cluster size:</b> ${p.cluster_size} routes<br>` +
            `<b>Mean Frechet distance:</b> ${formatNumber(p.medoid_mean_frechet_km, 3)} km`
          );
        }
      }).addTo(map);
      layerControlItems[`${clusterId} corridor (${properties.cluster_size} routes)`] = layer;
    });

    routeData.features.forEach((feature) => {
      const properties = feature.properties || {};
      const clusterId = properties.dbscan_cluster || "unclustered";
      const routeColor = clusterColors[clusterId] || "#666666";
      L.geoJSON(feature, {
        style: {
          color: routeColor,
          weight: 4,
          opacity: 0.95
        }
      }).addTo(map);
    });

    L.control.layers(null, layerControlItems, { collapsed: false }).addTo(map);

    const bounds = L.geoJSON(corridorData).getBounds();
    if (bounds.isValid()) {
      map.fitBounds(bounds.pad(0.12));
    }

    const legend = L.control({ position: "bottomright" });
    legend.onAdd = function () {
      const div = L.DomUtil.create("div", "panel");
      const rows = corridorData.features.map((feature) => {
        const p = feature.properties || {};
        const color = clusterColors[p.dbscan_cluster] || "#666666";
        return `
          <div class="legend-row">
            <span class="legend-line" style="border-top-color:${color};"></span>
            <span class="legend-label">
              <strong>${p.dbscan_cluster} (${p.cluster_size} routes)</strong>
              ${formatNumber(corridorWidthM, 0)} m wide,
              representative ${p.weight_id}
            </span>
          </div>
        `;
      }).join("");
      div.innerHTML = `
        <h1>Representative Route Corridors</h1>
        <div>${rows}</div>
      `;
      return div;
    };
    legend.addTo(map);
  </script>
</body>
</html>
"""


def build_corridors(routes):
    """Build route corridor polygons using a meter-based projected CRS."""
    projected = routes.to_crs(PROJECTED_CRS)
    corridors = projected.copy()
    corridors["geometry"] = projected.geometry.buffer(CORRIDOR_BUFFER_M)
    corridors["corridor_width_m"] = CORRIDOR_WIDTH_M
    corridors["corridor_buffer_m"] = CORRIDOR_BUFFER_M
    return corridors.to_crs(routes.crs)


def build_cluster_colors(routes):
    """Assign colors using the same size-ordered palette as the tetrahedron."""
    cluster_ids = routes["dbscan_cluster"].tolist()
    return {
        cluster_id: CLUSTER_COLORS[index % len(CLUSTER_COLORS)]
        for index, cluster_id in enumerate(cluster_ids)
    }


def write_html(routes, corridors, cluster_colors):
    """Write the corridor map HTML."""
    html = (
        HTML_TEMPLATE.replace("__CORRIDOR_DATA__", corridors.to_json())
        .replace("__ROUTE_DATA__", routes.to_json())
        .replace("__CLUSTER_COLORS__", json.dumps(cluster_colors))
        .replace("__CORRIDOR_WIDTH_M__", str(CORRIDOR_WIDTH_M))
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
    """Save the representative route corridor GeoJSON and HTML map."""
    if not os.path.exists(REPRESENTATIVE_ROUTES_GEOJSON):
        raise FileNotFoundError(
            f"Representative route GeoJSON not found: {REPRESENTATIVE_ROUTES_GEOJSON}"
        )

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes = gpd.read_file(REPRESENTATIVE_ROUTES_GEOJSON)
    if "path_node_ids" in routes.columns:
        routes = routes.drop(columns=["path_node_ids"])
    corridors = build_corridors(routes)
    cluster_colors = build_cluster_colors(routes)

    corridors.to_file(OUTPUT_GEOJSON, driver="GeoJSON")
    write_html(routes, corridors, cluster_colors)

    print(f"Saved representative route corridor GeoJSON: {OUTPUT_GEOJSON}")
    print(f"Saved representative route corridor HTML: {OUTPUT_HTML}")
    print(
        routes[
            [
                "dbscan_cluster",
                "cluster_size",
                "representative_route_run_id",
                "weight_id",
            ]
        ].to_string(index=False)
    )
    print(f"Corridor width: {CORRIDOR_WIDTH_M:.0f} m")
    print(f"Buffer on each side: {CORRIDOR_BUFFER_M:.0f} m")


if __name__ == "__main__":
    main()
