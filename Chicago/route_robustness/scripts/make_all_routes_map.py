"""Create an HTML map for the full 287-route robustness run."""

import json
import os

import geopandas as gpd


ALL_ROUTES_GEOJSON = "route_robustness/output/all_routes.geojson"
OUTPUT_FOLDER = "route_robustness/maps"
OUTPUT_HTML = os.path.join(OUTPUT_FOLDER, "all_routes.html")

ORIGIN_LABEL = "Clow International Airport"
ORIGIN_LAT = 41.695923717435235
ORIGIN_LON = -88.12876224517822
DESTINATION_LABEL = "Chicago Union Station"
DESTINATION_LAT = 41.87838051825937
DESTINATION_LON = -87.63905525207521

MAP_CENTER_LAT = 41.79
MAP_CENTER_LON = -87.88
MAP_ZOOM = 10

CATEGORY_COLORS = {
    "distance_dominant": "#1b9e77",
    "population_dominant": "#d95f02",
    "traffic_dominant": "#7570b3",
    "airspace_dominant": "#e7298a",
    "equal_weight": "#000000",
    "tied_dominant": "#666666",
}


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Route Robustness All Routes</title>
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
      max-width: 340px;
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
      margin: 5px 0;
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
    const categoryColors = {category_colors};
    const categoryCounts = {category_counts};
    const categoryLayers = {{}};
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
      const category = properties.dominant_weight_category || "tied_dominant";
      const routeColor = categoryColors[category] || "#666666";
      if (!categoryLayers[category]) {{
        categoryLayers[category] = L.layerGroup().addTo(map);
        const label = `${{category.replaceAll("_", " ")}} (${{categoryCounts[category] || 0}})`;
        layerControlItems[label] = categoryLayers[category];
      }}

      const layer = L.geoJSON(feature, {{
        style: {{
          color: routeColor,
          weight: properties.is_equal_weight === true || properties.is_equal_weight === "True" ? 7 : 4,
          opacity: properties.is_equal_weight === true || properties.is_equal_weight === "True" ? 1.0 : 0.28
        }},
        onEachFeature: (routeFeature, routeLayer) => {{
          const p = routeFeature.properties || {{}};
          routeLayer.bindPopup(
            `<b>${{p.weight_id}}</b><br>` +
            `<b>Category:</b> ${{p.dominant_weight_category}}<br>` +
            `<b>Distance:</b> ${{formatNumber(p.route_distance_km, 2)}} km<br>` +
            `<b>Score:</b> ${{formatNumber(p.total_weighted_score, 4)}}<br>` +
            `<b>Path nodes:</b> ${{p.path_nodes}}<br>` +
            `<b>Weights:</b> D=${{p.distance_weight}}, P=${{p.population_weight}}, T=${{p.traffic_weight}}, A=${{p.airspace_weight}}`
          );
        }}
      }});
      layer.addTo(categoryLayers[category]);
    }});

    L.control.layers(null, layerControlItems, {{ collapsed: false }}).addTo(map);

    const bounds = L.geoJSON(routeData).getBounds();
    if (bounds.isValid()) {{
      map.fitBounds(bounds.pad(0.12));
    }}

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "panel");
      const rows = Object.keys(categoryCounts).sort().map((category) => {{
        const color = categoryColors[category] || "#666666";
        const label = category.replaceAll("_", " ");
        return `
          <div class="legend-row">
            <span class="legend-line" style="border-top-color:${{color}};"></span>
            <span>${{label}} (${{categoryCounts[category]}})</span>
          </div>
        `;
      }}).join("");
      div.innerHTML = `
        <h1>All Weighted Routes</h1>
        <div>${{rows}}</div>
      `;
      return div;
    }};
    legend.addTo(map);
  </script>
</body>
</html>
"""


def dominant_weight_category(row):
    """Classify a route by the largest weight value."""
    if bool(row.get("is_equal_weight", False)):
        return "equal_weight"

    weights = {
        "distance_dominant": float(row["distance_weight"]),
        "population_dominant": float(row["population_weight"]),
        "traffic_dominant": float(row["traffic_weight"]),
        "airspace_dominant": float(row["airspace_weight"]),
    }
    maximum_weight = max(weights.values())
    matching_categories = [
        category
        for category, value in weights.items()
        if value == maximum_weight
    ]
    if len(matching_categories) != 1:
        return "tied_dominant"
    return matching_categories[0]


def main():
    """Write the full route map HTML."""
    if not os.path.exists(ALL_ROUTES_GEOJSON):
        raise FileNotFoundError(f"All routes GeoJSON not found: {ALL_ROUTES_GEOJSON}")

    routes = gpd.read_file(ALL_ROUTES_GEOJSON)
    if "path_node_ids" in routes.columns:
        routes = routes.drop(columns=["path_node_ids"])

    routes["dominant_weight_category"] = routes.apply(
        dominant_weight_category,
        axis=1,
    )
    category_counts = (
        routes["dominant_weight_category"]
        .value_counts()
        .sort_index()
        .to_dict()
    )
    route_data = json.loads(routes.to_json())

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    html = HTML_TEMPLATE.format(
        route_data=json.dumps(route_data),
        category_colors=json.dumps(CATEGORY_COLORS),
        category_counts=json.dumps(category_counts),
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

    print(f"Saved all-routes map: {OUTPUT_HTML}")
    print("Route counts by dominant category:")
    for category, count in category_counts.items():
        print(f"  {category}: {count}")


if __name__ == "__main__":
    main()
