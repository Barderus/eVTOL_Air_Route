"""Build HTML maps for direct-route Hierarchical + Jaccard representatives."""

import json
from pathlib import Path


ROUTE_ROBUSTNESS = Path("St Louis") / "route_robustness"
DIRECT_ROUTES_FOLDER = ROUTE_ROBUSTNESS / "route_clusters" / "direct_routes"
METHOD_SUFFIX = "hierarchical_jaccard"
ROUTE_PAIRS = [
    "midamerica_to_st_louis_lambert",
    "midamerica_to_st_louis_union_station",
    "st_louis_downtown_airport_to_st_louis_lambert",
]
COLORS = [
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
    "#e6ab02", "#1f78b4", "#a6761d", "#6a3d9a", "#33a02c",
]


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <link rel="stylesheet" href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin="">
  <style>
    html, body {{ height: 100%; margin: 0; color: #172026; font-family: Arial, Helvetica, sans-serif; }}
    #map {{ width: 100%; height: 100%; }}
    .panel {{ background: #fff; padding: 10px 12px; border-radius: 7px; box-shadow: 0 6px 18px rgba(23,32,38,.16); font-size: 13px; line-height: 1.35; max-width: 390px; max-height: 52vh; overflow: auto; }}
    .panel h1 {{ margin: 0 0 6px; font-size: 16px; line-height: 1.2; }}
    .row {{ display: flex; align-items: flex-start; gap: 8px; margin: 5px 0; }}
    .line {{ flex: 0 0 auto; width: 30px; margin-top: 7px; border-top: 4px solid #000; }}
    .label {{ min-width: 0; }}
    .label strong {{ display: block; }}
  </style>
</head>
<body>
  <main id="map" aria-label="Direct route representative map"></main>
  <script src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""></script>
  <script>
    const routeData = {route_data};
    const colors = {colors};
    const map = L.map("map", {{ preferCanvas: true }});
    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19, attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);
    const clusters = [...new Set(routeData.features.map((feature) => feature.properties.hierarchical_jaccard_cluster))];
    const largest = Math.max(...routeData.features.map((feature) => Number(feature.properties.cluster_size)), 1);
    const layers = {{}};
    function number(value, digits) {{
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toFixed(digits) : "n/a";
    }}
    routeData.features.forEach((feature, index) => {{
      const p = feature.properties;
      const cluster = p.hierarchical_jaccard_cluster;
      const layer = L.geoJSON(feature, {{
        style: {{ color: colors[index % colors.length], weight: 3 + 5 * Math.sqrt(Number(p.cluster_size) / largest), opacity: .9 }},
        onEachFeature: (_, routeLayer) => routeLayer.bindPopup(
          `<b>${{cluster}}</b><br>` +
          `<b>Representative:</b> ${{p.weight_id}}<br>` +
          `<b>Cluster size:</b> ${{p.cluster_size}} routes<br>` +
          `<b>Weight share:</b> ${{number(p.cluster_weight_space_percent, 1)}}%<br>` +
          `<b>Mean Jaccard similarity:</b> ${{number(p.representative_mean_jaccard_similarity, 3)}}<br>` +
          `<b>Distance:</b> ${{number(p.route_distance_km, 2)}} km<br>` +
          `<b>Score:</b> ${{number(p.total_weighted_score, 3)}}<br>` +
          `<b>Weights:</b> D=${{p.distance_weight}}, P=${{p.population_weight}}, T=${{p.traffic_weight}}, A=${{p.airspace_weight}}`
        )
      }}).addTo(map);
      layers[`${{cluster}} (${{p.cluster_size}} routes)`] = layer;
    }});
    L.control.layers(null, layers, {{ collapsed: false }}).addTo(map);
    const bounds = L.geoJSON(routeData).getBounds();
    if (bounds.isValid()) map.fitBounds(bounds.pad(.12));
    const origin = routeData.features[0].properties;
    L.marker([origin.origin_lat, origin.origin_lon]).bindPopup(`<b>Origin</b><br>${{origin.origin_label}}`).addTo(map);
    L.marker([origin.destination_lat, origin.destination_lon]).bindPopup(`<b>Destination</b><br>${{origin.destination_label}}`).addTo(map);
    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = () => {{
      const div = L.DomUtil.create("div", "panel");
      div.innerHTML = `<h1>Hierarchical + Jaccard Representatives</h1>` + routeData.features.map((feature, index) => {{
        const p = feature.properties;
        return `<div class="row"><span class="line" style="border-top-color:${{colors[index % colors.length]}}"></span><span class="label"><strong>${{p.hierarchical_jaccard_cluster}} (${{p.cluster_size}} routes)</strong>representative ${{p.weight_id}}, mean similarity ${{number(p.representative_mean_jaccard_similarity, 3)}}</span></div>`;
      }}).join("");
      return div;
    }};
    legend.addTo(map);
  </script>
</body>
</html>
"""


def main():
    """Write one representative map for each OD pair."""
    for route_pair in ROUTE_PAIRS:
        path = DIRECT_ROUTES_FOLDER / f"{route_pair}_{METHOD_SUFFIX}_representatives.geojson"
        with path.open("r", encoding="utf-8") as file_handle:
            route_data = json.load(file_handle)
        title = f"St. Louis Direct Route Representatives - {route_data['features'][0]['properties']['route_pair_label']}"
        html = HTML_TEMPLATE.format(
            title=title,
            route_data=json.dumps(route_data),
            colors=json.dumps(COLORS),
        )
        output_path = DIRECT_ROUTES_FOLDER / f"{route_pair}_hierarchical_jaccard_representatives.html"
        output_path.write_text(html, encoding="utf-8")
        print(f"Saved representative map: {output_path}")


if __name__ == "__main__":
    main()
