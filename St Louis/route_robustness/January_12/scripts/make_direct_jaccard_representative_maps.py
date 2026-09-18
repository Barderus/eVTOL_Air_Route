"""Build HTML maps for direct-route Hierarchical + Jaccard representatives."""

import json
from pathlib import Path

import geopandas as gpd

DATE_ROOT = Path(__file__).resolve().parents[1]
DIRECT_ROUTES_FOLDER = DATE_ROOT / "output" / "direct_routes"
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
CORRIDOR_WIDTH_M = 4828.0
CORRIDOR_BUFFER_M = CORRIDOR_WIDTH_M / 2.0
PROJECTED_CRS = "EPSG:32615"


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
    const corridorData = {corridor_data};
    const colors = {colors};
    const map = L.map("map", {{ preferCanvas: true }});
    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19, attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);
    const clusters = [...new Set(routeData.features.map((feature) => feature.properties.merged_cluster_id))];
    const largest = Math.max(...routeData.features.map((feature) => Number(feature.properties.cluster_size)), 1);
    const layers = {{}};
    function number(value, digits) {{
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toFixed(digits) : "n/a";
    }}
    routeData.features.forEach((feature, index) => {{
      const p = feature.properties;
      const cluster = p.merged_cluster_id;
      const layer = L.geoJSON(feature, {{
        style: {{ color: colors[index % colors.length], weight: 3 + 5 * Math.sqrt(Number(p.cluster_size) / largest), opacity: .9 }},
        onEachFeature: (_, routeLayer) => routeLayer.bindPopup(
          `<b>${{cluster}}</b><br>` +
          `<b>Representative:</b> ${{p.weight_id}}<br>` +
          `<b>Cluster size:</b> ${{p.cluster_size}} routes<br>` +
          `<b>Weight share:</b> ${{number(p.cluster_weight_space_percent, 1)}}%<br>` +
          `<b>Mean Frechet distance:</b> ${{number(p.representative_mean_frechet_miles, 3)}} miles<br>` +
          `<b>Distance:</b> ${{number(p.route_distance_km, 2)}} km<br>` +
          `<b>Score:</b> ${{number(p.total_weighted_score, 3)}}<br>` +
          `<b>Weights:</b> D=${{p.distance_weight}}, P=${{p.population_weight}}, T=${{p.traffic_weight}}, A=${{p.airspace_weight}}`
        )
      }});
      const corridorFeature = corridorData.features[index];
      const corridorLayer = L.geoJSON(corridorFeature, {{
        style: {{ color: colors[index % colors.length], weight: 1, opacity: .55, fillColor: colors[index % colors.length], fillOpacity: .18 }},
        onEachFeature: (_, corridorLayer) => corridorLayer.bindPopup(
          `<b>${{cluster}} corridor</b><br>` +
          `<b>Width:</b> ${{number(p.corridor_width_m, 0)}} m (3 miles)<br>` +
          `<b>Buffer each side:</b> ${{number(p.corridor_buffer_m, 0)}} m<br>` +
          `<b>Cluster size:</b> ${{p.cluster_size}} routes`
        )
      }});
      const routeGroup = L.layerGroup([corridorLayer, layer]).addTo(map);
      layers[`${{cluster}} (${{p.cluster_size}} routes)`] = routeGroup;
    }});
    L.control.layers(null, layers, {{ collapsed: false }}).addTo(map);
    const bounds = L.geoJSON(corridorData).getBounds();
    if (bounds.isValid()) map.fitBounds(bounds.pad(.12));
    const origin = routeData.features[0].properties;
    L.marker([origin.origin_lat, origin.origin_lon]).bindPopup(`<b>Origin</b><br>${{origin.origin_label}}`).addTo(map);
    L.marker([origin.destination_lat, origin.destination_lon]).bindPopup(`<b>Destination</b><br>${{origin.destination_label}}`).addTo(map);
    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = () => {{
      const div = L.DomUtil.create("div", "panel");
      const sample = routeData.features[0].properties;
      div.innerHTML = `<h1>Hierarchical + Jaccard Representatives<br><span style="font-size:12px;font-weight:400;">${{sample.traffic_dataset}} (${{sample.traffic_date}})</span></h1>` + routeData.features.map((feature, index) => {{
        const p = feature.properties;
        return `<div class="row"><span class="line" style="border-top-color:${{colors[index % colors.length]}}"></span><span class="label"><strong>${{p.merged_cluster_id}} (${{p.cluster_size}} routes)</strong>representative ${{p.weight_id}}, mean Frechet ${{number(p.representative_mean_frechet_miles, 3)}} miles</span></div>`;
      }}).join("");
      return div;
    }};
    legend.addTo(map);
  </script>
</body>
</html>
"""


def build_corridors(routes):
    """Build one independent 3-mile corridor around each representative route."""
    projected = routes.to_crs(PROJECTED_CRS)
    corridors = routes.copy()
    corridors["geometry"] = projected.geometry.buffer(CORRIDOR_BUFFER_M)
    corridors = corridors.to_crs("EPSG:4326")
    corridors["corridor_width_m"] = CORRIDOR_WIDTH_M
    corridors["corridor_buffer_m"] = CORRIDOR_BUFFER_M
    return corridors


def main():
    """Write one representative map for each OD pair."""
    for route_pair in ROUTE_PAIRS:
        path = DIRECT_ROUTES_FOLDER / f"{route_pair}_{METHOD_SUFFIX}_representatives.geojson"
        with path.open("r", encoding="utf-8") as file_handle:
            route_data = json.load(file_handle)
        routes = gpd.GeoDataFrame.from_features(route_data["features"], crs="EPSG:4326")
        corridors = build_corridors(routes)
        corridor_data = json.loads(corridors.to_json())
        properties = route_data["features"][0]["properties"]
        date_label = f"{properties['traffic_dataset']} ({properties['traffic_date']})"
        title = (
            "St. Louis Direct Route Representatives - "
            f"{properties['route_pair_label']} - {date_label}"
        )
        html = HTML_TEMPLATE.format(
            title=title,
            route_data=json.dumps(route_data),
            corridor_data=json.dumps(corridor_data),
            colors=json.dumps(COLORS),
        )
        corridor_path = DIRECT_ROUTES_FOLDER / f"{route_pair}_{METHOD_SUFFIX}_representatives_corridor.geojson"
        corridor_path.write_text(json.dumps(corridor_data), encoding="utf-8")
        output_path = DATE_ROOT / "maps" / "representatives" / f"{route_pair}_hierarchical_jaccard_representatives.html"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(html, encoding="utf-8")
        print(f"Saved representative map: {output_path}")


if __name__ == "__main__":
    main()
