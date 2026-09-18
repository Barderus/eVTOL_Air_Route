"""Build final all-date direct-route representative maps for St. Louis."""

import json
from pathlib import Path

import geopandas as gpd


ROUTE_ROBUSTNESS = Path(__file__).resolve().parents[1]
OUTPUT_FOLDER = ROUTE_ROBUSTNESS / "final_routes"
METHOD_SUFFIX = "hierarchical_jaccard"
PROJECTED_CRS = "EPSG:32615"
CORRIDOR_WIDTH_M = 4828.0
CORRIDOR_BUFFER_M = CORRIDOR_WIDTH_M / 2.0

DATE_FOLDERS = [
    ("July_12", "July 12th - Saturday", "2025-07-12"),
    ("July_14", "July 14th - Monday", "2025-07-14"),
    ("January_10", "January 10th - Saturday", "2026-01-10"),
    ("January_12", "January 12th - Monday", "2026-01-12"),
    ("March_7", "March 7th - Saturday", "2026-03-07"),
    ("March_9", "March 9th - Monday", "2026-03-09"),
]

ROUTE_PAIRS = [
    "midamerica_to_st_louis_lambert",
    "midamerica_to_st_louis_union_station",
    "st_louis_downtown_airport_to_st_louis_lambert",
]

REPRESENTATIVE_COLORS = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#1f78b4",
    "#a6761d",
    "#6a3d9a",
    "#33a02c",
    "#b15928",
    "#fb9a99",
]


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin=""
  >
  <style>
    :root {{
      --surface: #ffffff;
      --line: #d1d5db;
      --text: #172026;
      --muted: #64748b;
      --active: #165a72;
    }}
    * {{ box-sizing: border-box; }}
    html, body {{
      height: 100%;
      margin: 0;
      color: var(--text);
      font-family: Arial, Helvetica, sans-serif;
      overflow: hidden;
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      min-height: 62px;
      padding: 10px 14px;
      border-bottom: 1px solid var(--line);
      background: var(--surface);
      z-index: 1000;
    }}
    h1 {{
      margin: 0;
      font-size: 17px;
      line-height: 1.2;
    }}
    .status {{
      margin-top: 3px;
      color: var(--muted);
      font-size: 12px;
    }}
    .date-toggle {{
      display: grid;
      grid-template-columns: repeat(6, minmax(96px, 1fr));
      min-height: 38px;
      border: 1px solid var(--line);
      border-radius: 7px;
      overflow: hidden;
      background: var(--surface);
    }}
    .date-toggle button {{
      padding: 8px 9px;
      border: 0;
      border-right: 1px solid var(--line);
      color: var(--text);
      background: transparent;
      font: inherit;
      font-size: 12px;
      cursor: pointer;
      white-space: nowrap;
    }}
    .date-toggle button:last-child {{ border-right: 0; }}
    .date-toggle button:hover {{ background: #eef3f5; }}
    .date-toggle button.active {{
      color: #ffffff;
      background: var(--active);
    }}
    #map {{
      width: 100%;
      height: calc(100vh - 62px);
      min-height: 420px;
    }}
    .panel {{
      width: 310px;
      max-height: calc(100vh - 96px);
      overflow: auto;
      padding: 11px 12px;
      border-radius: 7px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.16);
      font-size: 12px;
      line-height: 1.35;
    }}
    .panel strong {{
      display: block;
      margin-bottom: 7px;
      font-size: 13px;
    }}
    .date-block {{
      margin: 9px 0 11px;
      padding-top: 8px;
      border-top: 1px solid rgba(209, 213, 219, 0.8);
    }}
    .date-heading {{
      display: flex;
      align-items: center;
      gap: 7px;
      margin-bottom: 5px;
      font-weight: 700;
    }}
    .swatch {{
      width: 12px;
      height: 12px;
      border-radius: 999px;
      border: 1px solid rgba(17, 24, 39, 0.25);
    }}
    .rep-row {{
      appearance: none;
      width: 100%;
      border: 0;
      background: transparent;
      display: grid;
      grid-template-columns: 16px 1fr auto;
      gap: 7px;
      align-items: center;
      padding: 4px 0;
      color: inherit;
      cursor: pointer;
      font: inherit;
      text-align: left;
    }}
    .rep-row:hover {{ background: rgba(229, 231, 235, 0.55); }}
    .rep-row.is-hidden {{ opacity: 0.42; }}
    .line {{
      width: 16px;
      border-top: 4px solid #000;
    }}
    @media (max-width: 1050px) {{
      .toolbar {{ align-items: stretch; flex-direction: column; }}
      .date-toggle {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
      .date-toggle button {{ white-space: normal; }}
      #map {{ height: calc(100vh - 120px); }}
    }}
  </style>
</head>
<body>
  <header class="toolbar">
    <div>
      <h1>{title}</h1>
      <div id="status" class="status">Select one traffic date at a time.</div>
    </div>
    <div class="date-toggle" id="dateToggle" role="group" aria-label="Traffic date filter"></div>
  </header>

  <main id="map" aria-label="Final St. Louis direct route representatives map"></main>

  <script
    src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""
  ></script>
  <script>
    const representativeData = {representative_data};
    const corridorData = {corridor_data};
    const dateOrder = {date_order};

    const map = L.map("map", {{ preferCanvas: true }});
    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);

    const routeLayer = L.layerGroup().addTo(map);
    const endpointLayer = L.layerGroup().addTo(map);
    const panelControl = L.control({{ position: "topright" }});
    const dateToggle = document.getElementById("dateToggle");
    const statusEl = document.getElementById("status");
    let panelDiv = null;
    let activeDate = dateOrder[0].key;
    let hiddenRepresentatives = new Set();
    let representativeLayers = new Map();

    function number(value, digits) {{
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed.toFixed(digits) : "n/a";
    }}

    function popupHtml(properties) {{
      return `
        <b>${{properties.date_label}}</b><br>
        <b>${{properties.merged_cluster_id}}</b><br>
        <b>Representative:</b> ${{properties.weight_id}}<br>
        <b>Cluster size:</b> ${{properties.cluster_size}} routes<br>
        <b>Weight share:</b> ${{number(properties.cluster_weight_space_percent, 1)}}%<br>
        <b>Mean Frechet:</b> ${{number(properties.representative_mean_frechet_miles, 3)}} miles<br>
        <b>Distance:</b> ${{number(properties.route_distance_km, 2)}} km<br>
        <b>Weights:</b>
        D=${{properties.distance_weight}},
        P=${{properties.population_weight}},
        T=${{properties.traffic_weight}},
        A=${{properties.airspace_weight}}
      `;
    }}

    function corridorPopupHtml(properties) {{
      return `
        <b>${{properties.date_label}}</b><br>
        <b>${{properties.merged_cluster_id}} corridor</b><br>
        <b>Width:</b> ${{number(properties.corridor_width_m, 0)}} m<br>
        <b>Buffer each side:</b> ${{number(properties.corridor_buffer_m, 0)}} m
      `;
    }}

    function featureKey(properties) {{
      return `${{properties.date_key}}::${{properties.merged_cluster_id}}`;
    }}

    function drawEndpoints() {{
      endpointLayer.clearLayers();
      const feature = representativeData.features[0];
      if (!feature) return;
      const p = feature.properties;
      [
        {{ label: "Origin", name: p.origin_label, point: [p.origin_lat, p.origin_lon] }},
        {{ label: "Destination", name: p.destination_label, point: [p.destination_lat, p.destination_lon] }}
      ].forEach((endpoint) => {{
        L.circleMarker(endpoint.point, {{
          radius: 7,
          color: "#172026",
          weight: 2,
          fillColor: endpoint.label === "Origin" ? "#ffffff" : "#facc15",
          fillOpacity: 1
        }}).bindPopup(`<b>${{endpoint.label}}</b><br>${{endpoint.name}}`).addTo(endpointLayer);
      }});
    }}

    function drawRoutes() {{
      routeLayer.clearLayers();
      representativeLayers = new Map();

      representativeData.features.forEach((feature, index) => {{
        const p = feature.properties;
        const key = featureKey(p);
        if (p.date_key !== activeDate || hiddenRepresentatives.has(key)) {{
          return;
        }}

        const corridorFeature = corridorData.features[index];
        const corridorLayer = L.geoJSON(corridorFeature, {{
          style: {{
            color: p.representative_color,
            weight: 1,
            opacity: 0.5,
            fillColor: p.representative_color,
            fillOpacity: 0.14
          }},
          onEachFeature: (_, layer) => layer.bindPopup(corridorPopupHtml(p))
        }});
        const routeLine = L.geoJSON(feature, {{
          style: {{
            color: p.representative_color,
            weight: 4.5,
            opacity: 0.92,
            dashArray: null
          }},
          onEachFeature: (_, layer) => layer.bindPopup(popupHtml(p))
        }});
        const group = L.layerGroup([corridorLayer, routeLine]).addTo(routeLayer);
        representativeLayers.set(key, group);
      }});

      renderPanel();
      updateStatus();
    }}

    function renderPanel() {{
      if (!panelDiv) return;
      const grouped = new Map(dateOrder.map((date) => [date.key, []]));
      representativeData.features.forEach((feature) => {{
        const p = feature.properties;
        if (grouped.has(p.date_key)) grouped.get(p.date_key).push(p);
      }});

      panelDiv.innerHTML = `
        <strong>Representatives</strong>
        ${{dateOrder.map((date) => {{
          const rows = grouped.get(date.key) || [];
          const hiddenDate = date.key !== activeDate;
          return `
            <div class="date-block" style="${{hiddenDate ? "opacity:.42" : ""}}">
              <div class="date-heading">
                <span class="swatch" style="background:${{date.key === activeDate ? "#165a72" : "#cbd5e1"}}"></span>
                <span>${{date.label}} (${{date.date}})</span>
              </div>
              ${{rows.map((p) => {{
                const key = featureKey(p);
                const isHidden = hiddenRepresentatives.has(key);
                return `
                  <button class="rep-row ${{isHidden ? "is-hidden" : ""}}" type="button" data-key="${{key}}">
                    <span class="line" style="border-top-color:${{p.representative_color}}"></span>
                    <span>${{p.merged_cluster_id}} · ${{p.weight_id}}</span>
                    <span>${{p.cluster_size}}</span>
                  </button>
                `;
              }}).join("")}}
            </div>
          `;
        }}).join("")}}
      `;

      panelDiv.querySelectorAll(".rep-row").forEach((row) => {{
        row.addEventListener("click", () => {{
          const key = row.dataset.key;
          if (hiddenRepresentatives.has(key)) {{
            hiddenRepresentatives.delete(key);
          }} else {{
            hiddenRepresentatives.add(key);
          }}
          drawRoutes();
        }});
      }});
    }}

    function updateStatus() {{
      const visible = representativeData.features.filter((feature) => {{
        const p = feature.properties;
        return p.date_key === activeDate && !hiddenRepresentatives.has(featureKey(p));
      }}).length;
      const active = dateOrder.find((date) => date.key === activeDate);
      statusEl.textContent = `${{visible}} representatives visible for ${{active.label}} (${{active.date}}).`;
    }}

    function setActiveDate(dateKey) {{
      activeDate = dateKey;
      dateToggle.querySelectorAll("button").forEach((button) => {{
        button.classList.toggle("active", button.dataset.dateKey === activeDate);
        button.setAttribute("aria-pressed", button.dataset.dateKey === activeDate ? "true" : "false");
      }});
      drawRoutes();
    }}

    dateOrder.forEach((date) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.dateKey = date.key;
      button.className = date.key === activeDate ? "active" : "";
      button.textContent = date.shortLabel;
      button.setAttribute("aria-pressed", date.key === activeDate ? "true" : "false");
      button.addEventListener("click", () => setActiveDate(date.key));
      dateToggle.appendChild(button);
    }});

    panelControl.onAdd = function() {{
      panelDiv = L.DomUtil.create("div", "panel");
      L.DomEvent.disableClickPropagation(panelDiv);
      L.DomEvent.disableScrollPropagation(panelDiv);
      return panelDiv;
    }};
    panelControl.addTo(map);

    drawEndpoints();
    drawRoutes();

    const bounds = L.geoJSON(corridorData).getBounds();
    if (bounds.isValid()) map.fitBounds(bounds.pad(0.12));
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


def load_route_pair_features(route_pair):
    """Load all dated representative features for one route pair."""
    features = []
    for folder, date_label, traffic_date in DATE_FOLDERS:
        path = (
            ROUTE_ROBUSTNESS
            / folder
            / "output"
            / "direct_routes"
            / f"{route_pair}_{METHOD_SUFFIX}_representatives.geojson"
        )
        with path.open("r", encoding="utf-8") as file_handle:
            geojson = json.load(file_handle)

        for feature_index, feature in enumerate(geojson["features"]):
            properties = feature.setdefault("properties", {})
            properties["date_key"] = folder
            properties["date_label"] = date_label
            properties["traffic_date"] = traffic_date
            properties["representative_color"] = REPRESENTATIVE_COLORS[
                feature_index % len(REPRESENTATIVE_COLORS)
            ]
            properties["display_cluster_id"] = (
                f"{date_label} - {properties.get('merged_cluster_id', '')}"
            )
            features.append(feature)
    return {"type": "FeatureCollection", "features": features}


def write_route_pair_map(route_pair):
    """Write the final all-date representative map for one route pair."""
    route_data = load_route_pair_features(route_pair)
    routes = gpd.GeoDataFrame.from_features(route_data["features"], crs="EPSG:4326")
    corridor_data = json.loads(build_corridors(routes).to_json())
    route_label = route_data["features"][0]["properties"]["route_pair_label"]
    title = f"St. Louis Direct Route Representatives - {route_label}"
    date_order = [
        {
            "key": folder,
            "label": label,
            "shortLabel": folder.replace("_", " "),
            "date": traffic_date,
        }
        for folder, label, traffic_date in DATE_FOLDERS
    ]

    html = HTML_TEMPLATE.format(
        title=title,
        representative_data=json.dumps(route_data),
        corridor_data=json.dumps(corridor_data),
        date_order=json.dumps(date_order),
    )
    output_path = OUTPUT_FOLDER / f"{route_pair}_all_days_hierarchical_jaccard_representatives.html"
    output_path.write_text(html, encoding="utf-8")
    print(f"Saved final representative map: {output_path}")


def main():
    """Build one final all-date representative map for each route pair."""
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
    for route_pair in ROUTE_PAIRS:
        write_route_pair_map(route_pair)


if __name__ == "__main__":
    main()
