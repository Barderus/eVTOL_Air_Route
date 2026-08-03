"""Create standalone Leaflet density maps from API grid outputs."""

import json
import os

import geopandas as gpd
import pandas as pd

from api_tool.select_study_area import select_study_area


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{title}</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
    integrity="sha256-p4NxAoJBhIINfQPDQo4O4FwBdxuV6pS2j3u2Bf8e6a4="
    crossorigin=""
  />
  <style>
    :root {{
      color-scheme: light;
      --surface: #ffffff;
      --surface-soft: #f5f7f8;
      --line: #cfd7dd;
      --text: #172026;
      --muted: #5f6f7a;
      --active: #165a72;
      --low: #2ca25f;
      --medium-air: #f1c40f;
      --high-air: #b30000;
      --medium-pop: #6baed6;
      --high-pop: #08306b;
      --traffic-low: #84cc16;
      --traffic-medium: #f97316;
      --traffic-high: #7f1d1d;
    }}

    * {{
      box-sizing: border-box;
    }}

    html, body {{
      height: 100%;
      margin: 0;
      overflow: hidden;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--text);
      background: var(--surface-soft);
    }}

    body {{
      display: flex;
      flex-direction: column;
    }}

    .toolbar {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      padding: 10px 14px;
      border-bottom: 1px solid var(--line);
      background: var(--surface);
      z-index: 1000;
    }}

    h1 {{
      margin: 0;
      font-size: 17px;
      font-weight: 700;
      line-height: 1.2;
    }}

    .status {{
      margin-top: 2px;
      color: var(--muted);
      font-size: 12px;
      line-height: 1.25;
      white-space: nowrap;
    }}

    #map {{
      flex: 1 1 auto;
      width: 100%;
      height: calc(100vh - 59px);
      min-height: 360px;
    }}

    .leaflet-container {{
      background: #edf1f2;
    }}

    .legend-box {{
      padding: 10px 12px;
      border: 1px solid rgba(23, 32, 38, 0.15);
      border-radius: 7px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.14);
      color: var(--text);
      font-size: 12px;
      line-height: 1.35;
    }}

    .legend-title {{
      margin-bottom: 7px;
      font-weight: 700;
    }}

    .legend-row {{
      display: grid;
      grid-template-columns: 14px auto;
      align-items: center;
      gap: 7px;
      margin: 5px 0;
      white-space: nowrap;
    }}

    .swatch {{
      width: 14px;
      height: 14px;
      display: inline-block;
      border: 1px solid rgba(23, 32, 38, 0.2);
    }}

    .popup-grid {{
      display: grid;
      grid-template-columns: auto auto;
      gap: 3px 12px;
      min-width: 190px;
    }}

    .popup-grid dt {{
      color: var(--muted);
      font-weight: 600;
    }}

    .popup-grid dd {{
      margin: 0;
      text-align: right;
    }}
  </style>
</head>
<body>
  <header class="toolbar">
    <div>
      <h1>{title}</h1>
      <div id="status" class="status">Loading density grid...</div>
    </div>
    <div class="status">Field: {value_field} | Cells: {cell_count}</div>
  </header>
  <div id="map"></div>

  <script
    src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"
    integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo="
    crossorigin=""
  ></script>
  <script>
    const gridData = {grid_data};
    const valueField = {value_field_json};
    const minValue = {minimum_json};
    const maxValue = {maximum_json};
    const studyBounds = {study_bounds_json};
    const statusEl = document.getElementById("status");

    const map = L.map("map", {{
      preferCanvas: true,
      zoomControl: true
    }});
    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);

    function numberOrZero(value) {{
      const numberValue = Number(value);
      return Number.isFinite(numberValue) ? numberValue : 0;
    }}

    function colorForValue(value) {{
      const numericValue = numberOrZero(value);
      if (!maxValue || maxValue <= minValue) return "#edf1f2";
      if (valueField === "airport_risk_combined") {{
        if (numericValue > 70) return "#b30000";
        if (numericValue > 30) return "#f1c40f";
        return "#2ca25f";
      }}
      if (valueField === "density_p_km2" || valueField === "city_risk") {{
        if (numericValue > 90) return "#08306b";
        if (numericValue > 40) return "#6baed6";
        return "#2ca25f";
      }}
      const ratio = Math.max(0, Math.min(1, (numericValue - minValue) / (maxValue - minValue)));
      if (ratio >= 0.8) return "#7f1d1d";
      if (ratio >= 0.6) return "#dc2626";
      if (ratio >= 0.4) return "#f97316";
      if (ratio >= 0.2) return "#facc15";
      if (ratio > 0) return "#84cc16";
      return "#e5e7eb";
    }}

    const gridLayer = L.geoJSON(gridData, {{
      renderer: L.canvas({{ padding: 0.35 }}),
      style: (feature) => {{
        const value = (feature.properties || {{}})[valueField];
        const color = colorForValue(value);
        return {{
          color,
          weight: 0.35,
          opacity: 0.45,
          fillColor: color,
          fillOpacity: 0.58
        }};
      }},
      onEachFeature: (feature, layer) => {{
        const p = feature.properties || {{}};
        layer.bindPopup(
          `<dl class="popup-grid">` +
          `<dt>${{valueField}}</dt><dd>${{numberOrZero(p[valueField]).toFixed(4)}}</dd>` +
          `<dt>Population</dt><dd>${{p.pop_class || "n/a"}}</dd>` +
          `<dt>Airspace</dt><dd>${{p.air_class || "n/a"}}</dd>` +
          `<dt>Pop risk</dt><dd>${{numberOrZero(p.city_risk).toFixed(2)}}</dd>` +
          `<dt>Air risk</dt><dd>${{numberOrZero(p.airport_risk_combined).toFixed(2)}}</dd>` +
          `<dt>Traffic count</dt><dd>${{numberOrZero(p.traffic_count).toFixed(0)}}</dd>` +
          `<dt>Traffic risk</dt><dd>${{numberOrZero(p.traffic_risk).toFixed(4)}}</dd>` +
          `</dl>`
        );
      }}
    }}).addTo(map);

    const bounds = gridLayer.getBounds();
    if (bounds.isValid()) {{
      map.fitBounds(bounds.pad(0.04));
      if (studyBounds) {{
        map.setMaxBounds(L.latLngBounds(studyBounds).pad(0.04));
      }}
      statusEl.textContent = `${{gridData.features ? gridData.features.length.toLocaleString() : 0}} cells loaded`;
    }} else {{
      map.setView([39.5, -98.35], 4);
      statusEl.textContent = "Density grid has no valid bounds";
    }}

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function () {{
      const div = L.DomUtil.create("div", "legend-box");
      div.innerHTML = `
        <div class="legend-title">Density Layer</div>
        <div class="legend-row"><span class="swatch" style="background:#2ca25f"></span><span>Low</span></div>
        <div class="legend-row"><span class="swatch" style="background:#6baed6"></span><span>Population medium</span></div>
        <div class="legend-row"><span class="swatch" style="background:#08306b"></span><span>Population high</span></div>
        <div class="legend-row"><span class="swatch" style="background:#f1c40f"></span><span>Airspace medium</span></div>
        <div class="legend-row"><span class="swatch" style="background:#b30000"></span><span>Airspace high</span></div>
        <div class="legend-row"><span class="swatch" style="background:#f97316"></span><span>Traffic medium</span></div>
        <div class="legend-row"><span class="swatch" style="background:#7f1d1d"></span><span>Traffic high</span></div>
      `;
      return div;
    }};
    legend.addTo(map);

    L.control.layers(null, {{
      "Density grid": gridLayer
    }}).addTo(map);

    map.whenReady(() => {{
      requestAnimationFrame(() => map.invalidateSize());
    }});
  </script>
</body>
</html>
"""


def summarize_field(grid, value_field):
    """Return display-ready summary values for one numeric grid field."""
    if value_field not in grid.columns:
        raise ValueError(f"Grid is missing density field: {value_field}")

    values = pd.to_numeric(grid[value_field], errors="coerce").fillna(0.0)
    return {
        "minimum": float(values.min()) if len(values) else 0.0,
        "maximum": float(values.max()) if len(values) else 0.0,
        "cell_count": int(len(grid)),
    }


def save_density_map_html(html_text, output_html_path):
    """Save density map HTML."""
    output_folder = os.path.dirname(output_html_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    with open(output_html_path, "w", encoding="utf-8") as output_file:
        output_file.write(html_text)


def make_density_map(
    city_name=None,
    grid=None,
    grid_path=None,
    value_field=None,
    output_html_path=None,
    title=None,
):
    """Create a standalone Leaflet HTML density map from a grid GeoJSON."""
    if city_name is None:
        city_name = input("City or study area name: ").strip()
        grid_path = input("Grid or scored-grid path: ").strip()
        value_field = input("Value field to map: ").strip()
        output_html_path = input("Map HTML output path: ").strip()
        title = input("Map title, blank for default: ").strip() or None

    if not value_field:
        raise ValueError("value_field is required.")
    if not output_html_path:
        raise ValueError("output_html_path is required.")

    if grid is None:
        if not grid_path or not os.path.exists(grid_path):
            raise FileNotFoundError(f"Grid not found: {grid_path}")
        grid = gpd.read_file(grid_path)

    grid = grid.to_crs("EPSG:4326")
    summary = summarize_field(grid, value_field)
    grid_data = json.loads(grid.to_json(default=str))
    west, south, east, north = grid.total_bounds.tolist()
    study_bounds = [[south, west], [north, east]]
    resolved_title = title or f"{city_name} {value_field} Density"

    html_text = HTML_TEMPLATE.format(
        title=resolved_title,
        value_field=value_field,
        value_field_json=json.dumps(value_field),
        minimum=round(summary["minimum"], 4),
        maximum=round(summary["maximum"], 4),
        minimum_json=json.dumps(summary["minimum"]),
        maximum_json=json.dumps(summary["maximum"]),
        study_bounds_json=json.dumps(study_bounds),
        cell_count=summary["cell_count"],
        grid_data=json.dumps(grid_data),
    )
    save_density_map_html(html_text, output_html_path)

    return {
        "study_area": select_study_area(city_name=city_name),
        "grid_path": grid_path,
        "value_field": value_field,
        "output_html_path": output_html_path,
        "output_exists": os.path.exists(output_html_path),
        "cell_count": summary["cell_count"],
        "minimum": summary["minimum"],
        "maximum": summary["maximum"],
    }


if __name__ == "__main__":
    result = make_density_map()
    print(result)
