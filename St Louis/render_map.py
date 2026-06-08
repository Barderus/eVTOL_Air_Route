"""Render the interactive St. Louis population-risk map."""

import json
import os

import settings


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>St. Louis eVTOL Risk Map</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
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
      --medium-combined: #f4b183;
      --high-combined: #c84c5a;
      --medium-air: #f1c40f;
      --high-air: #b30000;
      --medium-pop: #6baed6;
      --high-pop: #08306b;
      --no-data: #9ca3af;
    }}

    * {{ box-sizing: border-box; }}

    html,
    body {{
      height: 100%;
      margin: 0;
      color: var(--text);
      font-family: Arial, Helvetica, sans-serif;
      overflow: hidden;
    }}

    body {{
      display: flex;
      flex-direction: column;
      background: var(--surface-soft);
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
      line-height: 1.2;
    }}

    .status {{
      margin-top: 2px;
      color: var(--muted);
      font-size: 12px;
    }}

    .segmented {{
      display: inline-grid;
      grid-template-columns: repeat(3, minmax(108px, 1fr));
      min-height: 38px;
      border: 1px solid var(--line);
      border-radius: 7px;
      overflow: hidden;
      background: var(--surface);
    }}

    .segmented button {{
      padding: 8px 12px;
      border: 0;
      border-right: 1px solid var(--line);
      color: var(--text);
      background: transparent;
      font: inherit;
      font-size: 13px;
    }}

    .segmented button:last-child {{ border-right: 0; }}
    .segmented button:hover {{ background: #eef3f5; cursor: pointer; }}
    .segmented button.active {{ color: #ffffff; background: var(--active); }}

    #map {{
      flex: 1 1 auto;
      width: 100%;
      height: calc(100vh - 59px);
      min-height: 360px;
    }}

    .legend-box {{
      padding: 10px 12px;
      border: 1px solid rgba(23, 32, 38, 0.15);
      border-radius: 7px;
      background: rgba(255, 255, 255, 0.96);
      box-shadow: 0 6px 18px rgba(23, 32, 38, 0.14);
      font-size: 12px;
      line-height: 1.35;
    }}

    .legend-title {{ margin-bottom: 7px; font-weight: 700; }}
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
      border: 1px solid rgba(23, 32, 38, 0.28);
    }}

    .popup-grid {{
      display: grid;
      grid-template-columns: auto auto;
      gap: 3px 12px;
      min-width: 190px;
    }}

    .popup-grid dt {{ color: var(--muted); font-weight: 600; }}
    .popup-grid dd {{ margin: 0; text-align: right; }}

    .ring-label {{
      padding: 3px 7px;
      border: 1px solid rgba(23, 32, 38, 0.16);
      border-radius: 999px;
      background: rgba(255, 255, 255, 0.92);
      box-shadow: 0 4px 10px rgba(23, 32, 38, 0.12);
      color: var(--text);
      font-size: 11px;
      font-weight: 700;
      line-height: 1;
      white-space: nowrap;
    }}

    @media (max-width: 720px) {{
      .toolbar {{ align-items: stretch; flex-direction: column; gap: 9px; }}
      .segmented {{ width: 100%; grid-template-columns: repeat(3, 1fr); }}
      #map {{ height: calc(100vh - 105px); }}
    }}
  </style>
</head>
<body>
  <header class="toolbar">
    <div>
      <h1>St. Louis eVTOL Risk Map</h1>
      <div id="status" class="status">Loading risk grid...</div>
    </div>
    <div class="segmented" role="group" aria-label="Risk layers">
      <button id="btnCombined" class="active" type="button" aria-pressed="true">Combined</button>
      <button id="btnAir" type="button" aria-pressed="false">Airspace</button>
      <button id="btnPop" type="button" aria-pressed="false">Population</button>
    </div>
  </header>

  <main id="map" aria-label="St. Louis eVTOL population-risk map"></main>

  <script
    src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js"
    integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew=="
    crossorigin=""
  ></script>
  <script>
    const GRID_URL = "st_louis_risk_grid.geojson";
    const STUDY_BOUNDS = L.latLngBounds(
      [{south}, {west}],
      [{north}, {east}]
    );
    const locations = {locations};
    const labels = {labels};
    const airspaceSites = {airspace_sites};
    const CELL_SIZE_M = {cell_size_m};
    const NM_TO_M = 1852;
    const colors = {{
      low: "#2ca25f",
      mediumCombined: "#f4b183",
      highCombined: "#c84c5a",
      mediumAir: "#f1c40f",
      highAir: "#b30000",
      mediumPop: "#6baed6",
      highPop: "#08306b",
      noData: "#9ca3af"
    }};

    const statusEl = document.getElementById("status");
    const map = L.map("map", {{
      center: [{center_lat}, {center_lon}],
      zoom: {zoom},
      preferCanvas: true,
      zoomControl: true
    }});

    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      maxZoom: 19,
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO"
    }}).addTo(map);

    const airportIcon = L.icon({{
      iconUrl: "../../Chicago/images/icons/airplane.svg",
      iconSize: [26, 26],
      iconAnchor: [13, 13],
      popupAnchor: [0, -12]
    }});

    const downtownIcon = L.icon({{
      iconUrl: "icons/city.svg",
      iconSize: [28, 28],
      iconAnchor: [14, 14],
      popupAnchor: [0, -12]
    }});

    const landmarkLayer = L.layerGroup(
      Object.entries(locations).map(([name, coordinates]) => {{
        const isAirport = name.endsWith("airport");
        const icon = isAirport ? airportIcon : downtownIcon;
        return L.marker(coordinates, {{ icon }}).bindPopup(`<b>${{labels[name]}}</b>`);
      }})
    ).addTo(map);

    function radiusMeters(radius, unit) {{
      if (unit !== "NM") throw new Error(`Unsupported airspace unit: ${{unit}}`);
      return radius * NM_TO_M;
    }}

    function ringLabelLatLng(latitude, longitude, radius) {{
      return [latitude + (radius / 111320), longitude];
    }}

    function radiusLabel(radius) {{
      return Number.isInteger(radius) ? radius.toFixed(0) : radius.toFixed(1);
    }}

    const airspaceRingLayer = L.layerGroup(
      airspaceSites.flatMap((site) => {{
        const [latitude, longitude] = site.location;
        const rings = [{{
          label:
            `${{site.code}} ${{radiusLabel(site.inner_radius)}} ` +
            `${{site.inner_unit}}`,
          radius: radiusMeters(site.inner_radius, site.inner_unit),
          color: "#b30000",
          dashArray: "8 6"
        }}];

        if (site.outer_radius !== null) {{
          rings.push({{
            label:
              `${{site.code}} ${{radiusLabel(site.outer_radius)}} ` +
              `${{site.outer_unit}}`,
            radius: radiusMeters(site.outer_radius, site.outer_unit),
            color: "#f1c40f",
            dashArray: "14 7"
          }});
        }}

        return rings.flatMap((ring) => {{
          const circle = L.circle([latitude, longitude], {{
            radius: ring.radius,
            color: ring.color,
            weight: 2.2,
            opacity: 0.95,
            dashArray: ring.dashArray,
            fill: false
          }}).bindPopup(
            `<b>${{site.airport}}</b><br>` +
            `${{site.airspace_type}}<br>` +
            `${{ring.label}}<br>` +
            `${{site.vertical_range}}`
          );
          const label = L.marker(
            ringLabelLatLng(latitude, longitude, ring.radius),
            {{
              interactive: false,
              icon: L.divIcon({{
                className: "",
                html: `<span class="ring-label">${{ring.label}}</span>`,
                iconSize: null
              }})
            }}
          );
          return [circle, label];
        }});
      }})
    ).addTo(map);

    L.control.layers(
      null,
      {{
        "Landmarks": landmarkLayer,
        "Airspace rings": airspaceRingLayer
      }},
      {{ collapsed: false, position: "topright" }}
    ).addTo(map);

    function fmt(value, digits = 1) {{
      const number = Number(value);
      return Number.isFinite(number) ? number.toFixed(digits) : "n/a";
    }}

    function combinedColor(feature) {{
      const riskClass = feature.properties?.risk_class || "Low";
      if (riskClass === "High") return colors.highCombined;
      if (riskClass === "Medium") return colors.mediumCombined;
      return colors.low;
    }}

    function airspaceColor(feature) {{
      const airClass = feature.properties?.air_class || "Low";
      if (airClass === "High") return colors.highAir;
      if (airClass === "Medium") return colors.mediumAir;
      return colors.low;
    }}

    function populationColor(feature) {{
      const populationClass = feature.properties?.pop_class || "Low";
      if (populationClass === "No Data") return colors.noData;
      if (populationClass === "High") return colors.highPop;
      if (populationClass === "Medium") return colors.mediumPop;
      return colors.low;
    }}

    function gridStyle(colorFunction) {{
      return function(feature) {{
        const color = colorFunction(feature);
        return {{
          color,
          weight: 0.35,
          opacity: 0.45,
          fillColor: color,
          fillOpacity: 0.58
        }};
      }};
    }}

    function popupHtml(feature) {{
      const properties = feature.properties || {{}};
      return {{
        combined: `
          <dl class="popup-grid">
            <dt>Risk class</dt><dd>${{properties.risk_class ?? "n/a"}}</dd>
            <dt>Risk cost</dt><dd>${{fmt(properties.risk_cost)}}</dd>
            <dt>Dominant</dt><dd>${{properties.density_type ?? "n/a"}}</dd>
            <dt>Population</dt><dd>${{properties.pop_class ?? "n/a"}}</dd>
            <dt>Airspace</dt><dd>${{properties.air_class ?? "n/a"}}</dd>
            <dt>Pop risk</dt><dd>${{fmt(properties.city_risk)}}</dd>
            <dt>Air risk</dt><dd>${{fmt(properties.airport_risk_combined)}}</dd>
          </dl>
        `,
        air: `
          <dl class="popup-grid">
            <dt>Airspace risk</dt><dd>${{properties.air_class ?? "n/a"}}</dd>
            <dt>Airspace cost</dt><dd>${{fmt(properties.airport_risk_combined)}}</dd>
            <dt>Source</dt><dd>${{properties.airspace_source ?? "None"}}</dd>
            <dt>Vertical range</dt><dd>${{properties.airspace_vertical_range ?? "None"}}</dd>
            <dt>Cell size</dt><dd>${{CELL_SIZE_M}} m</dd>
          </dl>
        `,
        pop: `
          <dl class="popup-grid">
            <dt>Population risk</dt><dd>${{properties.pop_class ?? "n/a"}}</dd>
            <dt>Population density</dt><dd>${{fmt(properties.density_p_km2)}} people/km²</dd>
            <dt>Population cost</dt><dd>${{fmt(properties.city_risk)}}</dd>
            <dt>Cell size</dt><dd>${{CELL_SIZE_M}} m</dd>
          </dl>
        `
      }};
    }}

    function makeGridLayer(data, mode, colorFunction) {{
      let gridLayer;
      gridLayer = L.geoJSON(data, {{
        renderer: L.canvas({{ padding: 0.35 }}),
        style: gridStyle(colorFunction),
        onEachFeature: (feature, layer) => {{
          layer.bindPopup(() => popupHtml(feature)[mode], {{ maxWidth: 320 }});
          layer.on({{
            mouseover: () => layer.setStyle({{ weight: 1.5, color: "#172026", fillOpacity: 0.72 }}),
            mouseout: () => gridLayer.resetStyle(layer)
          }});
        }}
      }});
      return gridLayer;
    }}

    const legend = L.control({{ position: "bottomright" }});
    legend.onAdd = function() {{
      const div = L.DomUtil.create("div", "legend-box");
      div.id = "legend";
      return div;
    }};
    legend.addTo(map);

    const buttons = {{
      combined: document.getElementById("btnCombined"),
      air: document.getElementById("btnAir"),
      pop: document.getElementById("btnPop")
    }};
    const legendRows = {{
      combined: [
        ["Low combined risk", colors.low],
        ["Medium combined risk", colors.mediumCombined],
        ["High combined risk", colors.highCombined]
      ],
      air: [
        ["Low airspace risk", colors.low],
        ["Medium airspace risk", colors.mediumAir],
        ["High airspace risk", colors.highAir]
      ],
      pop: [
        ["Low population risk", colors.low],
        ["Medium population risk", colors.mediumPop],
        ["High population risk", colors.highPop],
        ["Population data unavailable", colors.noData]
      ]
    }};
    const layerLabels = {{
      combined: "Combined",
      air: "Airspace",
      pop: "Population"
    }};
    let activeMode = "combined";
    let gridLayers = {{}};

    function updateLegend() {{
      const rows = legendRows[activeMode].map(([label, color]) => (
        `<div class="legend-row"><span class="swatch" style="background:${{color}}"></span><span>${{label}}</span></div>`
      )).join("");
      document.getElementById("legend").innerHTML =
        `<div class="legend-title">${{layerLabels[activeMode]}} Layer</div>${{rows}}`;
    }}

    function setMode(mode) {{
      if (!gridLayers[mode]) return;
      map.removeLayer(gridLayers[activeMode]);
      activeMode = mode;
      gridLayers[activeMode].addTo(map);
      Object.entries(buttons).forEach(([buttonMode, button]) => {{
        const isActive = buttonMode === activeMode;
        button.classList.toggle("active", isActive);
        button.setAttribute("aria-pressed", String(isActive));
      }});
      updateLegend();
    }}

    Object.entries(buttons).forEach(([mode, button]) => {{
      button.addEventListener("click", () => setMode(mode));
    }});

    fetch(GRID_URL)
      .then((response) => {{
        if (!response.ok) throw new Error(`HTTP ${{response.status}}`);
        return response.json();
      }})
      .then((data) => {{
        gridLayers = {{
          combined: makeGridLayer(data, "combined", combinedColor),
          air: makeGridLayer(data, "air", airspaceColor),
          pop: makeGridLayer(data, "pop", populationColor)
        }};
        gridLayers.combined.addTo(map);
        updateLegend();
        statusEl.textContent = `${{data.features.length.toLocaleString()}} risk grid cells loaded`;
      }})
      .catch((error) => {{
        statusEl.textContent = `Unable to load population grid: ${{error.message}}`;
      }});
  </script>
</body>
</html>
"""


LOCATION_LABELS = {
    "st_louis_downtown": "Downtown St. Louis",
    "st_louis_downtown_airport": "St. Louis Downtown Airport",
    "st_louis_lambert_airport": "St. Louis Lambert International Airport",
    "midamerica_st_louis_airport": "MidAmerica St. Louis Airport",
}


def main():
    if not os.path.exists(settings.RISK_GRID_GEOJSON):
        raise FileNotFoundError(
            f"St. Louis risk grid not found: {settings.RISK_GRID_GEOJSON}. "
            "Run make_map.py first."
        )

    html = HTML_TEMPLATE.format(
        south=settings.SOUTH,
        west=settings.WEST,
        north=settings.NORTH,
        east=settings.EAST,
        center_lat=settings.MAP_CENTER_LAT,
        center_lon=settings.MAP_CENTER_LON,
        zoom=settings.MAP_ZOOM,
        cell_size_m=settings.CELL_SIZE_M,
        locations=json.dumps(settings.LOCATIONS),
        labels=json.dumps(LOCATION_LABELS),
        airspace_sites=json.dumps(settings.AIRSPACE_SITES),
    )

    with open(settings.MAP_HTML, "w", encoding="utf-8") as output_file:
        output_file.write(html)

    print(f"Saved St. Louis map: {settings.MAP_HTML}")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as error:
        raise SystemExit(str(error)) from error
