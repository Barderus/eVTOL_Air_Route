import json
import os

import numpy as np
import pandas as pd

import settings

TIME_WINDOWS = [1, 3, 6, 9, 12, 24]
EARTH_RADIUS_M = 6378137.0
HEATMAP_MIN_CELL_COUNT = 2
METERS_TO_FEET = 3.28084
GROUND_ELEVATION_FT_MSL = 466.0
TRACK_STRIDE = 3
ARROW_STRIDE = 18
OUTPUT_PATH = settings.TRAFFIC_MAP_HTML

LOCATION_LABELS = {
    "st_louis_downtown": "Downtown St. Louis",
    "st_louis_downtown_airport": "St. Louis Downtown Airport",
    "st_louis_lambert_airport": "St. Louis Lambert International Airport",
    "midamerica_st_louis_airport": "MidAmerica St. Louis Airport",
}


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>St. Louis Flight Density</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.6.0/dist/leaflet.css"
    integrity="sha512-xwE/Az9zrjBIphAcBb3F6JVqxf46+CDLwfLMHloNu6KEQCAWi6HcDUbeOfBIptF7tcCzusKFjFw2yuvEpDL9wQ=="
    crossorigin=""
  />
  <style>
    html, body {{
      height: 100%;
      margin: 0;
      font-family: Arial, sans-serif;
      background: #0d1117;
    }}

    #map {{
      width: 100vw;
      height: 100vh;
    }}

    .panel {{
      position: absolute;
      top: 12px;
      left: 12px;
      z-index: 1000;
      width: 320px;
      padding: 14px 16px;
      border-radius: 10px;
      background: rgba(255, 255, 255, 0.95);
      box-shadow: 0 10px 24px rgba(0, 0, 0, 0.18);
      line-height: 1.35;
    }}

    .panel h1 {{
      margin: 0 0 8px;
      font-size: 18px;
    }}

    .panel p {{
      margin: 0 0 8px;
      font-size: 13px;
      color: #334155;
    }}

    .metric {{
      display: inline-block;
      margin: 2px 6px 0 0;
      padding: 4px 8px;
      border-radius: 999px;
      background: #e2e8f0;
      font-size: 12px;
      color: #0f172a;
    }}

    .toggle-label {{
      margin-top: 10px;
      color: #334155;
      font-size: 12px;
      font-weight: 700;
    }}

    .date-toggle,
    .time-toggle {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 6px;
    }}

    .date-toggle button,
    .time-toggle button {{
      border: 0;
      border-radius: 999px;
      padding: 6px 10px;
      background: #cbd5e1;
      color: #0f172a;
      font-size: 12px;
      cursor: pointer;
    }}

    .date-toggle button.active,
    .time-toggle button.active {{
      background: #0f172a;
      color: #f8fafc;
    }}
  </style>
</head>
<body>
  <div id="map"></div>
  <div class="panel">
    <h1>St. Louis Flight Density</h1>
    <div class="metric" id="datasetMetric"></div>
    <div class="metric" id="rowMetric"></div>
    <div class="metric" id="flightMetric"></div>
    <div class="metric">Track stride: {track_stride}</div>
    <div class="metric">Arrow stride: {arrow_stride}</div>
    <div class="toggle-label">Traffic date</div>
    <div class="date-toggle" id="dateToggle"></div>
    <div class="toggle-label">Hours from local midnight</div>
    <div class="time-toggle" id="timeToggle"></div>
  </div>

  <script src="https://unpkg.com/leaflet@1.6.0/dist/leaflet.js" integrity="sha512-gZwIG9x3wUXg2hdXF6+rVkLF/0Vi9U8D2Ntg4Ga5I5BZpVkVxlJWbSQtXPSiUTtC0TjtGOmxa1AJPuV0CPthew==" crossorigin=""></script>
  <script src="https://unpkg.com/leaflet.heat/dist/leaflet-heat.js"></script>
  <script src="https://unpkg.com/leaflet-polylinedecorator@1.6.0/dist/leaflet.polylineDecorator.js"></script>
  <script>
    const trafficDatasets = {traffic_datasets};
    const bounds = {bounds};
    const timeButtons = {time_buttons};
    const locations = {locations};
    const locationLabels = {location_labels};

    const STUDY_BOUNDS = L.latLngBounds(
      [{south}, {west}],
      [{north}, {east}]
    );
    const map = L.map("map", {{
      center: [{center_lat}, {center_lon}],
      zoom: {zoom},
      preferCanvas: true,
      zoomSnap: 0.25
    }});

    L.tileLayer("https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png", {{
      attribution: "&copy; OpenStreetMap contributors &copy; CARTO",
      maxZoom: 19
    }}).addTo(map);

    const heatLayer = L.heatLayer([], {{
      radius: 18,
      blur: 14,
      maxZoom: 12,
      minOpacity: 0.35,
      gradient: {{
        0.2: "#1d4ed8",
        0.4: "#06b6d4",
        0.6: "#facc15",
        0.8: "#f97316",
        1.0: "#dc2626"
      }}
    }});

    const trackLayer = L.layerGroup();
    const originLayer = L.layerGroup();
    const dateToggle = document.getElementById("dateToggle");
    const timeToggle = document.getElementById("timeToggle");
    const datasetMetric = document.getElementById("datasetMetric");
    const rowMetric = document.getElementById("rowMetric");
    const flightMetric = document.getElementById("flightMetric");
    let activeDate = Object.keys(trafficDatasets)[0];
    let activeHours = 1;

    function pickTrackColor(meanAltitudeAglFt) {{
      if (meanAltitudeAglFt < 1500) return "#2563eb";
      if (meanAltitudeAglFt < 3000) return "#0891b2";
      if (meanAltitudeAglFt < 5000) return "#f59e0b";
      return "#dc2626";
    }}

    function drawSelection() {{
      const dataset = trafficDatasets[activeDate];
      const cutoff = dataset.start_time + activeHours * 3600;
      const heatRows = dataset.heatmap[String(activeHours)] || [];
      const rows = dataset.observations.filter((row) => row.time < cutoff);
      const heatMax = heatRows.reduce((maxValue, row) => Math.max(maxValue, row.count), 0) || 1;
      const heatPoints = heatRows.map((row) => [row.lat, row.lon, row.count / heatMax]);
      heatLayer.setLatLngs(heatPoints);

      trackLayer.clearLayers();
      originLayer.clearLayers();
      const flightMap = new Map();

      rows.forEach((row) => {{
        if (!flightMap.has(row.icao24)) {{
          flightMap.set(row.icao24, []);
        }}
        flightMap.get(row.icao24).push(row);
      }});

      // Rebuild tracks client-side from the sampled rows so a single HTML file
      // can support both the heat layer and the track/origin overlays.
      flightMap.forEach((flightRows, icao24) => {{
        if (flightRows.length < 2) {{
          return;
        }}

        flightRows.sort((a, b) => a.time - b.time);
        const points = flightRows.map((row) => [row.lat, row.lon]);
        const firstRow = flightRows[0];
        const meanAltitudeAglFt = flightRows.reduce((total, row) => total + row.altitude_agl_ft, 0) / flightRows.length;
        const color = pickTrackColor(meanAltitudeAglFt);

        const line = L.polyline(points, {{
          color: color,
          weight: 2,
          opacity: 0.55
        }}).bindPopup(
          `<b>ICAO24:</b> ${{icao24}}<br>` +
          `<b>Samples:</b> ${{points.length}}<br>` +
          `<b>Mean altitude:</b> ${{meanAltitudeAglFt.toFixed(0)}} ft AGL`
        );

        const arrows = L.polylineDecorator(line, {{
          patterns: [{{
            offset: "50%",
            repeat: "{arrow_repeat}px",
            symbol: L.Symbol.arrowHead({{
              pixelSize: 8,
              polygon: true,
              pathOptions: {{
                color: color,
                fillOpacity: 0.9,
                weight: 1
              }}
            }})
          }}]
        }});

        const origin = L.circleMarker([firstRow.lat, firstRow.lon], {{
          radius: 5,
          color: "#111827",
          weight: 1.5,
          fillColor: "#000000",
          fillOpacity: 0.95
        }}).bindPopup(
          `<b>Origin point</b><br>` +
          `<b>ICAO24:</b> ${{icao24}}<br>` +
          `<b>First sample:</b> ${{new Date(firstRow.time * 1000).toISOString()}}<br>` +
          `<b>Altitude:</b> ${{firstRow.altitude_agl_ft.toFixed(0)}} ft AGL`
        );

        line.addTo(trackLayer);
        arrows.addTo(trackLayer);
        origin.addTo(originLayer);
      }});

      datasetMetric.textContent = `Dataset: ${{dataset.label}}`;
      rowMetric.textContent = `Rows: ${{dataset.row_count.toLocaleString()}}`;
      flightMetric.textContent = `Flights: ${{dataset.flight_count.toLocaleString()}}`;
    }}

    function updateActiveButtons() {{
      Array.from(dateToggle.querySelectorAll("button")).forEach((button) => {{
        button.classList.toggle("active", button.dataset.date === activeDate);
      }});
      Array.from(timeToggle.querySelectorAll("button")).forEach((button) => {{
        button.classList.toggle("active", Number(button.dataset.hours) === activeHours);
      }});
    }}

    Object.entries(trafficDatasets).forEach(([date, dataset]) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.date = date;
      button.textContent = dataset.label;
      button.addEventListener("click", () => {{
        activeDate = date;
        drawSelection();
        updateActiveButtons();
      }});
      dateToggle.appendChild(button);
    }});

    timeButtons.forEach((hours) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.hours = String(hours);
      button.textContent = `${{hours}}h`;
      button.addEventListener("click", () => {{
        activeHours = hours;
        drawSelection();
        updateActiveButtons();
      }});
      timeToggle.appendChild(button);
    }});

    const landmarkLayer = L.layerGroup(
      Object.entries(locations).map(([name, coordinates]) => {{
        const isAirport = name.endsWith("airport");
        return L.circleMarker(coordinates, {{
          radius: isAirport ? 7 : 6,
          color: "#111827",
          weight: 2,
          fillColor: isAirport ? "#22c55e" : "#8b5cf6",
          fillOpacity: 0.95
        }}).bindPopup(`<b>${{locationLabels[name]}}</b>`);
      }})
    );

    heatLayer.addTo(map);
    trackLayer.addTo(map);
    originLayer.addTo(map);
    landmarkLayer.addTo(map);

    L.control.layers(null, {{
      "Density Heatmap": heatLayer,
      "Flight Tracks": trackLayer,
      "Flight Origins": originLayer,
      "Landmarks": landmarkLayer
    }}, {{
      collapsed: false
    }}).addTo(map);

    drawSelection();
    updateActiveButtons();
  </script>
</body>
</html>
"""


def detect_csv_encoding(csv_path):
    with open(csv_path, "rb") as file_handle:
        first_bytes = file_handle.read(4)

    if first_bytes.startswith(b"\xff\xfe") or first_bytes.startswith(b"\xfe\xff"):
        return "utf-16"
    if first_bytes.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8"


def load_flight_data(csv_path):
    data = pd.read_csv(csv_path, encoding=detect_csv_encoding(csv_path))
    if data.empty:
        raise SystemExit(
            f"No rows found in {csv_path}. The export only contains headers, so there is nothing to plot."
        )
    data = data.dropna(subset=["time", "lat", "lon", "baroaltitude", "icao24"]).copy()
    data["time"] = pd.to_numeric(data["time"], errors="coerce")
    data["lat"] = pd.to_numeric(data["lat"], errors="coerce")
    data["lon"] = pd.to_numeric(data["lon"], errors="coerce")
    data["baroaltitude"] = pd.to_numeric(data["baroaltitude"], errors="coerce")
    data = data.dropna(subset=["time", "lat", "lon", "baroaltitude"])
    if data.empty:
        raise SystemExit(
            f"No plottable rows found in {csv_path}. Check that the export returned valid time, position, and altitude values."
        )
    data["time"] = data["time"].astype(int)
    data["lat"] = data["lat"].astype(float)
    data["lon"] = data["lon"].astype(float)
    data["baroaltitude"] = data["baroaltitude"].astype(float)
    # Approximate AGL altitude from barometric altitude so track coloring stays
    # meaningful near the airport surface environment.
    altitude_ft = data["baroaltitude"] * METERS_TO_FEET
    data["altitude_agl_ft"] = (
        altitude_ft - GROUND_ELEVATION_FT_MSL
    ).clip(lower=0.0)
    data["icao24"] = data["icao24"].astype(str)
    return data


def lonlat_to_web_mercator(lon, lat):
    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)
    x = EARTH_RADIUS_M * lon_rad
    y = EARTH_RADIUS_M * np.log(np.tan(np.pi / 4.0 + lat_rad / 2.0))
    return x, y


def build_windowed_heatmap(data):
    start_epoch = int(data["time"].min())
    windowed = {}
    south, west = settings.SOUTH, settings.WEST
    xmin, ymin = lonlat_to_web_mercator(np.array([west]), np.array([south]))
    xmin = float(xmin[0])
    ymin = float(ymin[0])
    # Align to the same projected grid used by make_grid.py.

    for hours in TIME_WINDOWS:
        cutoff = start_epoch + hours * 3600
        rows = data[data["time"] < cutoff][["lat", "lon", "icao24"]].copy()
        if rows.empty:
            windowed[str(hours)] = []
            continue

        x, y = lonlat_to_web_mercator(
            rows["lon"].to_numpy(),
            rows["lat"].to_numpy(),
        )
        rows["cell_x"] = np.floor((x - xmin) / settings.CELL_SIZE_M).astype(int)
        rows["cell_y"] = np.floor((y - ymin) / settings.CELL_SIZE_M).astype(int)

        cell_counts = (
            rows.drop_duplicates(subset=["cell_x", "cell_y", "icao24"])
            .groupby(["cell_x", "cell_y"], as_index=False)
            .agg(count=("icao24", "size"))
        )
        cell_counts = cell_counts[cell_counts["count"] >= HEATMAP_MIN_CELL_COUNT]

        if cell_counts.empty:
            windowed[str(hours)] = []
            continue

        cell_x = cell_counts["cell_x"].to_numpy(dtype=float)
        cell_y = cell_counts["cell_y"].to_numpy(dtype=float)
        center_x = xmin + (cell_x + 0.5) * settings.CELL_SIZE_M
        center_y = ymin + (cell_y + 0.5) * settings.CELL_SIZE_M
        lon = np.rad2deg(center_x / EARTH_RADIUS_M)
        lat = np.rad2deg(2.0 * np.arctan(np.exp(center_y / EARTH_RADIUS_M)) - np.pi / 2.0)

        windowed[str(hours)] = [
            {
                "lat": float(lat[idx]),
                "lon": float(lon[idx]),
                "count": int(cell_counts.iloc[idx]["count"]),
            }
            for idx in range(len(cell_counts))
        ]

    return windowed


def build_traffic_datasets():
    traffic_datasets = {}
    observation_columns = [
        "time",
        "lat",
        "lon",
        "baroaltitude",
        "altitude_agl_ft",
        "icao24",
    ]

    for dataset in settings.TRAFFIC_DATASETS:
        csv_path = os.path.join(settings.TRAFFIC_FOLDER, dataset["filename"])
        if not os.path.isfile(csv_path):
            raise FileNotFoundError(f"Traffic file not found: {csv_path}")

        data = load_flight_data(csv_path)
        data = data.sort_values(["icao24", "time"]).reset_index(drop=True)
        sampled_data = data[
            data.groupby("icao24").cumcount() % TRACK_STRIDE == 0
        ].reset_index(drop=True)

        traffic_datasets[dataset["date"]] = {
            "label": dataset["label"],
            "row_count": len(data),
            "flight_count": int(data["icao24"].nunique()),
            "start_time": int(data["time"].min()),
            "observations": sampled_data[observation_columns].to_dict(
                orient="records"
            ),
            "heatmap": build_windowed_heatmap(data),
        }

        print(
            f"Loaded {dataset['label']}: "
            f"{len(data)} rows, {data['icao24'].nunique()} flights"
        )

    return traffic_datasets


def main():
    traffic_datasets = build_traffic_datasets()
    arrow_repeat = max(40, ARROW_STRIDE * 12)

    html = HTML_TEMPLATE.format(
        track_stride=TRACK_STRIDE,
        arrow_stride=ARROW_STRIDE,
        traffic_datasets=json.dumps(traffic_datasets),
        time_buttons=json.dumps(TIME_WINDOWS),
        bounds=json.dumps(
            [
                [settings.SOUTH, settings.WEST],
                [settings.NORTH, settings.EAST],
            ]
        ),
        arrow_repeat=arrow_repeat,
        south=settings.SOUTH,
        west=settings.WEST,
        north=settings.NORTH,
        east=settings.EAST,
        center_lat=settings.TRAFFIC_MAP_CENTER_LAT,
        center_lon=settings.TRAFFIC_MAP_CENTER_LON,
        zoom=settings.MAP_ZOOM,
        locations=json.dumps(settings.LOCATIONS),
        location_labels=json.dumps(LOCATION_LABELS),
    )
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

    with open(OUTPUT_PATH, "w", encoding="utf-8") as output_file:
        output_file.write(html)

    print(f"Saved Leaflet map: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
