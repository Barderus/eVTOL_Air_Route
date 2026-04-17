from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


OHARE_LAT = 41.97807408541273
OHARE_LON = -87.90902412382081
MIDWAY_LAT = 41.7856116663475
MIDWAY_LON = -87.75331135429448
TIME_WINDOWS = [1, 3, 6, 9, 12, 24]
METERS_TO_FEET = 3.28084
GROUND_ELEVATION_FT_MSL = 680.0


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Greater Chicago Flight Density</title>
  <link
    rel="stylesheet"
    href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"
    integrity="sha256-p4NxAoJBhIIN+hmNHrzRCf9tD/miZyoHS5obTRR9BMY="
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

    .time-toggle {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-top: 10px;
    }}

    .time-toggle button {{
      border: 0;
      border-radius: 999px;
      padding: 6px 10px;
      background: #cbd5e1;
      color: #0f172a;
      font-size: 12px;
      cursor: pointer;
    }}

    .time-toggle button.active {{
      background: #0f172a;
      color: #f8fafc;
    }}
  </style>
</head>
<body>
  <div id="map"></div>
  <div class="panel">
    <h1>Greater Chicago Flight Density</h1>
    <div class="metric">Rows: {row_count}</div>
    <div class="metric">Flights: {flight_count}</div>
    <div class="metric">Track stride: {track_stride}</div>
    <div class="metric">Arrow stride: {arrow_stride}</div>
    <div class="time-toggle" id="timeToggle"></div>
  </div>

  <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js" integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo=" crossorigin=""></script>
  <script src="https://unpkg.com/leaflet.heat/dist/leaflet-heat.js"></script>
  <script src="https://unpkg.com/leaflet-polylinedecorator@1.6.0/dist/leaflet.polylineDecorator.js"></script>
  <script>
    // Pre-bucket observations by time window so the browser only swaps arrays
    // when the user changes the horizon.
    const windowedObservations = {observation_data};
    const bounds = {bounds};
    const timeButtons = {time_buttons};

    const map = L.map("map", {{
      preferCanvas: true,
      zoomSnap: 0.25
    }});

    map.fitBounds(bounds, {{ padding: [20, 20] }});

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
    const timeToggle = document.getElementById("timeToggle");

    function pickTrackColor(meanAltitudeAglFt) {{
      if (meanAltitudeAglFt < 1500) return "#2563eb";
      if (meanAltitudeAglFt < 3000) return "#0891b2";
      if (meanAltitudeAglFt < 5000) return "#f59e0b";
      return "#dc2626";
    }}

    function drawForHours(hours) {{
      const rows = windowedObservations[String(hours)] || [];
      const heatPoints = rows.map((row) => [row.lat, row.lon, 1.0]);
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
    }}

    function setActiveButton(hours) {{
      Array.from(timeToggle.querySelectorAll("button")).forEach((button) => {{
        button.classList.toggle("active", Number(button.dataset.hours) === hours);
      }});
    }}

    timeButtons.forEach((hours) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.hours = String(hours);
      button.textContent = `${{hours}}h`;
      button.addEventListener("click", () => {{
        drawForHours(hours);
        setActiveButton(hours);
      }});
      timeToggle.appendChild(button);
    }});

    const ohare = L.circleMarker([{ohare_lat}, {ohare_lon}], {{
      radius: 7,
      color: "#111827",
      weight: 2,
      fillColor: "#22c55e",
      fillOpacity: 0.95
    }}).bindPopup("<b>O'Hare International Airport</b>");

    const midway = L.circleMarker([{midway_lat}, {midway_lon}], {{
      radius: 7,
      color: "#111827",
      weight: 2,
      fillColor: "#22c55e",
      fillOpacity: 0.95
    }}).bindPopup("<b>Midway International Airport</b>");

    heatLayer.addTo(map);
    trackLayer.addTo(map);
    originLayer.addTo(map);
    ohare.addTo(map);
    midway.addTo(map);

    L.control.layers(null, {{
      "Density Heatmap": heatLayer,
      "Flight Tracks": trackLayer,
      "Flight Origins": originLayer,
      "O'Hare": ohare,
      "Midway": midway
    }}, {{
      collapsed: false
    }}).addTo(map);

    drawForHours(1);
    setActiveButton(1);
  </script>
</body>
</html>
"""


def detect_csv_encoding(csv_path: Path) -> str:
    with csv_path.open("rb") as file_handle:
        first_bytes = file_handle.read(4)

    if first_bytes.startswith(b"\xff\xfe") or first_bytes.startswith(b"\xfe\xff"):
        return "utf-16"
    if first_bytes.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8"


def load_flight_data(csv_path: Path) -> pd.DataFrame:
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
    data["altitude_agl_ft"] = (data["baroaltitude"] * METERS_TO_FEET - GROUND_ELEVATION_FT_MSL).clip(lower=0.0)
    data["icao24"] = data["icao24"].astype(str)
    return data


def build_windowed_observations(data: pd.DataFrame) -> dict[str, list[dict[str, int | float | str]]]:
    start_epoch = int(data["time"].min())
    windowed: dict[str, list[dict[str, int | float | str]]] = {}

    for hours in TIME_WINDOWS:
        cutoff = start_epoch + hours * 3600
        rows = data[data["time"] < cutoff][["time", "lat", "lon", "baroaltitude", "altitude_agl_ft", "icao24"]]
        windowed[str(hours)] = rows.to_dict(orient="records")

    return windowed


def parse_args() -> argparse.Namespace:
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Create a Leaflet map from an OpenSky CSV file.")
    parser.add_argument(
        "csv_path",
        nargs="?",
        default=folder.parent / "opensky" / "output" / "ohare_2019-03-09_local_30s_15nm_bbox.csv",
        type=Path,
        help="CSV file to read.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=folder.parent / "html" / "ohare_density_leaflet.html",
        help="Output HTML file.",
    )
    parser.add_argument("--track-stride", type=int, default=3, help="Keep every Nth point for each flight.")
    parser.add_argument("--arrow-stride", type=int, default=18, help="Arrow spacing for the track layer.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_flight_data(args.csv_path)
    data = data.sort_values(["icao24", "time"]).reset_index(drop=True)

    track_stride = max(1, args.track_stride)
    arrow_repeat = max(40, args.arrow_stride * 12)

    # Sample each flight independently so dense tracks stay legible without
    # erasing shorter flights from the visualization.
    sampled_data = data[data.groupby("icao24").cumcount() % track_stride == 0].reset_index(drop=True)
    windowed_observations = build_windowed_observations(sampled_data)

    html = HTML_TEMPLATE.format(
        row_count=len(data),
        flight_count=data["icao24"].nunique(),
        track_stride=track_stride,
        arrow_stride=args.arrow_stride,
        observation_data=json.dumps(windowed_observations),
        time_buttons=json.dumps(TIME_WINDOWS),
        bounds=json.dumps(
            [
                [float(data["lat"].min()), float(data["lon"].min())],
                [float(data["lat"].max()), float(data["lon"].max())],
            ]
        ),
        arrow_repeat=arrow_repeat,
        ohare_lat=OHARE_LAT,
        ohare_lon=OHARE_LON,
        midway_lat=MIDWAY_LAT,
        midway_lon=MIDWAY_LON,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")

    print(f"Saved Leaflet map: {args.output}")
    print(f"Rows plotted: {len(data)}")
    print(f"Flights rendered: {data['icao24'].nunique()}")


if __name__ == "__main__":
    main()
