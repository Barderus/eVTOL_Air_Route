from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


OHARE_LAT = 41.97807408541273
OHARE_LON = -87.90902412382081


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>O'Hare Flight Density</title>
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
    <h1>O'Hare Flight Density</h1>
    <p>{summary}</p>
    <div class="metric">Rows: {row_count}</div>
    <div class="metric">Flights: {flight_count}</div>
    <div class="metric">Track stride: {track_stride}</div>
    <div class="metric">Arrow stride: {arrow_stride}</div>
    <div class="time-toggle" id="timeToggle"></div>
  </div>

  <script
    src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"
    integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo="
    crossorigin=""
  ></script>
  <script src="https://unpkg.com/leaflet.heat/dist/leaflet-heat.js"></script>
  <script src="https://unpkg.com/leaflet-polylinedecorator@1.6.0/dist/leaflet.polylineDecorator.js"></script>
  <script>
    const observationData = {observation_data};
    const bounds = {bounds};
    const durationOptions = [1, 3, 6, 9, 12, 24];
    const startEpoch = {start_epoch};
    const endEpoch = {end_epoch};

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
    const timeToggle = document.getElementById("timeToggle");

    function altitudeColor(meanAltitudeM) {{
      if (meanAltitudeM < 1500) return "#2563eb";
      if (meanAltitudeM < 3000) return "#0891b2";
      if (meanAltitudeM < 5000) return "#f59e0b";
      return "#dc2626";
    }}

    function rebuildLayers(hours) {{
      const cutoff = Math.min(endEpoch, startEpoch + hours * 3600);
      const filtered = observationData.filter((row) => row.time < cutoff);
      const heatPoints = filtered.map((row) => [row.lat, row.lon, 1.0]);
      heatLayer.setLatLngs(heatPoints);

      trackLayer.clearLayers();
      const grouped = new Map();
      filtered.forEach((row) => {{
        if (!grouped.has(row.icao24)) grouped.set(row.icao24, []);
        grouped.get(row.icao24).push(row);
      }});

      grouped.forEach((rows, icao24) => {{
        if (rows.length < 2) return;
        rows.sort((a, b) => a.time - b.time);
        const points = rows.map((row) => [row.lat, row.lon]);
        const meanAltitudeM = rows.reduce((acc, row) => acc + row.baroaltitude, 0) / rows.length;
        const color = altitudeColor(meanAltitudeM);

        const polyline = L.polyline(points, {{
          color,
          weight: 2,
          opacity: 0.55
        }}).bindPopup(
          `<b>ICAO24:</b> ${{icao24}}<br>` +
          `<b>Samples:</b> ${{points.length}}<br>` +
          `<b>Mean altitude:</b> ${{meanAltitudeM.toFixed(0)}} m`
        );

        const arrows = L.polylineDecorator(polyline, {{
          patterns: [{{
            offset: "50%",
            repeat: "{arrow_repeat}px",
            symbol: L.Symbol.arrowHead({{
              pixelSize: 8,
              polygon: true,
              pathOptions: {{
                color,
                fillOpacity: 0.9,
                weight: 1
              }}
            }})
          }}]
        }});

        polyline.addTo(trackLayer);
        arrows.addTo(trackLayer);
      }});
    }}

    function setActiveButton(hours) {{
      Array.from(timeToggle.querySelectorAll("button")).forEach((button) => {{
        button.classList.toggle("active", Number(button.dataset.hours) === hours);
      }});
    }}

    durationOptions.forEach((hours) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.hours = String(hours);
      button.textContent = `${{hours}}h`;
      button.addEventListener("click", () => {{
        rebuildLayers(hours);
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

    const overlays = {{
      "Density Heatmap": heatLayer,
      "Flight Tracks": trackLayer,
      "O'Hare": ohare
    }};

    heatLayer.addTo(map);
    trackLayer.addTo(map);
    ohare.addTo(map);
    L.control.layers(null, overlays, {{ collapsed: false }}).addTo(map);
    rebuildLayers(24);
    setActiveButton(24);
  </script>
</body>
</html>
"""


def detect_csv_encoding(csv_path: Path) -> str:
    with csv_path.open("rb") as handle:
        start = handle.read(4)

    if start.startswith(b"\xff\xfe") or start.startswith(b"\xfe\xff"):
        return "utf-16"
    if start.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8"


def altitude_color(mean_altitude_m: float) -> str:
    if mean_altitude_m < 1500:
        return "#2563eb"
    if mean_altitude_m < 3000:
        return "#0891b2"
    if mean_altitude_m < 5000:
        return "#f59e0b"
    return "#dc2626"


def parse_args() -> argparse.Namespace:
    project_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Generate a Leaflet flight-density map from an OpenSky CSV export.")
    parser.add_argument(
        "csv_path",
        nargs="?",
        default=project_dir / "ohare_2019-03-09_local_1min_15nm_bbox.csv",
        type=Path,
        help="Path to the OpenSky CSV export.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=project_dir / "ohare_density_leaflet.html",
        help="Output HTML map path.",
    )
    parser.add_argument(
        "--track-stride",
        type=int,
        default=3,
        help="Keep every Nth point per flight track to control HTML size.",
    )
    parser.add_argument(
        "--arrow-stride",
        type=int,
        default=18,
        help="Arrow spacing in track points; larger values render fewer arrows.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    encoding = detect_csv_encoding(args.csv_path)

    df = pd.read_csv(args.csv_path, encoding=encoding)
    df = df.dropna(subset=["time", "lat", "lon", "baroaltitude", "icao24"]).copy()
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
    df["baroaltitude"] = pd.to_numeric(df["baroaltitude"], errors="coerce")
    df = df.dropna(subset=["time", "lat", "lon", "baroaltitude"])
    df = df.sort_values(["icao24", "time"]).reset_index(drop=True)

    track_stride = max(1, args.track_stride)
    arrow_repeat = max(40, args.arrow_stride * 12)

    sampled_df = df[df.groupby("icao24").cumcount() % track_stride == 0].reset_index(drop=True)
    observation_frame = sampled_df[["time", "lat", "lon", "baroaltitude", "icao24"]].copy()
    observation_frame["time"] = observation_frame["time"].astype(int)
    observation_frame["lat"] = observation_frame["lat"].astype(float)
    observation_frame["lon"] = observation_frame["lon"].astype(float)
    observation_frame["baroaltitude"] = observation_frame["baroaltitude"].astype(float)
    observation_frame["icao24"] = observation_frame["icao24"].astype(str)
    observation_data = observation_frame.to_dict(orient="records")

    bounds = [
        [float(df["lat"].min()), float(df["lon"].min())],
        [float(df["lat"].max()), float(df["lon"].max())],
    ]
    start_epoch = int(df["time"].min())
    end_epoch = int(df["time"].max()) + 1

    summary = (
        f"Heat points come from every sampled observation in the CSV. "
        f"Tracks are grouped by ICAO24, downsampled every {track_stride} points, "
        f"and drawn with arrowheads so direction stays visible."
    )

    html = HTML_TEMPLATE.format(
        summary=summary,
        row_count=len(df),
        flight_count=df["icao24"].nunique(),
        track_stride=track_stride,
        arrow_stride=args.arrow_stride,
        observation_data=json.dumps(observation_data),
        bounds=json.dumps(bounds),
        arrow_repeat=arrow_repeat,
        ohare_lat=OHARE_LAT,
        ohare_lon=OHARE_LON,
        start_epoch=start_epoch,
        end_epoch=end_epoch,
    )
    args.output.write_text(html, encoding="utf-8")

    print(f"Saved Leaflet map: {args.output}")
    print(f"Rows plotted: {len(df)}")
    print(f"Flights rendered: {df['icao24'].nunique()}")
    print(f"CSV encoding: {encoding}")


if __name__ == "__main__":
    main()
