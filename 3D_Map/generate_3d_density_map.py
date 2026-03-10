from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


OHARE_LAT = 41.97807408541273
OHARE_LON = -87.90902412382081
TIME_WINDOWS_HOURS = [1, 3, 6, 9, 12, 24]


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>O'Hare 3D Traffic Density</title>
  <style>
    html, body, #app {
      margin: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      font-family: Arial, sans-serif;
      background: #020617;
    }

    .panel {
      position: absolute;
      top: 12px;
      left: 12px;
      z-index: 10;
      width: 360px;
      padding: 14px 16px;
      border-radius: 12px;
      background: rgba(255, 255, 255, 0.94);
      box-shadow: 0 10px 30px rgba(0, 0, 0, 0.22);
    }

    .panel h1 {
      margin: 0 0 6px;
      font-size: 19px;
      color: #0f172a;
    }

    .panel p {
      margin: 0 0 10px;
      font-size: 13px;
      line-height: 1.35;
      color: #334155;
    }

    .metrics {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-bottom: 10px;
    }

    .metric {
      padding: 4px 8px;
      border-radius: 999px;
      background: #e2e8f0;
      font-size: 12px;
      color: #0f172a;
    }

    .toggle-row {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-bottom: 10px;
    }

    .toggle-row button {
      border: 0;
      border-radius: 999px;
      padding: 7px 10px;
      background: #cbd5e1;
      color: #0f172a;
      font-size: 12px;
      cursor: pointer;
    }

    .toggle-row button.active {
      background: #0f172a;
      color: #f8fafc;
    }

    .slider-row label {
      display: block;
      margin-bottom: 6px;
      font-size: 12px;
      color: #334155;
    }

    .slider-row input {
      width: 100%;
    }
  </style>
</head>
<body>
  <div id="app"></div>
  <div class="panel">
    <h1>O'Hare 3D Traffic Density</h1>
    <p>{summary}</p>
    <div class="metrics">
      <div class="metric">Rows: {row_count}</div>
      <div class="metric">Flights: {flight_count}</div>
      <div class="metric">Horizontal bin: {lat_step_deg}° / {lon_step_deg}°</div>
      <div class="metric">Vertical bin: {altitude_step_m} m</div>
    </div>
    <div class="toggle-row" id="timeToggle"></div>
    <div class="slider-row">
      <label for="altitudeRange">Altitude ceiling: <span id="altitudeValue"></span></label>
      <input id="altitudeRange" type="range" min="{slider_min}" max="{slider_max}" step="{slider_step}" value="{slider_max}" />
    </div>
  </div>

  <script src="https://unpkg.com/deck.gl@9.0.13/dist.min.js"></script>
  <script>
    const densityByWindow = {density_by_window};
    const timeWindows = {time_windows};
    const maxAltitudeMeters = {slider_max};
    const ohare = [{ohare_lon}, {ohare_lat}, 0];
    const bounds = {bounds};

    const altitudeRange = document.getElementById("altitudeRange");
    const altitudeValue = document.getElementById("altitudeValue");
    const timeToggle = document.getElementById("timeToggle");

    let activeHours = 24;
    let activeAltitude = Number(altitudeRange.value);

    function updateAltitudeLabel() {{
      altitudeValue.textContent = `${{activeAltitude}} m`;
    }}

    function colorForCount(count) {{
      if (count >= 12) return [220, 38, 38, 220];
      if (count >= 8) return [249, 115, 22, 210];
      if (count >= 5) return [250, 204, 21, 195];
      if (count >= 3) return [14, 165, 233, 175];
      return [37, 99, 235, 160];
    }}

    function getFilteredData() {{
      const rows = densityByWindow[String(activeHours)] || [];
      return rows.filter((row) => row.altitude_top_m <= activeAltitude);
    }}

    function render() {{
      const data = getFilteredData();

      const tileLayer = new deck.TileLayer({{
        id: "osm-base",
        data: "https://c.tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png",
        minZoom: 0,
        maxZoom: 19,
        renderSubLayers: (props) => {{
          const {{
            bbox: {{west, south, east, north}}
          }} = props.tile;
          return new deck.BitmapLayer(props, {{
            data: null,
            image: props.data,
            bounds: [west, south, east, north]
          }});
        }}
      }});

      const densityLayer = new deck.ScatterplotLayer({{
        id: "density-points",
        data,
        pickable: true,
        stroked: false,
        filled: true,
        radiusUnits: "meters",
        radiusScale: 1,
        radiusMinPixels: 2,
        radiusMaxPixels: 40,
        getPosition: (d) => [d.lon, d.lat, d.altitude_mid_m],
        getRadius: (d) => d.radius_m,
        getFillColor: (d) => colorForCount(d.count),
        getLineColor: [15, 23, 42, 220],
        getElevation: (d) => d.altitude_mid_m
      }});

      const airportLayer = new deck.ScatterplotLayer({{
        id: "ohare-marker",
        data: [{{position: ohare}}],
        radiusUnits: "meters",
        getPosition: (d) => d.position,
        getRadius: 500,
        getFillColor: [34, 197, 94, 255],
        getLineColor: [15, 23, 42, 255],
        lineWidthMinPixels: 2
      }});

      deckgl.setProps({{
        layers: [tileLayer, densityLayer, airportLayer]
      }});
    }}

    function setActiveButton(hours) {{
      Array.from(timeToggle.querySelectorAll("button")).forEach((button) => {{
        button.classList.toggle("active", Number(button.dataset.hours) === hours);
      }});
    }}

    timeWindows.forEach((hours) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.dataset.hours = String(hours);
      button.textContent = `${{hours}}h`;
      button.addEventListener("click", () => {{
        activeHours = hours;
        setActiveButton(hours);
        render();
      }});
      timeToggle.appendChild(button);
    }});

    altitudeRange.addEventListener("input", () => {{
      activeAltitude = Number(altitudeRange.value);
      updateAltitudeLabel();
      render();
    }});

    const deckgl = new deck.DeckGL({{
      container: "app",
      controller: true,
      map: false,
      getTooltip: ({object}) => object && (
        `Count: ${{object.count}}\\n` +
        `Alt band: ${{object.altitude_bottom_m}}-${{object.altitude_top_m}} m\\n` +
        `Cell center: ${{object.lat.toFixed(4)}}, ${{object.lon.toFixed(4)}}`
      ),
      initialViewState: {{
        longitude: {center_lon},
        latitude: {center_lat},
        zoom: 9.8,
        pitch: 58,
        bearing: 20,
        maxPitch: 85
      }}
    }});

    updateAltitudeLabel();
    setActiveButton(activeHours);
    render();
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


def parse_args() -> argparse.Namespace:
    root_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description="Generate a 3D traffic density HTML map.")
    parser.add_argument(
        "csv_path",
        nargs="?",
        default=root_dir / "OpenSky-Trino" / "ohare_2019-03-09_local_1min_15nm_bbox.csv",
        type=Path,
        help="Path to the OpenSky CSV export.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "ohare_3d_density_map.html",
        help="Output HTML path.",
    )
    parser.add_argument("--lat-step", type=float, default=0.01, help="Latitude bin size in degrees.")
    parser.add_argument("--lon-step", type=float, default=0.01, help="Longitude bin size in degrees.")
    parser.add_argument("--altitude-step-m", type=int, default=250, help="Altitude bin size in meters.")
    return parser.parse_args()


def prepare_dataframe(csv_path: Path) -> pd.DataFrame:
    encoding = detect_csv_encoding(csv_path)
    df = pd.read_csv(csv_path, encoding=encoding)
    df = df.dropna(subset=["time", "lat", "lon", "baroaltitude", "icao24"]).copy()
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
    df["baroaltitude"] = pd.to_numeric(df["baroaltitude"], errors="coerce")
    df = df.dropna(subset=["time", "lat", "lon", "baroaltitude"])
    df["time"] = df["time"].astype(int)
    return df.sort_values("time").reset_index(drop=True)


def build_density_windows(
    df: pd.DataFrame,
    lat_step: float,
    lon_step: float,
    altitude_step_m: int,
) -> dict[str, list[dict[str, float | int]]]:
    start_epoch = int(df["time"].min())
    density_by_window: dict[str, list[dict[str, float | int]]] = {}

    for hours in TIME_WINDOWS_HOURS:
        cutoff = start_epoch + hours * 3600
        window_df = df[df["time"] < cutoff].copy()
        if window_df.empty:
            density_by_window[str(hours)] = []
            continue

        window_df["lat_bin"] = ((window_df["lat"] / lat_step).round().astype(int)) * lat_step
        window_df["lon_bin"] = ((window_df["lon"] / lon_step).round().astype(int)) * lon_step
        window_df["alt_bin"] = (window_df["baroaltitude"] // altitude_step_m).astype(int)

        grouped = (
            window_df.groupby(["lat_bin", "lon_bin", "alt_bin"], as_index=False)
            .agg(count=("time", "size"))
            .sort_values("count", ascending=False)
        )

        radius_m = max(150, int(min(lat_step, lon_step) * 111_000 * 0.45))
        density_by_window[str(hours)] = [
            {
                "lat": float(row.lat_bin),
                "lon": float(row.lon_bin),
                "altitude_bottom_m": int(row.alt_bin * altitude_step_m),
                "altitude_top_m": int((row.alt_bin + 1) * altitude_step_m),
                "altitude_mid_m": int((row.alt_bin * altitude_step_m) + altitude_step_m / 2),
                "count": int(row.count),
                "radius_m": radius_m,
            }
            for row in grouped.itertuples()
        ]

    return density_by_window


def main() -> None:
    args = parse_args()
    df = prepare_dataframe(args.csv_path)
    density_by_window = build_density_windows(df, args.lat_step, args.lon_step, args.altitude_step_m)

    bounds = [
        [float(df["lon"].min()), float(df["lat"].min())],
        [float(df["lon"].max()), float(df["lat"].max())],
    ]
    slider_max = int(((df["baroaltitude"].max() // args.altitude_step_m) + 1) * args.altitude_step_m)
    slider_min = 0
    slider_step = args.altitude_step_m

    summary = (
        "This view bins OpenSky state vectors into longitude, latitude, and altitude voxels. "
        "Each marker represents one occupied voxel, colored by observation count and positioned at the middle of its altitude band."
    )

    replacements = {
        "summary": summary,
        "row_count": len(df),
        "flight_count": df["icao24"].nunique(),
        "lat_step_deg": args.lat_step,
        "lon_step_deg": args.lon_step,
        "altitude_step_m": args.altitude_step_m,
        "density_by_window": json.dumps(density_by_window),
        "time_windows": json.dumps(TIME_WINDOWS_HOURS),
        "ohare_lon": OHARE_LON,
        "ohare_lat": OHARE_LAT,
        "bounds": json.dumps(bounds),
        "slider_min": slider_min,
        "slider_max": slider_max,
        "slider_step": slider_step,
        "center_lon": OHARE_LON,
        "center_lat": OHARE_LAT,
    }
    html = HTML_TEMPLATE
    for key, value in replacements.items():
        html = html.replace(f"{{{key}}}", str(value))
    html = html.replace("{{", "{").replace("}}", "}")
    args.output.write_text(html, encoding="utf-8")

    print(f"Saved 3D HTML map: {args.output}")
    print(f"Rows processed: {len(df)}")
    print(f"Flights processed: {df['icao24'].nunique()}")


if __name__ == "__main__":
    main()
