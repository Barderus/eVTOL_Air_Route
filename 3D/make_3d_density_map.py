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
TRACK_COLOR_BANDS = [
    (1500, [37, 99, 235, 180]),
    (3000, [8, 145, 178, 180]),
    (5000, [245, 158, 11, 180]),
]
HIGH_ALTITUDE_TRACK_COLOR = [220, 38, 38, 180]


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>O'Hare 3D Traffic Density</title>
  <style>
    html, body, #app {{
      margin: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      font-family: Arial, sans-serif;
      background: #020617;
    }}

    .panel {{
      position: absolute;
      top: 12px;
      left: 12px;
      z-index: 10;
      width: 360px;
      padding: 14px 16px;
      border-radius: 12px;
      background: rgba(255, 255, 255, 0.94);
      box-shadow: 0 10px 30px rgba(0, 0, 0, 0.22);
    }}

    .panel h1 {{
      margin: 0 0 6px;
      font-size: 19px;
      color: #0f172a;
    }}

    .panel p {{
      margin: 0 0 10px;
      font-size: 13px;
      line-height: 1.35;
      color: #334155;
    }}

    .metrics {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-bottom: 10px;
    }}

    .metric {{
      padding: 4px 8px;
      border-radius: 999px;
      background: #e2e8f0;
      font-size: 12px;
      color: #0f172a;
    }}

    .toggle-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      margin-bottom: 10px;
    }}

    .toggle-row button {{
      border: 0;
      border-radius: 999px;
      padding: 7px 10px;
      background: #cbd5e1;
      color: #0f172a;
      font-size: 12px;
      cursor: pointer;
    }}

    .toggle-row button.active {{
      background: #0f172a;
      color: #f8fafc;
    }}

    .slider-row label {{
      display: block;
      margin-bottom: 6px;
      font-size: 12px;
      color: #334155;
    }}

    .slider-row input {{
      width: 100%;
    }}
  </style>
</head>
<body>
  <div id="app"></div>
  <div class="panel">
    <h1>O'Hare 3D Traffic Density</h1>
    <p>{summary}</p>
    <div class="metrics">
      <div class="metric">Rows: {row_count}</div>
      <div class="metric">Density rows: {interpolated_row_count}</div>
      <div class="metric">Flights: {flight_count}</div>
      <div class="metric">Horizontal bin: {lat_step} deg / {lon_step} deg</div>
      <div class="metric">Vertical bin: {altitude_step} m</div>
    </div>
    <div class="toggle-row" id="timeToggle"></div>
    <div class="slider-row">
      <label for="altitudeRange">Altitude ceiling: <span id="altitudeValue"></span></label>
      <input id="altitudeRange" type="range" min="0" max="{slider_max}" step="{altitude_step}" value="{slider_max}" />
    </div>
  </div>

  <script src="https://unpkg.com/deck.gl@9.0.13/dist.min.js"></script>
  <script>
    const densityData = {density_data};
    const trackData = {track_data};
    const arrowData = {arrow_data};
    const originData = {origin_data};
    const timeButtons = {time_buttons};
    const airportPoints = [
      {{ position: [{ohare_lon}, {ohare_lat}, 0], color: [34, 197, 94, 255], radius: 500 }},
      {{ position: [{midway_lon}, {midway_lat}, 0], color: [245, 158, 11, 255], radius: 500 }}
    ];
    const altitudeSlider = document.getElementById("altitudeRange");
    const altitudeValue = document.getElementById("altitudeValue");
    const timeToggle = document.getElementById("timeToggle");

    let currentHours = 24;
    let currentAltitude = Number(altitudeSlider.value);
    let showDensity = true;
    let showTracks = true;
    let showArrows = true;
    let showOrigins = true;

    function getColor(cost) {{
      if (cost >= 12) return [220, 38, 38, 220];
      if (cost >= 8) return [249, 115, 22, 210];
      if (cost >= 5) return [250, 204, 21, 195];
      if (cost >= 3) return [14, 165, 233, 175];
      return [37, 99, 235, 160];
    }}

    function updateAltitudeText() {{
      altitudeValue.textContent = `${{currentAltitude}} m`;
    }}

    function getVisibleRows() {{
      const rows = densityData[String(currentHours)] || [];
      return rows.filter((row) => row.altitude_top_m <= currentAltitude);
    }}

    function getVisibleTracks() {{
      const rows = trackData[String(currentHours)] || [];
      return rows.filter((row) => row.max_altitude_m <= currentAltitude);
    }}

    function getVisibleArrows() {{
      const rows = arrowData[String(currentHours)] || [];
      return rows.filter((row) => row.target[2] <= currentAltitude);
    }}

    function getVisibleOrigins() {{
      const rows = originData[String(currentHours)] || [];
      return rows.filter((row) => row.position[2] <= currentAltitude);
    }}

    function renderMap() {{
      const rows = getVisibleRows();
      const visibleTracks = getVisibleTracks();
      const visibleArrows = getVisibleArrows();
      const visibleOrigins = getVisibleOrigins();

      const tileLayer = new deck.TileLayer({{
        id: "base-map",
        data: "https://c.tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png",
        minZoom: 0,
        maxZoom: 19,
        renderSubLayers: (props) => {{
          const west = props.tile.bbox.west;
          const south = props.tile.bbox.south;
          const east = props.tile.bbox.east;
          const north = props.tile.bbox.north;
          return new deck.BitmapLayer(props, {{
            data: null,
            image: props.data,
            bounds: [west, south, east, north]
          }});
        }}
      }});

      const pointLayer = new deck.ScatterplotLayer({{
        id: "density-points",
        data: rows,
        pickable: true,
        filled: true,
        stroked: false,
        radiusUnits: "meters",
        radiusMinPixels: 2,
        radiusMaxPixels: 40,
        getPosition: (row) => [row.lon, row.lat, row.altitude_mid_m],
        getRadius: (row) => row.radius_m,
        getFillColor: (row) => getColor(row.cost)
      }});

      const trackLayer = new deck.PathLayer({{
        id: "flight-tracks",
        data: visibleTracks,
        pickable: true,
        widthUnits: "meters",
        widthMinPixels: 2,
        getPath: (row) => row.path,
        getColor: (row) => row.color,
        getWidth: 140
      }});

      const arrowLayer = new deck.LineLayer({{
        id: "flight-arrows",
        data: visibleArrows,
        pickable: true,
        widthUnits: "meters",
        widthMinPixels: 2,
        getSourcePosition: (row) => row.source,
        getTargetPosition: (row) => row.target,
        getColor: (row) => row.color,
        getWidth: 120
      }});

      const airportLayer = new deck.ScatterplotLayer({{
        id: "airport-point",
        data: airportPoints,
        radiusUnits: "meters",
        getPosition: (row) => row.position,
        getRadius: (row) => row.radius,
        getFillColor: (row) => row.color,
        getLineColor: [15, 23, 42, 255],
        lineWidthMinPixels: 2
      }});

      const originLayer = new deck.ScatterplotLayer({{
        id: "flight-origins",
        data: visibleOrigins,
        pickable: true,
        filled: true,
        stroked: true,
        radiusUnits: "meters",
        radiusMinPixels: 3,
        radiusMaxPixels: 20,
        getPosition: (row) => row.position,
        getRadius: (row) => row.radius_m,
        getFillColor: [0, 0, 0, 230],
        getLineColor: [17, 24, 39, 255],
        lineWidthMinPixels: 1
      }});

      const layers = [tileLayer];
      if (showDensity) {{
        layers.push(pointLayer);
      }}
      if (showTracks) {{
        layers.push(trackLayer);
      }}
      if (showArrows) {{
        layers.push(arrowLayer);
      }}
      if (showOrigins) {{
        layers.push(originLayer);
      }}
      layers.push(airportLayer);

      deckMap.setProps({{
        layers: layers
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
        currentHours = hours;
        setActiveButton(hours);
        renderMap();
      }});
      timeToggle.appendChild(button);
    }});

    altitudeSlider.addEventListener("input", () => {{
      currentAltitude = Number(altitudeSlider.value);
      updateAltitudeText();
      renderMap();
    }});

    const viewToggles = [
      {{ label: "Density", value: "density" }},
      {{ label: "Tracks", value: "tracks" }},
      {{ label: "Arrows", value: "arrows" }},
      {{ label: "Origins", value: "origins" }}
    ];

    const toggleRow = document.createElement("div");
    toggleRow.className = "toggle-row";
    viewToggles.forEach((toggle) => {{
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = toggle.label;
      button.classList.add("active");
      button.addEventListener("click", () => {{
        if (toggle.value === "density") {{
          showDensity = !showDensity;
          button.classList.toggle("active", showDensity);
        }} else if (toggle.value === "tracks") {{
          showTracks = !showTracks;
          button.classList.toggle("active", showTracks);
        }} else if (toggle.value === "arrows") {{
          showArrows = !showArrows;
          button.classList.toggle("active", showArrows);
        }} else {{
          showOrigins = !showOrigins;
          button.classList.toggle("active", showOrigins);
        }}
        renderMap();
      }});
      toggleRow.appendChild(button);
    }});
    document.querySelector(".panel").appendChild(toggleRow);

    const deckMap = new deck.DeckGL({{
      container: "app",
      controller: true,
      map: false,
      getTooltip: ({object}) => object && (
        object.count !== undefined
          ? (
            `Cost: ${{object.cost}}\\n` +
            `Count: ${{object.count}}\\n` +
            `Altitude band: ${{object.altitude_bottom_m}}-${{object.altitude_top_m}} m\\n` +
            `Cell center: ${{object.lat.toFixed(4)}}, ${{object.lon.toFixed(4)}}`
          )
          : object.position
            ? `Origin point\\nICAO24: ${{object.icao24}}\\nFirst sample altitude: ${{object.altitude_m.toFixed(0)}} m`
            : object.icao24
            ? `ICAO24: ${{object.icao24}}\\nSamples: ${{object.samples}}\\nMean altitude: ${{object.mean_altitude_m.toFixed(0)}} m`
            : null
      ),
      initialViewState: {{
        longitude: {ohare_lon},
        latitude: {ohare_lat},
        zoom: 9.8,
        pitch: 58,
        bearing: 20,
        maxPitch: 85
      }}
    }});

    updateAltitudeText();
    setActiveButton(24);
    renderMap();
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
    data = data.dropna(subset=["time", "lat", "lon", "baroaltitude", "icao24"]).copy()
    data["time"] = pd.to_numeric(data["time"], errors="coerce")
    data["lat"] = pd.to_numeric(data["lat"], errors="coerce")
    data["lon"] = pd.to_numeric(data["lon"], errors="coerce")
    data["baroaltitude"] = pd.to_numeric(data["baroaltitude"], errors="coerce")
    data = data.dropna(subset=["time", "lat", "lon", "baroaltitude"])
    data["time"] = data["time"].astype(int)
    data["lat"] = data["lat"].astype(float)
    data["lon"] = data["lon"].astype(float)
    data["baroaltitude"] = data["baroaltitude"].astype(float)
    data["icao24"] = data["icao24"].astype(str)
    return data.sort_values("time").reset_index(drop=True)


def interpolate_flight_paths(
    data: pd.DataFrame,
    step_seconds: int,
    max_gap_seconds: int,
) -> pd.DataFrame:
    if step_seconds <= 0:
        return data

    interpolated_groups: list[pd.DataFrame] = []
    columns = ["time", "lat", "lon", "baroaltitude", "icao24"]

    for icao24, flight_rows in data.groupby("icao24", sort=False):
        ordered = flight_rows.sort_values("time").reset_index(drop=True)
        if len(ordered) < 2:
            interpolated_groups.append(ordered[columns].copy())
            continue

        flight_records: list[dict[str, float | int | str]] = []
        for index in range(len(ordered) - 1):
            start = ordered.iloc[index]
            end = ordered.iloc[index + 1]
            start_time = int(start["time"])
            end_time = int(end["time"])
            gap_seconds = end_time - start_time

            flight_records.append(
                {
                    "time": start_time,
                    "lat": float(start["lat"]),
                    "lon": float(start["lon"]),
                    "baroaltitude": float(start["baroaltitude"]),
                    "icao24": str(icao24),
                }
            )

            if gap_seconds <= step_seconds or gap_seconds > max_gap_seconds:
                continue

            steps = range(step_seconds, gap_seconds, step_seconds)
            for elapsed in steps:
                ratio = elapsed / gap_seconds
                flight_records.append(
                    {
                        "time": start_time + elapsed,
                        "lat": float(start["lat"] + (end["lat"] - start["lat"]) * ratio),
                        "lon": float(start["lon"] + (end["lon"] - start["lon"]) * ratio),
                        "baroaltitude": float(start["baroaltitude"] + (end["baroaltitude"] - start["baroaltitude"]) * ratio),
                        "icao24": str(icao24),
                    }
                )

        last_row = ordered.iloc[-1]
        flight_records.append(
            {
                "time": int(last_row["time"]),
                "lat": float(last_row["lat"]),
                "lon": float(last_row["lon"]),
                "baroaltitude": float(last_row["baroaltitude"]),
                "icao24": str(icao24),
            }
        )
        interpolated_groups.append(pd.DataFrame.from_records(flight_records, columns=columns))

    return pd.concat(interpolated_groups, ignore_index=True).sort_values("time").reset_index(drop=True)


def pick_track_color(mean_altitude: float) -> list[int]:
    for max_altitude, color in TRACK_COLOR_BANDS:
        if mean_altitude < max_altitude:
            return color
    return HIGH_ALTITUDE_TRACK_COLOR


def build_density_rows(data: pd.DataFrame, lat_step: float, lon_step: float, altitude_step: int) -> dict[str, list[dict[str, float | int]]]:
    start_time = int(data["time"].min())
    all_windows: dict[str, list[dict[str, float | int]]] = {}

    for hours in TIME_WINDOWS:
        cutoff = start_time + hours * 3600
        window_data = data[data["time"] < cutoff].copy()
        if window_data.empty:
            all_windows[str(hours)] = []
            continue

        window_data["lat_bin"] = ((window_data["lat"] / lat_step).round().astype(int)) * lat_step
        window_data["lon_bin"] = ((window_data["lon"] / lon_step).round().astype(int)) * lon_step
        window_data["altitude_bin"] = (window_data["baroaltitude"] // altitude_step).astype(int)

        grouped = (
            window_data.groupby(["lat_bin", "lon_bin", "altitude_bin"], as_index=False)
            .agg(count=("time", "size"))
            .sort_values("count", ascending=False)
        )
        horizontal_costs = (
            window_data.groupby(["lat_bin", "lon_bin"], as_index=False)
            .agg(cost=("time", "size"))
        )
        grouped = grouped.merge(horizontal_costs, on=["lat_bin", "lon_bin"], how="left")

        radius_meters = max(150, int(min(lat_step, lon_step) * 111_000 * 0.45))
        rows = []
        for row in grouped.itertuples():
            bottom = int(row.altitude_bin * altitude_step)
            top = int((row.altitude_bin + 1) * altitude_step)
            middle = int(bottom + altitude_step / 2)
            rows.append(
                {
                    "lat": float(row.lat_bin),
                    "lon": float(row.lon_bin),
                    "altitude_bottom_m": bottom,
                    "altitude_top_m": top,
                    "altitude_mid_m": middle,
                    "count": int(row.count),
                    "cost": int(row.cost),
                    "radius_m": radius_meters,
                }
            )

        all_windows[str(hours)] = rows

    return all_windows


def flatten_cost_rows(density_rows: dict[str, list[dict[str, float | int]]]) -> pd.DataFrame:
    flattened: list[dict[str, float | int | str]] = []
    for hours, rows in density_rows.items():
        for row in rows:
            flattened.append(
                {
                    "hours": str(hours),
                    "lat": row["lat"],
                    "lon": row["lon"],
                    "altitude_bottom_m": row["altitude_bottom_m"],
                    "altitude_top_m": row["altitude_top_m"],
                    "altitude_mid_m": row["altitude_mid_m"],
                    "observation_count": row["count"],
                    "cost": row["cost"],
                }
            )

    return pd.DataFrame(flattened)


def build_track_rows(
    data: pd.DataFrame,
    track_stride: int,
    arrow_stride: int,
) -> tuple[
    dict[str, list[dict[str, object]]],
    dict[str, list[dict[str, object]]],
    dict[str, list[dict[str, object]]],
]:
    start_time = int(data["time"].min())
    all_tracks: dict[str, list[dict[str, object]]] = {}
    all_arrows: dict[str, list[dict[str, object]]] = {}
    all_origins: dict[str, list[dict[str, object]]] = {}

    for hours in TIME_WINDOWS:
        cutoff = start_time + hours * 3600
        window_data = data[data["time"] < cutoff].copy()
        if window_data.empty:
            all_tracks[str(hours)] = []
            all_arrows[str(hours)] = []
            all_origins[str(hours)] = []
            continue

        sampled = window_data[window_data.groupby("icao24").cumcount() % track_stride == 0]
        track_rows: list[dict[str, object]] = []
        arrow_rows: list[dict[str, object]] = []
        origin_rows: list[dict[str, object]] = []

        for icao24, flight_rows in sampled.groupby("icao24", sort=False):
            ordered = flight_rows.sort_values("time")
            if len(ordered) < 2:
                continue

            mean_altitude = float(ordered["baroaltitude"].mean())
            color = pick_track_color(mean_altitude)
            first_row = ordered.iloc[0]
            path = [
                [float(row.lon), float(row.lat), float(row.baroaltitude)]
                for row in ordered.itertuples()
            ]
            track_rows.append(
                {
                    "icao24": str(icao24),
                    "path": path,
                    "color": color,
                    "samples": len(path),
                    "mean_altitude_m": mean_altitude,
                    "max_altitude_m": float(ordered["baroaltitude"].max()),
                }
            )
            origin_rows.append(
                {
                    "icao24": str(icao24),
                    "position": [float(first_row["lon"]), float(first_row["lat"]), float(first_row["baroaltitude"])],
                    "altitude_m": float(first_row["baroaltitude"]),
                    "radius_m": 250.0,
                }
            )

            for index in range(1, len(path), arrow_stride):
                arrow_rows.append(
                    {
                        "icao24": str(icao24),
                        "source": path[index - 1],
                        "target": path[index],
                        "color": color,
                        "samples": len(path),
                        "mean_altitude_m": mean_altitude,
                    }
                )

        all_tracks[str(hours)] = track_rows
        all_arrows[str(hours)] = arrow_rows
        all_origins[str(hours)] = origin_rows

    return all_tracks, all_arrows, all_origins


def parse_args() -> argparse.Namespace:
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Create a 3D density map from an OpenSky CSV file.")
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
        default=folder.parent / "html" / "ohare_3d_density_map.html",
        help="Output HTML file.",
    )
    parser.add_argument(
        "--cost-output",
        type=Path,
        default=folder / "output" / "ohare_3d_cell_costs.csv",
        help="Output CSV file for 3D cell costs.",
    )
    parser.add_argument("--lat-step", type=float, default=0.01, help="Latitude bin size in degrees.")
    parser.add_argument("--lon-step", type=float, default=0.01, help="Longitude bin size in degrees.")
    parser.add_argument("--altitude-step-m", type=int, default=250, help="Altitude bin size in meters.")
    parser.add_argument(
        "--interpolate-step-seconds",
        type=int,
        default=5,
        help="Insert interpolated flight positions every N seconds before density binning.",
    )
    parser.add_argument(
        "--max-interpolate-gap-seconds",
        type=int,
        default=180,
        help="Only interpolate gaps up to this size to avoid bridging missing segments.",
    )
    parser.add_argument("--track-stride", type=int, default=3, help="Keep every Nth point for each flight track.")
    parser.add_argument("--arrow-stride", type=int, default=6, help="Keep one directional segment every N track points.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_flight_data(args.csv_path)
    density_data = interpolate_flight_paths(
        data,
        step_seconds=max(1, args.interpolate_step_seconds),
        max_gap_seconds=max(1, args.max_interpolate_gap_seconds),
    )
    density_rows = build_density_rows(density_data, args.lat_step, args.lon_step, args.altitude_step_m)
    cost_rows = flatten_cost_rows(density_rows)
    track_rows, arrow_rows, origin_rows = build_track_rows(
        data,
        track_stride=max(1, args.track_stride),
        arrow_stride=max(1, args.arrow_stride),
    )

    slider_max = int(((data["baroaltitude"].max() // args.altitude_step_m) + 1) * args.altitude_step_m)
    summary = (
        "This view groups OpenSky state vectors into longitude, latitude, and altitude boxes. "
        "Each point shows one occupied box and the color gets stronger when more observations fall inside it."
    )

    replacements = {
        "summary": summary,
        "row_count": len(data),
        "interpolated_row_count": len(density_data),
        "flight_count": data["icao24"].nunique(),
        "lat_step": args.lat_step,
        "lon_step": args.lon_step,
        "altitude_step": args.altitude_step_m,
        "slider_max": slider_max,
        "density_data": json.dumps(density_rows),
        "track_data": json.dumps(track_rows),
        "arrow_data": json.dumps(arrow_rows),
        "origin_data": json.dumps(origin_rows),
        "time_buttons": json.dumps(TIME_WINDOWS),
        "ohare_lon": OHARE_LON,
        "ohare_lat": OHARE_LAT,
        "midway_lon": MIDWAY_LON,
        "midway_lat": MIDWAY_LAT,
    }
    html = HTML_TEMPLATE
    for key, value in replacements.items():
        html = html.replace(f"{{{key}}}", str(value))
    html = html.replace("{{", "{").replace("}}", "}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    args.cost_output.parent.mkdir(parents=True, exist_ok=True)
    cost_rows.to_csv(args.cost_output, index=False)

    print(f"Saved 3D HTML map: {args.output}")
    print(f"Saved 3D cell costs: {args.cost_output}")
    print(f"Rows processed: {len(data)}")
    print(f"Flights processed: {data['icao24'].nunique()}")


if __name__ == "__main__":
    main()
