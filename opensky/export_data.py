from __future__ import annotations

import argparse
import math
import os
import subprocess
from datetime import UTC, datetime, timedelta
from pathlib import Path
from string import Template
from zoneinfo import ZoneInfo


OHARE_LAT = 41.97807408541273
OHARE_LON = -87.90902412382081
NM_TO_KM = 1.852
FEET_TO_METERS = 0.3048
JOLIET_LAT = 41.52519
ARLINGTON_HEIGHTS_LAT = 42.08836
PINGREE_GROVE_LON = -88.413416


def resolve_trino_path(trino_path: Path) -> Path:
    # Accept either a direct JAR path or a containing folder so local setup can
    # stay flexible across machines.
    if trino_path.is_file():
        return trino_path

    if trino_path.is_dir():
        jar_candidates = sorted(
            candidate
            for candidate in trino_path.iterdir()
            if candidate.is_file() and candidate.suffix.lower() == ".jar"
        )
        if len(jar_candidates) == 1:
            return jar_candidates[0]
        if len(jar_candidates) > 1:
            preferred = [candidate for candidate in jar_candidates if "trino" in candidate.name.lower()]
            if len(preferred) == 1:
                return preferred[0]
            names = ", ".join(candidate.name for candidate in jar_candidates)
            raise SystemExit(
                f"Multiple JAR files found under {trino_path}. Pass --trino-path explicitly. Candidates: {names}"
            )

    if not trino_path.suffix and trino_path.parent.is_dir():
        jar_candidates = sorted(
            candidate
            for candidate in trino_path.parent.iterdir()
            if candidate.is_file() and candidate.suffix.lower() == ".jar" and "trino" in candidate.name.lower()
        )
        if len(jar_candidates) == 1:
            return jar_candidates[0]

    raise SystemExit(
        "Trino CLI JAR not found. "
        f"Checked {trino_path}. Pass --trino-path to a local Trino CLI JAR such as C:\\tools\\trino-cli.jar."
    )


def get_day_window(date_text: str, timezone_name: str) -> tuple[int, int]:
    # Interpret the requested date in local time first, then hand UTC epochs to
    # Trino so the SQL stays unambiguous.
    timezone = ZoneInfo(timezone_name)
    day_start = datetime.strptime(date_text, "%Y-%m-%d").replace(tzinfo=timezone)
    day_end = day_start + timedelta(days=1)
    return int(day_start.timestamp()), int(day_end.timestamp())


def get_hour_window_utc(start_epoch: int, end_epoch: int) -> tuple[str, str]:
    # Trino partition pruning works at the UTC hour level, so round the day
    # window outward to matching hour buckets.
    start_dt = datetime.fromtimestamp(start_epoch, tz=UTC).replace(minute=0, second=0, microsecond=0)
    end_dt = datetime.fromtimestamp(end_epoch, tz=UTC).replace(minute=0, second=0, microsecond=0)
    start_hour = int(start_dt.timestamp())
    end_hour = int(end_dt.timestamp())
    return start_hour, end_hour


def get_bounding_box(center_lat: float, center_lon: float, radius_nm: float) -> tuple[float, float, float, float]:
    radius_km = radius_nm * NM_TO_KM
    lat_change = radius_km / 110.574
    lon_change = radius_km / (111.320 * math.cos(math.radians(center_lat)))
    south = center_lat - lat_change
    north = center_lat + lat_change
    west = center_lon - lon_change
    east = center_lon + lon_change
    return south, north, west, east


def build_query(template_path: Path, values: dict[str, object]) -> str:
    template_text = template_path.read_text(encoding="utf-8")
    return Template(template_text).substitute(values)


def get_unique_output_path(output_path: Path) -> Path:
    if not output_path.exists():
        return output_path

    for index in range(1, 1000):
        candidate = output_path.with_name(f"{output_path.stem}_{index}{output_path.suffix}")
        if not candidate.exists():
            return candidate

    raise SystemExit(f"Could not find an unused output filename near {output_path}.")


def run_trino_query(
    trino_path: Path,
    user_name: str,
    server: str,
    catalog: str,
    schema: str,
    query_text: str,
    output_path: Path,
) -> None:
    resolved_trino_path = resolve_trino_path(trino_path)
    # Shell out to the Trino CLI rather than bundling a client dependency into
    # the repo.
    command = [
        "java",
        "-jar",
        str(resolved_trino_path),
        "--server",
        server,
        "--catalog",
        catalog,
        "--schema",
        schema,
        "--external-authentication",
        "--user",
        user_name,
        "--output-format",
        "CSV_HEADER",
        "--execute",
        query_text,
    ]

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise SystemExit("Java was not found on PATH. Install Java or make sure `java` is available in this shell.") from exc
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else ""
        stdout = exc.stdout.strip() if exc.stdout else ""
        details = stderr or stdout or "Trino CLI exited without output."
        raise SystemExit(f"Trino query failed.\nCommand: {' '.join(command)}\n\n{details}") from exc

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(result.stdout, encoding="utf-8")

    if result.stderr:
        print(result.stderr.strip())


def parse_args() -> argparse.Namespace:
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Export OpenSky flight data with Trino.")
    parser.add_argument("--date", default="2019-03-09", help="Local date in YYYY-MM-DD format.")
    parser.add_argument("--timezone", default="America/Chicago", help="Timezone for the local day.")
    parser.add_argument("--user", default=os.environ.get("TRINO_USER"), help="Trino user name.")
    parser.add_argument(
        "--trino-path",
        type=Path,
        default=Path.home() / "Downloads" / "Tools" / "trino",
        help="Path to the Trino CLI JAR file or a folder containing it.",
    )
    parser.add_argument(
        "--query-file",
        dest="query_file",
        type=Path,
        default=folder / "query.sql",
        help="SQL template file.",
    )
    parser.add_argument("--server", default="https://trino.opensky-network.org", help="Trino server URL.")
    parser.add_argument("--catalog", default="minio", help="Trino catalog name.")
    parser.add_argument("--schema", default="osky", help="Trino schema name.")
    parser.add_argument("--radius-nm", type=float, default=15.0, help="Bounding box radius around O'Hare in nautical miles.")
    parser.add_argument(
        "--south",
        type=float,
        default=JOLIET_LAT,
        help="South latitude bound. Default is Joliet, IL.",
    )
    parser.add_argument(
        "--north",
        type=float,
        default=ARLINGTON_HEIGHTS_LAT,
        help="North latitude bound. Default is Arlington Heights, IL.",
    )
    parser.add_argument(
        "--west",
        type=float,
        default=PINGREE_GROVE_LON,
        help="West longitude bound. Default is Pingree Grove, IL.",
    )
    parser.add_argument("--east", type=float, default=None, help="Override east longitude bound.")
    parser.add_argument("--altitude-max-ft", type=float, default=32808.4, help="Maximum barometric altitude in feet MSL.")
    parser.add_argument("--sample-seconds", type=int, default=1, help="Keep one row every N seconds.")
    parser.add_argument("--output", type=Path, default=None, help="Output CSV file path.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.user:
        raise SystemExit("Set TRINO_USER or pass --user before running this script.")

    start_epoch, end_epoch = get_day_window(args.date, args.timezone)
    start_hour_utc, end_hour_utc = get_hour_window_utc(start_epoch, end_epoch)
    south, north, west, east = get_bounding_box(OHARE_LAT, OHARE_LON, args.radius_nm)
    south = args.south if args.south is not None else south
    north = args.north if args.north is not None else north
    west = args.west if args.west is not None else west
    east = args.east if args.east is not None else east

    # Materialize one concrete SQL query after all time, altitude, and bbox
    # parameters have been resolved from CLI input.
    query_text = build_query(
        args.query_file,
        {
            "start_epoch": start_epoch,
            "end_epoch": end_epoch,
            "start_hour_utc": start_hour_utc,
            "end_hour_utc": end_hour_utc,
            "south": south,
            "north": north,
            "west": west,
            "east": east,
            "altitude_max_m": args.altitude_max_ft * FEET_TO_METERS,
            "sample_seconds": args.sample_seconds,
        },
    )

    output_file = args.output
    if output_file is None:
        file_name = f"ohare_{args.date}_local_{args.sample_seconds}s_{int(args.radius_nm)}nm_bbox.csv"
        output_file = Path(__file__).resolve().parent / "output" / file_name
    output_file = get_unique_output_path(output_file)

    run_trino_query(
        trino_path=args.trino_path,
        user_name=args.user,
        server=args.server,
        catalog=args.catalog,
        schema=args.schema,
        query_text=query_text,
        output_path=output_file,
    )

    print(f"Saved CSV: {output_file}")
    print(f"UTC epoch window: {start_epoch} to {end_epoch}")
    print(f"Bounding box: lat {south} to {north}, lon {west} to {east}")


if __name__ == "__main__":
    main()
