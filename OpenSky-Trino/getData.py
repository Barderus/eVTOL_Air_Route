from __future__ import annotations

import argparse
import math
import os
import subprocess
from datetime import datetime, timedelta
from pathlib import Path
from string import Template
from zoneinfo import ZoneInfo


OHARE_LAT = 41.97807408541273
OHARE_LON = -87.90902412382081
NM_TO_KM = 1.852


def local_day_to_utc_epoch(date_str: str, timezone_name: str) -> tuple[int, int]:
    tz = ZoneInfo(timezone_name)
    day_start = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=tz)
    day_end = day_start + timedelta(days=1)
    return int(day_start.timestamp()), int(day_end.timestamp())


def bounding_box(lat: float, lon: float, radius_nm: float) -> tuple[float, float, float, float]:
    radius_km = radius_nm * NM_TO_KM
    lat_delta = radius_km / 110.574
    lon_delta = radius_km / (111.320 * math.cos(math.radians(lat)))
    south = lat - lat_delta
    north = lat + lat_delta
    west = lon - lon_delta
    east = lon + lon_delta
    return south, north, west, east


def render_query(template_path: Path, params: dict[str, object]) -> str:
    template = Template(template_path.read_text(encoding="utf-8"))
    return template.substitute(params)


def run_export(
    trino_path: Path,
    user: str,
    query: str,
    output_path: Path,
    server: str,
    catalog: str,
    schema: str,
) -> None:
    command = [
        "java",
        "-jar",
        str(trino_path),
        "--server",
        server,
        "--catalog",
        catalog,
        "--schema",
        schema,
        "--external-authentication",
        "--user",
        user,
        "--output-format",
        "CSV_HEADER",
        "--execute",
        query,
    ]

    result = subprocess.run(command, check=True, capture_output=True, text=True)
    output_path.write_text(result.stdout, encoding="utf-8")
    if result.stderr:
        print(result.stderr.strip())


def parse_args() -> argparse.Namespace:
    project_dir = Path(__file__).resolve().parent
    default_trino = Path.home() / "Downloads" / "Tools" / "trino"
    default_user = os.environ.get("TRINO_USER")

    parser = argparse.ArgumentParser(description="Export OpenSky flight data with Trino.")
    parser.add_argument("--date", default="2019-03-09", help="Local date in YYYY-MM-DD.")
    parser.add_argument("--timezone", default="America/Chicago", help="Timezone for the day window.")
    parser.add_argument(
        "--user",
        default=default_user,
        help="Trino user name. Defaults to the TRINO_USER environment variable if set.",
    )
    parser.add_argument("--trino-path", type=Path, default=default_trino, help="Path to the Trino CLI JAR.")
    parser.add_argument("--query-template", type=Path, default=project_dir / "query.sql", help="SQL template file.")
    parser.add_argument("--output", type=Path, default=None, help="Output CSV path.")
    parser.add_argument("--server", default="https://trino.opensky-network.org", help="Trino server URL.")
    parser.add_argument("--catalog", default="minio", help="Trino catalog.")
    parser.add_argument("--schema", default="osky", help="Trino schema.")
    parser.add_argument("--radius-nm", type=float, default=15.0, help="Bounding-box radius around O'Hare in nautical miles.")
    parser.add_argument("--altitude-max-m", type=float, default=10000.0, help="Upper barometric altitude filter in meters.")
    parser.add_argument("--sample-seconds", type=int, default=60, help="Time sampling interval in seconds.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.user:
        raise SystemExit("Provide --user or set the TRINO_USER environment variable.")

    start_epoch, end_epoch = local_day_to_utc_epoch(args.date, args.timezone)
    south, north, west, east = bounding_box(OHARE_LAT, OHARE_LON, args.radius_nm)

    params = {
        "start_epoch": start_epoch,
        "end_epoch": end_epoch,
        "south": south,
        "north": north,
        "west": west,
        "east": east,
        "altitude_max_m": args.altitude_max_m,
        "sample_seconds": args.sample_seconds,
    }
    query = render_query(args.query_template, params)

    output_path = args.output
    if output_path is None:
        output_path = Path(__file__).resolve().parent / f"ohare_{args.date}_local_{args.sample_seconds}s_{int(args.radius_nm)}nm_bbox.csv"

    run_export(
        trino_path=args.trino_path,
        user=args.user,
        query=query,
        output_path=output_path,
        server=args.server,
        catalog=args.catalog,
        schema=args.schema,
    )

    print(f"Saved CSV: {output_path}")
    print(f"UTC epoch window: {start_epoch} to {end_epoch}")
    print(f"Bounding box: lat {south} to {north}, lon {west} to {east}")


if __name__ == "__main__":
    main()
