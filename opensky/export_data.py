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


def get_day_window(date_text: str, timezone_name: str) -> tuple[int, int]:
    timezone = ZoneInfo(timezone_name)
    day_start = datetime.strptime(date_text, "%Y-%m-%d").replace(tzinfo=timezone)
    day_end = day_start + timedelta(days=1)
    return int(day_start.timestamp()), int(day_end.timestamp())


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


def run_trino_query(
    trino_path: Path,
    user_name: str,
    server: str,
    catalog: str,
    schema: str,
    query_text: str,
    output_path: Path,
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
        user_name,
        "--output-format",
        "CSV_HEADER",
        "--execute",
        query_text,
    ]

    result = subprocess.run(command, check=True, capture_output=True, text=True)
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
        help="Path to the Trino CLI jar file.",
    )
    parser.add_argument("--query-file", type=Path, default=folder / "query.sql", help="SQL template file.")
    parser.add_argument("--server", default="https://trino.opensky-network.org", help="Trino server URL.")
    parser.add_argument("--catalog", default="minio", help="Trino catalog name.")
    parser.add_argument("--schema", default="osky", help="Trino schema name.")
    parser.add_argument("--radius-nm", type=float, default=15.0, help="Bounding box radius around O'Hare in nautical miles.")
    parser.add_argument("--altitude-max-m", type=float, default=10000.0, help="Maximum barometric altitude in meters.")
    parser.add_argument("--sample-seconds", type=int, default=60, help="Keep one row every N seconds.")
    parser.add_argument("--output", type=Path, default=None, help="Output CSV file path.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.user:
        raise SystemExit("Set TRINO_USER or pass --user before running this script.")

    start_epoch, end_epoch = get_day_window(args.date, args.timezone)
    south, north, west, east = get_bounding_box(OHARE_LAT, OHARE_LON, args.radius_nm)

    query_text = build_query(
        args.query_file,
        {
            "start_epoch": start_epoch,
            "end_epoch": end_epoch,
            "south": south,
            "north": north,
            "west": west,
            "east": east,
            "altitude_max_m": args.altitude_max_m,
            "sample_seconds": args.sample_seconds,
        },
    )

    output_file = args.output
    if output_file is None:
        file_name = f"ohare_{args.date}_local_{args.sample_seconds}s_{int(args.radius_nm)}nm_bbox.csv"
        output_file = Path(__file__).resolve().parent / "output" / file_name

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
