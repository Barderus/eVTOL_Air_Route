"""Export St. Louis OpenSky traffic data through Trino."""

import math
import os
import subprocess
from datetime import UTC, datetime, timedelta
from string import Template
from zoneinfo import ZoneInfo

import settings


EXPORT_DATES = [
    "2026-03-07",
    "2026-03-09",
    "2026-01-10",
    "2026-01-12",
    "2025-07-12",
    "2025-07-14",
]
SAMPLE_SECONDS = 1
ALTITUDE_MAX_FT = 5000
OVERWRITE_EXISTING = True
TRINO_USER = os.environ.get("TRINO_USER")
TRINO_JAR = r"C:\Users\Barderus_Legion\Downloads\trino-cli-480.jar"
TRINO_SERVER = "https://trino.opensky-network.org"
TRINO_CATALOG = "minio"
TRINO_SCHEMA = "osky"
QUERY_FILE = "traffic/query.sql"

NM_TO_KM = 1.852
FEET_TO_METERS = 0.3048


def get_day_window(export_date):
    """Return the configured local day as UTC epoch seconds."""
    timezone = ZoneInfo(settings.LOCAL_TIMEZONE)
    day_start = datetime.strptime(export_date, "%Y-%m-%d").replace(tzinfo=timezone)
    day_end = day_start + timedelta(days=1)
    return int(day_start.timestamp()), int(day_end.timestamp())


def get_hour_window(start_epoch, end_epoch):
    """Return UTC hour buckets for Trino partition filtering."""
    start_time = datetime.fromtimestamp(start_epoch, tz=UTC)
    end_time = datetime.fromtimestamp(end_epoch, tz=UTC)
    start_hour = start_time.replace(minute=0, second=0, microsecond=0)
    end_hour = end_time.replace(minute=0, second=0, microsecond=0)
    return int(start_hour.timestamp()), int(end_hour.timestamp())


def get_bounding_box():
    """Return the configured or radius-based traffic bounding box."""
    radius_km = settings.TRAFFIC_RADIUS_NM * NM_TO_KM
    latitude_change = radius_km / 110.574
    longitude_change = radius_km / (
        111.320 * math.cos(math.radians(settings.TRAFFIC_CENTER_LAT))
    )

    south = settings.TRAFFIC_CENTER_LAT - latitude_change
    north = settings.TRAFFIC_CENTER_LAT + latitude_change
    west = settings.TRAFFIC_CENTER_LON - longitude_change
    east = settings.TRAFFIC_CENTER_LON + longitude_change

    if settings.TRAFFIC_SOUTH is not None:
        south = settings.TRAFFIC_SOUTH
    if settings.TRAFFIC_NORTH is not None:
        north = settings.TRAFFIC_NORTH
    if settings.TRAFFIC_WEST is not None:
        west = settings.TRAFFIC_WEST
    if settings.TRAFFIC_EAST is not None:
        east = settings.TRAFFIC_EAST

    return south, north, west, east


def build_query(export_date):
    """Fill the St. Louis SQL template with the configured values."""
    start_epoch, end_epoch = get_day_window(export_date)
    start_hour, end_hour = get_hour_window(start_epoch, end_epoch)
    south, north, west, east = get_bounding_box()

    with open(QUERY_FILE, "r", encoding="utf-8") as query_file:
        query_template = Template(query_file.read())

    return query_template.substitute(
        start_epoch=start_epoch,
        end_epoch=end_epoch,
        start_hour_utc=start_hour,
        end_hour_utc=end_hour,
        south=south,
        north=north,
        west=west,
        east=east,
        altitude_max_m=ALTITUDE_MAX_FT * FEET_TO_METERS,
        sample_seconds=SAMPLE_SECONDS,
    )


def main():
    settings.validate_traffic_settings()

    if not TRINO_USER:
        raise ValueError("Set the TRINO_USER environment variable before exporting data.")
    if not os.path.isfile(TRINO_JAR):
        raise FileNotFoundError(f"Trino CLI JAR not found: {TRINO_JAR}")

    os.makedirs(settings.TRAFFIC_FOLDER, exist_ok=True)

    for export_date in EXPORT_DATES:
        try:
            datetime.strptime(export_date, "%Y-%m-%d")
        except ValueError as error:
            raise ValueError(
                f"Export date must use YYYY-MM-DD format: {export_date}"
            ) from error

        output_name = (
            f"{settings.TRAFFIC_CENTER_NAME}_{export_date}_"
            f"{SAMPLE_SECONDS}s_{int(settings.TRAFFIC_RADIUS_NM)}nm_bbox.csv"
        )
        output_path = os.path.join(settings.TRAFFIC_FOLDER, output_name)

        if os.path.isfile(output_path) and not OVERWRITE_EXISTING:
            print(f"Skipping existing traffic data: {output_path}")
            continue

        query = build_query(export_date)
        command = [
            "java",
            "-jar",
            TRINO_JAR,
            "--server",
            TRINO_SERVER,
            "--catalog",
            TRINO_CATALOG,
            "--schema",
            TRINO_SCHEMA,
            "--external-authentication",
            "--user",
            TRINO_USER,
            "--output-format",
            "CSV_HEADER",
            "--execute",
            query,
        ]

        print(f"Exporting OpenSky traffic for {export_date}...")
        result = subprocess.run(command, check=True, capture_output=True, text=True)

        with open(output_path, "w", encoding="utf-8") as output_file:
            output_file.write(result.stdout)

        print(f"Saved traffic data: {output_path}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error
