from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


OHARE_LON = -87.90902412382081
OHARE_LAT = 41.97807408541273
MIDWAY_LON = -87.75331135429448
MIDWAY_LAT = 41.7856116663475


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
    data = data.dropna(subset=["lat", "lon", "baroaltitude"]).copy()
    data["time"] = pd.to_numeric(data["time"], errors="coerce")
    data["lat"] = pd.to_numeric(data["lat"], errors="coerce")
    data["lon"] = pd.to_numeric(data["lon"], errors="coerce")
    data["baroaltitude"] = pd.to_numeric(data["baroaltitude"], errors="coerce")
    data = data.dropna(subset=["time", "lat", "lon", "baroaltitude"])
    return data


def parse_args() -> argparse.Namespace:
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Create a 2D density plot from an OpenSky CSV file.")
    parser.add_argument(
        "csv_path",
        nargs="?",
        default=folder.parent / "opensky" / "output" / "ohare_2019-03-09_local_1s_15nm_bbox.csv",
        type=Path,
        help="CSV file to read.",
    )
    parser.add_argument("--bins", type=int, default=180, help="Number of bins for the histogram.")
    parser.add_argument(
        "--output",
        type=Path,
        default=folder.parent / "images" / "visualizations" / "ohare_density_map.png",
        help="Output image file.",
    )
    parser.add_argument("--show", action="store_true", help="Show the figure after saving it.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_flight_data(args.csv_path)

    figure, axis = plt.subplots(figsize=(9, 9))
    histogram = axis.hist2d(data["lon"], data["lat"], bins=args.bins, cmap="hot")

    axis.scatter(OHARE_LON, OHARE_LAT, color="cyan", edgecolors="black", linewidths=0.8, s=60, label="O'Hare")
    axis.scatter(MIDWAY_LON, MIDWAY_LAT, color="orange", edgecolors="black", linewidths=0.8, s=60, label="Midway")
    axis.set_xlabel("Longitude")
    axis.set_ylabel("Latitude")
    axis.set_title("OpenSky 2D Density Map Around O'Hare")
    axis.legend(loc="upper right")
    figure.colorbar(histogram[3], ax=axis, label="Aircraft observations")
    figure.tight_layout()
    figure.savefig(args.output, dpi=200)

    print(f"Saved density map: {args.output}")
    print(f"Rows plotted: {len(data)}")

    if args.show:
        plt.show()
    else:
        plt.close(figure)


if __name__ == "__main__":
    main()
