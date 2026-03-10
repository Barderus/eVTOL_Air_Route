from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


OHARE = (-87.90902412382081, 41.97807408541273)


def detect_csv_encoding(csv_path: Path) -> str:
    with csv_path.open("rb") as handle:
        start = handle.read(4)

    if start.startswith(b"\xff\xfe") or start.startswith(b"\xfe\xff"):
        return "utf-16"
    if start.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8"


def parse_args() -> argparse.Namespace:
    project_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Render a 2D density map from an OpenSky CSV export.")
    parser.add_argument(
        "csv_path",
        nargs="?",
        default=project_dir / "ohare_2019-03-09_local_1min_15nm_bbox.csv",
        type=Path,
        help="Path to the OpenSky CSV export.",
    )
    parser.add_argument("--bins", type=int, default=180, help="Number of bins per axis.")
    parser.add_argument("--output", type=Path, default=project_dir / "ohare_density_map.png", help="Output image path.")
    parser.add_argument("--show", action="store_true", help="Display the figure interactively.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    encoding = detect_csv_encoding(args.csv_path)
    df = pd.read_csv(args.csv_path, encoding=encoding)
    df = df.dropna(subset=["lat", "lon", "baroaltitude"]).copy()
    df["datetime_utc"] = pd.to_datetime(df["time"], unit="s", utc=True)

    fig, ax = plt.subplots(figsize=(9, 9))
    hist = ax.hist2d(df["lon"], df["lat"], bins=args.bins, cmap="hot")

    ax.scatter(*OHARE, color="cyan", edgecolors="black", linewidths=0.8, s=60, label="O'Hare")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("OpenSky 2D Density Map Around O'Hare")
    ax.legend(loc="upper right")
    fig.colorbar(hist[3], ax=ax, label="Aircraft observations")
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)

    print(f"Saved density map: {args.output}")
    print(f"Rows plotted: {len(df)}")
    print(f"CSV encoding: {encoding}")

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
