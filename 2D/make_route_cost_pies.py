from __future__ import annotations

import math
import re
import runpy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
GENERATE_SCRIPT = ROOT / "2D" / "generate_astar_toggle_pages.py"
OPENSKY_OUTPUT = ROOT / "opensky" / "output"
OUTPUT_DIR = ROOT / "images" / "visualizations"

COMPONENT_COLORS = {
    "Distance": "#3b82f6",
    "Population": "#10b981",
    "Airspace": "#f59e0b",
    "Traffic": "#8b5cf6",
}


def detect_csv_encoding(csv_path: Path) -> str:
    with csv_path.open("rb") as file_handle:
        first_bytes = file_handle.read(4)

    if first_bytes.startswith(b"\xff\xfe") or first_bytes.startswith(b"\xfe\xff"):
        return "utf-16"
    if first_bytes.startswith(b"\xef\xbb\xbf"):
        return "utf-8-sig"
    return "utf-8"


def parse_dataset_date(csv_path: Path) -> pd.Timestamp | None:
    match = re.search(r"ohare_(\d{4}-\d{2}-\d{2})_", csv_path.name)
    if not match:
        return None
    return pd.Timestamp(match.group(1))


def latest_opensky_csv() -> Path:
    candidates = list(OPENSKY_OUTPUT.glob("ohare_*_1s_15nm_bbox.csv"))
    if not candidates:
        raise FileNotFoundError(f"No OpenSky CSV files found in {OPENSKY_OUTPUT}")

    dated_candidates: list[tuple[pd.Timestamp, Path]] = []
    fallback_candidates: list[tuple[float, Path]] = []
    for candidate in candidates:
        parsed = parse_dataset_date(candidate)
        if parsed is not None:
            dated_candidates.append((parsed, candidate))
        else:
            fallback_candidates.append((candidate.stat().st_mtime, candidate))

    if dated_candidates:
        return max(dated_candidates, key=lambda item: item[0])[1]
    return max(fallback_candidates, key=lambda item: item[0])[1]


def load_helpers() -> dict[str, object]:
    return runpy.run_path(str(GENERATE_SCRIPT))


def compute_combined_route_breakdown(
    graph: nx.Graph,
    path: list[int],
    cent_x: np.ndarray,
    cent_y: np.ndarray,
    city_risk_norm: np.ndarray,
    airspace_risk_norm: np.ndarray,
    traffic_risk_norm: np.ndarray,
    max_edge_distance: float,
    weights: dict[str, float],
) -> dict[str, float]:
    contributions = {name: 0.0 for name in COMPONENT_COLORS}

    for u, v in zip(path, path[1:]):
        dx = cent_x[u] - cent_x[v]
        dy = cent_y[u] - cent_y[v]
        distance_norm = math.hypot(dx, dy) / max_edge_distance
        population_norm = 0.5 * (city_risk_norm[u] + city_risk_norm[v])
        airspace_norm = 0.5 * (airspace_risk_norm[u] + airspace_risk_norm[v])
        traffic_norm = 0.5 * (traffic_risk_norm[u] + traffic_risk_norm[v])

        contributions["Distance"] += weights["distance"] * distance_norm
        contributions["Population"] += weights["population"] * population_norm
        contributions["Airspace"] += weights["airspace"] * airspace_norm
        contributions["Traffic"] += weights["traffic"] * traffic_norm

    return contributions


def save_pie_chart(
    destination_label: str,
    dataset_label: str,
    contributions: dict[str, float],
    total_cost: float,
    output_path: Path,
) -> None:
    labels = list(contributions.keys())
    values = np.array([contributions[label] for label in labels], dtype=float)
    colors = [COMPONENT_COLORS[label] for label in labels]

    figure, axis = plt.subplots(figsize=(7.5, 7.5))
    axis.pie(
        values,
        labels=labels,
        colors=colors,
        startangle=90,
        counterclock=False,
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 1 else "",
        textprops={"fontsize": 11},
        wedgeprops={"linewidth": 1, "edgecolor": "white"},
    )
    axis.set_aspect("equal")
    axis.set_title(
        f"{destination_label} combined cost breakdown\n{dataset_label}",
        fontsize=15,
        pad=18,
    )
    figure.text(
        0.5,
        0.04,
        f"Total combined path cost: {total_cost:.4f}",
        ha="center",
        va="center",
        fontsize=11,
        color="#374151",
    )
    figure.tight_layout(rect=(0, 0.06, 1, 1))
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    helpers = load_helpers()
    GRID_PATH = helpers["GRID_PATH"]
    DESTINATIONS = helpers["DESTINATIONS"]
    START = helpers["START"]
    normalize = helpers["normalize"]
    build_graph = helpers["build_graph"]
    node_for_point = helpers["node_for_point"]
    load_traffic_counts = helpers["load_traffic_counts"]
    assign_route_edge_weights = helpers["assign_route_edge_weights"]
    make_heuristic = helpers["make_heuristic"]
    path_cost = helpers["path_cost"]
    ROUTE_SPECS = helpers["ROUTE_SPECS"]

    combined_spec = next(spec for spec in ROUTE_SPECS if spec["name"] == "combined")
    latest_csv = latest_opensky_csv()
    dataset_date = parse_dataset_date(latest_csv)
    dataset_label = dataset_date.strftime("%Y-%m-%d") if dataset_date is not None else latest_csv.stem

    grid = gpd.read_file(GRID_PATH)
    graph, centroids_m, cent_x, cent_y = build_graph(grid)
    sindex = grid.sindex
    geoms = grid.geometry

    city_risk = grid["city_risk"].to_numpy(dtype=float)
    airspace_risk = grid["airport_risk_combined"].to_numpy(dtype=float)
    city_risk_norm, _ = normalize(city_risk)
    airspace_risk_norm, _ = normalize(airspace_risk, clip_percentile=99.0, power=0.75)

    traffic_counts = load_traffic_counts(latest_csv, grid)
    traffic_risk_norm, _ = normalize(
        traffic_counts,
        clip_percentile=95.0,
        transform="log1p",
        power=0.75,
    )
    max_edge_distance = assign_route_edge_weights(
        graph,
        cent_x,
        cent_y,
        city_risk_norm,
        airspace_risk_norm,
        traffic_risk_norm,
    )

    start_node = node_for_point(START["lat"], START["lon"], sindex, geoms, centroids_m)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for destination in DESTINATIONS:
        end_node = node_for_point(destination["lat"], destination["lon"], sindex, geoms, centroids_m)
        heuristic = make_heuristic(combined_spec["distance_weight"], cent_x, cent_y, max_edge_distance)
        path = nx.astar_path(
            graph,
            start_node,
            end_node,
            heuristic=heuristic,
            weight="weight_combined",
        )
        total_cost = path_cost(graph, path, "weight_combined")
        contributions = compute_combined_route_breakdown(
            graph,
            path,
            cent_x,
            cent_y,
            city_risk_norm,
            airspace_risk_norm,
            traffic_risk_norm,
            max_edge_distance,
            {
                "distance": combined_spec["distance_weight"],
                "population": combined_spec["population_weight"],
                "airspace": combined_spec["airspace_weight"],
                "traffic": combined_spec["traffic_weight"],
            },
        )
        contribution_total = sum(contributions.values())
        if not math.isclose(total_cost, contribution_total, rel_tol=1e-9, abs_tol=1e-9):
            raise RuntimeError(
                f"Contribution sum mismatch for {destination['slug']}: "
                f"{contribution_total:.8f} vs {total_cost:.8f}"
            )

        output_path = OUTPUT_DIR / f"clow_to_{destination['slug']}_cost_breakdown_{dataset_label}.png"
        save_pie_chart(
            destination_label=destination["label"],
            dataset_label=f"Most recent OpenSky dataset: {dataset_label}",
            contributions=contributions,
            total_cost=total_cost,
            output_path=output_path,
        )
        print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
