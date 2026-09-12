"""Build St. Louis barycentric tetrahedrons for direct Jaccard routes."""

import sys
from pathlib import Path

import pandas as pd


SCRIPT_FOLDER = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_FOLDER))

import make_weight_barycentric_plots as barycentric


ROUTE_ROBUSTNESS = Path("St Louis") / "route_robustness"
OUTPUT_FOLDER = ROUTE_ROBUSTNESS / "output"
WEIGHTED_ROUTES_CSV = OUTPUT_FOLDER / "weighted_routes" / "st_louis_weighted_route_runs.csv"
DIRECT_ROUTES_FOLDER = OUTPUT_FOLDER / "direct_routes"
BARYCENTRIC_OUTPUT_FOLDER = OUTPUT_FOLDER / "barycentric"
MAPS_FOLDER = ROUTE_ROBUSTNESS / "maps" / "barycentric"
CLUSTER_COLUMN = "hierarchical_jaccard_cluster"
METHOD = {
    "label": "Hierarchical On Jaccard (Direct Routes)",
    "cluster_column": CLUSTER_COLUMN,
    "suffix": "hierarchical_jaccard_direct_routes",
}
ROUTE_PAIRS = [
    "midamerica_to_st_louis_lambert",
    "midamerica_to_st_louis_union_station",
    "st_louis_downtown_airport_to_st_louis_lambert",
]


def load_direct_route_runs(route_pair):
    """Load full route records selected by the direct Jaccard subset."""
    direct_path = DIRECT_ROUTES_FOLDER / (
        f"{route_pair}_hierarchical_jaccard.csv"
    )
    direct_routes = pd.read_csv(direct_path)
    all_routes = pd.read_csv(WEIGHTED_ROUTES_CSV)
    selected = all_routes[all_routes["route_run_id"].isin(direct_routes["route_run_id"])].copy()
    if len(selected) != len(direct_routes):
        raise ValueError(f"Could not match every direct route for {route_pair}")
    selected = selected.merge(
        direct_routes[["route_run_id", CLUSTER_COLUMN]],
        on="route_run_id",
        how="left",
    )
    return selected


def main():
    """Save direct-route barycentric coordinates and tetrahedron pages."""
    if not WEIGHTED_ROUTES_CSV.exists():
        raise FileNotFoundError(f"Weighted routes not found: {WEIGHTED_ROUTES_CSV}")
    MAPS_FOLDER.mkdir(parents=True, exist_ok=True)
    BARYCENTRIC_OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    route_frames = [load_direct_route_runs(route_pair) for route_pair in ROUTE_PAIRS]
    route_runs = pd.concat(route_frames, ignore_index=True)
    route_runs = route_runs[route_runs["status"] == "success"].copy()
    route_runs = barycentric.add_barycentric_coordinates(route_runs)
    route_runs = barycentric.add_route_variant_ids(route_runs)

    output_columns = [
        "route_run_id",
        "route_pair",
        "route_pair_label",
        "route_variant",
        "traffic_date",
        "weight_id",
        *barycentric.WEIGHT_COLUMNS,
        "barycentric_x",
        "barycentric_y",
        "barycentric_z",
        CLUSTER_COLUMN,
    ]
    output_csv = BARYCENTRIC_OUTPUT_FOLDER / "st_louis_direct_jaccard_barycentric_coordinates.csv"
    route_runs[output_columns].to_csv(output_csv, index=False)
    print(f"Saved direct-route barycentric coordinates: {output_csv}")

    for route_pair in ROUTE_PAIRS:
        direct_frame = route_runs[route_runs["route_pair"] == route_pair]
        variant_output = MAPS_FOLDER / f"{route_pair}_direct_routes_weight_barycentric.html"
        barycentric.make_route_pair_plot(direct_frame, route_pair, variant_output)
        cluster_output = MAPS_FOLDER / (
            f"{route_pair}_direct_routes_weight_barycentric_hierarchical_jaccard.html"
        )
        barycentric.make_cluster_plot(
            route_runs,
            route_pair,
            METHOD,
            cluster_output,
        )
        cluster_count = direct_frame[CLUSTER_COLUMN].nunique()
        print(
            f"Saved direct tetrahedrons for {route_pair}: "
            f"{len(direct_frame)} weights, {cluster_count} Jaccard clusters"
        )


if __name__ == "__main__":
    main()
