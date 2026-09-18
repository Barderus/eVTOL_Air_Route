"""Build Leaflet maps from the selected direct-route cluster subsets."""

import sys
from pathlib import Path

import pandas as pd


SCRIPT_FOLDER = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_FOLDER))

import make_route_cluster_maps as route_maps


DATE_ROOT = Path(__file__).resolve().parents[1]
DIRECT_ROUTES_FOLDER = DATE_ROOT / "output" / "direct_routes"

METHODS = [
    {
        "key": "hierarchical_jaccard",
        "label": "Hierarchical + Jaccard",
        "file_suffix": "hierarchical_jaccard",
        "cluster_column": "hierarchical_jaccard_cluster",
    },
]

ROUTE_PAIRS = route_maps.ROUTE_PAIR_ORDER


def load_direct_cluster_tables(route_pair):
    """Load the four direct-route tables for one route pair."""
    tables = {}
    for method in METHODS:
        path = DIRECT_ROUTES_FOLDER / f"{route_pair}_{method['file_suffix']}.csv"
        table = pd.read_csv(path)
        route_maps.require_columns(
            table,
            ["route_run_id", "route_pair", method["cluster_column"]],
            str(path),
        )
        tables[method["key"]] = table
    return tables


def direct_route_ids(tables):
    """Return every selected direct-route id across the active direct methods."""
    route_ids = set()
    for table in tables.values():
        route_ids.update(table["route_run_id"])
    return route_ids


def main():
    """Create one direct-route map for each St. Louis route pair."""
    if not (DIRECT_ROUTES_FOLDER / "direct_routes_manifest.csv").exists():
        raise FileNotFoundError(
            "Run extract_direct_routes.py before creating direct-route maps."
        )

    features_by_id = route_maps.load_route_features()
    DIRECT_ROUTES_FOLDER.mkdir(parents=True, exist_ok=True)

    original_methods = route_maps.METHODS
    route_maps.METHODS = METHODS
    try:
        for route_pair in ROUTE_PAIRS:
            tables = load_direct_cluster_tables(route_pair)
            payload = route_maps.build_route_pair_payload(
                route_pair, tables, features_by_id
            )
            payload = route_maps.add_background_routes(
                payload,
                features_by_id,
                route_pair,
                direct_route_ids(tables),
            )
            output_path = DATE_ROOT / "maps" / "direct_routes" / f"{route_pair}_direct_routes.html"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(
                route_maps.html_template(
                    payload, title_prefix="St. Louis Direct Route Clusters"
                ),
                encoding="utf-8",
            )
            print(f"Saved direct-route map: {output_path}")
    finally:
        route_maps.METHODS = original_methods


if __name__ == "__main__":
    main()
