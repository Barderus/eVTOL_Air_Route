"""Extract user-selected direct route clusters for each St. Louis OD pair."""

import json
from pathlib import Path

import pandas as pd


DATE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_FOLDER = DATE_ROOT / "output"
DIRECT_ROUTES_FOLDER = OUTPUT_FOLDER / "direct_routes"
ROUTES_GEOJSON = OUTPUT_FOLDER / "weighted_routes" / "st_louis_weighted_routes.geojson"

METHODS = {
    "hierarchical_jaccard": {
        "input": OUTPUT_FOLDER / "cluster_assignments" / "st_louis_route_clusters_hierarchical_jaccard.csv",
        "cluster_column": "hierarchical_jaccard_cluster",
        "cluster_prefix": "hier_jaccard",
        "distance_method": "jaccard",
    },
}

DIRECT_CLUSTER_IDS = {
    "midamerica_to_st_louis_lambert": {
        "hierarchical_jaccard": [1, 2, 4, 5, 6, 8, 18],
    },
    "midamerica_to_st_louis_union_station": {
        "hierarchical_jaccard": [2, 3, 6, 7, 8, 12, 17],
    },
    "st_louis_downtown_airport_to_st_louis_lambert": {
        "hierarchical_jaccard": [1, 2, 3, 4, 5],
    },
}


def cluster_label(prefix, cluster_number):
    """Build the stored cluster label from its numeric suffix."""
    return f"{prefix}_{cluster_number:03d}"


def load_features():
    """Load route features indexed by their route run ID."""
    with ROUTES_GEOJSON.open("r", encoding="utf-8") as file_handle:
        geojson = json.load(file_handle)
    return {
        feature["properties"]["route_run_id"]: feature
        for feature in geojson["features"]
    }


def extract_method(route_pair, method_key, method, features):
    """Write one selected-cluster CSV and GeoJSON subset."""
    data = pd.read_csv(method["input"])
    pair_data = data[data["route_pair"] == route_pair].copy()
    selected_numbers = DIRECT_CLUSTER_IDS[route_pair][method_key]
    selected_labels = {
        cluster_label(method["cluster_prefix"], number)
        for number in selected_numbers
    }

    available_labels = set(pair_data[method["cluster_column"]])
    missing_labels = selected_labels - available_labels
    if missing_labels:
        missing_text = ", ".join(sorted(missing_labels))
        raise ValueError(f"{method_key} / {route_pair} has missing clusters: {missing_text}")

    selected = pair_data[pair_data[method["cluster_column"]].isin(selected_labels)].copy()
    selected.insert(
        selected.columns.get_loc(method["cluster_column"]),
        "direct_cluster_number",
        selected[method["cluster_column"]].str.extract(r"_(\d+)$")[0].astype(int),
    )
    selected.insert(0, "clustering_method", method_key)
    selected.insert(1, "distance_method", method["distance_method"])

    stem = f"{route_pair}_{method_key}"
    selected.to_csv(DIRECT_ROUTES_FOLDER / f"{stem}.csv", index=False)

    selected_ids = set(selected["route_run_id"])
    selected_features = [features[route_id] for route_id in selected_ids]
    geojson = {"type": "FeatureCollection", "features": selected_features}
    with (DIRECT_ROUTES_FOLDER / f"{stem}.geojson").open("w", encoding="utf-8") as file_handle:
        json.dump(geojson, file_handle)

    cluster_counts = selected.groupby(method["cluster_column"]).size()
    return {
        "route_pair": route_pair,
        "clustering_method": method_key,
        "distance_method": method["distance_method"],
        "selected_cluster_numbers": ", ".join(f"{number:02d}" for number in selected_numbers),
        "selected_cluster_count": len(selected_numbers),
        "selected_route_count": len(selected),
        "cluster_sizes": "; ".join(
            f"{label.removeprefix(method['cluster_prefix'] + '_')}={int(size)}"
            for label, size in cluster_counts.items()
        ),
        "csv_file": f"{stem}.csv",
        "geojson_file": f"{stem}.geojson",
    }


def main():
    """Extract every configured direct-route subset."""
    DIRECT_ROUTES_FOLDER.mkdir(parents=True, exist_ok=True)
    features = load_features()
    manifest_rows = []
    for route_pair, method_clusters in DIRECT_CLUSTER_IDS.items():
        for method_key, method in METHODS.items():
            manifest_rows.append(
                extract_method(route_pair, method_key, method, features)
            )

    pd.DataFrame(manifest_rows).to_csv(
        DIRECT_ROUTES_FOLDER / "direct_routes_manifest.csv", index=False
    )
    print(f"Saved {len(manifest_rows)} direct-route subsets to {DIRECT_ROUTES_FOLDER}")


if __name__ == "__main__":
    main()
