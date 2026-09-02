"""Select one Jaccard-medoid route for each direct-route cluster."""

import json
from pathlib import Path

import pandas as pd


ROUTE_ROBUSTNESS = Path("St Louis") / "route_robustness"
DIRECT_ROUTES_FOLDER = ROUTE_ROBUSTNESS / "route_clusters" / "direct_routes"
ROUTES_GEOJSON = ROUTE_ROBUSTNESS / "output" / "st_louis_weighted_routes.geojson"

CLUSTER_COLUMN = "hierarchical_jaccard_cluster"
METHOD_SUFFIX = "hierarchical_jaccard"
ROUTE_PAIRS = [
    "midamerica_to_st_louis_lambert",
    "midamerica_to_st_louis_union_station",
    "st_louis_downtown_airport_to_st_louis_lambert",
]


def load_route_features():
    """Load full route features by route run ID."""
    with ROUTES_GEOJSON.open("r", encoding="utf-8") as file_handle:
        geojson = json.load(file_handle)
    return {
        feature["properties"]["route_run_id"]: feature
        for feature in geojson["features"]
    }


def jaccard_similarity(sequence_a, sequence_b):
    """Return shared-node Jaccard similarity for two route sequences."""
    set_a = set(sequence_a)
    set_b = set(sequence_b)
    union_size = len(set_a | set_b)
    return len(set_a & set_b) / union_size if union_size else 1.0


def select_cluster_medoid(cluster_rows, route_features):
    """Select the route with the highest mean within-cluster Jaccard similarity."""
    route_ids = cluster_rows["route_run_id"].tolist()
    if len(route_ids) == 1:
        return route_ids[0], 1.0

    sequences = {
        route_id: route_features[route_id]["properties"]["path_node_ids"]
        for route_id in route_ids
    }
    mean_similarities = {}
    for route_id in route_ids:
        similarities = [
            jaccard_similarity(sequences[route_id], sequences[other_id])
            for other_id in route_ids
        ]
        mean_similarities[route_id] = sum(similarities) / len(similarities)

    representative_id = max(
        route_ids,
        key=lambda route_id: (mean_similarities[route_id], route_id),
    )
    return representative_id, mean_similarities[representative_id]


def select_representatives(route_pair, route_features):
    """Select one representative row for each direct Jaccard cluster."""
    input_path = DIRECT_ROUTES_FOLDER / f"{route_pair}_{METHOD_SUFFIX}.csv"
    routes = pd.read_csv(input_path)
    total_routes = len(routes)
    representative_rows = []
    representative_features = []

    for cluster_id, cluster_rows in routes.groupby(CLUSTER_COLUMN, sort=True):
        representative_id, mean_similarity = select_cluster_medoid(
            cluster_rows, route_features
        )
        representative_row = cluster_rows[
            cluster_rows["route_run_id"] == representative_id
        ].iloc[0].to_dict()
        cluster_size = len(cluster_rows)
        representative_row.update(
            {
                "representative_route_run_id": representative_id,
                "representative_mean_jaccard_similarity": mean_similarity,
                "cluster_size": cluster_size,
                "cluster_weight_space_percent": 100.0 * cluster_size / total_routes,
            }
        )
        representative_rows.append(representative_row)

        feature = json.loads(json.dumps(route_features[representative_id]))
        feature["properties"].update(representative_row)
        representative_features.append(feature)

    representatives = pd.DataFrame(representative_rows).sort_values(
        ["cluster_size", CLUSTER_COLUMN], ascending=[False, True]
    )
    return representatives, {"type": "FeatureCollection", "features": representative_features}


def main():
    """Save direct-route Jaccard representatives for every OD pair."""
    DIRECT_ROUTES_FOLDER.mkdir(parents=True, exist_ok=True)
    route_features = load_route_features()
    for route_pair in ROUTE_PAIRS:
        representatives, geojson = select_representatives(route_pair, route_features)
        stem = f"{route_pair}_{METHOD_SUFFIX}_representatives"
        representatives.to_csv(DIRECT_ROUTES_FOLDER / f"{stem}.csv", index=False)
        with (DIRECT_ROUTES_FOLDER / f"{stem}.geojson").open("w", encoding="utf-8") as file_handle:
            json.dump(geojson, file_handle)
        print(
            f"Saved {len(representatives)} representatives for {route_pair}: "
            f"{DIRECT_ROUTES_FOLDER / stem}"
        )


if __name__ == "__main__":
    main()
