"""Select one Frechet-medoid route within each direct-route Jaccard cluster."""

import json
from pathlib import Path

import pandas as pd

from cluster_weighted_routes import build_frechet_distance_matrix


DATE_ROOT = Path(__file__).resolve().parents[1]
DIRECT_ROUTES_FOLDER = DATE_ROOT / "output" / "direct_routes"
ROUTES_GEOJSON = DATE_ROOT / "output" / "weighted_routes" / "st_louis_weighted_routes.geojson"

CLUSTER_COLUMN = "hierarchical_jaccard_cluster"
METHOD_SUFFIX = "hierarchical_jaccard"
FRECHET_MERGE_THRESHOLD_KM = 2.5
KM_TO_MILES = 0.621371
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


def select_cluster_medoid(cluster_rows, route_features):
    """Select the route with the lowest mean within-cluster Frechet distance."""
    route_ids = cluster_rows["route_run_id"].tolist()
    distance_matrix = build_frechet_distance_matrix(cluster_rows, route_features)
    total_distances = [sum(row) for row in distance_matrix]
    representative_index = min(
        range(len(route_ids)),
        key=lambda index: (total_distances[index], route_ids[index]),
    )
    denominator = max(len(route_ids) - 1, 1)
    mean_distance = total_distances[representative_index] / denominator
    return route_ids[representative_index], mean_distance


def find_similar_cluster_groups(cluster_representatives, route_features):
    """Group Jaccard clusters whose Frechet medoids are within the threshold."""
    cluster_ids = cluster_representatives[CLUSTER_COLUMN].tolist()
    representative_ids = cluster_representatives["representative_route_run_id"].tolist()
    representative_routes = pd.DataFrame({"route_run_id": representative_ids})
    distance_matrix = build_frechet_distance_matrix(representative_routes, route_features)
    remaining = set(range(len(cluster_ids)))
    groups = []

    while remaining:
        group = {remaining.pop()}
        changed = True
        while changed:
            changed = False
            for index in list(remaining):
                if any(
                    distance_matrix[index][other] <= FRECHET_MERGE_THRESHOLD_KM
                    for other in group
                ):
                    group.add(index)
                    remaining.remove(index)
                    changed = True
        groups.append(sorted(group))
    return groups


def select_representatives(route_pair, route_features):
    """Select one representative row for each direct Jaccard cluster."""
    input_path = DIRECT_ROUTES_FOLDER / f"{route_pair}_{METHOD_SUFFIX}.csv"
    routes = pd.read_csv(input_path)
    total_routes = len(routes)
    initial_representatives = []
    grouped_routes = {}
    for cluster_id, cluster_rows in routes.groupby(CLUSTER_COLUMN, sort=True):
        representative_id, mean_distance = select_cluster_medoid(cluster_rows, route_features)
        initial_representatives.append({
            CLUSTER_COLUMN: cluster_id,
            "representative_route_run_id": representative_id,
            "initial_mean_frechet_km": mean_distance,
        })
        grouped_routes[cluster_id] = cluster_rows

    initial_representatives = pd.DataFrame(initial_representatives)
    merged_groups = find_similar_cluster_groups(initial_representatives, route_features)
    representative_rows = []
    representative_features = []

    for merged_index, group in enumerate(merged_groups, start=1):
        source_cluster_ids = initial_representatives.iloc[group][CLUSTER_COLUMN].tolist()
        cluster_rows = pd.concat(
            [grouped_routes[cluster_id] for cluster_id in source_cluster_ids],
            ignore_index=True,
        )
        representative_id, mean_distance = select_cluster_medoid(cluster_rows, route_features)
        representative_row = cluster_rows[
            cluster_rows["route_run_id"] == representative_id
        ].iloc[0].to_dict()
        cluster_size = len(cluster_rows)
        merged_cluster_id = f"frechet_merged_{merged_index:03d}"
        representative_row.update(
            {
                "representative_route_run_id": representative_id,
                "representative_mean_frechet_miles": mean_distance * KM_TO_MILES,
                "cluster_size": cluster_size,
                "cluster_weight_space_percent": 100.0 * cluster_size / total_routes,
                "merged_cluster_id": merged_cluster_id,
                "source_jaccard_clusters": ", ".join(source_cluster_ids),
            }
        )
        representative_rows.append(representative_row)

        feature = json.loads(json.dumps(route_features[representative_id]))
        feature["properties"].update(representative_row)
        representative_features.append(feature)

    representatives = pd.DataFrame(representative_rows).sort_values(
        ["cluster_size", "merged_cluster_id"], ascending=[False, True]
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
