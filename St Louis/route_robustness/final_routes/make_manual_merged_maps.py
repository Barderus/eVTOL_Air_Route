"""Build manually consolidated maps from the 70% buffered-overlap representatives."""

import importlib.util
import json
from pathlib import Path

import geopandas as gpd
import pandas as pd


HERE = Path(__file__).resolve().parent
BASE_GENERATOR_PATH = HERE / "make_final_representative_maps.py"
BASE_OVERLAP_THRESHOLD = 0.70

MANUAL_MERGES = {
    "midamerica_to_st_louis_lambert": {
        "July_12": [
            ["buffer_merged_003", "buffer_merged_004"],
            ["buffer_merged_002", "buffer_merged_005"],
        ],
        "January_10": [["buffer_merged_001", "buffer_merged_004"]],
        "January_12": [
            ["buffer_merged_001", "buffer_merged_002", "buffer_merged_003"],
            ["buffer_merged_004", "buffer_merged_005"],
        ],
        "March_9": [
            ["buffer_merged_003", "buffer_merged_004"],
            ["buffer_merged_001", "buffer_merged_002"],
        ],
    },
    "midamerica_to_st_louis_union_station": {
        "January_10": [["buffer_merged_002", "buffer_merged_004"]],
    },
    "st_louis_downtown_airport_to_st_louis_lambert": {
        "July_14": [
            ["frechet_merged_002", "frechet_merged_004", "frechet_merged_003"]
        ],
        "March_9": [["frechet_merged_002", "frechet_merged_003"]],
    },
}


def load_base_generator():
    """Load the original map module so both map sets share one HTML structure."""
    spec = importlib.util.spec_from_file_location("base_final_maps", BASE_GENERATOR_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BASE = load_base_generator()


def mutual_overlap(geometry_a, geometry_b):
    """Measure shared area as a fraction of the larger corridor."""
    intersection_area = geometry_a.intersection(geometry_b).area
    larger_area = max(geometry_a.area, geometry_b.area)
    return intersection_area / larger_area if larger_area else 0.0


def overlap_matrix(corridors):
    """Build a symmetric mutual-overlap matrix."""
    count = len(corridors)
    matrix = [[1.0 for _ in range(count)] for _ in range(count)]
    for index_a in range(count):
        for index_b in range(index_a + 1, count):
            overlap = mutual_overlap(
                corridors.geometry.iloc[index_a], corridors.geometry.iloc[index_b]
            )
            matrix[index_a][index_b] = overlap
            matrix[index_b][index_a] = overlap
    return matrix


def complete_link_groups(matrix, threshold):
    """Reconstruct the prior all-pairs buffered-overlap groups."""
    groups = [{index} for index in range(len(matrix))]
    while True:
        best_pair = None
        best_mean = -1.0
        for left in range(len(groups)):
            for right in range(left + 1, len(groups)):
                values = [
                    matrix[index_a][index_b]
                    for index_a in groups[left]
                    for index_b in groups[right]
                ]
                mean = sum(values) / len(values)
                if min(values) >= threshold and mean > best_mean:
                    best_pair = (left, right)
                    best_mean = mean
        if best_pair is None:
            break
        left, right = best_pair
        groups[left] = groups[left].union(groups[right])
        del groups[right]
    return sorted((sorted(group) for group in groups), key=lambda group: group[0])


def medoid_index(group, matrix):
    """Select an existing center route with the greatest mean corridor overlap."""
    if len(group) == 1:
        return group[0]
    return max(
        group,
        key=lambda index: (
            sum(matrix[index][other] for other in group if other != index)
            / (len(group) - 1),
            -index,
        ),
    )


def build_70_percent_baseline(dated_routes):
    """Rebuild the named buffer_merged groups used in the user's merge list."""
    dated_routes = dated_routes.reset_index(drop=True)
    corridors = BASE.build_corridors(dated_routes).to_crs(BASE.PROJECTED_CRS)
    matrix = overlap_matrix(corridors)
    baseline = []
    for group_number, group in enumerate(
        complete_link_groups(matrix, BASE_OVERLAP_THRESHOLD), start=1
    ):
        selected = medoid_index(group, matrix)
        row = dated_routes.iloc[selected].copy()
        sources = dated_routes.iloc[group]
        row["merged_cluster_id"] = f"buffer_merged_{group_number:03d}"
        row["source_cluster_ids"] = list(sources["merged_cluster_id"].astype(str))
        row["source_weight_ids"] = list(sources["weight_id"].astype(str))
        row["cluster_size"] = int(sources["cluster_size"].sum())
        row["cluster_weight_space_percent"] = float(
            sources["cluster_weight_space_percent"].sum()
        )
        baseline.append(row)
    return gpd.GeoDataFrame(baseline, crs=dated_routes.crs).reset_index(drop=True)


def resolve_manual_group(baseline, identifiers):
    """Resolve buffer IDs or original Frechet IDs to baseline row indexes."""
    resolved = set()
    for identifier in identifiers:
        matches = [
            index
            for index, row in baseline.iterrows()
            if row["merged_cluster_id"] == identifier
            or identifier in row["source_cluster_ids"]
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one match for {identifier}, found {len(matches)}"
            )
        resolved.add(matches[0])
    return resolved


def apply_manual_merges(baseline, requested_merges):
    """Apply the explicitly requested unions and keep all other rows separate."""
    groups = [{index} for index in range(len(baseline))]
    claimed = set()
    for identifiers in requested_merges:
        indexes = resolve_manual_group(baseline, identifiers)
        if claimed.intersection(indexes):
            raise ValueError(f"Overlapping manual merge request: {identifiers}")
        claimed.update(indexes)
        groups = [group for group in groups if group.isdisjoint(indexes)]
        groups.append(indexes)
    groups = sorted(groups, key=lambda group: min(group))

    corridors = BASE.build_corridors(baseline).to_crs(BASE.PROJECTED_CRS)
    matrix = overlap_matrix(corridors)
    merged = []
    for group_number, group_set in enumerate(groups, start=1):
        group = sorted(group_set)
        selected = medoid_index(group, matrix)
        row = baseline.iloc[selected].copy()
        sources = baseline.iloc[group]
        row["merged_cluster_id"] = f"manual_merged_{group_number:03d}"
        row["manual_source_buffer_ids"] = ", ".join(
            sources["merged_cluster_id"].astype(str)
        )
        row["manual_source_weight_ids"] = ", ".join(
            sources["weight_id"].astype(str)
        )
        row["source_cluster_ids"] = [
            cluster_id
            for values in sources["source_cluster_ids"]
            for cluster_id in values
        ]
        row["source_weight_ids"] = [
            weight_id
            for values in sources["source_weight_ids"]
            for weight_id in values
        ]
        row["cluster_size"] = int(sources["cluster_size"].sum())
        row["cluster_weight_space_percent"] = float(
            sources["cluster_weight_space_percent"].sum()
        )
        row["representative_color"] = BASE.REPRESENTATIVE_COLORS[
            (group_number - 1) % len(BASE.REPRESENTATIVE_COLORS)
        ]
        row["display_cluster_id"] = (
            f"{row['date_label']} - {row['merged_cluster_id']}"
        )
        merged.append(row)
    return gpd.GeoDataFrame(merged, crs=baseline.crs).reset_index(drop=True)


def build_manual_routes(route_pair):
    """Create the manually merged representatives for every displayed date."""
    route_data = BASE.load_route_pair_features(route_pair)
    routes = gpd.GeoDataFrame.from_features(route_data["features"], crs="EPSG:4326")
    outputs = []
    pair_merges = MANUAL_MERGES.get(route_pair, {})
    for date_key, dated_routes in routes.groupby("date_key", sort=False):
        baseline = build_70_percent_baseline(dated_routes)
        requested = pair_merges.get(date_key, [])
        merged = apply_manual_merges(baseline, requested)
        outputs.append(merged)
        print(
            f"{route_pair} {date_key}: {len(dated_routes)} original -> "
            f"{len(baseline)} buffer groups -> {len(merged)} manual groups"
        )
    return gpd.GeoDataFrame(
        pd.concat(outputs, ignore_index=True), crs="EPSG:4326"
    )


def write_manual_map(route_pair):
    """Write one manual-merge HTML map without changing the reference map."""
    routes = build_manual_routes(route_pair)
    corridors = BASE.build_corridors(routes)
    route_data = json.loads(routes.to_json(drop_id=True))
    corridor_data = json.loads(corridors.to_json(drop_id=True))
    route_label = route_data["features"][0]["properties"]["route_pair_label"]
    title = f"St. Louis Manually Merged Route Representatives - {route_label}"
    date_order = [
        {
            "key": folder,
            "label": label,
            "shortLabel": folder.replace("_", " "),
            "date": traffic_date,
        }
        for folder, label, traffic_date in BASE.DATE_FOLDERS
    ]
    template = BASE.HTML_TEMPLATE.replace(
        "<b>Representative:</b> ${{properties.weight_id}}<br>",
        "<b>Representative:</b> ${{properties.weight_id}}<br>\n"
        "        <b>Merged buffer groups:</b> "
        "${{properties.manual_source_buffer_ids}}<br>\n"
        "        <b>Source weights:</b> ${{properties.manual_source_weight_ids}}<br>",
    )
    html = template.format(
        title=title,
        representative_data=json.dumps(route_data),
        corridor_data=json.dumps(corridor_data),
        date_order=json.dumps(date_order),
    )
    output_path = HERE / f"{route_pair}_all_days_manual_merged_representatives.html"
    output_path.write_text(html, encoding="utf-8")
    print(f"Saved manual representative map: {output_path}")


def main():
    """Build all three manually merged maps."""
    for route_pair in BASE.ROUTE_PAIRS:
        write_manual_map(route_pair)


if __name__ == "__main__":
    main()
