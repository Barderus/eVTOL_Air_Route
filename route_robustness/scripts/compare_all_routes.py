"""Compare all weighted routes using shared path-node similarity."""

import json
import os

import numpy as np
import pandas as pd


ALL_ROUTE_RUNS_CSV = "route_robustness/output/all_route_runs.csv"
OUTPUT_FOLDER = "route_robustness/output"
SIMILARITY_MATRIX_CSV = os.path.join(
    OUTPUT_FOLDER,
    "all_route_similarity_matrix.csv",
)
SIMILARITY_PAIRS_CSV = os.path.join(
    OUTPUT_FOLDER,
    "all_route_similarity_pairs.csv",
)


def load_route_runs():
    """Load successful full route runs with saved path node IDs."""
    if not os.path.exists(ALL_ROUTE_RUNS_CSV):
        raise FileNotFoundError(f"All route runs not found: {ALL_ROUTE_RUNS_CSV}")

    route_runs = pd.read_csv(ALL_ROUTE_RUNS_CSV)
    required_columns = {"weight_id", "route_run_id", "path_node_ids", "status"}
    missing_columns = required_columns.difference(route_runs.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"All-route table is missing columns: {missing_text}")

    route_runs = route_runs[route_runs["status"] == "success"].copy()
    if route_runs.empty:
        raise ValueError("No successful routes found.")

    route_runs["path_node_set"] = route_runs["path_node_ids"].apply(
        lambda value: set(json.loads(value))
    )
    return route_runs


def calculate_node_similarity(first_nodes, second_nodes):
    """Calculate shared-node Jaccard similarity for two routes."""
    union_nodes = first_nodes.union(second_nodes)
    if not union_nodes:
        return 0.0, 0, 0

    shared_nodes = first_nodes.intersection(second_nodes)
    similarity = len(shared_nodes) / len(union_nodes)
    return similarity, len(shared_nodes), len(union_nodes)


def build_similarity_outputs(route_runs):
    """Create matrix and long-form pairwise similarity tables."""
    route_ids = route_runs["route_run_id"].tolist()
    weight_ids = dict(zip(route_runs["route_run_id"], route_runs["weight_id"]))
    route_lookup = {
        row.route_run_id: row.path_node_set
        for row in route_runs.itertuples(index=False)
    }

    matrix = pd.DataFrame(
        np.eye(len(route_ids)),
        index=route_ids,
        columns=route_ids,
    )
    pair_rows = []

    for first_index, first_route_id in enumerate(route_ids):
        for second_index in range(first_index + 1, len(route_ids)):
            second_route_id = route_ids[second_index]
            similarity, shared_node_count, union_node_count = calculate_node_similarity(
                route_lookup[first_route_id],
                route_lookup[second_route_id],
            )
            matrix.loc[first_route_id, second_route_id] = similarity
            matrix.loc[second_route_id, first_route_id] = similarity
            pair_rows.append(
                {
                    "route_run_a": first_route_id,
                    "route_run_b": second_route_id,
                    "weight_id_a": weight_ids[first_route_id],
                    "weight_id_b": weight_ids[second_route_id],
                    "similarity": similarity,
                    "shared_node_count": shared_node_count,
                    "union_node_count": union_node_count,
                }
            )

    pairs = pd.DataFrame(pair_rows)
    pairs = pairs.sort_values(
        by="similarity",
        ascending=False,
    ).reset_index(drop=True)
    return matrix, pairs


def main():
    """Compare every successful full-route pair."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    route_runs = load_route_runs()
    matrix, pairs = build_similarity_outputs(route_runs)

    matrix.to_csv(SIMILARITY_MATRIX_CSV)
    pairs.to_csv(SIMILARITY_PAIRS_CSV, index=False)

    print(f"Saved similarity matrix: {SIMILARITY_MATRIX_CSV}")
    print(f"Saved similarity pairs: {SIMILARITY_PAIRS_CSV}")
    print("Route count:", len(route_runs))
    print("Pair count:", len(pairs))
    print("Most similar route pairs:")
    print(pairs.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
