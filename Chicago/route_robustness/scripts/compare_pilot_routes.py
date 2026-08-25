"""Compare pilot routes using shared path-node similarity."""

import json
import os

import numpy as np
import pandas as pd


PILOT_ROUTE_RUNS_CSV = "Chicago/route_robustness/output/pilot_route_runs.csv"
OUTPUT_FOLDER = "Chicago/route_robustness/output"
SIMILARITY_MATRIX_CSV = os.path.join(
    OUTPUT_FOLDER,
    "pilot_route_similarity_matrix.csv",
)
SIMILARITY_PAIRS_CSV = os.path.join(
    OUTPUT_FOLDER,
    "pilot_route_similarity_pairs.csv",
)


def load_route_runs():
    """Load successful pilot route runs with saved path node IDs."""
    if not os.path.exists(PILOT_ROUTE_RUNS_CSV):
        raise FileNotFoundError(f"Pilot route runs not found: {PILOT_ROUTE_RUNS_CSV}")

    route_runs = pd.read_csv(PILOT_ROUTE_RUNS_CSV)
    required_columns = {"pilot_name", "route_run_id", "path_node_ids", "status"}
    missing_columns = required_columns.difference(route_runs.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Pilot route table is missing columns: {missing_text}")

    route_runs = route_runs[route_runs["status"] == "success"].copy()
    if route_runs.empty:
        raise ValueError("No successful pilot routes found.")

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
    pilot_names = route_runs["pilot_name"].tolist()
    route_lookup = {
        row.pilot_name: row.path_node_set
        for row in route_runs.itertuples(index=False)
    }

    matrix = pd.DataFrame(
        np.eye(len(pilot_names)),
        index=pilot_names,
        columns=pilot_names,
    )
    pair_rows = []

    for first_index, first_name in enumerate(pilot_names):
        for second_index, second_name in enumerate(pilot_names):
            if second_index <= first_index:
                continue

            similarity, shared_node_count, union_node_count = calculate_node_similarity(
                route_lookup[first_name],
                route_lookup[second_name],
            )
            matrix.loc[first_name, second_name] = similarity
            matrix.loc[second_name, first_name] = similarity
            pair_rows.append(
                {
                    "route_a": first_name,
                    "route_b": second_name,
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
    """Compare every successful pilot route pair."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    route_runs = load_route_runs()
    matrix, pairs = build_similarity_outputs(route_runs)

    matrix.to_csv(SIMILARITY_MATRIX_CSV)
    pairs.to_csv(SIMILARITY_PAIRS_CSV, index=False)

    print(f"Saved similarity matrix: {SIMILARITY_MATRIX_CSV}")
    print(f"Saved similarity pairs: {SIMILARITY_PAIRS_CSV}")
    print("Most similar route pairs:")
    print(pairs.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
