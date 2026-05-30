from __future__ import annotations

"""Approach 1: smooth routes by penalizing turns during A* search.

This approach changes the search itself. The returned path still follows the
underlying grid exactly, but sharp heading changes become more expensive, which
nudges the solver toward smoother-looking routes.
"""

from pathlib import Path

import json

from route_smoothing_common import (
    HTML_OUTPUT,
    RouteResult,
    astar_with_turn_penalty,
    base_routes,
    build_dataset_contexts,
    destination_display_name,
    evaluate_combined_cost_for_line,
    route_result_to_row,
    sort_rows,
    write_cost_csv,
    write_html,
)


OUTPUT_CSV_PATH = Path(__file__).resolve().parent / "approach1_turn_penalty_costs.csv"


def build_results() -> list[RouteResult]:
    """Build one turn-penalty route result for every dataset and destination."""
    contexts = build_dataset_contexts()
    results: list[RouteResult] = []

    for dataset in base_routes.DATASETS:
        context = contexts[dataset["slug"]]
        for destination in base_routes.DESTINATIONS:
            # The search cost includes the turn penalty. The recomputed final
            # cost uses the baseline combined weights so the three smoothing
            # approaches can still be compared on the same scale.
            path, search_cost = astar_with_turn_penalty(context, destination["slug"])
            display_line = base_routes.route_line(context.grid, path)
            final_cost, final_path_nodes, _ = evaluate_combined_cost_for_line(display_line, context)
            results.append(
                RouteResult(
                    destination_slug=destination["slug"],
                    destination_label=destination["label"],
                    dataset_slug=dataset["slug"],
                    dataset_label=dataset["label"],
                    total_cost=float(search_cost),
                    path_nodes=len(path),
                    final_combined_cost=final_cost,
                    final_path_nodes=final_path_nodes,
                    display_line=display_line,
                    raw_path=path,
                    metadata={
                        "technique": "approach1_turn_penalty",
                        "turn_penalty_weight": 1.5,
                        "destination_lat": destination["lat"],
                        "destination_lon": destination["lon"],
                    },
                )
            )
    return results


def main() -> None:
    """Run the turn-penalty experiment and write its CSV and HTML outputs."""
    results = build_results()
    rows = [route_result_to_row(result) for result in results]
    sort_rows(rows)
    write_cost_csv(OUTPUT_CSV_PATH, rows)
    html_path = write_html(
        output_filename="combined_route_turn_penalty.html",
        page_title="Combined Route Smoothing: Turn-Penalty A*",
        technique_name="Approach 1: Turn-penalty A*",
        results=results,
    )

    print("Approach 1 route costs")
    print("======================")
    for row in rows:
        print(
            f"{row['destination']} | {row['dataset_label']} | "
            f"search={row['total_cost']:.4f} | final={row['final_combined_cost']:.4f}"
        )
    print()
    print(f"Saved CSV: {OUTPUT_CSV_PATH}")
    print(f"Saved HTML: {html_path}")


if __name__ == "__main__":
    main()
