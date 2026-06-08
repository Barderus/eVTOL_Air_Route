from __future__ import annotations

"""Approach 3: combine turn-penalty A* with cubic-spline display smoothing.

This approach first generates the route with the turn-penalty search from
Approach 1, then applies the same spline display smoothing used in Approach 2.
It is the most aggressive smoothing variant in the repository.
"""

from pathlib import Path

from route_smoothing_common import (
    RouteResult,
    astar_with_turn_penalty,
    base_routes,
    build_dataset_contexts,
    catmull_rom_spline,
    evaluate_combined_cost_for_line,
    route_result_to_row,
    sort_rows,
    write_cost_csv,
    write_html,
)


OUTPUT_CSV_PATH = Path(__file__).resolve().parent / "approach3_turn_penalty_plus_spline_costs.csv"


def build_results() -> list[RouteResult]:
    """Build one turn-penalty-plus-spline result per dataset and destination."""
    contexts = build_dataset_contexts()
    results: list[RouteResult] = []

    for dataset in base_routes.DATASETS:
        context = contexts[dataset["slug"]]
        for destination in base_routes.DESTINATIONS:
            path, search_cost = astar_with_turn_penalty(context, destination["slug"])
            base_line = base_routes.route_line(context.grid, path)
            # Search is smoothed once by the turn penalty and then a second time
            # visually by the spline. Re-score the rendered line so comparisons
            # stay on the same combined-weight scale as the other approaches.
            display_line = catmull_rom_spline(base_line)
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
                        "technique": "approach3_turn_penalty_plus_spline",
                        "turn_penalty_weight": 1.5,
                        "spline": "catmull_rom",
                        "destination_lat": destination["lat"],
                        "destination_lon": destination["lon"],
                    },
                )
            )
    return results


def main() -> None:
    """Run the combined smoothing experiment and write its outputs."""
    results = build_results()
    rows = [route_result_to_row(result) for result in results]
    sort_rows(rows)
    write_cost_csv(OUTPUT_CSV_PATH, rows)
    html_path = write_html(
        output_filename="combined_route_turn_penalty_plus_spline.html",
        page_title="Combined Route Smoothing: Turn-Penalty + Cubic Spline",
        technique_name="Approach 3: Turn-penalty A* plus cubic spline",
        results=results,
    )

    print("Approach 3 route costs")
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
