from __future__ import annotations

"""Approach 2: smooth the displayed route geometry with a cubic spline.

This approach does not change A* search. It starts from the already-generated
combined route, smooths the rendered line, and then re-scores that displayed
line against the base combined route weights.
"""

from pathlib import Path

from route_smoothing_common import (
    RouteResult,
    base_combined_routes_by_destination,
    base_routes,
    build_dataset_contexts,
    catmull_rom_spline,
    evaluate_combined_cost_for_line,
    route_result_to_row,
    sort_rows,
    write_cost_csv,
    write_html,
)


OUTPUT_CSV_PATH = Path(__file__).resolve().parent / "approach2_cubic_spline_costs.csv"


def build_results() -> list[RouteResult]:
    """Build one spline-smoothed result for every baseline combined route."""
    contexts = build_dataset_contexts()
    base_routes_by_destination = base_combined_routes_by_destination()
    results: list[RouteResult] = []

    for destination in base_routes.DESTINATIONS:
        destination_gdf = base_routes_by_destination[destination["slug"]]
        context_lookup = {dataset["slug"]: contexts[dataset["slug"]] for dataset in base_routes.DATASETS}

        for row in destination_gdf.itertuples(index=False):
            context = context_lookup[row.dataset_slug]
            # `row.total_cost` is the original combined A* result from the base
            # route generator. After smoothing, we re-sample the displayed line
            # back onto the grid to estimate its comparable final cost.
            display_line = catmull_rom_spline(row.geometry)
            final_cost, final_path_nodes, _ = evaluate_combined_cost_for_line(display_line, context)
            results.append(
                RouteResult(
                    destination_slug=row.destination_slug,
                    destination_label=row.destination_label,
                    dataset_slug=row.dataset_slug,
                    dataset_label=row.dataset_label,
                    total_cost=float(row.total_cost),
                    path_nodes=int(row.path_nodes),
                    final_combined_cost=final_cost,
                    final_path_nodes=final_path_nodes,
                    display_line=display_line,
                    raw_path=[],
                    metadata={
                        "technique": "approach2_cubic_spline",
                        "spline": "catmull_rom",
                        "destination_lat": destination["lat"],
                        "destination_lon": destination["lon"],
                    },
                )
            )
    return results


def main() -> None:
    """Run the cubic-spline experiment and write its CSV and HTML outputs."""
    results = build_results()
    rows = [route_result_to_row(result) for result in results]
    sort_rows(rows)
    write_cost_csv(OUTPUT_CSV_PATH, rows)
    html_path = write_html(
        output_filename="combined_route_cubic_spline.html",
        page_title="Combined Route Smoothing: Cubic Spline",
        technique_name="Approach 2: Cubic spline smoothing",
        results=results,
    )

    print("Approach 2 route costs")
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
