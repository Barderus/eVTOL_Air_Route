"""Create average route cost breakdown pie charts from the analysis CSV."""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import settings


COMPONENTS = [
    ("Distance", "distance_cost", "#3b82f6"),
    ("Population", "population_cost", "#10b981"),
    ("Airspace", "airspace_cost", "#f59e0b"),
    ("Traffic", "traffic_cost", "#8b5cf6"),
]


def percentage_label(percent):
    """Show percentages large enough to read on the pie chart."""
    if percent < 1.0:
        return ""
    return f"{percent:.1f}%"


def average_route_costs(analysis):
    """Average every cost component across the available dates by route."""
    cost_columns = [column for _, column, _ in COMPONENTS]
    average_columns = ["route_distance_km", "total_cost", *cost_columns]
    averages = (
        analysis.groupby("route", sort=False)[average_columns]
        .mean()
        .reset_index()
    )
    averages["calculated_total_cost"] = averages[cost_columns].sum(axis=1)
    return averages


def make_average_chart(route_averages, date_count):
    """Create one figure containing the average breakdown for every route."""
    figure, axes = plt.subplots(
        1,
        len(route_averages),
        figsize=(17, 6),
    )

    if len(route_averages) == 1:
        axes = [axes]

    for axis, row in zip(axes, route_averages.itertuples(index=False)):
        values = [getattr(row, column) for _, column, _ in COMPONENTS]
        labels = [name for name, _, _ in COMPONENTS]
        colors = [color for _, _, color in COMPONENTS]

        axis.pie(
            values,
            labels=labels,
            colors=colors,
            startangle=90,
            counterclock=False,
            autopct=percentage_label,
            textprops={"fontsize": 9},
            wedgeprops={"linewidth": 1, "edgecolor": "white"},
        )
        axis.set_title(
            f"{row.route}\n"
            f"Average total cost: {row.calculated_total_cost:.3f} | "
            f"Average distance: {row.route_distance_km:.2f} km",
            fontsize=11,
            pad=14,
        )
        axis.set_aspect("equal")

    figure.suptitle(
        f"St. Louis Average Route Cost Breakdown\n"
        f"Average of {date_count} traffic dates",
        fontsize=16,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.9))

    output_name = "route_cost_breakdown_average.png"
    output_path = os.path.join(settings.ROUTE_ANALYSIS_FOLDER, output_name)
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved average route cost chart: {output_path}")


def main():
    if not os.path.isfile(settings.ROUTE_ANALYSIS_CSV):
        raise FileNotFoundError(
            f"Route analysis CSV not found: {settings.ROUTE_ANALYSIS_CSV}. "
            "Run astar_routes.py first."
        )

    analysis = pd.read_csv(settings.ROUTE_ANALYSIS_CSV)
    os.makedirs(settings.ROUTE_ANALYSIS_FOLDER, exist_ok=True)
    date_count = analysis["date"].nunique()
    route_averages = average_route_costs(analysis)
    make_average_chart(route_averages, date_count)


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as error:
        raise SystemExit(str(error)) from error
