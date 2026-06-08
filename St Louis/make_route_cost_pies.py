"""Create route cost breakdown pie charts from the analysis CSV."""

import os
import re

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


def clean_filename(value):
    """Convert a label into a simple filename."""
    filename = value.lower()
    filename = re.sub(r"[^a-z0-9]+", "_", filename)
    return filename.strip("_")


def percentage_label(percent):
    """Show percentages large enough to read on the pie chart."""
    if percent < 1.0:
        return ""
    return f"{percent:.1f}%"


def make_date_chart(date_rows):
    """Create one figure containing the three route breakdowns for a date."""
    dataset_label = date_rows.iloc[0]["dataset"]
    figure, axes = plt.subplots(1, len(date_rows), figsize=(17, 6))

    if len(date_rows) == 1:
        axes = [axes]

    for axis, row in zip(axes, date_rows.itertuples(index=False)):
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
            f"Total cost: {row.total_cost:.3f} | "
            f"Distance: {row.route_distance_km:.2f} km",
            fontsize=11,
            pad=14,
        )
        axis.set_aspect("equal")

    figure.suptitle(
        f"St. Louis Route Cost Breakdown\n{dataset_label}",
        fontsize=16,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.9))

    output_name = f"route_cost_breakdown_{clean_filename(dataset_label)}.png"
    output_path = os.path.join(settings.ROUTE_ANALYSIS_FOLDER, output_name)
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved route cost chart: {output_path}")


def main():
    if not os.path.isfile(settings.ROUTE_ANALYSIS_CSV):
        raise FileNotFoundError(
            f"Route analysis CSV not found: {settings.ROUTE_ANALYSIS_CSV}. "
            "Run astar_routes.py first."
        )

    analysis = pd.read_csv(settings.ROUTE_ANALYSIS_CSV)
    os.makedirs(settings.ROUTE_ANALYSIS_FOLDER, exist_ok=True)

    for dataset in settings.TRAFFIC_DATASETS:
        date_rows = analysis[analysis["date"] == dataset["date"]].copy()
        date_rows = date_rows.reset_index(drop=True)
        make_date_chart(date_rows)


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as error:
        raise SystemExit(str(error)) from error
