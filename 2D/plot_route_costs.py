from pathlib import Path
import json

import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parent.parent
GEOJSON_FOLDER = PROJECT_ROOT / "geojson"
OUTPUT_FOLDER = PROJECT_ROOT / "images"

DESTINATIONS = [
    {"slug": "midway", "label": "Midway International Airport"},
    {"slug": "ohare", "label": "O'Hare International Airport"},
    {"slug": "union_station", "label": "Chicago"},
]

DATASETS = [
    {"slug": "2026-03-07", "label": "Sat 2026-03-07"},
    {"slug": "2026-03-09", "label": "Mon 2026-03-09"},
    {"slug": "2026-01-10", "label": "Sat 2026-01-10"},
    {"slug": "2026-01-12", "label": "Mon 2026-01-12"},
    {"slug": "2025-07-14", "label": "Mon 2025-07-14"},
    {"slug": "2025-07-12", "label": "Sat 2025-07-12"},
]

COMBINED_ROUTE_NAME = "combined"
COMBINED_ROUTE_COLOR = "#dc2626"


def load_geojson(destination_slug):
    geojson_path = GEOJSON_FOLDER / f"clow_to_{destination_slug}_astar_routes.geojson"

    with geojson_path.open("r", encoding="utf-8") as file_handle:
        data = json.load(file_handle)

    return data


def get_cost_map(destination_slug):
    geojson_data = load_geojson(destination_slug)
    cost_map = {}

    for feature in geojson_data["features"]:
        properties = feature["properties"]
        dataset_slug = properties["dataset_slug"]
        route_name = properties["route_name"]
        total_cost = float(properties["total_cost"])

        if dataset_slug not in cost_map:
            cost_map[dataset_slug] = {}

        cost_map[dataset_slug][route_name] = total_cost

    return cost_map


def make_bar_chart(destination):
    destination_slug = destination["slug"]
    destination_label = destination["label"]
    cost_map = get_cost_map(destination_slug)

    figure, axis = plt.subplots(figsize=(11, 6))

    x_positions = []
    x_labels = []
    bar_heights = []

    for i, dataset in enumerate(DATASETS):
        x_positions.append(i)
        x_labels.append(dataset["label"])
        dataset_slug = dataset["slug"]
        bar_heights.append(cost_map[dataset_slug][COMBINED_ROUTE_NAME])

    axis.bar(
        x_positions,
        bar_heights,
        width=0.6,
        color=COMBINED_ROUTE_COLOR,
    )

    axis.set_xticks(x_positions)
    axis.set_xticklabels(x_labels, rotation=25, ha="right")
    axis.set_xlabel("Traffic Dataset")
    axis.set_ylabel("Combined Route Cost")
    axis.set_title(f"Clow to {destination_label} Combined Route Cost")
    axis.grid(axis="y", alpha=0.25)

    figure.tight_layout()

    output_path = OUTPUT_FOLDER / f"clow_to_{destination_slug}_route_costs.png"
    figure.savefig(output_path, dpi=200)
    plt.close(figure)

    return output_path


def main():
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    for destination in DESTINATIONS:
        output_path = make_bar_chart(destination)
        print(f"Saved route cost plot: {output_path}")


if __name__ == "__main__":
    main()
