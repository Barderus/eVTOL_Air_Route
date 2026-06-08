from pathlib import Path
import csv
import io
from contextlib import redirect_stdout

from generate_astar_toggle_pages import DATASETS
from generate_astar_toggle_pages import DESTINATIONS
from generate_astar_toggle_pages import build_route_features


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_CSV_PATH = PROJECT_ROOT / "2D" / "combined_route_costs.csv"


def get_destination_name(destination_slug, original_label):
    if destination_slug == "union_station":
        return "Chicago"
    return original_label


def build_rows():
    rows = []

    with redirect_stdout(io.StringIO()):
        all_routes = build_route_features()

    for destination in DESTINATIONS:
        destination_slug = destination["slug"]
        destination_routes = all_routes[destination_slug]

        combined_only = destination_routes[destination_routes["route_name"] == "combined"]

        for row in combined_only.itertuples(index=False):
            row_data = {
                "destination_slug": row.destination_slug,
                "destination": get_destination_name(row.destination_slug, row.destination_label),
                "dataset_slug": row.dataset_slug,
                "dataset_label": row.dataset_label,
                "total_cost": float(row.total_cost),
                "path_nodes": int(row.path_nodes),
            }
            rows.append(row_data)

    return rows


def sort_rows(rows):
    destination_positions = {}
    dataset_positions = {}

    for i, destination in enumerate(DESTINATIONS):
        destination_positions[destination["slug"]] = i

    for i, dataset in enumerate(DATASETS):
        dataset_positions[dataset["slug"]] = i

    rows.sort(
        key=lambda row: (
            destination_positions[row["destination_slug"]],
            dataset_positions[row["dataset_slug"]],
        )
    )


def save_csv(rows):
    field_names = [
        "destination_slug",
        "destination",
        "dataset_slug",
        "dataset_label",
        "total_cost",
        "path_nodes",
    ]

    with OUTPUT_CSV_PATH.open("w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=field_names)
        writer.writeheader()
        writer.writerows(rows)


def print_rows(rows):
    print("Combined route costs")
    print("====================")

    for row in rows:
        destination_name = row["destination"]
        dataset_label = row["dataset_label"]
        total_cost = row["total_cost"]
        path_nodes = row["path_nodes"]

        print(f"{destination_name} | {dataset_label} | cost={total_cost:.4f} | nodes={path_nodes}")

    print()
    print(f"Saved CSV: {OUTPUT_CSV_PATH}")


def main():
    rows = build_rows()
    sort_rows(rows)
    save_csv(rows)
    print_rows(rows)


if __name__ == "__main__":
    main()
