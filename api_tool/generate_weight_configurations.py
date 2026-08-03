"""Generate route-weight configurations for robustness workflows."""

import csv
import os


def save_weight_configurations(rows, output_path):
    """Save weight configurations as CSV."""
    output_folder = os.path.dirname(output_path)
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def generate_weight_configurations(step=0.10, output_path=None, include_equal_weight=True):
    """
    Generate valid four-factor weights that sum to 1.0.

    The step is converted to integer units to avoid floating-point duplicates.
    """
    if step <= 0 or step > 1:
        raise ValueError("step must be greater than 0 and less than or equal to 1.")

    unit_count = round(1 / step)
    resolved_step = 1 / unit_count
    rows = []
    seen_weights = set()

    for distance_units in range(unit_count + 1):
        for population_units in range(unit_count + 1):
            for traffic_units in range(unit_count + 1):
                airspace_units = unit_count - distance_units - population_units - traffic_units
                if airspace_units < 0:
                    continue

                weights = (
                    round(distance_units * resolved_step, 6),
                    round(population_units * resolved_step, 6),
                    round(traffic_units * resolved_step, 6),
                    round(airspace_units * resolved_step, 6),
                )
                seen_weights.add(weights)

    if include_equal_weight:
        seen_weights.add((0.25, 0.25, 0.25, 0.25))

    for index, weights in enumerate(sorted(seen_weights), start=1):
        distance_weight, population_weight, traffic_weight, airspace_weight = weights
        rows.append(
            {
                "weight_id": f"w_{index:04d}",
                "distance_weight": distance_weight,
                "population_weight": population_weight,
                "traffic_weight": traffic_weight,
                "airspace_weight": airspace_weight,
                "weight_sum": round(sum(weights), 6),
                "weight_step": round(resolved_step, 6),
            }
        )

    if output_path:
        save_weight_configurations(rows, output_path)

    return {
        "weight_configurations": rows,
        "count": len(rows),
        "step": round(resolved_step, 6),
        "output_path": output_path,
    }


if __name__ == "__main__":
    result = generate_weight_configurations()
    print(result)

