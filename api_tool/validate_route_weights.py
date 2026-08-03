"""Validate route weights for an API-style workflow."""


WEIGHT_FIELDS = [
    "distance_weight",
    "population_weight",
    "traffic_weight",
    "airspace_weight",
]


def prompt_float(prompt_text):
    """Ask for a required floating-point value."""
    while True:
        value = input(prompt_text).strip()
        try:
            return float(value)
        except ValueError:
            print("Enter a valid number.")


def validate_route_weights(
    distance_weight=None,
    population_weight=None,
    traffic_weight=None,
    airspace_weight=None,
    require_sum_to_one=True,
    tolerance=0.000001,
):
    """Validate route weights and return a structured result."""
    if distance_weight is None:
        distance_weight = prompt_float("Distance weight: ")
        population_weight = prompt_float("Population weight: ")
        traffic_weight = prompt_float("Traffic weight: ")
        airspace_weight = prompt_float("Airspace weight: ")

    weights = {
        "distance_weight": float(distance_weight),
        "population_weight": float(population_weight),
        "traffic_weight": float(traffic_weight),
        "airspace_weight": float(airspace_weight),
    }

    errors = []
    for field_name, value in weights.items():
        if value < 0:
            errors.append(f"{field_name} must be greater than or equal to 0.")

    weight_sum = sum(weights.values())
    if require_sum_to_one and abs(weight_sum - 1.0) > tolerance:
        errors.append(f"Weights must sum to 1.0. Current sum is {weight_sum}.")

    return {
        "weights": weights,
        "weight_sum": weight_sum,
        "require_sum_to_one": require_sum_to_one,
        "tolerance": tolerance,
        "is_valid": len(errors) == 0,
        "errors": errors,
    }


def request_route_weights_from_user(require_sum_to_one=True):
    """Ask for route weights until the configuration is valid."""
    while True:
        result = validate_route_weights(require_sum_to_one=require_sum_to_one)
        if result["is_valid"]:
            return result["weights"]

        print("Invalid weight configuration:")
        for error in result["errors"]:
            print("-", error)
        print("Enter the weights again.")


if __name__ == "__main__":
    result = validate_route_weights()
    print(result)
