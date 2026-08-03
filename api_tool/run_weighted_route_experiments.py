"""Run weighted route experiments for robustness workflows."""

from api_tool.generate_weight_configurations import generate_weight_configurations
from api_tool.run_astar_route import run_astar_route


def run_weighted_route_experiments(
    city_name=None,
    route_pair=None,
    weight_configurations=None,
    graph_package=None,
    graph_path=None,
    grid_path=None,
    traffic_dataset=None,
    route_runner=None,
    run_routes=False,
):
    """
    Prepare or run one route experiment per weight configuration.

    If run_routes is True, the API A* runner is called once per weight row.
    A custom route_runner can also be supplied.
    """
    if city_name is None:
        city_name = input("City or study area name: ").strip()
        grid_path = input("Scored grid path (optional): ").strip() or None

    if weight_configurations is None:
        weight_configurations = generate_weight_configurations()["weight_configurations"]

    experiment_rows = []
    for weight_row in weight_configurations:
        experiment = {
            "city": city_name,
            "route_pair": route_pair,
            "traffic_dataset": traffic_dataset,
            "graph_path": graph_path,
            "grid_path": grid_path,
            "weight_id": weight_row["weight_id"],
            "distance_weight": weight_row["distance_weight"],
            "population_weight": weight_row["population_weight"],
            "traffic_weight": weight_row["traffic_weight"],
            "airspace_weight": weight_row["airspace_weight"],
            "status": "prepared",
        }

        if route_runner:
            route_result = route_runner(experiment)
            experiment["status"] = "complete"
            experiment["route_result"] = route_result
        elif run_routes:
            if not route_pair:
                raise ValueError("route_pair is required when run_routes is True.")
            route_result = run_astar_route(
                city_name=city_name,
                graph_package=graph_package,
                graph_path=graph_path,
                grid_path=grid_path,
                origin=route_pair["origin"],
                destination=route_pair["destination"],
                weights={
                    "distance_weight": weight_row["distance_weight"],
                    "population_weight": weight_row["population_weight"],
                    "traffic_weight": weight_row["traffic_weight"],
                    "airspace_weight": weight_row["airspace_weight"],
                },
                route_id=weight_row["weight_id"],
            )
            experiment["status"] = "complete"
            experiment["route_result"] = route_result

        experiment_rows.append(experiment)

    return {
        "city": city_name,
        "route_pair": route_pair,
        "experiment_count": len(experiment_rows),
        "experiments": experiment_rows,
    }


if __name__ == "__main__":
    result = run_weighted_route_experiments()
    print(result)
