"""Compare route similarity for robustness workflows."""


def parse_node_sequence(value):
    """Convert a route node sequence into a set."""
    if value is None:
        return set()
    if isinstance(value, str):
        cleaned = value.replace("[", "").replace("]", "").replace("'", "").replace('"', "")
        return {item.strip() for item in cleaned.split(",") if item.strip()}
    return {str(item) for item in value}


def jaccard_similarity(route_a_nodes, route_b_nodes):
    """Calculate shared-node Jaccard similarity."""
    nodes_a = parse_node_sequence(route_a_nodes)
    nodes_b = parse_node_sequence(route_b_nodes)
    union_nodes = nodes_a | nodes_b
    if not union_nodes:
        return 0.0
    return len(nodes_a & nodes_b) / len(union_nodes)


def compare_route_similarity(routes=None, node_field="path_node_ids"):
    """Return route-pair Jaccard similarity rows."""
    if routes is None:
        return {
            "similarities": [],
            "message": "Pass routes to compare route similarity.",
        }

    similarities = []
    for first_index, route_a in enumerate(routes):
        for second_index in range(first_index + 1, len(routes)):
            route_b = routes[second_index]
            similarities.append(
                {
                    "route_a": route_a.get("route_run_id") or route_a.get("route_id") or first_index,
                    "route_b": route_b.get("route_run_id") or route_b.get("route_id") or second_index,
                    "similarity": jaccard_similarity(
                        route_a.get(node_field),
                        route_b.get(node_field),
                    ),
                    "method": "shared_node_jaccard",
                }
            )

    return {"similarities": similarities, "method": "shared_node_jaccard"}


if __name__ == "__main__":
    result = compare_route_similarity()
    print(result)

