"""Select representative routes for robustness clusters."""


def select_representative_routes(routes=None, clusters=None, similarities=None):
    """Select the route with highest average within-cluster similarity."""
    if routes is None or clusters is None:
        return {
            "representatives": [],
            "message": "Pass routes, clusters, and similarities to select representatives.",
        }

    route_ids = [
        route.get("route_run_id") or route.get("route_id") or f"route_{index + 1:04d}"
        for index, route in enumerate(routes)
    ]
    cluster_by_route = {row["route_id"]: row["cluster_id"] for row in clusters}
    similarity_lookup = {}
    for row in similarities or []:
        similarity_lookup[(row["route_a"], row["route_b"])] = row["similarity"]
        similarity_lookup[(row["route_b"], row["route_a"])] = row["similarity"]

    cluster_ids = sorted(set(cluster_by_route.values()))
    representatives = []
    for cluster_id in cluster_ids:
        members = [route_id for route_id in route_ids if cluster_by_route.get(route_id) == cluster_id]
        best_route_id = members[0]
        best_average = -1

        for route_id in members:
            if len(members) == 1:
                average_similarity = 1.0
            else:
                scores = [
                    similarity_lookup.get((route_id, other_route_id), 0.0)
                    for other_route_id in members
                    if other_route_id != route_id
                ]
                average_similarity = sum(scores) / len(scores)

            if average_similarity > best_average:
                best_average = average_similarity
                best_route_id = route_id

        representatives.append(
            {
                "cluster_id": cluster_id,
                "representative_route_id": best_route_id,
                "member_count": len(members),
                "average_similarity": best_average,
            }
        )

    return {"representatives": representatives}


if __name__ == "__main__":
    result = select_representative_routes()
    print(result)

