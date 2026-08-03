"""Cluster similar routes for robustness workflows."""


def cluster_routes(routes=None, similarities=None, threshold=0.25):
    """Cluster routes by connecting route pairs above a similarity threshold."""
    if routes is None:
        return {"clusters": [], "message": "Pass routes and similarities to cluster routes."}

    route_ids = [
        route.get("route_run_id") or route.get("route_id") or f"route_{index + 1:04d}"
        for index, route in enumerate(routes)
    ]
    parent = {route_id: route_id for route_id in route_ids}

    def find(route_id):
        while parent[route_id] != route_id:
            parent[route_id] = parent[parent[route_id]]
            route_id = parent[route_id]
        return route_id

    def union(route_a, route_b):
        root_a = find(route_a)
        root_b = find(route_b)
        if root_a != root_b:
            parent[root_b] = root_a

    for similarity_row in similarities or []:
        if similarity_row["similarity"] >= threshold:
            union(similarity_row["route_a"], similarity_row["route_b"])

    root_to_cluster = {}
    clustered_rows = []
    for route_id in route_ids:
        root = find(route_id)
        if root not in root_to_cluster:
            root_to_cluster[root] = f"cluster_{len(root_to_cluster) + 1:03d}"
        clustered_rows.append(
            {
                "route_id": route_id,
                "cluster_id": root_to_cluster[root],
            }
        )

    return {
        "clusters": clustered_rows,
        "threshold": threshold,
        "cluster_count": len(set(row["cluster_id"] for row in clustered_rows)),
    }


if __name__ == "__main__":
    result = cluster_routes()
    print(result)

