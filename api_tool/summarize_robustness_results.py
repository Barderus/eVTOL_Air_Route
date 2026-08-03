"""Summarize route-robustness cluster results."""


def summarize_robustness_results(clusters=None, representatives=None):
    """Return route cluster counts and representative route summary rows."""
    if clusters is None:
        return {
            "cluster_summary": [],
            "representatives": representatives or [],
            "message": "Pass clusters to summarize robustness results.",
        }

    counts = {}
    for row in clusters:
        cluster_id = row["cluster_id"]
        counts[cluster_id] = counts.get(cluster_id, 0) + 1

    total_routes = sum(counts.values())
    representative_by_cluster = {
        row["cluster_id"]: row.get("representative_route_id")
        for row in representatives or []
    }

    summary_rows = []
    for cluster_id, route_count in sorted(counts.items()):
        summary_rows.append(
            {
                "cluster_id": cluster_id,
                "route_count": route_count,
                "weight_space_percent": round((route_count / total_routes) * 100, 2)
                if total_routes
                else 0,
                "representative_route_id": representative_by_cluster.get(cluster_id),
            }
        )

    return {
        "cluster_summary": summary_rows,
        "representatives": representatives or [],
        "total_routes": total_routes,
    }


if __name__ == "__main__":
    result = summarize_robustness_results()
    print(result)

