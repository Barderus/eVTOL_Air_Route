"""Select representative routes for the 0.15-threshold clusters."""

import os

import geopandas as gpd
import pandas as pd


ALL_ROUTES_GEOJSON = "route_robustness/output/all_routes.geojson"
ALL_ROUTE_CLUSTERS_CSV = "route_robustness/output/all_route_clusters_threshold_015.csv"
SIMILARITY_MATRIX_CSV = "route_robustness/output/all_route_similarity_matrix.csv"
OUTPUT_FOLDER = "route_robustness/output"
REPRESENTATIVES_CSV = os.path.join(
    OUTPUT_FOLDER,
    "all_cluster_representative_routes_threshold_015.csv",
)
REPRESENTATIVES_GEOJSON = os.path.join(
    OUTPUT_FOLDER,
    "all_cluster_representative_routes_threshold_015.geojson",
)


def load_inputs():
    """Load routes, 0.15 clusters, and the similarity matrix."""
    for file_path in [
        ALL_ROUTES_GEOJSON,
        ALL_ROUTE_CLUSTERS_CSV,
        SIMILARITY_MATRIX_CSV,
    ]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Required input file not found: {file_path}")

    routes = gpd.read_file(ALL_ROUTES_GEOJSON)
    clusters = pd.read_csv(ALL_ROUTE_CLUSTERS_CSV)
    similarity_matrix = pd.read_csv(SIMILARITY_MATRIX_CSV, index_col=0)

    if "path_node_ids" in routes.columns:
        routes = routes.drop(columns=["path_node_ids"])

    cluster_columns = [
        "route_run_id",
        "cluster_id",
        "cluster_size",
        "similarity_threshold",
        "dominant_weight_category",
    ]
    routes = routes.drop(
        columns=[
            column
            for column in [
                "cluster_id",
                "cluster_size",
                "similarity_threshold",
                "dominant_weight_category",
            ]
            if column in routes.columns
        ]
    )
    routes = routes.merge(
        clusters[cluster_columns],
        on="route_run_id",
        how="left",
    )
    return routes, similarity_matrix


def select_representatives(routes, similarity_matrix):
    """Select the medoid route for every 0.15-threshold cluster."""
    representative_rows = []
    representative_geometries = []
    total_routes = len(routes)

    for cluster_id, cluster_routes in routes.groupby("cluster_id", sort=True):
        route_ids = cluster_routes["route_run_id"].tolist()
        if len(route_ids) == 1:
            representative_route_id = route_ids[0]
            mean_similarity = 1.0
        else:
            cluster_similarity = similarity_matrix.loc[route_ids, route_ids]
            mean_similarities = cluster_similarity.mean(axis=1)
            representative_route_id = mean_similarities.sort_values(
                ascending=False,
            ).index[0]
            mean_similarity = float(mean_similarities.loc[representative_route_id])

        representative = cluster_routes[
            cluster_routes["route_run_id"] == representative_route_id
        ].iloc[0]
        representative_row = representative.drop(labels=["geometry"]).to_dict()
        representative_row["representative_route_run_id"] = representative_route_id
        representative_row["representative_mean_similarity"] = mean_similarity
        representative_row["cluster_weight_space_percent"] = (
            100.0 * int(representative_row["cluster_size"]) / total_routes
        )
        representative_rows.append(representative_row)
        representative_geometries.append(representative.geometry)

    representatives = gpd.GeoDataFrame(
        representative_rows,
        geometry=representative_geometries,
        crs=routes.crs,
    )
    representatives = representatives.sort_values(
        by=["cluster_size", "cluster_id"],
        ascending=[False, True],
    )
    return representatives


def main():
    """Save one representative route per 0.15-threshold cluster."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    routes, similarity_matrix = load_inputs()
    representatives = select_representatives(routes, similarity_matrix)

    representatives.drop(columns=["geometry"]).to_csv(REPRESENTATIVES_CSV, index=False)
    representatives.to_file(REPRESENTATIVES_GEOJSON, driver="GeoJSON")

    print(f"Saved 0.15 representative CSV: {REPRESENTATIVES_CSV}")
    print(f"Saved 0.15 representative GeoJSON: {REPRESENTATIVES_GEOJSON}")
    print(
        representatives[
            [
                "cluster_id",
                "cluster_size",
                "weight_id",
                "route_distance_km",
                "representative_mean_similarity",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
