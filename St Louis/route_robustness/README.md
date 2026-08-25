# St. Louis Route Robustness

This folder extends the St. Louis routing workflow with normalized weight-grid
experiments and barycentric weight-space visualizations.

The first pass uses the Chicago-style four-factor normalized weight table:

- distance
- population
- traffic
- airspace

The default grid step is `0.10`, plus the explicit equal-weight case. That
creates `287` weight configurations.

## Current Barycentric Pass

The initial St. Louis pass runs the `287` weights for each configured
origin-destination pair in `St Louis/settings.py`, using the first configured
traffic date by default: `2026-03-07`.

Run from the repository root:

```powershell
$env:UV_CACHE_DIR = (Join-Path (Get-Location) '.uv-cache')
uv run python "St Louis/route_robustness/scripts/run_weighted_routes.py"
uv run python "St Louis/route_robustness/scripts/cluster_weighted_routes.py"
uv run python "St Louis/route_robustness/scripts/make_weight_barycentric_plots.py"
uv run python "St Louis/route_robustness/scripts/make_route_cluster_maps.py"
```

Generated data:

```text
St Louis/route_robustness/output/st_louis_weight_configurations.csv
St Louis/route_robustness/output/st_louis_weighted_route_runs.csv
St Louis/route_robustness/output/st_louis_weighted_routes.geojson
St Louis/route_robustness/output/st_louis_route_clusters_dbscan.csv
St Louis/route_robustness/output/st_louis_route_clusters_edit_distance.csv
St Louis/route_robustness/output/st_louis_route_clusters_frechet.csv
St Louis/route_robustness/output/st_louis_route_clusters_hierarchical_jaccard.csv
St Louis/route_robustness/output/st_louis_route_weight_barycentric_coordinates.csv
```

Generated pre-clustering visualizations:

```text
St Louis/route_robustness/tetrahedrons_barycentric/midamerica_to_st_louis_union_station_weight_barycentric.html
St Louis/route_robustness/tetrahedrons_barycentric/midamerica_to_st_louis_lambert_weight_barycentric.html
St Louis/route_robustness/tetrahedrons_barycentric/st_louis_downtown_airport_to_st_louis_lambert_weight_barycentric.html
```

Generated method-specific barycentric visualizations:

```text
St Louis/route_robustness/tetrahedrons_barycentric/*_weight_barycentric_dbscan.html
St Louis/route_robustness/tetrahedrons_barycentric/*_weight_barycentric_edit_distance.html
St Louis/route_robustness/tetrahedrons_barycentric/*_weight_barycentric_frechet.html
St Louis/route_robustness/tetrahedrons_barycentric/*_weight_barycentric_hierarchical_jaccard.html
```

The barycentric plots color points by exact route variant, using the route's
grid-cell path as the temporary variant signature. This is a pre-clustering
view that shows how the tested weight space maps to route changes before
formal clustering is applied.

The method-specific barycentric plots color points by route clusters from:

- DBSCAN using Frechet distance
- average-linkage hierarchical clustering using edit distance
- average-linkage hierarchical clustering using Frechet distance
- average-linkage hierarchical clustering using Jaccard distance

Generated route-cluster Leaflet maps:

```text
St Louis/route_robustness/route_clusters/midamerica_to_st_louis_lambert_route_clusters.html
St Louis/route_robustness/route_clusters/midamerica_to_st_louis_union_station_route_clusters.html
St Louis/route_robustness/route_clusters/st_louis_downtown_airport_to_st_louis_lambert_route_clusters.html
```

Each route-cluster map fixes one origin-destination pair and lets the viewer
toggle between DBSCAN, edit-distance, Frechet, and hierarchical-Jaccard
clusters. The cluster list on the right toggles individual route clusters on
and off.
