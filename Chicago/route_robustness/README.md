# Chicago Route Robustness

The baseline route pair is:

```text
Clow International Airport to Chicago Union Station
```

The additional airport route pairs are:

```text
Clow International Airport to O'Hare International Airport
Clow International Airport to Midway International Airport
```

## Weight Grid

The normalized weight grid uses `0.10` increments for four factors:

- distance
- population
- traffic
- airspace

The systematic grid has `286` valid combinations. The exact equal-weight case
is added separately:

```text
Distance = 0.25
Population = 0.25
Traffic = 0.25
Airspace = 0.25
```

This makes `287` weight configurations for each route pair.

## Baseline Union Station Run

The original full robustness run covers Clow to Chicago Union Station.

Main outputs:

```text
output/all_route_runs.csv
output/all_routes.geojson
output/all_route_clusters.csv
output/all_route_similarity_pairs.csv
output/all_route_similarity_matrix.csv
maps/all_routes.html
maps/all_route_clusters.html
maps/all_cluster_representative_routes.html
```

The first pass used shared path-node Jaccard similarity:

```text
similarity = shared_nodes / union_nodes
```

The `0.25` similarity threshold produced `12` clusters. A broader `0.15`
threshold was also tested and produced `7` clusters.

Representative routes use a medoid-style selection: the representative is the
real route with the highest average similarity to the other routes in the same
cluster.

## Direct-Route Clustering

The direct-route subset starts from selected route families in the `0.25`
full-route clustering result:

```text
cluster_001
cluster_003
cluster_004
cluster_006
cluster_012
```

In this context, direct routes are the route families that stay closest to the
practical Clow-to-Union Station corridor instead of taking large detours around
the study area. They were used as the main follow-up subset for three reasons:

- They were the largest route family in the full clustering result. The selected
  direct-route clusters produce `188` weight configurations, with `cluster_001`
  alone containing `150`.
- They were robust across different weight settings, so the analysis could focus
  on feasible corridor variation instead of route families that were indirect or
  less practical.
- They were the most similar to the initial route concept prepared for IDOT.

The direct-route clustering methods are:

| Method | Distance measure | Current result |
| --- | --- | ---: |
| Jaccard threshold | Shared route-node set distance | 5 clusters |
| DBSCAN | Discrete Frechet distance | 3 clusters |
| OPTICS | Discrete Frechet distance | 3 clusters |
| Hierarchical Jaccard | Average Jaccard distance | 10 clusters |
| Hierarchical Frechet | Average Frechet distance | 8 clusters |
| Hierarchical edit distance | Average normalized route-node edit distance | 27 clusters |

DBSCAN and OPTICS use:

```text
eps = 2.5 km
min_samples = 4
```

Frechet geometry is converted from longitude/latitude to approximate local
kilometers and resampled to at most `80` points before distance calculation.

Direct-route maps:

```text
maps/direct_route_clustering.html
maps/direct_route_clusters_dbscan.html
maps/direct_route_clusters_edit_distance.html
maps/direct_route_clusters_frechet.html
maps/direct_route_clusters_hierarchical_jaccard.html
maps/direct_route_clusters_jaccard.html
maps/direct_route_clusters_optics.html
```

## O'Hare And Midway Robustness Runs

The O'Hare and Midway runs use the same Chicago risk grid, traffic file,
normalization settings, weight grid, and A* cost model as the Union Station
run.

The additional-route outputs are kept separate so they do not replace the
existing Union Station artifacts:

```text
output/additional_route_runs.csv
output/additional_routes.geojson
```

Each route pair has `287` successful weighted routes:

| Route pair | Successful runs |
| --- | ---: |
| Clow to O'Hare | 287 |
| Clow to Midway | 287 |

Additional-route cluster outputs:

```text
output/additional_route_clusters_dbscan.csv
output/additional_route_clusters_edit_distance.csv
output/additional_route_clusters_frechet.csv
output/additional_route_clusters_hierarchical_jaccard.csv
```

Cluster counts:

| Route pair | DBSCAN | Edit distance | Frechet | Hierarchical Jaccard |
| --- | ---: | ---: | ---: | ---: |
| Clow to O'Hare | 30 | 36 | 25 | 15 |
| Clow to Midway | 20 | 45 | 17 | 16 |

Generated route-cluster maps:

```text
additional_route_clusters/clow_to_ohare_route_clusters.html
additional_route_clusters/clow_to_midway_route_clusters.html
```

## Barycentric Tetrahedron Views

The four route weights are projected into a 3D tetrahedron. Each point is one
weight configuration.

For the Union Station route, the direct-route barycentric output is:

```text
output/route_weight_barycentric_plot.html
```

For the O'Hare and Midway routes, the route-pair-specific tetrahedron outputs
are:

```text
additional_tetrahedrons_barycentric/clow_to_ohare_weight_barycentric.html
additional_tetrahedrons_barycentric/clow_to_midway_weight_barycentric.html
additional_tetrahedrons_barycentric/*_weight_barycentric_dbscan.html
additional_tetrahedrons_barycentric/*_weight_barycentric_edit_distance.html
additional_tetrahedrons_barycentric/*_weight_barycentric_frechet.html
additional_tetrahedrons_barycentric/*_weight_barycentric_hierarchical_jaccard.html
```

Exact route variants before clustering:

| Route pair | Exact route variants |
| --- | ---: |
| Clow to O'Hare | 148 |
| Clow to Midway | 224 |
