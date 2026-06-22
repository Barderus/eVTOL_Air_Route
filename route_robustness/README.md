# Route Robustness

## What Is This Project

This folder contains the route robustness and weight sensitivity analysis.

The goal is to test how much the selected eVTOL route changes when the weights
for distance, population exposure, flight density, and airspace cost are
adjusted.

## Current Status

The workflow currently includes:

- normalized factor summaries
- normalized weight configurations
- a 10-route pilot
- all 287 weighted routes
- route similarity tables
- route clusters at multiple thresholds
- representative routes for each cluster
- HTML maps for pilot routes, all routes, clustered routes, and representative
  cluster routes

## Weight Normalization Note

Weight combinations are generated from integer units before being converted to
decimal weights. For example, the initial `0.10` grid uses units that add up to
`10`, then divides each unit by `10`.

The saved `weight_sum` value is rounded for clean CSV output. This avoids
floating-point display noise such as `0.9999999999999999`, while still keeping
the actual validation based on numeric closeness to `1.0`.

The initial systematic weight grid uses `0.10` increments and produces `286`
valid combinations. The exact equal-weight case is added as a separate special
row:

```text
Distance = 0.25
Population = 0.25
Traffic = 0.25
Airspace = 0.25
```

This makes the first weight table `287` rows total: `286` grid rows plus `1`
equal-weight row.

## Route Similarity Method

Routes are compared using shared path-node Jaccard similarity:

```text
similarity = shared_nodes / union_nodes
```

This means two routes are considered similar when they use many of the same
routing-grid nodes.

The first version uses path-node similarity because it is simple, reproducible,
and works directly from the saved `path_node_ids` column in the route output
tables.

## Pilot Clustering

The pilot used 10 route configurations:

- distance only
- population only
- traffic only
- airspace only
- equal weight
- distance heavy
- population heavy
- traffic heavy
- airspace heavy
- mixed balanced

The pilot threshold is:

```text
0.25
```

At `0.25`, the main mixed-route pilot cluster became:

```text
distance_heavy
equal_weight
mixed_balanced
population_heavy
```

That threshold was used as the first full-route clustering baseline.

## Full 287-Route Run

The `all_routes.html` map displays all 287 routes grouped by dominant weight
category:

- distance dominant
- population dominant
- traffic dominant
- airspace dominant
- tied dominant
- equal weight

This map is useful for seeing how the tested weight space spreads across the
study area before cluster labels are applied.

## Full Clustering At Threshold 0.25

The first full clustering pass used:

```text
similarity threshold = 0.25
```

This produced:

```text
12 clusters
```

Cluster sizes:

| Cluster | Routes | Share of tested weights |
| --- | ---: | ---: |
| cluster_001 | 150 | 52.3% |
| cluster_002 | 89 | 31.0% |
| cluster_003 | 28 | 9.8% |
| cluster_004 | 7 | 2.4% |
| cluster_005 | 5 | 1.7% |
| cluster_006 | 2 | 0.7% |
| cluster_007 to cluster_012 | 1 each | 0.3% each |

This is the stricter full clustering style. It keeps more small route families
separate and is useful when the analysis needs stronger evidence that routes
share the same corridor.

## Representative Routes At Threshold 0.25

Each cluster is also reduced to one representative route.

The representative route is selected using a medoid-style approach:

1. For every route in a cluster, calculate its average similarity to every
   other route in the same cluster.
2. Select the route with the highest average within-cluster similarity.
3. Use that real route geometry as the cluster representative.

This avoids creating an artificial merged geometry that may not correspond to
a valid A* route.

## Broader Clustering At Threshold 0.15

After reviewing the `0.25` result, a broader threshold was tested:

```text
similarity threshold = 0.15
```

This produced:

```text
7 clusters
```

Cluster sizes:

| Cluster | Routes | Share of tested weights |
| --- | ---: | ---: |
| cluster_001 | 187 | 65.2% |
| cluster_002 | 89 | 31.0% |
| cluster_003 | 5 | 1.7% |
| cluster_004 | 2 | 0.7% |
| cluster_005 | 2 | 0.7% |
| cluster_006 | 1 | 0.3% |
| cluster_007 | 1 | 0.3% |

The `0.15` threshold is the broader clustering style. It combines more routes
into the main route family while still preserving a second major route family
of 89 routes.

The `0.15` threshold is currently the best broader option because it reduces
fragmentation without merging the two largest route families.

## Representative Routes At Threshold 0.15

The same medoid-style representative method was used for the `0.15` clusters.
