# Chicago Graph Output

## Routing setup summary

### Grid and graph model

- Source grid: `geojson/risk_grid_v5.geojson`
- Nodes: grid cells
- Edges: any two cells whose polygons touch, including diagonal neighbors
- Start point: Clow
- End point: Union Station

### Cost normalization

For the normalized routing runs, each edge cost was built from three normalized terms:

- `distance_norm = edge_distance_m / max_edge_distance_m`
- `population_norm = 0.5 * (city_risk_norm[u] + city_risk_norm[v])`
- `airspace_norm = 0.5 * (airspace_risk_norm[u] + airspace_risk_norm[v])`

Where:

- `city_risk_norm = city_risk / max(city_risk)`
- `airspace_risk_norm = airport_risk_combined / max(airport_risk_combined)`

### A* heuristic / distance formula

For the grid-based A* runs, the heuristic used Euclidean distance between cell centroids in Web Mercator meters (`EPSG:3857`):

```text
dx = cent_x[u] - cent_x[v]
dy = cent_y[u] - cent_y[v]
euclidean_distance_m = sqrt(dx^2 + dy^2)
```

The heuristic then scaled that Euclidean distance by the route's distance weight:

```text
h(u, v) = distance_weight * (euclidean_distance_m / max_edge_distance_m)
```

This means:

- if a route includes a distance term, A* uses a normalized Euclidean heuristic
- if a route has `distance_weight = 0`, the heuristic is `0`, which reduces A* to Dijkstra-style search for that mode while remaining admissible

### Historical weighting phases

The previous `chicago_graph.py` went through three practical phases:

1. Distance scaled by average total risk
   - Used in the earliest Dijkstra-only and early A* versions.
   - Edge weight:

```text
edge_distance_m = distance between neighboring cell centroids in EPSG:3857
risk_avg = 0.5 * (risk_cost[u] + risk_cost[v])
weight = edge_distance_m * (1.0 + risk_avg)
```

   - A* heuristic in that phase:

```text
h(u, v) = euclidean_distance_m
```

2. Explicit weighted distance plus risk
   - This is the phase that produced the large benchmark costs such as `3212055.1838673027`.
   - Constants:
     - `ALPHA = 1.0`
     - `BETA = 5.0`
   - Edge weight:

```text
edge_distance_m = sqrt(dx_edge^2 + dy_edge^2)
risk_avg = 0.5 * (risk_cost[u] + risk_cost[v])
weight = ALPHA * edge_distance_m + BETA * risk_avg
```

   - A* heuristic in that phase:

```text
h(u, v) = euclidean_distance_m
```

3. Normalized combined routing with equal weights
   - `DISTANCE_WEIGHT = 1.0`
   - `POPULATION_WEIGHT = 1.0`
   - `AIRSPACE_WEIGHT = 1.0`
   - Edge weight:

```text
weight = 1.0 * distance_norm
       + 1.0 * population_norm
       + 1.0 * airspace_norm
```

4. Current multi-route normalized routing
   - Combined route:
     - `DISTANCE_WEIGHT = 0.8`
     - `POPULATION_WEIGHT = 0.9`
     - `AIRSPACE_WEIGHT = 1.4`
   - Population-only route:
     - `distance_weight = 0.0`
     - `population_weight = 1.0`
     - `airspace_weight = 0.0`
   - Distance-only route:
     - `distance_weight = 1.0`
     - `population_weight = 0.0`
     - `airspace_weight = 0.0`
   - Airspace-only route:
     - `distance_weight = 0.0`
     - `population_weight = 0.0`
     - `airspace_weight = 1.0`

## Prior `chicago_graph.py` output

These were the historical outputs embedded in the previous version of `2D/chicago_graph.py`.

```text
Output:
Dijkstra seconds: 0.26059980000718497
A* seconds: 0.8576807000063127
Same path? False
Dijkstra cost: 3212055.1838673027
A* cost: 3212055.1838673023

After changing heuristic calculation:
Dijkstra seconds: 0.28847320000932086
A* seconds: 0.2662500000005821
Same path? False
Dijkstra cost: 3212055.1838673027
A* cost: 3212055.1838673023

# New performance metrics:
Nodes: 101286
Edges: 403235
Start node: 77708 End node: 87693
Path nodes: 114
A* path nodes: 114
Dijkstra seconds: 0.2512913000246044
A* seconds: 0.04089380000368692
Same path? True
Dijkstra cost: 103983.27083941635
A* cost: 103983.27083941635

After normalization:
Nodes: 101286
Edges: 403235
Normalization maxima:
  distance max edge (m): 707.1067811911573
  city_risk max: 120.0
  airport_risk_combined max: 186.6979748846326
Start node: 77708 End node: 87693
Path nodes: 114
A* path nodes: 114
Saved geojson/routes.geojson
Dijkstra seconds: 0.2288434000001871
A* seconds: 0.036822400000346533
Same path? True
Dijkstra cost: 147.42153282431707
A* cost: 147.42153282431707
```

## New A* costs

These are the outputs from the updated multi-route version of `2D/chicago_graph.py`.

### Weight formulas used in the current version

Combined route:

```text
weight_combined = 0.8 * distance_norm
                + 0.9 * population_norm
                + 1.4 * airspace_norm
```

Population-only route:

```text
weight_population_only = 1.0 * population_norm
```

Distance-only route:

```text
weight_distance_only = 1.0 * distance_norm
```

Airspace-only route:

```text
weight_airspace_only = 1.0 * airspace_norm
```

### Heuristic used in the current version

- Combined route:

```text
h(u, v) = 0.8 * (euclidean_distance_m / max_edge_distance_m)
```

- Population-only route:

```text
h(u, v) = 0
```

- Distance-only route:

```text
h(u, v) = 1.0 * (euclidean_distance_m / max_edge_distance_m)
```

- Airspace-only route:

```text
h(u, v) = 0
```

```text
Nodes: 101286
Edges: 403235
Normalization maxima:
  distance max edge (m): 707.1067811911573
  city_risk max: 120.0
  airport_risk_combined max: 186.6979748846326
Route weights:
  combined: 0.8 0.9 1.4
Start node: 77708 End node: 87693
combined: nodes=114 cost=130.5409 seconds=0.2810
population_only: nodes=124 cost=33.0211 seconds=0.4194
distance_only: nodes=110 cost=93.1838 seconds=0.0194
airspace_only: nodes=659 cost=0.0144 seconds=0.3632
Saved C:\Users\Owner\PycharmProjects\EVOTL-MAP\geojson\routes.geojson
```
