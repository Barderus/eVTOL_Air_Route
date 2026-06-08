# Chicago Graph Output

## Current routing setup

### Grid and graph model

- Source grid: `geojson/risk_grid_v6.geojson`
- Airport airspace source: `airport_risk_combined` from the grid
- Air traffic source: `opensky/output/ohare_2019-03-09_local_1s_15nm_bbox.csv`
- Nodes: grid cells
- Edges: any two cells whose polygons touch, including diagonal neighbors
- Start point: Clow
- End point: Union Station

### Risk layers used by the current script

The current `2D/chicago_graph.py` uses four independent route terms:

1. Distance
   - centroid-to-centroid edge distance in `EPSG:3857`
2. Population
   - `city_risk` from `geojson/risk_grid_v6.geojson`
3. Airport airspace
   - `airport_risk_combined` from `geojson/risk_grid_v6.geojson`
4. Air traffic density
   - derived from the `1s` OpenSky CSV by spatially joining flight observations into the grid and computing:
     - `traffic_count`
     - `traffic_density` in observations per square kilometer

This distinction matters:

- `airspace` means the airport-related static cost from the risk grid
- `air traffic density` means the observed flight-density layer from the `1s` traffic data

### Cost normalization

For the current normalized routing runs, each edge cost is built from four normalized terms:

- `distance_norm = edge_distance_m / max_edge_distance_m`
- `population_norm = 0.5 * (city_risk_norm[u] + city_risk_norm[v])`
- `airspace_norm = 0.5 * (airport_risk_norm[u] + airport_risk_norm[v])`
- `traffic_norm = 0.5 * (traffic_density_norm[u] + traffic_density_norm[v])`

Where:

- `city_risk_norm = city_risk / max(city_risk)`
- `airport_risk_norm = airport_risk_combined / max(airport_risk_combined)`
- `traffic_density_norm = traffic_density / max(traffic_density)`

### A* heuristic

For the grid-based A* runs, the heuristic uses Euclidean distance between cell centroids in Web Mercator meters (`EPSG:3857`):

```text
dx = cent_x[u] - cent_x[v]
dy = cent_y[u] - cent_y[v]
euclidean_distance_m = sqrt(dx^2 + dy^2)
```

The heuristic is scaled only by the route's distance weight:

```text
h(u, v) = distance_weight * (euclidean_distance_m / max_edge_distance_m)
```

If a route has `distance_weight = 0`, the heuristic is `0`, which reduces A* to Dijkstra-style search for that route mode.

## Current route definitions

Combined route:

```text
weight_combined = 0.8 * distance_norm
                + 0.9 * population_norm
                + 1.1 * airspace_norm
                + 1.4 * traffic_norm
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

Air-traffic-only route:

```text
weight_air_traffic_only = 1.0 * traffic_norm
```

## Verified run output

These are the outputs from the current five-route version of `2D/chicago_graph.py`.

```text
Nodes: 51984
Edges: 206513
Traffic density source:
  csv: C:\Users\Owner\PycharmProjects\EVOTL-MAP\opensky\output\ohare_2019-03-09_local_1s_15nm_bbox.csv
  input rows used: 1152330
  rows matched to grid: 1152330
  nonzero traffic cells: 25387
Normalization maxima:
  distance max edge (m): 707.1067811924744
  city_risk max: 120.0
  airport_risk_combined max: 185.83015359113148
  traffic_density max (obs/km^2): 116059.99999956765
Route weights:
  combined: 0.8 0.9 1.1 1.4
Start node: 27975 End node: 45997
combined: nodes=114 cost=126.2033 seconds=0.1448
population_only: nodes=124 cost=33.0211 seconds=0.2376
distance_only: nodes=110 cost=93.1838 seconds=0.0297
airspace_only: nodes=411 cost=0.0252 seconds=0.1450
air_traffic_only: nodes=326 cost=0.0042 seconds=0.1351
Saved C:\Users\Owner\PycharmProjects\EVOTL-MAP\geojson\routes.geojson
```

## Map outputs

- `geojson/routes.geojson` now contains both:
  - `airspace_only`
  - `air_traffic_only`
- `html/routes_map.html` exposes both as separate toggleable layers:
  - `Airspace Only`
  - `Air Flight Density Only`
