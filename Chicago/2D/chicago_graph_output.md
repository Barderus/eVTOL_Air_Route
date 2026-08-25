# Chicago 2D Route Output

## Current Status

The active Chicago 2D route workflow is `generate_astar_toggle_pages.py`. It builds A* route comparisons from Clow International Airport to three destinations:

- Union Station
- O'Hare International Airport
- Midway International Airport

The older `chicago_graph.py` script is still useful as a smaller one-file route builder. It currently runs the same cost model for Clow to Union Station and writes a single GeoJSON file.

## Grid And Graph Model

- Source grid: `Chicago/geojson/risk_grid_v7.geojson`
- Nodes: grid cells
- Edges: touching cells, including diagonal neighbors
- Static airspace source: `airport_risk_combined`
- Population source: `city_risk`
- Air traffic source: OpenSky `1s` CSV observations joined into the grid

The graph connects adjacent grid cells and routes through centroid-to-centroid line segments. This is a planning and comparison model, not an operational flight path.

## Route Cost Terms

The current route scripts combine four normalized terms:

```text
distance_norm
population_norm
airspace_norm
traffic_norm
```

The fixed combined-route weights are:

```text
weight_combined = 0.6 * distance_norm
                + 0.9 * population_norm
                + 1.4 * airspace_norm
                + 1.0 * traffic_norm
```

These are fixed A* route-generation weights. They are separate from the route robustness experiments, where normalized weight configurations sum to `1.0` and are clustered after route generation.

## Route Modes

The 2D A* workflow exports five route modes for each destination and traffic date:

- `combined`
- `airspace_only`
- `population_only`
- `distance_only`
- `air_traffic_only`

The single-factor routes isolate one cost layer at a time. The combined route uses all four cost layers with the fixed weights above.

## Normalization

Distance is normalized by the maximum neighbor edge distance. Population, airspace, and traffic are normalized before being assigned to graph edges.

The current scripts use percentile clipping and power scaling where needed so dense airspace and traffic values do not dominate the route search:

- `city_risk`: max normalization
- `airport_risk_combined`: 99th percentile clipping with power scaling
- `traffic_density`: `log1p`, 95th percentile clipping, and power scaling

Each edge receives the average normalized cell cost from its two endpoint cells.

## A* Heuristic

The A* heuristic uses Euclidean distance between cell centroids in Web Mercator meters (`EPSG:3857`):

```text
h(u, v) = distance_weight * (euclidean_distance_m / max_edge_distance_m)
```

When a route mode has `distance_weight = 0`, the heuristic is `0`, so the route search behaves like Dijkstra-style search for that route mode.

## Traffic Datasets

The active toggle-page workflow uses these OpenSky traffic files:

- `Chicago/opensky/output/ohare_2026-03-07_1s_15nm_bbox.csv`
- `Chicago/opensky/output/ohare_2026-03-09_1s_15nm_bbox.csv`
- `Chicago/opensky/output/ohare_2026-01-10_1s_15nm_bbox.csv`
- `Chicago/opensky/output/ohare_2026-01-12_1s_15nm_bbox.csv`
- `Chicago/opensky/output/ohare_2025-07-14_1s_15nm_bbox.csv`
- `Chicago/opensky/output/ohare_2025-07-12_1s_15nm_bbox.csv`

## Outputs

The active multi-destination workflow writes one GeoJSON and one self-contained Leaflet page per destination:

- `Chicago/geojson/clow_to_union_station_astar_routes.geojson`
- `Chicago/geojson/clow_to_ohare_astar_routes.geojson`
- `Chicago/geojson/clow_to_midway_astar_routes.geojson`
- `Chicago/html/clow_to_union_station_astar.html`
- `Chicago/html/clow_to_ohare_astar.html`
- `Chicago/html/clow_to_midway_astar.html`

The standalone `chicago_graph.py` script writes:

- `Chicago/geojson/routes.geojson`

## Relationship To Route Robustness

This file documents the 2D A* route output layer. The route robustness work lives under `Chicago/route_robustness/` and studies how route geometry changes across many normalized weight configurations and clustering methods.
