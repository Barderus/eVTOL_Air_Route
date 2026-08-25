# St. Louis Route Results

## What Is This Project

This folder summarizes the St. Louis A* route outputs and cost analysis. The
workflow compares three origin-destination pairs across six OpenSky traffic
dates:

- MidAmerica St. Louis Airport to St. Louis Union Station
- MidAmerica St. Louis Airport to St. Louis Lambert International Airport
- St. Louis Downtown Airport to St. Louis Lambert International Airport

The routes use a `500 m` grid and a combined edge cost from four factors:

- distance
- population density
- airport airspace
- observed air traffic

The combined edge cost is:

```text
0.6 * distance
+ 0.9 * population
+ 1.4 * airspace
+ 1.0 * traffic
```

These are weighted model scores, not physical distance or monetary cost.

## How To Set It Up

Install dependencies from the repository root:

```powershell
uv sync
```

Route definitions, traffic datasets, and weights are configured in:

```text
St Louis/settings.py
```

## How To Regenerate Results

Run the route and cost-analysis scripts from the `St Louis/` folder:

```powershell
uv run python astar_routes.py
uv run python make_route_cost_pies.py
uv run python generate_astar_toggle_pages.py
uv run python make_route_traffic_map.py
```

## Route Cost Summary

Average route cost and distance across the six traffic dates:

| Route | Average cost | Average distance |
| --- | ---: | ---: |
| MidAmerica to St. Louis Union Station | 93.448 | 51.60 km |
| MidAmerica to St. Louis Lambert | 158.239 | 76.53 km |
| St. Louis Downtown Airport to St. Louis Lambert | 107.195 | 35.13 km |

Route cost by traffic date:

| Route | 2025-07-12 | 2025-07-14 | 2026-01-10 | 2026-01-12 | 2026-03-07 | 2026-03-09 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| MidAmerica to St. Louis Union Station | 97.475 | 84.416 | 94.159 | 98.586 | 99.271 | 86.783 |
| MidAmerica to St. Louis Lambert | 161.195 | 147.074 | 151.402 | 167.005 | 164.157 | 158.598 |
| St. Louis Downtown Airport to St. Louis Lambert | 111.294 | 104.007 | 106.072 | 109.604 | 106.189 | 106.002 |

Average weighted factor cost by route:

| Route | Distance | Population | Airspace | Traffic |
| --- | ---: | ---: | ---: | ---: |
| MidAmerica to St. Louis Union Station | 34.073 (36.5%) | 13.044 (14.0%) | 33.040 (35.4%) | 13.292 (14.2%) |
| MidAmerica to St. Louis Lambert | 50.537 (31.9%) | 19.719 (12.5%) | 70.233 (44.4%) | 17.749 (11.2%) |
| St. Louis Downtown Airport to St. Louis Lambert | 23.199 (21.6%) | 21.548 (20.1%) | 47.367 (44.2%) | 15.081 (14.1%) |

Across all route and date combinations, the average cost distribution is:

| Factor | Average cost | Share |
| --- | ---: | ---: |
| Distance | 35.936 | 30.0% |
| Population | 18.104 | 15.1% |
| Airspace | 50.213 | 42.0% |
| Traffic | 15.374 | 12.9% |

Airspace is the largest average factor for both routes ending at Lambert.

## Output Files

The route workflow creates:

```text
routes/
|-- analysis/
|   |-- route_cost_breakdown_*.png
|-- output/
|   |-- downtown_airport_to_lambert_astar_routes.geojson
|   |-- midamerica_to_lambert_astar_routes.geojson
|   |-- midamerica_to_union_station_astar_routes.geojson
|   |-- st_louis_astar_routes.geojson
|   |-- st_louis_route_cost_analysis.csv
`-- README.md
```

`st_louis_astar_routes.geojson` contains geometry and metadata for all three
routes across the six traffic dates.

`st_louis_route_cost_analysis.csv` contains total route cost, route distance,
raw weighted factor costs, and factor percentages.

The three route-specific GeoJSON files each contain `30` route variants:

- six traffic dates
- combined weighting
- distance-only weighting
- population-only weighting
- airspace-only weighting
- air-traffic-only weighting

The related interactive route pages are stored in `St Louis/maps/`:

```text
downtown_airport_to_lambert_astar.html
midamerica_to_lambert_astar.html
midamerica_to_union_station_astar.html
st_louis_routes_24h_traffic.html
```

## Limitations

- Airport airspace is modeled with simplified radial zones.
- The route model is two-dimensional and does not optimize altitude.
- OpenSky traffic is historical and depends on receiver coverage.
- The current weights are experimental and require sensitivity testing.

## Future Improvements

- Compare these fixed-weight route results with the St. Louis route robustness
  outputs.
- Add route-stability summaries after the St. Louis cluster maps are reviewed.
- Revisit traffic normalization if route results are sensitive to a small
  number of high-density cells.
