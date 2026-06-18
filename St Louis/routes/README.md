# St. Louis Route Analysis

## What Is This Project

This folder contains the generated A* routes and route cost analysis for the
St. Louis eVTOL study area. The workflow compares three origin-destination
pairs across six OpenSky traffic dates:

- MidAmerica St. Louis Airport to St. Louis Union Station
- MidAmerica St. Louis Airport to St. Louis Lambert International Airport
- St. Louis Downtown Airport to St. Louis Lambert International Airport

The routes use a 500-meter grid. Each graph edge receives a combined cost from
four factors:

- Distance
- Population density
- Airport airspace
- Observed air traffic

This is an experimental routing model. The results show how the current cost
assumptions affect route selection, but they are not operational flight plans.

## Cost Calculation

The combined edge cost is:

```text
0.6 * distance
+ 0.9 * population
+ 1.4 * airspace
+ 1.0 * traffic
```

Distance is charged on every edge. Population, airspace, and traffic use the
average normalized risk of the two cells connected by the edge.

Airspace risk uses the simplified airport zones stored in the grid:

- High-risk airport cores use a risk value of 100.
- The St. Louis Lambert outer shelf uses a medium risk value of 60.
- Cells outside the modeled zones use an airspace risk value of 0.

Traffic risk comes from 1-second OpenSky observations counted inside each grid
cell. The raw counts contain a small number of extreme hotspots, so traffic is
normalized with:

1. `log1p` transformation
2. 95th-percentile clipping across occupied cells
3. A power transformation of `0.75`

This keeps one unusually busy cell from reducing most traffic values to nearly
zero. The route generator and cost analysis use the same normalized values.

## Current Results

The total weighted cost for each route and traffic date is:

| Route | 2025-07-12 | 2025-07-14 | 2026-01-10 | 2026-01-12 | 2026-03-07 | 2026-03-09 | Average cost | Average distance |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MidAmerica to St. Louis Union Station | 97.475 | 84.416 | 94.159 | 98.586 | 99.271 | 86.783 | 93.448 | 51.60 km |
| MidAmerica to St. Louis Lambert | 161.195 | 147.074 | 151.402 | 167.005 | 164.157 | 158.598 | 158.239 | 76.53 km |
| St. Louis Downtown Airport to St. Louis Lambert | 111.294 | 104.007 | 106.072 | 109.604 | 106.189 | 106.002 | 107.195 | 35.13 km |

These values come from
`routes/output/st_louis_route_cost_analysis.csv`. They are weighted model
scores, not monetary costs.

The following table averages each weighted factor cost across the six traffic
dates. The percentage in parentheses is that factor's share of the route's
average total cost:

| Route | Distance | Population | Airspace | Traffic |
| --- | ---: | ---: | ---: | ---: |
| MidAmerica to St. Louis Union Station | 34.073 (36.5%) | 13.044 (14.0%) | 33.040 (35.4%) | 13.292 (14.2%) |
| MidAmerica to St. Louis Lambert | 50.537 (31.9%) | 19.719 (12.5%) | 70.233 (44.4%) | 17.749 (11.2%) |
| St. Louis Downtown Airport to St. Louis Lambert | 23.199 (21.6%) | 21.548 (20.1%) | 47.367 (44.2%) | 15.081 (14.1%) |

Across all 18 route and date combinations, the average distribution is:

| Factor | Average cost | Share of average total |
| --- | ---: | ---: |
| Distance | 35.936 | 30.0% |
| Population | 18.104 | 15.1% |
| Airspace | 50.213 | 42.0% |
| Traffic | 15.374 | 12.9% |

Airspace is the largest average factor for both routes ending at Lambert.

Traffic changes with the selected OpenSky date. Its contribution ranges from:

- 3.5% to 19.9% for MidAmerica to Union Station
- 6.3% to 17.3% for MidAmerica to Lambert
- 11.8% to 17.5% for Downtown Airport to Lambert

The total weighted cost is not a physical distance or monetary value. It is
the sum of the normalized, weighted edge penalties used by A*.

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
|-- README.md
```

`st_louis_astar_routes.geojson` contains the route geometry and basic route
metadata for all three routes and six dates.

`st_louis_route_cost_analysis.csv` contains:

- Total route cost
- Route distance in kilometers
- Raw weighted cost for each factor
- Percentage contribution for each factor

`route_cost_breakdown_average.png` contains one pie chart for each configured
route. Each pie uses the average distance, population, airspace, and traffic
cost across the six traffic dates. The percentages are calculated from those
average factor costs, so each pie shows the average contribution of every
factor to that route's total cost.

The three route-specific GeoJSON files each contain 30 route variants:

- Six traffic dates
- Combined weighting
- Distance-only weighting
- Population-only weighting
- Airspace-only weighting
- Air-traffic-only weighting

The related interactive pages are stored in `St Louis/maps/`:

```text
downtown_airport_to_lambert_astar.html
midamerica_to_lambert_astar.html
midamerica_to_union_station_astar.html
st_louis_routes_24h_traffic.html
```

## How To Set It Up

Use Python 3.12 or newer and install the project dependencies from the
repository root:

```powershell
uv sync
```

Before generating routes, confirm that these files exist:

```text
St Louis/maps/st_louis_risk_grid.geojson
St Louis/traffic/output/st_louis_*.csv
```

Route definitions, traffic datasets, and weights are configured in
`St Louis/settings.py`.

## How To Run It

Use `St Louis/` as the working directory.

Generate the routes and analysis CSV:

```powershell
uv run python astar_routes.py
```

Generate the average cost breakdown chart:

```powershell
uv run python make_route_cost_pies.py
```

Generate the interactive route map:

```powershell
uv run python make_route_map.py
```

Generate the single-factor comparison pages and route-specific GeoJSON:

```powershell
uv run python generate_astar_toggle_pages.py
```

Generate the route overlay on the full 24-hour traffic-density heatmap:

```powershell
uv run python make_route_traffic_map.py
```

The combined route map is a separate output. Cost percentages are stored in
the CSV and pie charts rather than added to the maps.

## Limitations

- Airport airspace is modeled with simplified radial zones rather than full
  FAA airspace geometry and altitude shelves.
- The route model is two-dimensional and does not optimize altitude.
- OpenSky traffic is historical and depends on receiver coverage.
- One-second observations can favor aircraft that remain inside a cell longer.
- The 24-hour overlay aggregates traffic for an entire day rather than a
  selected departure time.
- The current weights are experimental and require sensitivity testing.
- A lower-cost grid route is not automatically a safe or legal flight route.

## Future Improvements

- Test several distance, population, airspace, and traffic weight combinations.
- Add published Class B, Class C, Class D, and restricted-airspace geometry.
- Separate traffic by hour so routes can match a planned departure time.
- Account for aircraft altitude when calculating traffic conflicts.
- Compare raw observation counts with unique-aircraft counts per cell.
- Add route smoothing and turn penalties after the cost model is validated.
- Review route changes against straight-line distance and known flight
  corridors.
