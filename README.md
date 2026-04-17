# eVTOL Air Route

This project explores low-altitude routing around Chicago for eVTOL-style flights.

It combines five route views:

- `combined`
- `airspace_only`
- `flight_density_only`
- `population_only`
- `distance_only`

The idea is simple:

- population cost comes from census block-group density
- airspace cost comes from airport corridors and airport control radii
- flight-density cost comes from observed OpenSky traffic
- distance cost comes from Euclidean distance between neighboring grid cells

## Main Folders

- `2D/` builds the grid, routing graph, and 2D traffic maps
- `3D/` builds the 3D traffic-density map
- `opensky/` exports or stores OpenSky CSV data
- `geojson/` stores the grid and route outputs
- `html/` stores the map pages
- `Notes/` explains how the costs and routing work

## Main Workflow

### 1. Build the population layer

```powershell
python population\population_calculation.py
```

This creates:

- `geojson\il_blockgroups_population_density.geojson`

### 2. Build the risk grid

```powershell
python 2D\make_grid.py
```

This creates:

- `geojson\risk_grid_v5.geojson`

### 3. Generate the A* route pages

```powershell
uv run python 2D\generate_astar_toggle_pages.py
```

This creates GeoJSON and HTML outputs for:

- Clow to Union Station
- Clow to O'Hare
- Clow to Midway

The generated HTML pages are:

- `html\clow_to_union_station_astar.html`
- `html\clow_to_ohare_astar.html`
- `html\clow_to_midway_astar.html`

## Optional Traffic Visualizations

### 2D density plot

```powershell
uv run python 2D\make_density_plot.py opensky\output\ohare_2019-03-09_1s_15nm_bbox.csv
```

### 2D Leaflet traffic map

```powershell
uv run python 2D\make_leaflet_map.py opensky\output\ohare_2019-03-09_1s_15nm_bbox.csv
```

### 3D traffic-density map

```powershell
uv run python 3D\generate_3d_density_map.py opensky\output\ohare_2019-03-09_1s_15nm_bbox.csv
```

## Viewing The HTML Files

Run a local server from the repo root:

```powershell
python -m http.server 8080
```

Then open pages like:

```text
http://localhost:8080/html/clow_to_ohare_astar.html
http://localhost:8080/html/ohare_density_leaflet.html
http://localhost:8080/html/ohare_3d_density_map.html
```

## Notes

- `Notes\astar_cost_breakdown.md` explains how each routing cost is calculated.
- `Notes\distance_cost_explained.md` explains why the distance route is not a single straight line.
- The route pages are precomputed and embedded directly in the HTML for fast switching between dates and route types.
