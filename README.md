# eVTOL Air Route

`eVTOL Air Route` is a research-oriented Python project for modeling low-altitude airspace around Chicago for eVTOL and urban air mobility routing. It combines:

- population density
- airport-related airspace risk
- observed aircraft traffic from OpenSky

## Project Layout

- `2D/`: active 2D grid, routing, and 2D traffic visualization scripts
- `3D/`: active 3D traffic visualization scripts
- `opensky/`: OpenSky export/query utilities and CSV outputs
- `geojson/`: shared GeoJSON assets and routing/grid outputs
- `html/`: hand-authored and generated HTML map views
- `images/icons/`: SVG icons
- `images/visualizations/`: PNG visualizations
- `population/`: population source files and shapefiles

## Requirements

- Python `>=3.12,<3.15`
- PowerShell examples assume Windows
- Java installed if you want to export OpenSky data through Trino
- A local Trino CLI JAR if you want to query OpenSky

Install dependencies with `uv`:

```powershell
uv sync
```

If you plan to generate OpenSky outputs or 3D CSV outputs, create the output folders once:

```powershell
New-Item -ItemType Directory -Force opensky\output, 3D\output
```

## How To Run

### 1. Build the population layer

```powershell
python population\population_calculation.py
```

Output:

- `geojson\il_blockgroups_population_density.geojson`

### 2. Build the Chicago risk grid

```powershell
python 2D\make_grid.py
```

Output:

- `geojson\risk_grid_v6.geojson`

To view the risk map:

```powershell
python -m http.server 8080
```

Open:

```text
http://localhost:8080/html/map.html
```

### 3. Run the routing experiment

```powershell
python 2D\chicago_graph.py
```

Output:

- `geojson\routes.geojson`

### 4. Export OpenSky traffic data through Trino

Set your username:

```powershell
$env:TRINO_USER = "your_username"
```

Run the exporter:

```powershell
python opensky\export_data.py --date 2019-03-09
```

Default output:

- `opensky\output\ohare_2019-03-09_local_1s_15nm_bbox.csv`

### 5. Generate the 2D traffic visualizations

Static density image:

```powershell
python 2D\make_density_plot.py opensky\output\ohare_2019-03-09_local_1s_15nm_bbox.csv
```

Output:

- `images\visualizations\ohare_density_map.png`

Interactive Leaflet map:

```powershell
python 2D\make_leaflet_map.py opensky\output\ohare_2019-03-09_local_30s_15nm_bbox.csv
python -m http.server 8080
```

Open:

```text
http://localhost:8080/html/ohare_density_leaflet.html
```

### 6. Generate the 3D traffic-density map

```powershell
python 3D\make_3d_density_map.py opensky\output\ohare_2019-03-09_local_30s_15nm_bbox.csv
python -m http.server 8080
```

Open:

```text
http://localhost:8080/html/ohare_3d_density_map.html
```

## Notes

- `2D\chicago_graph.py` expects `geojson\risk_grid_v6.geojson` to exist first.
- `2D\make_grid.py` expects `geojson\il_blockgroups_population_density.geojson` to exist first.
- The OpenSky exporter expects CSV columns consistent with `opensky\query.sql`.
