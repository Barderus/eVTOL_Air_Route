# eVTOL Air Route

`eVTOL Air Route` is a research-oriented Python project for modeling low-altitude airspace around Chicago for eVTOL and urban air mobility routing. It combines three inputs:

- population density
- airport-related airspace risk
- observed aircraft traffic from OpenSky

The repository is organized as a script-driven workflow that:

- builds a population-density layer
- generates a risk-weighted Chicago grid
- runs graph-based routing experiments
- exports historical traffic data near O'Hare through Trino
- renders 2D and 3D traffic-density visualizations

## Requirements

- Python `>=3.12,<3.15`
- PowerShell examples assume Windows
- Java installed if you want to export OpenSky data through Trino
- A local Trino CLI JAR if you want to query OpenSky

Install dependencies with either `uv` or `pip`:

```powershell
uv sync
```

or

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

If you plan to generate OpenSky outputs, create the output folders once:

```powershell
New-Item -ItemType Directory -Force opensky\output, traffic_3d\output
```

## How To Run

### 1. Build the population layer

```powershell
python population\population_calculation.py
```

Output:

- `population\il_blockgroups_population_density.geojson`

### 2. Build the Chicago risk grid

```powershell
python map\make_grid.py
```

Output:

- `map\risk_grid_v5.geojson`

To view the risk map:

```powershell
python -m http.server 8080
```

Open:

```text
http://localhost:8080/map/map.html
```

### 3. Run the routing experiment

```powershell
python graph\chicago_graph.py
```

Output:

- `graph\routes.geojson`

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

The script accepts either a direct JAR path or a folder containing the Trino CLI JAR. If your local setup is elsewhere, pass it explicitly:

```powershell
python opensky\export_data.py --date 2019-03-09 --trino-path C:\tools\trino-cli.jar
```

### 5. Generate the 2D traffic visualizations

Static density image:

```powershell
python opensky\make_density_plot.py opensky\output\ohare_2019-03-09_local_1s_15nm_bbox.csv
```

Interactive Leaflet map:

```powershell
python opensky\make_leaflet_map.py opensky\output\ohare_2019-03-09_local_1s_15nm_bbox.csv
python -m http.server 8080
```

Open:

```text
http://localhost:8080/opensky/output/ohare_density_leaflet.html
```

### 6. Generate the 3D traffic-density map

```powershell
python traffic_3d\make_3d_density_map.py opensky\output\ohare_2019-03-09_local_1s_15nm_bbox.csv
python -m http.server 8080
```

Open:

```text
http://localhost:8080/traffic_3d/output/ohare_3d_density_map.html
```

## Notes

- `graph\chicago_graph.py` expects `map\risk_grid_v5.geojson` to exist first.
- The OpenSky visualization scripts expect CSV columns consistent with [`opensky/query.sql`](opensky/query.sql).
- Legacy wrapper folders such as `OpenSky-Trino\` and `3D_Map\` still exist for older command compatibility.
