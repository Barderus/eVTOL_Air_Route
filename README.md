# eVTOL Air Route

## What Is This Project

This repository studies low-altitude eVTOL routing with population, airspace,
flight-traffic, and distance costs.

The project is split into two independent city folders:

```text
eVTOL_Air_Route/
|-- Chicago/
|-- St Louis/
|-- docs/
|-- tests/
|-- pyproject.toml
`-- README.md
```

- `Chicago/` contains the original and more experimental research workflow.
- `St Louis/` contains a smaller standalone research workflow.

The city folders do not depend on each other's scripts. This allows Chicago
and St. Louis work to continue separately.

This is a research prototype, not an FAA compliance or operational planning
tool. Airspace costs are modeling assumptions used to compare routes.

### Chicago

Chicago currently includes:

- population preprocessing
- a `500 m` risk grid
- modeled airport and controlled-airspace costs
- OpenSky traffic exports
- A* routes from Clow International Airport
- 2D and 3D maps
- route smoothing experiments

The active Chicago grid is:

```text
Chicago/geojson/risk_grid_v7.geojson
```

The combined A* weights are:

- Distance: `0.6`
- Population: `0.9`
- Airspace: `1.4`
- Air traffic: `1.0`

### St. Louis

St. Louis has separate files for:

- population preprocessing
- risk-grid and map development
- OpenSky traffic export
- A* route generation

## How To Set It Up

The repository uses Python `3.12+` and `uv`.

1. Install Python `3.12` or newer.
2. Install `uv`.
3. Create the environment and install dependencies:

```powershell
uv sync
```

Core dependencies are listed in `pyproject.toml`:

- `geopandas`
- `matplotlib`
- `networkx`
- `pandas`

OpenSky exports also require:

- Java on `PATH`
- a Trino CLI JAR
- OpenSky/Trino access
- `TRINO_USER`

## How To Run It

Run commands from the repository root.

### Chicago

Build the risk grid:

```powershell
uv run python Chicago/2D/make_grid.py
```

Generate the main A* route GeoJSON and HTML pages:

```powershell
uv run python Chicago/2D/generate_astar_toggle_pages.py
```

Export Chicago OpenSky traffic:

```powershell
uv run python Chicago/opensky/export_data.py --date 2026-03-07 --user YOUR_TRINO_USER --trino-path C:\path\to\trino-cli.jar
```

Run the smoothing experiments:

```powershell
uv run python Chicago/Smoothing/approach1_turn_penalty.py
uv run python Chicago/Smoothing/approach2_cubic_spline.py
uv run python Chicago/Smoothing/approach3_turn_penalty_plus_spline.py
```

Tracked Chicago HTML maps are stored in `Chicago/html/`.

### St. Louis

- TBD


## Future Improvements

- Complete the St. Louis population and airspace risk layers.
- Add verified St. Louis endpoints and OpenSky query bounds.
- Add a tracked St. Louis HTML route map.
- Add more focused tests for Chicago routing and smoothing calculations.
- Remove old Chicago experiments when they are no longer useful.
