# eVTOL Air Route

This repository models low-altitude eVTOL routing in the Chicago region.
It combines a static risk grid, date-specific flight-traffic data, A* route
generation, and optional route-smoothing experiments.

The project is a research/prototyping codebase, not an FAA compliance tool.
Its airspace handling is an explicit modeling assumption used to compare route
behavior under different cost layers.

## What This Project Produces

- A static Chicago-area risk grid in GeoJSON form
- Route comparisons from Clow International Airport to:
  - Union Station
  - O'Hare
  - Midway
- Prebuilt Leaflet HTML pages for route review
- OpenSky traffic exports used to build traffic-density layers
- Smoothing experiment outputs for post-processing route geometry

## Current Source Of Truth

- Active routing grid: `geojson/risk_grid_v7.geojson`
- Grid cell size: `500 m`
- Combined A* route weights:
  - distance: `0.6`
  - population: `0.9`
  - airspace: `1.4`
  - air traffic: `1.0`
- Static map combined cell cost:
  - `risk_cost = 0.9 * city_risk + 1.4 * airport_risk_combined`

## Airspace Assumptions

The current airspace model is a project assumption for low-altitude routing.
It is documented in more detail in `PROJECT_AIRSPACE_SUMMARY.md` and
`docs/airspace.md`.

Current assumptions in code:

- modeled low-altitude envelope: `3000 ft AGL / 4000 ft MSL`
- O'Hare:
  - coded shelf geometry with a `5 NM` core and `10 NM` outer shelf
- Midway:
  - `5 NM` core and `10 NM` outer shelf
- Lewis University Airport:
  - Class D polygon with a Clow-facing cutout
- shelf-edge radial decay:
  - up to `1500 m` beyond the modeled footprint
- runway-aligned corridor risk:
  - added on top of shelf/polygon airspace risk

## Repository Layout

- `2D/`
  - static grid creation, A* routing, route-cost summaries, and HTML map generation
- `3D/`
  - 3D traffic-density visualization outputs
- `Smoothing/`
  - route-smoothing experiments and shared helpers
- `opensky/`
  - Trino/OpenSky export workflow and raw CSV outputs
- `geojson/`
  - main GIS artifacts used by the scripts
- `html/`
  - generated Leaflet pages
- `docs/`
  - focused notes on airspace, routing, traffic, and smoothing

## Environment Setup

This repo is configured with `uv` and Python `3.12+`.

1. Install Python `3.12` or newer.
2. Install `uv`.
3. Create the environment and install dependencies:

```powershell
uv sync
```

Core Python dependencies are declared in `pyproject.toml`:

- `geopandas`
- `matplotlib`
- `networkx`
- `pandas`

You will also need a working local geospatial stack compatible with
`geopandas` on your machine.

## Typical Workflow

Most handoff work should follow this order:

1. Build or verify the population layer.
2. Build the static risk grid.
3. Export OpenSky traffic CSVs if new dates are needed.
4. Generate the multi-date A* route pages.
5. Run smoothing experiments if you need the alternate post-processed outputs.

## How To Run

### 1. Build the static risk grid

This script rebuilds the main GeoJSON used by routing:

```powershell
uv run python 2D/make_grid.py
```

Primary output:

- `geojson/risk_grid_v7.geojson`

### 2. Generate the multi-date A* route outputs

This script reads the active grid plus the configured OpenSky CSV datasets and
generates destination-specific GeoJSON and HTML pages.

```powershell
uv run python 2D/generate_astar_toggle_pages.py
```

Primary outputs:

- `geojson/clow_to_union_station_astar_routes.geojson`
- `geojson/clow_to_ohare_astar_routes.geojson`
- `geojson/clow_to_midway_astar_routes.geojson`
- matching HTML files in `html/`

### 3. Run a single route build directly

This is the simplest one-off route script for the Clow to Union Station case.

```powershell
uv run python 2D/chicago_graph.py
```

Primary output:

- `geojson/routes.geojson`

### 4. Export OpenSky traffic data

This requires:

- Java on `PATH`
- a local Trino CLI JAR
- OpenSky/Trino access
- `TRINO_USER` set or passed on the command line

Example:

```powershell
uv run python opensky/export_data.py --date 2026-03-07 --user YOUR_TRINO_USER --trino-path C:\path\to\trino-cli.jar
```

Default output location:

- `opensky/output/`

### 5. Run smoothing experiments

```powershell
uv run python Smoothing/approach1_turn_penalty.py
uv run python Smoothing/approach2_cubic_spline.py
uv run python Smoothing/approach3_turn_penalty_plus_spline.py
```

Outputs:

- HTML files in `html/`
- CSV summaries in `Smoothing/`

## Key Files For A New Maintainer

Read these first:

- `2D/make_grid.py`
  - builds the static grid and airspace-related risk fields
- `2D/generate_astar_toggle_pages.py`
  - main route-generation pipeline used for handoff artifacts
- `2D/chicago_graph.py`
  - smallest end-to-end routing script
- `opensky/export_data.py`
  - exports raw traffic data from OpenSky via Trino
- `Smoothing/route_smoothing_common.py`
  - shared machinery behind the smoothing experiments

Then use the focused docs:

- `PROJECT_AIRSPACE_SUMMARY.md`
- `docs/airspace.md`
- `docs/astar-routing.md`
- `docs/flight-traffic.md`
- `docs/smoothing.md`

## Important Project Notes

- The static map and the route solver are not the same product view.
- `html/map.html` reflects population, airspace, and combined cell cost only.
- Traffic is route-specific and date-specific, so it is handled during route
  generation rather than in the static map.
- Distance is an edge cost, not a cell cost.
- Several generated artifacts already exist in the repo. Rebuild only the ones
  you intend to replace.
- The worktree may contain experimental outputs or local-only files. Check
  `git status` before assuming an artifact is canonical.

