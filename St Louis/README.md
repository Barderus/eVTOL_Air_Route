# St. Louis eVTOL Routing

## What Is This Project

This folder contains the St. Louis-specific version of the eVTOL routing
project. It is separate from the Chicago scripts so both city workflows can be
developed without editing the same files.

The St. Louis workflow is currently a scaffold. It does not contain guessed
study bounds, endpoints, population data, or airspace rules.

## Folder Layout

- `settings.py`
  - editable study area, paths, endpoints, traffic settings, and route weights
- `population/`
  - source population data and generated population-density GeoJSON
- `maps/`
  - generated risk-grid GeoJSON and tracked HTML maps
- `traffic/`
  - OpenSky SQL and ignored CSV exports
- `routes/`
  - ignored route GeoJSON outputs
- `build_population.py`
  - population preprocessing entry point
- `make_map.py`
  - risk-grid and map-data entry point
- `export_traffic.py`
  - OpenSky export entry point
- `astar_routes.py`
  - standalone A* route generator

## How To Set It Up

1. Add the St. Louis population source under `St Louis/population/`.
2. Update `POPULATION_SOURCE` in `settings.py`.
3. Fill in the study area, endpoints, datasets, and traffic settings.
4. Add verified St. Louis airport and controlled-airspace assumptions to
   `make_map.py`.
5. Set `TRINO_USER` and update `TRINO_JAR` in `export_traffic.py`.

## How To Run It

Run each script from the project root:

```powershell
uv run python "St Louis/build_population.py"
uv run python "St Louis/make_map.py"
uv run python "St Louis/export_traffic.py"
uv run python "St Louis/astar_routes.py"
```

Each script stops with a clear message until its required St. Louis inputs are
configured.

Generated traffic CSV, route GeoJSON, risk-grid GeoJSON, and processed
population files are ignored. Source population data and HTML maps can be
tracked.

## Future Improvements

- Complete the population preprocessing for the selected source.
- Implement the St. Louis airspace risk model.
- Add route endpoints and traffic datasets.
- Add a simple tracked HTML route map.
