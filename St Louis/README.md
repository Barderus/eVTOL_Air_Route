# St. Louis eVTOL Routing

## What Is This Project

This folder contains the St. Louis-specific version of the eVTOL routing
project. It is separate from the Chicago scripts so both city workflows can be
developed without editing the same files.

The St. Louis workflow now includes Missouri block-group population
preprocessing and configured study bounds. Route endpoints and airspace rules
still need to be configured.

## Study Area

The rectangular study bounds use these reference locations:

| Boundary | Place | State | Latitude | Longitude |
| --- | --- | --- | ---: | ---: |
| West | Wentzville | MO | 38.812962 | -90.839767 |
| South | Arnold | MO | 38.432831 | -90.377619 |
| North | Bethalto | IL | 38.902497 | -90.041333 |
| East | New Memphis | IL | 38.479170 | -89.678330 |

The map center is the midpoint of the rectangular bounds at latitude
38.667664 and longitude -90.259049.

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
- `render_map.py`
  - shared HTML renderer that always includes downtown and airport landmarks
- `export_traffic.py`
  - OpenSky export entry point
- `astar_routes.py`
  - standalone A* route generator

## How To Set It Up

1. Run `build_population.py` to prepare the population-density GeoJSON.
2. Fill in the study area, endpoints, datasets, and traffic settings.
3. Add verified St. Louis airport and controlled-airspace assumptions to
   `make_map.py`.
4. Set `TRINO_USER` and update `TRINO_JAR` in `export_traffic.py`.

## How To Run It

Run each script from the project root:

```powershell
uv run python "St Louis/build_population.py"
uv run python "St Louis/make_map.py"
uv run python "St Louis/render_map.py"
uv run python "St Louis/export_traffic.py"
uv run python "St Louis/astar_routes.py"
```

Scripts that still require St. Louis configuration stop with a clear message.
The HTML renderer adds these locations to every St. Louis map:

- Downtown St. Louis
- St. Louis Downtown Airport
- St. Louis Lambert International Airport
- MidAmerica St. Louis Airport

Generated traffic CSV, route GeoJSON, risk-grid GeoJSON, and processed
population files are ignored. Source population data and HTML maps can be
tracked.

## Future Improvements

- Implement the St. Louis airspace risk model.
- Add route endpoints and traffic datasets.
- Add a simple tracked HTML route map.
