# Chicago eVTOL Routing

## What Is This Project

This folder contains the original Chicago eVTOL routing research. It includes
the larger experimental workflow for population risk, airspace modeling,
OpenSky traffic, A* routes, smoothing methods, and 2D and 3D maps.

## Folder Layout

- `population/`
  - Chicago population preprocessing and source data
- `2D/`
  - risk-grid creation, A* routing, cost analysis, and 2D maps
- `3D/`
  - 3D traffic-density maps
- `Smoothing/`
  - route-smoothing experiments
- `opensky/`
  - Chicago OpenSky SQL, export script, and ignored CSV outputs
- `geojson/`
  - ignored generated grids and routes
- `html/`
  - tracked Chicago HTML maps
- `images/`
  - ignored generated charts and local image assets

## Future Improvements

- Continue improving Chicago airspace and traffic assumptions.
- Reduce older experimental scripts when they are no longer useful.
- Add focused tests for shared Chicago route calculations.
