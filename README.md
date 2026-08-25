# eVTOL Air Route

## What Is This Project

This repository studies low-altitude eVTOL route planning with four cost
factors:

- distance
- population exposure
- observed flight traffic
- airport and controlled-airspace risk

The current work compares route formation across different weight settings and
route clustering methods. The main research question is how much the selected
route changes when the cost weights or route-distance method changes.

This is a research prototype. The routes and airspace costs are modeling
outputs, not FAA compliance guidance or operational flight plans.

## Project Layout

```text
eVTOL_Air_Route/
|-- Chicago/
|   |-- 2D/
|   |-- 3D/
|   |-- Smoothing/
|   |-- route_robustness/
|-- St Louis/
|   |-- route_robustness/
|   |-- routes/
|-- docs/
|-- pyproject.toml
`-- README.md
```

`Chicago/` contains the original and broader workflow. It includes the risk
grid, A* routes, route smoothing experiments, and the route robustness work.

`St Louis/` contains a separate city workflow using the same general routing
ideas on St. Louis route pairs.

## Current Route Robustness Work

The route robustness workflow generates normalized weight configurations for
the four route-cost factors. The baseline weight grid uses `0.10` increments
plus an explicit equal-weight case, producing `287` weight configurations.

The Chicago robustness work includes:

- Clow International Airport to Chicago Union Station
- Clow International Airport to O'Hare International Airport
- Clow International Airport to Midway International Airport

The St. Louis robustness work includes:

- MidAmerica St. Louis Airport to St. Louis Union Station
- MidAmerica St. Louis Airport to St. Louis Lambert International Airport
- St. Louis Downtown Airport to St. Louis Lambert International Airport

Routes are clustered with several methods so the analysis is not tied to one
definition of route similarity:

- shared-node Jaccard similarity
- Levenshtein edit distance on route-node sequences
- discrete Frechet distance on route geometry
- DBSCAN using Frechet distance
- OPTICS using Frechet distance for the direct-route subset
- average-linkage hierarchical clustering using Jaccard, edit distance, or
  Frechet distance

The route-weight visualizations use a 3D barycentric tetrahedron because the
route model has four weights. Each point represents one tested weight
configuration.

## Main Outputs

Chicago route robustness outputs are under:

```text
Chicago/route_robustness/
|-- maps/
|-- output/
|-- additional_route_clusters/
`-- additional_tetrahedrons_barycentric/
```

St. Louis route robustness outputs are under:

```text
St Louis/route_robustness/
|-- output/
|-- route_clusters/
`-- tetrahedrons_barycentric/
```

Route cost summaries for St. Louis are documented in:

```text
St Louis/routes/README.md
```

## How To Set It Up

Use Python `3.12+` and install dependencies from the repository root:

```powershell
uv sync
```

## How To Run It

Most scripts are run from the repository root:

```powershell
uv run python path/to/script.py
```

City-specific workflow notes are in:

```text
Chicago/route_robustness/README.md
St Louis/routes/README.md
St Louis/route_robustness/README.md
```
