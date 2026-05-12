# eVTOL Air Route

This project studies low-altitude eVTOL routing in the Chicago region. It
builds a risk-aware routing surface, compares alternate routes to major
destinations, and visualizes how population exposure, airport airspace, and
observed flight activity influence path selection.

The project combines GIS processing, graph search, and route post-processing.
Its main methods include grid-based risk modeling, A* search, traffic-density
aggregation from OpenSky observations, and optional smoothing experiments for
reducing sharp turns in routed paths.

The data used in the project comes from census-derived population density,
airport and controlled-airspace assumptions encoded as geometric risk surfaces,
and OpenSky state-vector traffic exports filtered to the Chicago study area.
The result is a set of map-based artifacts for exploring how different cost
factors change preferred low-altitude flight paths.
