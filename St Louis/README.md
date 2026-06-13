# St. Louis eVTOL Routing

## What Is This Project

This folder contains the St. Louis-specific version of the eVTOL routing
project. It is separate from the Chicago scripts so both city workflows can be
developed without editing the same files.

The St. Louis workflow now includes Missouri and Illinois block-group
population preprocessing, configured study bounds, and a simplified radial
airport airspace model. It also includes OpenSky traffic maps, three configured
route pairs, route-cost analysis, single-factor comparisons, and a 24-hour
traffic-density route overlay.

## Study Area

The rectangular study bounds use these reference locations:

| Boundary | Place | State | Latitude | Longitude |
| --- | --- | --- | ---: | ---: |
| West | Lindenwood University | MO | 38.787800 | -90.501400 |
| South | Arnold | MO | 38.432831 | -90.377619 |
| North | Wood River | IL | 38.863100 | -90.088600 |
| East | New Memphis | IL | 38.479170 | -89.678330 |

## Airspace Model

The current grid uses simplified radial airspace assumptions:

| Code | Airport | Type | Core | Second shelf | Vertical range |
| --- | --- | --- | ---: | ---: | --- |
| STL | St. Louis Lambert International Airport | Class B | 6 NM | 15 NM | Surface to 4000 ft MSL |
| CPS | St. Louis Downtown-Parks Airport | Class D | 3.475905 NM | None | Surface to 2900 ft MSL |
| BLV | Scott AFB / MidAmerica St. Louis Airport | Class D | 4.257984 NM | None | Surface to 3000 ft MSL |

The STL first layer is a high-risk circle from the airport center through
`6 NM`. Its second layer is a medium-risk ring greater than `6 NM` and through
`15 NM`, excluding the core. The CPS and BLV cores are high risk. All
configured airspace distances use nautical miles.

## How To Run It

Set the working directory to the `St Louis` folder. Run each script
from that folder so the editable paths in `settings.py` resolve consistently:

The configured route pairs are:

- MidAmerica St. Louis Airport to St. Louis Union Station
- MidAmerica St. Louis Airport to St. Louis Lambert International Airport
- St. Louis Downtown Airport to St. Louis Lambert International Airport

Generate each route for every configured traffic date.

Each page compares the configured combined route with distance-only,
population-only, airspace-only, and air-traffic-only routes for every traffic
date. The generated pages are:

```text
maps/midamerica_to_union_station_astar.html
maps/midamerica_to_lambert_astar.html
maps/downtown_airport_to_lambert_astar.html
```

In `maps/st_louis_routes_24h_traffic.html`. The page includes traffic
date, route-pair, and route-weighting controls, with no shorter-hour controls.
It embeds the 24-hour density layer only, keeping it substantially smaller than
the full traffic page with sampled flight tracks.

Route cost totals and distance, population, airspace, and traffic percentages
are saved to `routes/output/st_louis_route_cost_analysis.csv`. Pie charts for
each traffic date are saved in `routes/analysis/`.

Scripts that still require configuration stop with a clear message.
The HTML renderer adds these locations to every St. Louis map:

- Downtown St. Louis
- St. Louis Downtown Airport
- St. Louis Lambert International Airport
- MidAmerica St. Louis Airport

Generated traffic CSV, route GeoJSON, risk-grid GeoJSON, and processed
population files are ignored. Source population data and HTML maps can be
tracked.

Population risk thresholds are calculated from block groups that intersect the
configured study area. Low, medium, and high classes use the 70th and 90th
percentile density thresholds for that local subset.
