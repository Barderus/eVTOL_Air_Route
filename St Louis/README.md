# St. Louis eVTOL Routing

## What Is This Project

This folder contains the St. Louis-specific version of the eVTOL routing
project. It is separate from the Chicago scripts so both city workflows can be
developed without editing the same files.

The St. Louis workflow now includes Missouri and Illinois block-group
population preprocessing, configured study bounds, and a simplified radial
airport airspace model. Route endpoints still need to be configured.

## Study Area

The rectangular study bounds use these reference locations:

| Boundary | Place | State | Latitude | Longitude |
| --- | --- | --- | ---: | ---: |
| West | Wentzville | MO | 38.812962 | -90.839767 |
| South | Arnold | MO | 38.432831 | -90.377619 |
| North | Newbern | IL | 39.006940 | -90.336940 |
| East | New Memphis | IL | 38.479170 | -89.678330 |

The map center is the midpoint of the rectangular bounds at latitude
38.719886 and longitude -90.259049.

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

```powershell
cd "St Louis"
uv run python build_population.py
uv run python make_map.py
uv run python render_map.py
```

After traffic and route settings are filled in, run:

```powershell
uv run python export_traffic.py
uv run python astar_routes.py
```

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
