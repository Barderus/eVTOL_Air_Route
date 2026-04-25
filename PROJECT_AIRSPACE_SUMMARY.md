# Project Airspace Summary

This file summarizes the main study-area and airspace constants currently encoded in the project source.

Primary source files:

- `2D/make_grid.py`
- `3D/make_3d_density_map.py`

## 1. Study Bounds

Source: `2D/make_grid.py`

- West: `-88.29`
- South: `41.104190944576466`
- East: `-87.52534901500378`
- North: `42.12361111111111`

Bounding box order used by the project:

```text
WEST, SOUTH, EAST, NORTH
```

## 2. No-Fly Zones

Source: `2D/make_grid.py`

The project defines two no-fly circles:

| Site | Center Latitude | Center Longitude | Radius (m) | Radius (km) | Radius (mi) | Area (km^2) | Area (mi^2) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Argonne National Lab | 41.7150 | -87.9830 | 3000 | 3.000 | 1.864 | 28.274 | 10.917 |
| Fermilab | 41.8407 | -88.2620 | 3000 | 3.000 | 1.864 | 28.274 | 10.917 |

Area formula used for the derived values:

```text
area = pi * r^2
```

## 3. Airport Airspace / Simplified Class B Cylinders

Source: `2D/make_grid.py`

The current model does not encode full FAA shelf geometry. It uses simplified concentric 2D radii around airport centers.

### O'Hare

- Airport center: `41.97807408541273, -87.90902412382081`
- Inner cylinder radius: `10.0 NM` = `18,520 m`
- Outer cylinder radius: `15.0 NM` = `27,780 m`
- Inner disk area: `1077.536 km^2` (`416.039 mi^2`)
- Outer disk area: `2424.456 km^2` (`936.088 mi^2`)
- Outer ring only (`10 NM` to `15 NM`): `1346.920 km^2` (`520.049 mi^2`)

### Midway

- Airport center: `41.7856116663475, -87.75331135429448`
- Inner cylinder radius: `5.0 NM` = `9,260 m`
- Outer cylinder radius: `10.0 NM` = `18,520 m`
- Inner disk area: `269.384 km^2` (`104.010 mi^2`)
- Outer disk area: `1077.536 km^2` (`416.039 mi^2`)
- Outer ring only (`5 NM` to `10 NM`): `808.152 km^2` (`312.029 mi^2`)

### Lewis University Airport

- Airport center: `41.60690586957971, -88.09526487573747`
- Single cylinder radius: `5.1 mi` = `8,207.654 m`
- Disk area: `211.635 km^2` (`81.713 mi^2`)

## 4. Radius of Each Cylinder

These are the current cylinder radii in the simplified model:

| Airport | Radius 1 | Radius 2 |
| --- | --- | --- |
| O'Hare | `10 NM` (`18,520 m`) | `15 NM` (`27,780 m`) |
| Midway | `5 NM` (`9,260 m`) | `10 NM` (`18,520 m`) |
| Lewis University Airport | `5.1 mi` (`8,207.654 m`) | none |

## 5. Drone Flight Height / Altitude Handling

Source: `3D/make_3d_density_map.py`

The project does not currently define one fixed drone flight height for routing.

What is explicitly in the code:

- Altitudes are handled as `AGL` (`above ground level`) in the 3D traffic workflow.
- Ground elevation constant: `680 ft MSL`
- Default altitude bin size: `1000 ft AGL` via `--altitude-step-ft`
- The 3D map altitude ceiling is data-driven and is computed from the maximum observed flight altitude in the loaded CSV.

So, if this file is being used as a project reference, the safest wording is:

```text
No single drone cruise altitude is hard-coded in the current project.
The project works with altitude bands in feet AGL, with a default 3D bin size of 1000 ft.
```

Project assumption to use going forward:

```text
Drone max height: 4000 ft MSL
Equivalent height above ground: 3320 ft AGL
Assumed ground reference: 680 ft MSL
```

Calculation:

```text
4000 ft MSL - 680 ft ground elevation = 3320 ft AGL
```

## 6. Class B Information in This Project

This repo uses a simplified Class B approximation rather than full FAA shelf boundaries:

- O'Hare Class B approximation: `10 NM` inner disk and `15 NM` outer disk
- Midway Class B approximation: `5 NM` inner disk and `10 NM` outer disk
- Risk handling:
  - inside inner radius: treated as high airspace risk
  - between inner and outer radius: treated as a medium band with exponential radial decay

Related source constants from `2D/make_grid.py`:

- `AIRSPACE_HIGH_VAL = 1.0`
- `AIRSPACE_MED_PEAK = 0.6`
- `AIRSPACE_DECAY_S_M = 8000.0`

## 7. Notes

- All Class B values above are project-model values, not a full legal/operational FAA depiction.
- The no-fly zones are modeled as simple circles with a fixed `3000 m` radius.
- The airport airspace areas above are derived from the coded radii, not stored directly in the repo.
