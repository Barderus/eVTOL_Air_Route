# Pilot Results Summary

This note summarizes the first route robustness pilot for the Clow
International Airport to Chicago Union Station route.

The single-factor routes are already familiar from earlier work, so this note
focuses mostly on the mixed-weight and heavy-weight route cases.

## Normalization Check

All route factors were normalized to the expected `0` to `1` range.

| Factor | Normalized minimum | Normalized maximum | Method |
| --- | ---: | ---: | --- |
| Distance | 0.415120 | 1.000000 | max edge distance |
| Population | 0.000000 | 1.000000 | max value |
| Traffic | 0.000000 | 1.000000 | log1p, 95th percentile clipping, power 0.75 |
| Airspace | ~0.000000 | 1.000000 | 99th percentile clipping, power 0.75 |

Distance does not start at `0` because every graph edge has a positive travel
distance.

## Novel Route Cases

The most useful pilot routes for the next robustness step are the mixed and
heavy-weight cases.

| Pilot route | Weights D/P/T/A | Distance | Path nodes | Weighted score |
| --- | --- | ---: | ---: | ---: |
| Equal weight | 0.25 / 0.25 / 0.25 / 0.25 | 67.34 km | 110 | 58.5004 |
| Distance heavy | 0.70 / 0.10 / 0.10 / 0.10 | 65.68 km | 110 | 79.9895 |
| Population heavy | 0.10 / 0.70 / 0.10 / 0.10 | 69.75 km | 114 | 45.1760 |
| Traffic heavy | 0.10 / 0.10 / 0.70 / 0.10 | 177.29 km | 319 | 43.7613 |
| Airspace heavy | 0.10 / 0.10 / 0.10 / 0.70 | 125.12 km | 213 | 38.2640 |
| Mixed balanced | 0.20 / 0.20 / 0.30 / 0.30 | 68.17 km | 110 | 55.8830 |

The traffic-heavy and airspace-heavy routes change more than the other mixed
or heavy-weight cases.

## Initial Interpretation

The equal-weight, distance-heavy, population-heavy, and mixed-balanced routes
stay relatively close in distance.

This suggests that moderate mixed weighting may keep the route in a similar
general corridor for this test case.

The traffic-heavy and airspace-heavy routes move farther away from that
corridor. This suggests that traffic and airspace become more influential when
they receive dominant weights.

This is the first clear sign that route robustness depends strongly on whether
traffic or airspace is treated as the main cost factor.

## Reference Single-Factor Routes

The single-factor routes are useful as endpoints for understanding model
behavior, but they are not the main focus of this pilot summary.

| Pilot route | Distance | Path nodes | Weighted score |
| --- | ---: | ---: | ---: |
| Distance only | 65.68 km | 110 | 92.8909 |
| Population only | 74.90 km | 116 | 33.4750 |
| Traffic only | 195.94 km | 329 | 11.5214 |
| Airspace only | 239.35 km | 422 | 0.2066 |

This is expected behavior for pure single-factor optimization. The more useful
question is how quickly routes start moving toward those extreme corridors as
traffic or airspace weights increase.
