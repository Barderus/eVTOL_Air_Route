# Chicago and St. Louis Comparison

This note compares the Chicago and St. Louis routing models using the six
OpenSky dates that both workflows share.

These values describe the modeled routing areas in this repository. They do
not describe the legal city limits or all real-world flights in either metro
area.

## Study Area and Air Traffic

Chicago has the larger modeled routing-grid area. The Chicago grid covers
about `1.61` times more area than the St. Louis grid. Chicago also has more
observed aircraft in the OpenSky samples.

| Metric | Chicago | St. Louis |
| --- | ---: | ---: |
| Modeled routing-grid area | 5,518.23 km^2 | 3,425.01 km^2 |
| Average unique aircraft per date | 1,470.00 | 326.67 |
| Sum of daily unique-aircraft counts | 8,820 | 1,960 |
| Average air traffic density | 26.64 aircraft per 100 km^2 | 9.54 aircraft per 100 km^2 |

Across the six shared dates, Chicago has about `4.50` times more unique
aircraft in these daily exports. After adjusting for modeled area, Chicago
still has higher traffic density.

This suggests that the Chicago test area is both larger and busier in the
current model.

## Flights by Date

The Chicago aircraft counts are higher on every shared date.

| Date | Chicago unique aircraft | St. Louis unique aircraft | Chicago / St. Louis |
| --- | ---: | ---: | ---: |
| 2025-07-12 | 1,544 | 347 | 4.45 |
| 2025-07-14 | 1,702 | 342 | 4.98 |
| 2026-01-10 | 1,196 | 267 | 4.48 |
| 2026-01-12 | 1,379 | 347 | 3.97 |
| 2026-03-07 | 1,385 | 292 | 4.74 |
| 2026-03-09 | 1,614 | 365 | 4.42 |
| **Average** | **1,470.00** | **326.67** | **4.50** |

These counts use distinct `icao24` aircraft identifiers from each 24-hour
OpenSky export. They are a reproducible proxy for traffic volume, but they are
not a perfect flight count. One aircraft can make more than one flight in a
day, and receiver coverage can differ between cities.

## Route Costs

Both workflows use the same combined weights.

For Chicago, the lowest average modeled route cost is Clow International
Airport to Midway International Airport.

| Route from Clow International Airport | Average cost | Minimum | Maximum |
| --- | ---: | ---: | ---: |
| To Chicago Union Station | 131.541 | 131.489 | 131.647 |
| To O'Hare International Airport | 120.894 | 120.840 | 120.969 |
| To Midway International Airport | 88.038 | 87.958 | 88.331 |

For St. Louis, the lowest average modeled route cost is MidAmerica to St.
Louis Union Station.

| Route | Average cost | Minimum | Maximum | Average distance |
| --- | ---: | ---: | ---: | ---: |
| MidAmerica to St. Louis Union Station | 93.448 | 84.416 | 99.271 | 51.60 km |
| MidAmerica to St. Louis Lambert | 158.239 | 147.074 | 167.005 | 76.53 km |
| St. Louis Downtown Airport to St. Louis Lambert | 107.195 | 104.007 | 111.294 | 35.13 km |

The MidAmerica to Lambert route has the highest average cost in this
comparison. The longer route and modeled Lambert airspace both appear to
contribute to that higher score.

## Interpretation

The main takeaway is that Chicago is larger, has more observed aircraft, and
has higher area-normalized air traffic density in this model.

The route cost comparison should be read carefully. These costs are
dimensionless weighted A* scores. They are not dollars, travel time, or
physical distance.

Cost values are most useful within each workflow. Cross-city comparisons are
less direct because the grids, endpoints, normalized risk distributions, and
airspace models are different.
