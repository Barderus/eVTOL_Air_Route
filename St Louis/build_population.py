"""Prepare population density data for the St. Louis study area."""

import os

import settings


def main():
    if not os.path.exists(settings.POPULATION_SOURCE):
        raise FileNotFoundError(
            "Add the St. Louis population source and update POPULATION_SOURCE "
            "in St Louis/settings.py."
        )

    raise NotImplementedError(
        "Add the preprocessing steps required by the selected St. Louis "
        "population source, then save the result to POPULATION_GEOJSON."
    )


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, NotImplementedError) as error:
        raise SystemExit(str(error)) from error
