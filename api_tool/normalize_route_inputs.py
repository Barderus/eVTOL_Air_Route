"""Normalize route input fields for route robustness workflows."""


def min_max_normalize(value, minimum, maximum):
    """Normalize one value to the 0-1 range."""
    if maximum == minimum:
        return 0.0
    return (value - minimum) / (maximum - minimum)


def normalize_route_inputs(rows=None, fields=None):
    """
    Normalize selected numeric fields in a list of route dictionaries.

    Returns copied rows with *_norm columns added.
    """
    if rows is None:
        return {
            "rows": [],
            "fields": fields or [],
            "field_ranges": {},
            "message": "Pass rows and fields to normalize route inputs.",
        }

    fields = fields or ["distance", "population", "traffic", "airspace"]
    field_ranges = {}
    for field_name in fields:
        values = [float(row[field_name]) for row in rows if field_name in row]
        if values:
            field_ranges[field_name] = {"min": min(values), "max": max(values)}

    normalized_rows = []
    for row in rows:
        normalized_row = row.copy()
        for field_name, value_range in field_ranges.items():
            if field_name in row:
                normalized_row[f"{field_name}_norm"] = min_max_normalize(
                    float(row[field_name]),
                    value_range["min"],
                    value_range["max"],
                )
        normalized_rows.append(normalized_row)

    return {
        "rows": normalized_rows,
        "fields": fields,
        "field_ranges": field_ranges,
    }


if __name__ == "__main__":
    result = normalize_route_inputs()
    print(result)

