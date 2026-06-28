"""Build an interactive 3D Plotly plot for route weight configurations."""

import json
import os

import pandas as pd


ALL_ROUTE_RUNS_CSV = "route_robustness/output/all_route_runs.csv"
DIRECT_ROUTE_WEIGHTS_CSV = "route_robustness/output/direct_route_weight_configurations.csv"
OUTPUT_FOLDER = "route_robustness/output"
OUTPUT_CSV = os.path.join(OUTPUT_FOLDER, "route_weight_barycentric_coordinates.csv")
OUTPUT_HTML = os.path.join(OUTPUT_FOLDER, "route_weight_barycentric_plot.html")

WEIGHT_COLUMNS = [
    "distance_weight",
    "population_weight",
    "traffic_weight",
    "airspace_weight",
]

# Four simplex corners. Each plotted point is the weighted average of these
# tetrahedron anchors, using route weights that sum to 1.
BARYCENTRIC_ANCHORS = {
    "distance_weight": (1.0, 1.0, 1.0),
    "population_weight": (1.0, -1.0, -1.0),
    "traffic_weight": (-1.0, 1.0, -1.0),
    "airspace_weight": (-1.0, -1.0, 1.0),
}

LABELS = {
    "distance_weight": "Distance",
    "population_weight": "Population",
    "traffic_weight": "Traffic",
    "airspace_weight": "Airspace",
}


def require_columns(dataframe, required_columns, table_name):
    """Validate that a table has the expected columns."""
    missing_columns = set(required_columns).difference(dataframe.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"{table_name} is missing columns: {missing_text}")


def load_route_weight_table():
    """Load all successful routes and flag direct-route weight configurations."""
    route_runs = pd.read_csv(ALL_ROUTE_RUNS_CSV)
    require_columns(
        route_runs,
        ["route_run_id", "weight_id", "status", *WEIGHT_COLUMNS],
        ALL_ROUTE_RUNS_CSV,
    )

    route_runs = route_runs[route_runs["status"] == "success"].copy()
    if route_runs.empty:
        raise ValueError("No successful routes found in all-route run table.")

    direct_routes = pd.read_csv(DIRECT_ROUTE_WEIGHTS_CSV)
    require_columns(
        direct_routes,
        ["route_run_id", "weight_id"],
        DIRECT_ROUTE_WEIGHTS_CSV,
    )
    direct_route_ids = set(direct_routes["route_run_id"])

    route_runs["is_direct_route"] = route_runs["route_run_id"].isin(direct_route_ids)
    route_runs["route_type"] = route_runs["is_direct_route"].map(
        {True: "Direct route", False: "Non-direct route"}
    )
    return route_runs


def add_barycentric_coordinates(route_weights):
    """Convert route weights to 3D barycentric tetrahedron coordinates."""
    route_weights = route_weights.copy()
    weight_sum = route_weights[WEIGHT_COLUMNS].sum(axis=1)
    if not weight_sum.sub(1.0).abs().le(1e-6).all():
        raise ValueError("Route weights must sum to 1 before plotting.")

    route_weights["barycentric_x"] = 0.0
    route_weights["barycentric_y"] = 0.0
    route_weights["barycentric_z"] = 0.0
    for column_name, (anchor_x, anchor_y, anchor_z) in BARYCENTRIC_ANCHORS.items():
        route_weights["barycentric_x"] += route_weights[column_name] * anchor_x
        route_weights["barycentric_y"] += route_weights[column_name] * anchor_y
        route_weights["barycentric_z"] += route_weights[column_name] * anchor_z

    return route_weights


def build_route_trace(route_weights, label, is_direct, color, marker_size):
    """Build one Plotly 3D scatter trace for a route class."""
    subset = route_weights[route_weights["is_direct_route"] == is_direct]
    customdata = subset[
        [
            "route_run_id",
            "weight_id",
            "distance_weight",
            "population_weight",
            "traffic_weight",
            "airspace_weight",
        ]
    ].values.tolist()

    return {
        "type": "scatter3d",
        "mode": "markers",
        "name": f"{label} ({len(subset)})",
        "x": subset["barycentric_x"].round(6).tolist(),
        "y": subset["barycentric_y"].round(6).tolist(),
        "z": subset["barycentric_z"].round(6).tolist(),
        "customdata": customdata,
        "marker": {
            "size": marker_size,
            "color": color,
            "opacity": 0.86 if is_direct else 0.7,
            "line": {"color": "white", "width": 0.8},
        },
        "hovertemplate": (
            "D: %{customdata[2]:.2f}<br>"
            "A: %{customdata[5]:.2f}<br>"
            "P: %{customdata[3]:.2f}<br>"
            "T: %{customdata[4]:.2f}"
            "<extra></extra>"
        ),
    }


def build_frame_traces():
    """Build Plotly traces for the four-weight tetrahedron frame."""
    anchor_columns = list(BARYCENTRIC_ANCHORS)
    anchors = [BARYCENTRIC_ANCHORS[column] for column in anchor_columns]
    edge_pairs = [
        ("distance_weight", "population_weight"),
        ("distance_weight", "traffic_weight"),
        ("distance_weight", "airspace_weight"),
        ("population_weight", "traffic_weight"),
        ("population_weight", "airspace_weight"),
        ("traffic_weight", "airspace_weight"),
    ]

    edge_x = []
    edge_y = []
    edge_z = []
    for start_column, end_column in edge_pairs:
        start_anchor = BARYCENTRIC_ANCHORS[start_column]
        end_anchor = BARYCENTRIC_ANCHORS[end_column]
        edge_x.extend([start_anchor[0], end_anchor[0], None])
        edge_y.extend([start_anchor[1], end_anchor[1], None])
        edge_z.extend([start_anchor[2], end_anchor[2], None])

    face_trace = {
        "type": "mesh3d",
        "name": "Weight simplex",
        "x": [anchor[0] for anchor in anchors],
        "y": [anchor[1] for anchor in anchors],
        "z": [anchor[2] for anchor in anchors],
        "i": [0, 0, 0, 1],
        "j": [1, 1, 2, 2],
        "k": [2, 3, 3, 3],
        "color": "#d1d5db",
        "opacity": 0.16,
        "hoverinfo": "skip",
        "showlegend": False,
    }
    edge_trace = {
        "type": "scatter3d",
        "mode": "lines",
        "name": "Tetrahedron edges",
        "x": edge_x,
        "y": edge_y,
        "z": edge_z,
        "line": {"color": "#4b5563", "width": 5},
        "hoverinfo": "skip",
        "showlegend": False,
    }
    anchor_trace = {
        "type": "scatter3d",
        "mode": "markers+text",
        "name": "Weight anchors",
        "x": [anchor[0] for anchor in anchors],
        "y": [anchor[1] for anchor in anchors],
        "z": [anchor[2] for anchor in anchors],
        "text": [LABELS[column] for column in anchor_columns],
        "textposition": "top center",
        "marker": {
            "size": 5,
            "color": "#111827",
            "line": {"color": "white", "width": 1},
        },
        "hovertemplate": "<b>%{text}</b><extra></extra>",
        "showlegend": False,
    }

    return [face_trace, edge_trace, anchor_trace]


def make_plot(route_weights):
    """Create and save the interactive route weight barycentric plot."""
    traces = [
        *build_frame_traces(),
        build_route_trace(route_weights, "Non-direct route", False, "#2563eb", 4.8),
        build_route_trace(route_weights, "Direct route", True, "#16a34a", 5.8),
    ]
    layout = {
        "title": {
            "text": "3D Route Weight Barycentric Tetrahedron",
            "x": 0.5,
            "xanchor": "center",
        },
        "scene": {
            "xaxis": {
                "title": "X",
                "range": [-1.18, 1.18],
                "backgroundcolor": "white",
                "gridcolor": "#e5e7eb",
                "zerolinecolor": "#d1d5db",
            },
            "yaxis": {
                "title": "Y",
                "range": [-1.18, 1.18],
                "backgroundcolor": "white",
                "gridcolor": "#e5e7eb",
                "zerolinecolor": "#d1d5db",
            },
            "zaxis": {
                "title": "Z",
                "range": [-1.18, 1.18],
                "backgroundcolor": "white",
                "gridcolor": "#e5e7eb",
                "zerolinecolor": "#d1d5db",
            },
            "aspectmode": "cube",
            "camera": {
                "eye": {"x": 1.65, "y": 1.55, "z": 1.25},
                "center": {"x": 0.0, "y": 0.0, "z": 0.0},
            },
        },
        "legend": {"x": 0.73, "y": 0.97, "bgcolor": "rgba(255,255,255,0.88)"},
        "paper_bgcolor": "white",
        "margin": {"l": 0, "r": 0, "t": 64, "b": 0},
        "hovermode": "closest",
    }
    html = build_html(traces, layout)
    with open(OUTPUT_HTML, "w", encoding="utf-8") as file_handle:
        file_handle.write(html)


def build_html(traces, layout):
    """Build a Plotly HTML document with embedded route data."""
    traces_json = json.dumps(traces)
    layout_json = json.dumps(layout)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>3D Route Weight Barycentric Tetrahedron</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    html, body {{
      height: 100%;
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      background: #ffffff;
    }}
    #plot {{
      width: 100%;
      height: 100vh;
      min-height: 720px;
    }}
    #selection {{
      position: fixed;
      left: 16px;
      top: 16px;
      z-index: 10;
      min-width: 230px;
      padding: 12px 14px;
      border: 1px solid #d1d5db;
      border-radius: 6px;
      background: rgba(255, 255, 255, 0.94);
      box-shadow: 0 6px 20px rgba(15, 23, 42, 0.12);
      color: #111827;
      font-size: 14px;
      line-height: 1.45;
    }}
    #selection strong {{
      display: block;
      margin-bottom: 6px;
      font-size: 15px;
    }}
    #selection .weights {{
      margin-top: 6px;
      font-family: Consolas, Monaco, monospace;
    }}
  </style>
</head>
<body>
  <div id="selection">
    <strong>Selected Weights</strong>
    Click a point to show D, A, P, T.
  </div>
  <div id="plot"></div>
  <script>
    const traces = {traces_json};
    const layout = {layout_json};
    const config = {{
      responsive: true,
      displaylogo: false,
      scrollZoom: true,
      modeBarButtonsToRemove: ["lasso2d", "select2d"]
    }};
    const plot = document.getElementById("plot");
    const selection = document.getElementById("selection");

    Plotly.newPlot(plot, traces, layout, config).then(() => {{
      plot.on("plotly_click", (eventData) => {{
        const point = eventData.points[0];
        const data = point.customdata;
        if (!data) {{
          return;
        }}
        selection.innerHTML = `
          <strong>Selected Weights</strong>
          <div class="weights">
            D=${{Number(data[2]).toFixed(2)}} |
            A=${{Number(data[5]).toFixed(2)}} |
            P=${{Number(data[3]).toFixed(2)}} |
            T=${{Number(data[4]).toFixed(2)}}
          </div>
        `;
      }});
    }});
  </script>
</body>
</html>
"""


def main():
    """Save route barycentric coordinates and interactive 3D plot."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    route_weights = load_route_weight_table()
    route_weights = add_barycentric_coordinates(route_weights)

    output_columns = [
        "route_run_id",
        "weight_id",
        "route_type",
        "is_direct_route",
        *WEIGHT_COLUMNS,
        "barycentric_x",
        "barycentric_y",
        "barycentric_z",
    ]
    route_weights[output_columns].to_csv(OUTPUT_CSV, index=False)
    make_plot(route_weights)

    print(f"Saved barycentric coordinates: {OUTPUT_CSV}")
    print(f"Saved interactive 3D barycentric plot: {OUTPUT_HTML}")
    print("All successful route count:", len(route_weights))
    print("Direct route count:", int(route_weights["is_direct_route"].sum()))
    print("Non-direct route count:", int((~route_weights["is_direct_route"]).sum()))


if __name__ == "__main__":
    main()
