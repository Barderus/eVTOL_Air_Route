"""Build Chicago barycentric route-weight visualizations."""

import colorsys
import json
from pathlib import Path

import pandas as pd


CHICAGO_ROBUSTNESS = Path("Chicago") / "route_robustness"
OUTPUT_FOLDER = CHICAGO_ROBUSTNESS / "output"
MAPS_FOLDER = CHICAGO_ROBUSTNESS / "additional_tetrahedrons_barycentric"
ROUTE_RUNS_CSV = OUTPUT_FOLDER / "additional_route_runs.csv"
BARYCENTRIC_CSV = OUTPUT_FOLDER / "additional_route_weight_barycentric_coordinates.csv"

CLUSTER_METHODS = [
    {
        "label": "DBSCAN On Frechet",
        "file": OUTPUT_FOLDER / "additional_route_clusters_dbscan.csv",
        "cluster_column": "dbscan_cluster",
        "suffix": "dbscan",
    },
    {
        "label": "Hierarchical On Edit Distance",
        "file": OUTPUT_FOLDER / "additional_route_clusters_edit_distance.csv",
        "cluster_column": "edit_distance_cluster",
        "suffix": "edit_distance",
    },
    {
        "label": "Hierarchical On Frechet",
        "file": OUTPUT_FOLDER / "additional_route_clusters_frechet.csv",
        "cluster_column": "frechet_cluster",
        "suffix": "frechet",
    },
    {
        "label": "Hierarchical On Jaccard",
        "file": OUTPUT_FOLDER / "additional_route_clusters_hierarchical_jaccard.csv",
        "cluster_column": "hierarchical_jaccard_cluster",
        "suffix": "hierarchical_jaccard",
    },
]

WEIGHT_COLUMNS = [
    "distance_weight",
    "population_weight",
    "traffic_weight",
    "airspace_weight",
]

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

VARIANT_COLORS = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
    "#1f78b4",
    "#b15928",
    "#6a3d9a",
    "#fb9a99",
    "#33a02c",
]


def color_for_index(index, total_count):
    """Return a deterministic distinct color for large cluster counts."""
    if total_count <= len(VARIANT_COLORS):
        return VARIANT_COLORS[index]

    hue = (index * 137.508) % 360.0
    lightness = 43.0 + (index % 4) * 7.0
    saturation = 0.68
    red, green, blue = colorsys.hls_to_rgb(hue / 360.0, lightness / 100.0, saturation)
    return f"#{round(red * 255):02x}{round(green * 255):02x}{round(blue * 255):02x}"


def require_columns(dataframe, required_columns, table_name):
    """Validate that a table has the expected columns."""
    missing_columns = set(required_columns).difference(dataframe.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"{table_name} is missing columns: {missing_text}")


def add_barycentric_coordinates(route_runs):
    """Convert four route weights to 3D tetrahedron coordinates."""
    route_runs = route_runs.copy()
    weight_sum = route_runs[WEIGHT_COLUMNS].sum(axis=1)
    if not weight_sum.sub(1.0).abs().le(1e-6).all():
        raise ValueError("Route weights must sum to 1 before plotting.")

    route_runs["barycentric_x"] = 0.0
    route_runs["barycentric_y"] = 0.0
    route_runs["barycentric_z"] = 0.0

    for column_name, (anchor_x, anchor_y, anchor_z) in BARYCENTRIC_ANCHORS.items():
        route_runs["barycentric_x"] += route_runs[column_name] * anchor_x
        route_runs["barycentric_y"] += route_runs[column_name] * anchor_y
        route_runs["barycentric_z"] += route_runs[column_name] * anchor_z

    return route_runs


def add_route_variant_ids(route_runs):
    """Assign exact-path variant labels within each route pair."""
    route_runs = route_runs.copy()
    route_runs["route_variant"] = ""

    for route_pair in sorted(route_runs["route_pair"].unique()):
        route_pair_mask = route_runs["route_pair"] == route_pair
        route_pair_runs = route_runs[route_pair_mask]
        signature_counts = route_pair_runs["path_node_ids"].value_counts()
        signature_to_variant = {
            signature: f"variant_{index:03d}"
            for index, signature in enumerate(signature_counts.index, start=1)
        }
        route_runs.loc[route_pair_mask, "route_variant"] = route_pair_runs[
            "path_node_ids"
        ].map(signature_to_variant)

    return route_runs


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

    return [
        {
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
        },
        {
            "type": "scatter3d",
            "mode": "lines",
            "name": "Tetrahedron edges",
            "x": edge_x,
            "y": edge_y,
            "z": edge_z,
            "line": {"color": "#4b5563", "width": 5},
            "hoverinfo": "skip",
            "showlegend": False,
        },
        {
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
        },
    ]


def build_variant_trace(route_runs, route_variant, color):
    """Build one Plotly 3D scatter trace for an exact route variant."""
    subset = route_runs[route_runs["route_variant"] == route_variant]
    customdata = subset[
        [
            "weight_id",
            "distance_weight",
            "population_weight",
            "traffic_weight",
            "airspace_weight",
            "route_distance_km",
            "total_weighted_score",
            "path_nodes",
        ]
    ].values.tolist()

    return {
        "type": "scatter3d",
        "mode": "markers",
        "name": f"{route_variant} ({len(subset)})",
        "x": subset["barycentric_x"].round(6).tolist(),
        "y": subset["barycentric_y"].round(6).tolist(),
        "z": subset["barycentric_z"].round(6).tolist(),
        "customdata": customdata,
        "marker": {
            "size": 5.6,
            "color": color,
            "opacity": 0.84,
            "line": {"color": "white", "width": 0.8},
        },
        "hovertemplate": (
            "Weight: %{customdata[0]}<br>"
            "D: %{customdata[1]:.2f}<br>"
            "P: %{customdata[2]:.2f}<br>"
            "T: %{customdata[3]:.2f}<br>"
            "A: %{customdata[4]:.2f}<br>"
            "Distance: %{customdata[5]:.2f} km<br>"
            "Score: %{customdata[6]:.3f}<br>"
            "Nodes: %{customdata[7]}"
            "<extra></extra>"
        ),
    }


def build_cluster_trace(route_runs, cluster_column, cluster_id, color):
    """Build one Plotly 3D scatter trace for a route-distance cluster."""
    subset = route_runs[route_runs[cluster_column] == cluster_id]
    customdata = subset[
        [
            "weight_id",
            "distance_weight",
            "population_weight",
            "traffic_weight",
            "airspace_weight",
            "route_distance_km",
            "total_weighted_score",
            "path_nodes",
            cluster_column,
        ]
    ].values.tolist()

    return {
        "type": "scatter3d",
        "mode": "markers",
        "name": f"{cluster_id} ({len(subset)})",
        "x": subset["barycentric_x"].round(6).tolist(),
        "y": subset["barycentric_y"].round(6).tolist(),
        "z": subset["barycentric_z"].round(6).tolist(),
        "customdata": customdata,
        "marker": {
            "size": 5.8,
            "color": color,
            "opacity": 0.86,
            "line": {"color": "white", "width": 0.8},
        },
        "hovertemplate": (
            "Weight: %{customdata[0]}<br>"
            "D: %{customdata[1]:.2f}<br>"
            "P: %{customdata[2]:.2f}<br>"
            "T: %{customdata[3]:.2f}<br>"
            "A: %{customdata[4]:.2f}<br>"
            "Distance: %{customdata[5]:.2f} km<br>"
            "Score: %{customdata[6]:.3f}<br>"
            "Nodes: %{customdata[7]}<br>"
            "Cluster: %{customdata[8]}"
            "<extra></extra>"
        ),
    }


def build_html(traces, layout, cluster_summary=None):
    """Build a self-contained Plotly HTML document."""
    traces_json = json.dumps(traces)
    layout_json = json.dumps(layout)
    cluster_summary_json = json.dumps(cluster_summary or [])
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Chicago Route Weight Barycentric Tetrahedron</title>
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
      min-width: 250px;
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
    #clusterIndex {{
      position: fixed;
      right: 16px;
      top: 16px;
      z-index: 10;
      width: 248px;
      max-height: calc(100vh - 32px);
      overflow: auto;
      padding: 12px 12px;
      border: 1px solid #d1d5db;
      border-radius: 6px;
      background: rgba(255, 255, 255, 0.94);
      box-shadow: 0 6px 20px rgba(15, 23, 42, 0.12);
      color: #111827;
      font-size: 13px;
      line-height: 1.35;
    }}
    #clusterIndex:empty {{
      display: none;
    }}
    #clusterIndex strong {{
      display: block;
      margin-bottom: 8px;
      font-size: 14px;
    }}
    .cluster-row {{
      appearance: none;
      width: 100%;
      border: 0;
      background: transparent;
      display: grid;
      grid-template-columns: 12px 1fr auto;
      gap: 7px;
      align-items: center;
      padding: 3px 0;
      border-bottom: 1px solid rgba(209, 213, 219, 0.55);
      color: inherit;
      cursor: pointer;
      font: inherit;
      text-align: left;
    }}
    .cluster-row:hover {{
      background: rgba(229, 231, 235, 0.55);
    }}
    .cluster-row.is-hidden {{
      opacity: 0.42;
    }}
    .cluster-swatch {{
      width: 10px;
      height: 10px;
      border-radius: 999px;
      border: 1px solid rgba(17, 24, 39, 0.25);
    }}
  </style>
</head>
<body>
  <div id="selection">
    <strong>Selected Weights</strong>
    Click a point to show D, P, T, A.
  </div>
  <div id="clusterIndex"></div>
  <div id="plot"></div>
  <script>
    const traces = {traces_json};
    const layout = {layout_json};
    const clusterSummary = {cluster_summary_json};
    const config = {{
      responsive: true,
      displaylogo: false,
      scrollZoom: true,
      modeBarButtonsToRemove: ["lasso2d", "select2d"]
    }};
    const plot = document.getElementById("plot");
    const selection = document.getElementById("selection");
    const clusterIndex = document.getElementById("clusterIndex");

    if (clusterSummary.length > 0) {{
      clusterIndex.innerHTML = `
        <strong>Clusters (${{clusterSummary.length}})</strong>
        ${{clusterSummary.map((cluster) => `
          <button class="cluster-row" type="button" data-trace-index="${{cluster.traceIndex}}">
            <span class="cluster-swatch" style="background:${{cluster.color}}"></span>
            <span>${{cluster.label}}</span>
            <span>${{cluster.count}}</span>
          </button>
        `).join("")}}
      `;
    }}

    Plotly.newPlot(plot, traces, layout, config).then(() => {{
      plot.on("plotly_click", (eventData) => {{
        const point = eventData.points[0];
        const data = point.customdata;
        if (!data) {{
          return;
        }}
        selection.innerHTML = `
          <strong>Selected Weights</strong>
          <div>${{data[0]}}</div>
          <div class="weights">
            D=${{Number(data[1]).toFixed(2)}} |
            P=${{Number(data[2]).toFixed(2)}} |
            T=${{Number(data[3]).toFixed(2)}} |
            A=${{Number(data[4]).toFixed(2)}}
          </div>
          <div>Distance: ${{Number(data[5]).toFixed(2)}} km</div>
          <div>Score: ${{Number(data[6]).toFixed(3)}}</div>
        `;
      }});

      clusterIndex.querySelectorAll(".cluster-row").forEach((row) => {{
        row.addEventListener("click", () => {{
          const traceIndex = Number(row.dataset.traceIndex);
          const isHidden = row.classList.toggle("is-hidden");
          Plotly.restyle(
            plot,
            {{ visible: isHidden ? "legendonly" : true }},
            [traceIndex]
          );
        }});
      }});
    }});
  </script>
</body>
</html>
"""


def make_route_pair_plot(route_runs, route_pair, output_html):
    """Create one route-pair barycentric plot."""
    subset = route_runs[route_runs["route_pair"] == route_pair].copy()
    route_label = subset["route_pair_label"].iloc[0]
    traffic_dataset = subset["traffic_dataset"].iloc[0]
    variant_counts = subset["route_variant"].value_counts()

    traces = build_frame_traces()
    for color_index, route_variant in enumerate(variant_counts.index):
        traces.append(
            build_variant_trace(
                subset,
                route_variant,
                VARIANT_COLORS[color_index % len(VARIANT_COLORS)],
            )
        )

    layout = {
        "title": {
            "text": (
                "Chicago Route Weight Barycentric Tetrahedron - "
                f"{route_label} ({traffic_dataset})"
            ),
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

    output_html.write_text(build_html(traces, layout), encoding="utf-8")
    return output_html, len(subset), len(variant_counts)


def make_cluster_plot(route_runs, route_pair, method, output_html):
    """Create one route-pair barycentric plot colored by cluster assignment."""
    cluster_column = method["cluster_column"]
    subset = route_runs[route_runs["route_pair"] == route_pair].copy()
    route_label = subset["route_pair_label"].iloc[0]
    traffic_dataset = subset["traffic_dataset"].iloc[0]
    cluster_counts = subset[cluster_column].value_counts()

    traces = build_frame_traces()
    cluster_summary = []
    for color_index, cluster_id in enumerate(cluster_counts.index):
        color = color_for_index(color_index, len(cluster_counts))
        trace_index = len(traces)
        cluster_summary.append(
            {
                "label": cluster_id,
                "count": int(cluster_counts.loc[cluster_id]),
                "color": color,
                "traceIndex": trace_index,
            }
        )
        traces.append(
            build_cluster_trace(
                subset,
                cluster_column,
                cluster_id,
                color,
            )
        )

    layout = {
        "title": {
            "text": (
                "Chicago Route Weight Barycentric Tetrahedron - "
                f"{route_label} ({traffic_dataset}) - {method['label']}"
            ),
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
        "showlegend": False,
        "paper_bgcolor": "white",
        "margin": {"l": 0, "r": 0, "t": 64, "b": 0},
        "hovermode": "closest",
    }

    output_html.write_text(
        build_html(traces, layout, cluster_summary=cluster_summary),
        encoding="utf-8",
    )
    return output_html, len(subset), len(cluster_counts)


def main():
    """Save barycentric coordinates and one plot per Chicago route pair."""
    if not ROUTE_RUNS_CSV.exists():
        raise FileNotFoundError(
            f"Run Chicago weighted routes first: {ROUTE_RUNS_CSV}"
        )

    MAPS_FOLDER.mkdir(parents=True, exist_ok=True)
    route_runs = pd.read_csv(ROUTE_RUNS_CSV)
    require_columns(
        route_runs,
        [
            "route_run_id",
            "route_pair",
            "route_pair_label",
            "traffic_dataset",
            "weight_id",
            "status",
            "path_node_ids",
            *WEIGHT_COLUMNS,
        ],
        str(ROUTE_RUNS_CSV),
    )

    route_runs = route_runs[route_runs["status"] == "success"].copy()
    if route_runs.empty:
        raise ValueError("No successful Chicago route runs found.")

    route_runs = add_barycentric_coordinates(route_runs)
    route_runs = add_route_variant_ids(route_runs)

    output_columns = [
        "route_run_id",
        "route_pair",
        "route_pair_label",
        "route_variant",
        "traffic_dataset",
        "weight_id",
        *WEIGHT_COLUMNS,
        "barycentric_x",
        "barycentric_y",
        "barycentric_z",
    ]
    route_runs[output_columns].to_csv(BARYCENTRIC_CSV, index=False)
    print(f"Saved barycentric coordinates: {BARYCENTRIC_CSV}")

    for route_pair in sorted(route_runs["route_pair"].unique()):
        output_html = MAPS_FOLDER / f"{route_pair}_weight_barycentric.html"
        output_html, route_count, variant_count = make_route_pair_plot(
            route_runs,
            route_pair,
            output_html,
        )
        print(
            f"Saved {output_html}: {route_count} weights, "
            f"{variant_count} exact route variants"
        )

    for method in CLUSTER_METHODS:
        if not method["file"].exists():
            print(f"Skipping {method['label']}; missing {method['file']}")
            continue

        clusters = pd.read_csv(method["file"])
        require_columns(
            clusters,
            ["route_run_id", method["cluster_column"]],
            str(method["file"]),
        )
        cluster_weights = route_runs.merge(
            clusters[["route_run_id", method["cluster_column"]]],
            on="route_run_id",
            how="left",
        )
        if cluster_weights[method["cluster_column"]].isna().any():
            raise ValueError(f"Missing cluster assignments in {method['file']}")

        for route_pair in sorted(cluster_weights["route_pair"].unique()):
            output_html = (
                MAPS_FOLDER
                / f"{route_pair}_weight_barycentric_{method['suffix']}.html"
            )
            output_html, route_count, cluster_count = make_cluster_plot(
                cluster_weights,
                route_pair,
                method,
                output_html,
            )
            print(
                f"Saved {output_html}: {route_count} weights, "
                f"{cluster_count} {method['label']} clusters"
            )


if __name__ == "__main__":
    main()
