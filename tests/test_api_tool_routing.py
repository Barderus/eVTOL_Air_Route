import unittest
import json
import os
import pickle
import tempfile
from unittest.mock import patch

import geopandas as gpd
import networkx as nx
from shapely.geometry import box

from api_tool.astar_route_tool import ASTAR_ROUTE_TOOL, generate_astar_route
from api_tool.run_astar_route import make_route_geometry


class RouteGeometryTests(unittest.TestCase):
    def test_make_route_geometry_rejects_single_cell_path(self):
        projected_grid = gpd.GeoDataFrame(
            {"cell_id": [1]},
            geometry=[box(0, 0, 1, 1)],
            crs="EPSG:3857",
        )

        with self.assertRaisesRegex(
            ValueError,
            "Route path must contain at least two grid cells",
        ):
            make_route_geometry(projected_grid, [0])


class GenerateAstarRouteTests(unittest.TestCase):
    def make_graph_package(self):
        grid = gpd.GeoDataFrame(
            {"cell_id": [0, 1]},
            geometry=[
                box(-90.01, 38.00, -90.00, 38.01),
                box(-90.00, 38.00, -89.99, 38.01),
            ],
            crs="EPSG:4326",
        )
        projected_grid = grid.to_crs("EPSG:3857")
        centroids = projected_grid.geometry.centroid
        graph = nx.Graph()
        graph.add_node(0)
        graph.add_node(1)
        graph.add_edge(
            0,
            1,
            distance_norm=1.0,
            population_norm=0.2,
            traffic_norm=0.3,
            airspace_norm=0.4,
            distance_m=1000.0,
        )
        return {
            "graph": graph,
            "grid": grid,
            "projected_grid": projected_grid,
            "centroid_x": centroids.x.to_numpy(),
            "centroid_y": centroids.y.to_numpy(),
            "maximum_distance": 1000.0,
        }

    def test_generate_astar_route_returns_json_safe_success(self):
        weights = {
            "distance_weight": 0.4,
            "population_weight": 0.3,
            "traffic_weight": 0.2,
            "airspace_weight": 0.1,
        }

        with tempfile.TemporaryDirectory() as temp_dir:
            graph_path = os.path.join(temp_dir, "routing_graph.pkl")
            with open(graph_path, "wb") as graph_file:
                pickle.dump(self.make_graph_package(), graph_file)

            with patch("builtins.input", side_effect=AssertionError("unexpected prompt")):
                result = generate_astar_route(
                    graph_path=graph_path,
                    origin={"lat": 38.005, "lon": -90.005},
                    destination={"lat": 38.005, "lon": -89.995},
                    weights=weights,
                    route_id="test_route",
                )

        self.assertEqual(result["status"], "success")
        self.assertEqual(result["tool"], "generate_astar_route")
        self.assertEqual(result["route"]["route_id"], "test_route")
        self.assertEqual(result["route"]["path_node_ids"], [0, 1])
        self.assertEqual(result["route"]["path_nodes"], 2)
        self.assertEqual(result["route"]["route_geojson"]["type"], "LineString")
        self.assertAlmostEqual(result["route"]["costs"]["route_distance_km"], 1.0)
        json.dumps(result)

    def test_generate_astar_route_returns_failure_without_weights(self):
        with patch("builtins.input", side_effect=AssertionError("unexpected prompt")):
            result = generate_astar_route(
                graph_path="missing.pkl",
                origin={"lat": 38.005, "lon": -90.005},
                destination={"lat": 38.005, "lon": -89.995},
                weights=None,
                route_id="bad_route",
            )

        self.assertEqual(result["status"], "failure")
        self.assertIsNone(result["route"])
        self.assertEqual(result["error"]["type"], "ValueError")

    def test_astar_route_tool_contract_lists_only_supported_inputs(self):
        self.assertEqual(ASTAR_ROUTE_TOOL["name"], "generate_astar_route")
        self.assertEqual(
            sorted(ASTAR_ROUTE_TOOL["inputs"]),
            ["destination", "graph_path", "origin", "route_id", "weights"],
        )


if __name__ == "__main__":
    unittest.main()
