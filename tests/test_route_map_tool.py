import os
import json
import tempfile
import unittest

import geopandas as gpd
from shapely.geometry import LineString

from api_tool.route_map_tool import ROUTE_MAP_TOOL, create_route_map


class RouteMapToolTests(unittest.TestCase):
    def make_route_geojson(self, output_path):
        route = gpd.GeoDataFrame(
            [
                {
                    "route_id": "test_route",
                    "total_cost": 1.5,
                    "route_distance_km": 2.0,
                    "path_nodes": 2,
                }
            ],
            geometry=[LineString([(-88.1, 41.7), (-87.7, 41.9)])],
            crs="EPSG:4326",
        )
        route.to_file(output_path, driver="GeoJSON")

    def test_create_route_map_returns_structured_success(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            route_geojson_path = os.path.join(temp_dir, "routes.geojson")
            output_html_path = os.path.join(temp_dir, "route_map.html")
            self.make_route_geojson(route_geojson_path)

            result = create_route_map(
                route_geojson_path=route_geojson_path,
                output_html_path=output_html_path,
                locations={
                    "origin": {
                        "label": "Origin",
                        "lat": 41.7,
                        "lon": -88.1,
                        "type": "origin",
                    }
                },
                title="Test Route Map",
            )

            self.assertEqual(result["status"], "success")
            self.assertEqual(result["tool"], "make_route_map")
            self.assertEqual(result["map"]["route_ids"], ["test_route"])
            self.assertEqual(result["map"]["route_geojson_path"], route_geojson_path)
            self.assertTrue(result["map"]["output_exists"])
            self.assertEqual(result["map"]["output_html_path"], output_html_path)
            self.assertEqual(result["map"]["route_count"], 1)
            self.assertEqual(result["map"]["grid_cell_count"], 0)
            self.assertIsNone(result["error"])
            self.assertTrue(os.path.exists(output_html_path))
            json.dumps(result)

    def test_create_route_map_returns_failure_without_output_path(self):
        result = create_route_map(
            route_geojson_path="routes.geojson",
            output_html_path=None,
        )

        self.assertEqual(result["status"], "failure")
        self.assertEqual(result["error"]["type"], "ValueError")
        self.assertIsNone(result["map"])

    def test_route_map_tool_contract_lists_supported_inputs(self):
        self.assertEqual(ROUTE_MAP_TOOL["name"], "make_route_map")
        self.assertEqual(
            sorted(ROUTE_MAP_TOOL["inputs"]),
            [
                "locations",
                "output_html_path",
                "route_geojson_path",
                "scored_grid_path",
                "title",
            ],
        )


if __name__ == "__main__":
    unittest.main()
