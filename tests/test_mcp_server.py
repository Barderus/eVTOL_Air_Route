import json
import os
import sys
import tempfile
import unittest

from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client


class McpServerTests(unittest.IsolatedAsyncioTestCase):
    async def test_mcp_server_lists_and_calls_generate_astar_route(self):
        server_params = StdioServerParameters(
            command=sys.executable,
            args=["-m", "agent.mcp_server"],
            cwd=".",
        )

        async with stdio_client(server_params) as (read_stream, write_stream):
            async with ClientSession(read_stream, write_stream) as session:
                await session.initialize()

                tools_result = await session.list_tools()
                tool_names = [tool.name for tool in tools_result.tools]
                self.assertEqual(tool_names, ["generate_astar_route", "make_route_map"])
                self.assertIn(
                    "output_geojson_path",
                    tools_result.tools[0].input_schema["properties"],
                )
                self.assertIn(
                    "output_directory",
                    tools_result.tools[0].input_schema["properties"],
                )

                call_result = await session.call_tool(
                    "generate_astar_route",
                    arguments={
                        "graph_path": "api_tool/output/chicago_smoke/routing_graph.pkl",
                        "origin": {
                            "lat": 41.695923717435235,
                            "lon": -88.12876224517822,
                        },
                        "destination": {
                            "lat": 41.87838051825937,
                            "lon": -87.63905525207521,
                        },
                        "weights": {
                            "distance_weight": 0.4,
                            "population_weight": 0.3,
                            "traffic_weight": 0.2,
                            "airspace_weight": 0.1,
                        },
                        "route_id": "chicago_mcp_smoke_test",
                    },
                )

        self.assertFalse(call_result.is_error)
        self.assertEqual(call_result.structured_content["status"], "success")
        self.assertEqual(
            call_result.structured_content["route"]["route_id"],
            "chicago_mcp_smoke_test",
        )
        self.assertEqual(call_result.structured_content["route"]["path_nodes"], 113)
        self.assertEqual(
            call_result.structured_content["route"]["requested"]["origin"]["lat"],
            41.695923717435235,
        )
        self.assertEqual(
            call_result.structured_content["route"]["snapped"]["origin"]["node_id"],
            13722,
        )
        self.assertEqual(
            call_result.structured_content["route"]["geometry_summary"]["geometry_type"],
            "LineString",
        )
        self.assertAlmostEqual(
            call_result.structured_content["route"]["distance_km"],
            68.6335136523831,
        )
        self.assertAlmostEqual(
            call_result.structured_content["route"]["total_cost"],
            51.9725486283908,
        )
        self.assertNotIn("geometry_geojson", call_result.structured_content["route"])
        self.assertNotIn("path_coordinates", call_result.structured_content["route"])
        self.assertIsNone(
            call_result.structured_content["route"]["outputs"]["route_map_path"],
        )
        json.loads(call_result.content[0].text)

    async def test_mcp_server_calls_make_route_map(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_directory = os.path.join(temp_dir, "qwen_mcp_test")
            route_geojson_path = os.path.join(
                output_directory,
                "route.geojson",
            ).replace("\\", "/")
            output_html_path = os.path.join(temp_dir, "route_map.html")

            server_params = StdioServerParameters(
                command=sys.executable,
                args=["-m", "agent.mcp_server"],
                cwd=".",
            )

            async with stdio_client(server_params) as (read_stream, write_stream):
                async with ClientSession(read_stream, write_stream) as session:
                    await session.initialize()
                    route_result = await session.call_tool(
                        "generate_astar_route",
                        arguments={
                            "graph_path": "api_tool/output/chicago_smoke/routing_graph.pkl",
                            "origin": {
                                "lat": 41.695923717435235,
                                "lon": -88.12876224517822,
                            },
                            "destination": {
                                "lat": 41.87838051825937,
                                "lon": -87.63905525207521,
                            },
                            "weights": {
                                "distance_weight": 0.4,
                                "population_weight": 0.3,
                                "traffic_weight": 0.2,
                                "airspace_weight": 0.1,
                            },
                            "route_id": "mcp_map_test",
                            "output_directory": output_directory,
                        },
                    )
                    route_geojson_from_tool = route_result.structured_content["outputs"][
                        "route_geojson_path"
                    ]
                    call_result = await session.call_tool(
                        "make_route_map",
                        arguments={
                            "route_geojson_path": route_geojson_from_tool,
                            "output_html_path": output_html_path,
                            "title": "MCP Map Test",
                        },
                    )

            self.assertFalse(route_result.is_error)
            self.assertEqual(route_geojson_from_tool, route_geojson_path)
            self.assertTrue(os.path.exists(route_geojson_path))
            self.assertFalse(call_result.is_error)
            self.assertEqual(call_result.structured_content["status"], "success")
            self.assertEqual(call_result.structured_content["tool"], "make_route_map")
            self.assertEqual(
                call_result.structured_content["map"]["route_ids"],
                ["mcp_map_test"],
            )
            self.assertEqual(
                call_result.structured_content["map"]["route_geojson_path"],
                route_geojson_path,
            )
            self.assertTrue(call_result.structured_content["map"]["output_exists"])
            self.assertEqual(
                call_result.structured_content["map"]["output_html_path"],
                output_html_path,
            )
            self.assertEqual(call_result.structured_content["map"]["route_count"], 1)
            self.assertTrue(os.path.exists(output_html_path))


if __name__ == "__main__":
    unittest.main()
