import json
import sys
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
                self.assertEqual(tool_names, ["generate_astar_route"])

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
            call_result.structured_content["route"]["route_geojson"]["type"],
            "LineString",
        )
        json.loads(call_result.content[0].text)


if __name__ == "__main__":
    unittest.main()
