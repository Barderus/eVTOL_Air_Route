"""MCP server exposing eVTOL route tools."""

import json

import anyio
from mcp import types
from mcp.server.lowlevel import Server
from mcp.server.stdio import stdio_server

from api_tool.astar_route_tool import ASTAR_ROUTE_TOOL, generate_astar_route


GENERATE_ASTAR_ROUTE_INPUT_SCHEMA = {
    "type": "object",
    "properties": {
        "graph_path": {
            "type": "string",
            "description": "Path to a saved routing graph .pkl package.",
        },
        "origin": {
            "type": "object",
            "description": "Origin coordinate.",
            "properties": {
                "lat": {"type": "number"},
                "lon": {"type": "number"},
            },
            "required": ["lat", "lon"],
            "additionalProperties": False,
        },
        "destination": {
            "type": "object",
            "description": "Destination coordinate.",
            "properties": {
                "lat": {"type": "number"},
                "lon": {"type": "number"},
            },
            "required": ["lat", "lon"],
            "additionalProperties": False,
        },
        "weights": {
            "type": "object",
            "description": "Route cost weights. Values must sum to 1.0.",
            "properties": {
                "distance_weight": {"type": "number"},
                "population_weight": {"type": "number"},
                "traffic_weight": {"type": "number"},
                "airspace_weight": {"type": "number"},
            },
            "required": [
                "distance_weight",
                "population_weight",
                "traffic_weight",
                "airspace_weight",
            ],
            "additionalProperties": False,
        },
        "route_id": {
            "type": "string",
            "description": "Caller-supplied route identifier.",
        },
    },
    "required": ["graph_path", "origin", "destination", "weights", "route_id"],
    "additionalProperties": False,
}


GENERATE_ASTAR_ROUTE_OUTPUT_SCHEMA = {
    "type": "object",
    "properties": {
        "status": {"type": "string", "enum": ["success", "failure"]},
        "tool": {"type": "string"},
        "route": {"type": ["object", "null"]},
        "error": {"type": ["object", "null"]},
    },
    "required": ["status", "tool", "route", "error"],
}


async def list_tools(_context, _params):
    """Return the tools exposed by this MCP server."""
    return types.ListToolsResult(
        tools=[
            types.Tool(
                name=ASTAR_ROUTE_TOOL["name"],
                description=ASTAR_ROUTE_TOOL["description"],
                inputSchema=GENERATE_ASTAR_ROUTE_INPUT_SCHEMA,
                outputSchema=GENERATE_ASTAR_ROUTE_OUTPUT_SCHEMA,
            )
        ]
    )


async def call_tool(_context, params):
    """Dispatch MCP tool calls to local Python wrappers."""
    if params.name != ASTAR_ROUTE_TOOL["name"]:
        error_result = {
            "status": "failure",
            "tool": params.name,
            "route": None,
            "error": {
                "type": "UnknownToolError",
                "message": f"Unknown tool: {params.name}",
            },
        }
        return make_tool_result(error_result, is_error=True)

    arguments = params.arguments or {}
    result = generate_astar_route(
        graph_path=arguments.get("graph_path"),
        origin=arguments.get("origin"),
        destination=arguments.get("destination"),
        weights=arguments.get("weights"),
        route_id=arguments.get("route_id"),
    )
    return make_tool_result(result, is_error=result["status"] != "success")


def make_tool_result(result, is_error=False):
    """Create an MCP tool result with structured and text content."""
    return types.CallToolResult(
        content=[
            types.TextContent(
                text=json.dumps(result, indent=2),
            )
        ],
        structuredContent=result,
        isError=is_error,
    )


server = Server(
    "evtol-route-tools",
    version="0.1.0",
    description="MCP tools for the eVTOL route API prototype.",
    on_list_tools=list_tools,
    on_call_tool=call_tool,
)


async def run_stdio_server():
    """Run the MCP server over stdio."""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options(),
        )


def main():
    """CLI entry point."""
    anyio.run(run_stdio_server)


if __name__ == "__main__":
    main()
