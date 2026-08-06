import unittest

from api_tool.create_routing_grid import create_grid_from_bounds


class RoutingGridTests(unittest.TestCase):
    def test_chicago_grid_uses_local_meter_crs(self):
        bounds = {
            "west": -88.34,
            "south": 41.48,
            "east": -87.52,
            "north": 42.21,
        }

        grid = create_grid_from_bounds(bounds, cell_size_m=500)

        self.assertEqual(len(grid), 22363)
        self.assertEqual(str(grid.crs), "EPSG:4326")
        west, south, east, north = grid.total_bounds.tolist()
        self.assertAlmostEqual(west, bounds["west"], places=4)
        self.assertAlmostEqual(south, bounds["south"], places=4)
        self.assertAlmostEqual(east, bounds["east"], places=4)
        self.assertAlmostEqual(north, 42.210733, places=4)


if __name__ == "__main__":
    unittest.main()
