import networkx as nx
import math
import matplotlib.pyplot as plt


chicago = nx.MultiGraph()

chicago.add_nodes_from([
    ("Bolingbrook", {"pos": (41.6811, -88.1236)}),
    ("Romeoville", {"pos": (41.647812, -88.087128)}),
    ("Woodridge", {"pos": (41.750179, -88.043602)}),
    ("Downers Grove", {"pos": (41.808861, -88.011131)}),
    ("Lisle", {"pos": (41.801250, -88.074783)}),
    ("Naperville", {"pos": (41.750839, -88.153534)}),
    ("Darien", {"pos": (41.7456, -87.9811)}),
    ("Westmont", {"pos": (41.7947, -87.9742)}),
    ("Hinsdale", {"pos": (41.8006, -87.9273)}),
    ("Willowbrook", {"pos": (41.7589, -87.9492)}),
    ("Burr Ridge", {"pos": (41.7558, -87.9183)}),
    ("Countryside", {"pos": (41.7789, -87.8739)}),
    ("La Grange", {"pos": (41.8081, -87.8733)}),
    ("Brookfield", {"pos": (41.8228, -87.8475)}),
    ("Lyons", {"pos": (41.8133, -87.8220)}),
    ("Forest View", {"pos": (41.8075, -87.7861)}),
    ("Berwyn", {"pos": (41.8425, -87.7900)}),
    ("Cicero", {"pos": (41.8444, -87.7592)}),
    ("Oak Park", {"pos": (41.8883, -87.7894)}),
     ("Chicago",  {"pos": (41.8820, -87.6278)})
        ])


# Naperville
chicago.add_edge("Naperville", "Romeoville", weight=9)
chicago.add_edge("Naperville", "Lisle", weight=7)
chicago.add_edge("Naperville", "Bolingbrook", weight=8)

# Lisle
chicago.add_edge("Lisle", "Bolingbrook", weight=7)
chicago.add_edge("Lisle", "Downers Grove", weight=4)
chicago.add_edge("Lisle", "Woodridge", weight=4)

# Romeoville
chicago.add_edge("Romeoville", "Bolingbrook", weight=5)

# Bolingbrook
chicago.add_edge("Bolingbrook", "Downers Grove", weight=14)
chicago.add_edge("Bolingbrook", "Woodridge", weight=4)
chicago.add_edge("Bolingbrook", "Darien", weight=8)

# Downers Grove
chicago.add_edge("Downers Grove", "Darien", weight=6)
chicago.add_edge("Downers Grove", "Willowbrook", weight=7)
chicago.add_edge("Downers Grove", "Westmont", weight=2)
chicago.add_edge("Downers Grove", "Woodridge", weight=7)

# Woodridge
chicago.add_edge("Woodridge", "Darien", weight=4)

# Darien
chicago.add_edge("Darien", "Willowbrook", weight=2)
chicago.add_edge("Darien", "Burr Ridge", weight=5)
chicago.add_edge("Darien", "Westmont", weight=3)

# Westmont
chicago.add_edge("Westmont", "Willowbrook", weight=4)
chicago.add_edge("Westmont", "Hinsdale", weight=2)

# Hinsdale
chicago.add_edge("Hinsdale", "La Grange", weight=3)
chicago.add_edge("Hinsdale", "Willowbrook", weight=3)

# Willowbrook
chicago.add_edge("Willowbrook", "Burr Ridge", weight=3)

# Burr Ridge
chicago.add_edge("Burr Ridge", "La Grange", weight=1)
chicago.add_edge("Burr Ridge", "Countryside", weight=1)

# Countryside
chicago.add_edge("Countryside", "La Grange", weight=6)
chicago.add_edge("Countryside", "Lyons", weight=9)

# La Grange
chicago.add_edge("La Grange", "Brookfield", weight=2)
chicago.add_edge("La Grange", "Lyons", weight=3)

# Brookfield
chicago.add_edge("Brookfield", "Lyons", weight=2)
chicago.add_edge("Brookfield", "Berwyn", weight=4)

# Lyons
chicago.add_edge("Lyons", "Forest View", weight=2)
chicago.add_edge("Lyons", "Berwyn", weight=3)

# Forest View
chicago.add_edge("Forest View", "Cicero", weight=5)
chicago.add_edge("Forest View", "Berwyn", weight=3)

# Berwyn
chicago.add_edge("Berwyn", "Oak Park", weight=3)
chicago.add_edge("Berwyn", "Cicero", weight=2)

# Cicero
chicago.add_edge("Cicero", "Oak Park", weight=3)
chicago.add_edge("Cicero", "Chicago", weight=8)

# Oak Park
chicago.add_edge("Oak Park", "Chicago", weight=9)

#print(chicago.nodes)
nx.draw(chicago, with_labels=True, font_weight='bold')
plt.savefig("chicago.png")
plt.show()

# Harvesine formula allows us to transform our degree coord into miles
def dist(u, v):
    lat1, lon1 = chicago.nodes[u]["pos"]
    lat2, lon2 = chicago.nodes[v]["pos"]

    # Convert degrees → radians
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    # Haversine formula
    a = math.sin(dphi / 2) ** 2 + \
        math.cos(phi1) * math.cos(phi2) * \
        math.sin(dlambda / 2) ** 2

    c = 2 * math.asin(math.sqrt(a))

    R = 3958.8  # Earth radius in miles
    return R * c

path = nx.astar_path(chicago, "Bolingbrook", "Chicago", heuristic=dist,  weight="weight")

step_costs = []
total = 0
for u, v in zip(path, path[1:]):
    w = min(data.get("weight", 1) for data in chicago.get_edge_data(u, v).values())
    step_costs.append((u, v, w))
    total += w

print("Path:", path)
print("Steps:", step_costs)
print("Total:", total)

