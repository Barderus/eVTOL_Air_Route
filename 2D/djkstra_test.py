import networkx as nx
import matplotlib.pyplot as plt


# MultiGraph keeps room for multiple parallel road segments between the same
# towns, even though this toy example usually adds one edge per pair.
chicago = nx.MultiGraph()

chicago.add_nodes_from(["Bolingbrook", "Romeoville", "Woodridge", "Downers Grove", "Lisle", "Naperville", "Darien",
                            "Westmont", "Hinsdale", "Willowbrook", "Burr Ridge", "Countryside", "La Grange", "Brookfield",
                            "Lyons", "Forest View", "Berwyn", "Cicero", "Oak Park", "Chicago"])

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

path = nx.dijkstra_path(chicago, "Bolingbrook", "Chicago", weight="weight")

step_costs = []
total = 0
for u, v in zip(path, path[1:]):
    # Dijkstra returns node order; for a MultiGraph we still need the lightest
    # concrete edge between each pair when summarizing the path cost.
    w = min(data.get("weight", 1) for data in chicago.get_edge_data(u, v).values())
    step_costs.append((u, v, w))
    total += w

print("Path:", path)
print("Steps:", step_costs)
print("Total:", total)

