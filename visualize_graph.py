import networkx as nx
import matplotlib.pyplot as plt

# Load the graph
G = nx.read_graphml("real_output.graphml")

# Color nodes based on whether they appear in multiple genomes
node_colors = []
for node, data in G.nodes(data=True):
    genomes = data.get("genomes", "")
    count = len(genomes.split(",")) if isinstance(genomes, str) else 1
    node_colors.append("tomato" if count > 1 else "lightgray")

# Draw the graph
plt.figure(figsize=(14, 8))
nx.draw(
    G,
    with_labels=True,
    node_color=node_colors,
    node_size=500,
    font_size=8,
    arrows=True
)
plt.title("Pangenome Graph (shared = red, unique = gray)")
plt.tight_layout()
plt.show()
