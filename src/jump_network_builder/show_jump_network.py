import networkx as nx
from networkx.readwrite.multiline_adjlist import read_multiline_adjlist
import matplotlib.pyplot as plt

# adjlist = "../../tests/bcc_int.adjlist"
adjlist = "../../examples/hcp_vac/assets/jump_structure.adjlist"
G = read_multiline_adjlist(
    adjlist,
    comments="#", create_using=nx.MultiDiGraph
)
# G = nx.dodecahedral_graph()
labels = {0: 'A', 1: 'C', 2: 'D'}

for edge in G.edges():
    print(f"{labels[int(edge[0])]}-{labels[int(edge[1])]}")
nx.draw(G)
plt.draw()

