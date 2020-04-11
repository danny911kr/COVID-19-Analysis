import networkx as nx

ace2 = {'59272': 2}
furin = {}
tmprss2 = {}
cathepsin = {'1508': 1, '1509': 1, '1513': 1, '1514': 2}
mannose = {}

G = nx.Graph()
G.add_node(0, label='COVID-19', modularity_class=5)
G.nodes[0]['viz'] = {'size': '120.0'}

print(len(G))
if len(ace2) > 0:
    for a in ace2:
        G.add_node(len(G), label=a, modularity_class=0)
        G.nodes[len(G) - 1]['viz'] = {'size': '40.0'}
        G.add_edge(0, len(G)-1, weight=ace2[a], label='ace2')

if len(furin) > 0:
    for a in furin:
        G.add_node(len(G), label=a, modularity_class=1)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        G.add_edge(0, len(G)-1, weight=furin[a], label='ace2')

if len(tmprss2) > 0:
    for a in tmprss2:
        G.add_node(len(G), label=a, modularity_class=2)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        G.add_edge(0, len(G)-1, weight=tmprss2[a], label='ace2')

if len(cathepsin) > 0:
    for a in cathepsin:
        G.add_node(len(G), label=a, modularity_class=3)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        G.add_edge(0, len(G)-1, weight=cathepsin[a], label='ace2')

if len(mannose) > 0:
    for a in mannose:
        G.add_node(len(G), label=a, modularity_class=4)
        G.nodes[len(G)-1]['viz'] = {'size': '40.0'}
        G.add_edge(0, len(G)-1, weight=mannose[a], label='ace2')

nx.write_gexf(G, "test.gexf")