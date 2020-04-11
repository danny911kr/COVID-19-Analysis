import networkx as nx


# <node id="0" label="Myriel">
#         <attvalues>
#           <attvalue for="modularity_class" value="0"></attvalue>
#         </attvalues>
#         <viz:size value="28.685715"></viz:size>
#         <viz:position x="-266.82776" y="299.6904" z="0.0"></viz:position>
#         <viz:color r="235" g="81" b="72"></viz:color>
#       </node>


G = nx.Graph()

G.add_node(0, label='ACE2', modularity_class=0)
G.nodes[0]['viz'] = {'size': '60.0'}
G.add_node(1, label='furin', modularity_class=1)
G.nodes[1]['viz'] = {'size': '60.0'}
G.add_node(2, label='TMPRSS2', modularity_class=2)
G.nodes[2]['viz'] = {'size': '60.0'}
G.add_node(3, label='cathepsin', modularity_class=3)
G.nodes[3]['viz'] = {'size': '60.0'}
G.add_node(4, label='mannose', modularity_class=4)
G.nodes[4]['viz'] = {'size': '60.0'}
G.add_node(5, label='COVID-19', modularity_class=5)
G.nodes[5]['viz'] = {'size': '120.0'}

nx.write_gexf(G, "test.gexf")