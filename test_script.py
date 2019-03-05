from KGutils import (kegg_graph, get_organism_codes,kegg_url,
                     get_nodetype_nodes, connected_components, kegg_infos)

from networkx.algorithms.bipartite import projected_graph

import numpy as np
import pandas as pd

#%%
KEGG = kegg_graph("disease", "hsa")

genes = get_nodetype_nodes(KEGG, "hsa")
diseases = get_nodetype_nodes(KEGG, "disease")
#%%
gene_net = projected_graph(KEGG,genes)
disease_net = projected_graph(KEGG, diseases)

subgraphs = connected_components(disease_net)