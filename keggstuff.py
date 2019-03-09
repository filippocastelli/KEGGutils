import KEGGutils as kg

import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
#%%
dis_gen = kg.kegg_graph("disease", "hsa")

genes_enzymes = kg.kegg_graph("hsa", "enzyme")

enzyme_reaction = kg.kegg_graph("enzyme", "reaction")

#%%

disease ='ds:H00039'

kg.get_infos(disease,verbose = False)

bc_genes = list(dis_gen[disease])

kg.get_unique_nodetypes(genes_enzymes)

bc_gene_enzyme_graph = kg.descendant_graph(genes_enzymes, bc_genes)

bc_enzymes = kg.get_nodetype_nodes(bc_gene_enzyme_graph, 'enzyme')

bc_enzyme_reaction = kg.descendant_graph(enzyme_reaction, bc_enzymes)

bc_reactions = kg.get_nodetype_nodes(bc_enzyme_reaction, 'reaction')

bc_enzyme_enzyme = kg.projected_graph(enzyme_reaction, bc_enzymes)

kg.draw(bc_enzyme_enzyme)






