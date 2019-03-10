import KEGGutils as kg

import pandas as pd
#%%
dis_gen = kg.kegg_graph("disease", "hsa")

genes_enzymes = kg.kegg_graph("hsa", "enzyme")

enzyme_reaction = kg.kegg_graph("enzyme", "reaction")

    #%%
df = pd.DataFrame(columns =['disease','enzymes','enzymegraph', 'totnodes', 'totedges', "ncliques", "radius", "diameter"])


kegg_diseases = kg.get_list("disease")
disease = 'ds:H00773'

ds_genes = list(dis_gen[disease])

ds_gene_enzyme_graph = kg.neighbor_graph(genes_enzymes,ds_genes, name = "DSneig_graph_genes_enzymes")


ds_gene_enzyme_graph = kg.descendant_graph(genes_enzymes, ds_genes, name = 'dsgene_enzyme')

ds_enzymes = kg.get_nodetype_nodes(ds_gene_enzyme_graph, 'enzyme')

ds_enzyme_reaction = kg.descendant_graph(enzyme_reaction, ds_enzymes, name = 'DSenzyme_reaction')

ds_reactions = kg.get_nodetype_nodes(ds_enzyme_reaction, 'reaction')

ds_enzyme_enzyme = kg.projected_graph(ds_enzyme_reaction, ds_enzymes)

#    kg.draw(ds_enzyme_enzyme)

s_element = pd.Series(data = kg.graph_measures(ds_enzyme_enzyme))
s_element['enzymegraph'] = ds_enzyme_enzyme
s_element['disease'] = disease
s_element['enzymes'] = ds_enzymes

df.loc[df.shape[0]] = s_element


df.set_index('disease', inplace = True)

#%%
# malattia a maggiore numero di cliques
max_cliques_disease = df[df['ncliques']==df['ncliques'].max()].index

kg.draw(df.loc[max_cliques_disease]['enzymegraph'][0], layout = "random_layout")



