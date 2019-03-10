import KEGGutils as kg

import pandas as pd
import tqdm
#%%
dis_gen = kg.kegg_graph("disease", "hsa")

genes_enzymes = kg.kegg_graph("hsa", "enzyme")

enzyme_reaction = kg.kegg_graph("enzyme", "reaction")

    #%%
df = pd.DataFrame(columns =['disease','enzymes','enzymegraph', 'totnodes', 'totedges', "ncliques", "radius", "diameter"])


kegg_diseases = kg.get_list("disease")

disease_list = list(set(kegg_diseases) & set(list(kg.get_nodetype_nodes(dis_gen, 'disease'))))
#disease_list = ['ds:H00773']
print("Scanning over all KEGG diseases...")
for disease in tqdm.tqdm(disease_list):
    
#    kg.get_infos(disease,verbose = False)
    try:
        ds_genes = list(dis_gen[disease])
        
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
        
    except (kg.NoDescendantError, kg.NoProjectedError):
        pass

df.set_index('disease', inplace = True)

#%%
# malattia a maggiore numero di cliques
max_cliques_disease = df[df['ncliques']==df['ncliques'].max()].index

kg.draw(df.loc[max_cliques_disease]['enzymegraph'][0], layout = "random_layout")






