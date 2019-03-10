from KEGGutils import (kegg_graph, get_organism_codes,kegg_url,
                     get_nodetype_nodes, connected_components, kegg_infos)

#%%
print("The human genes vs diseases db is at", kegg_url("disease", "hsa"))

KEGG = kegg_graph("disease", "hsa")

subgraphs = connected_components(KEGG)


genes = get_nodetype_nodes(KEGG, "hsa")
diseases = get_nodetype_nodes(KEGG, "disease")

kegg_infos('hsa:64109', verbose = True)