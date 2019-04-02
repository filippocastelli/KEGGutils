#!/usr/bin/env python
# coding: utf-8

# # KEGGutils Tutorial 1 : Enzymatic Correlation Graphs
# ## ...using just KEGGutils basic functions
# ***

# ### First things first
# Let's import `KEGGutils` and `networkx`, on which `KEGGutils` is based, with aliases ( with maybe some module path scope adjusting first )

# In[1]:


import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)


# In[2]:


import KEGGutils as kg
import networkx as nx


# check that the latest version of `KEGGutils` is installed

# In[3]:


kg.__version__


# and we're ready to start our tutorial

# let's remove the cached files to use the freshest available data

# In[4]:


kg.delete_cached_files()


# ### The interesting stuff

# let's download the disease list from KEGG:

# In[5]:


kegg_diseases = kg.keggapi_list("disease")


# In[6]:


kegg_diseases[:10] 


# we can obtain some description on the `ds:H00773` disease or really any kind of KEGG entry using `get_infos()`

# In[7]:


print(kg.get_infos.__doc__)


# In[8]:


kg.get_infos("ds:H00773", verbose = False)


# enabling the `verbose` option will show the full description.

# let's try to download the bipartite graph linking each disease to a set of genes by using `kegg_graph()`, you must specify the source and the target categories: note that human genes are referred with the `hsa` key.
# 
# You can use the `force_download` option to download the file again overwriting previous copies that you may have already downloaded

# In[9]:


print(kg.kegg_link_graph.__doc__)


# In[10]:


dis_gene0 = kg.KEGGlinkgraph(source_db = "disease", target_db = "hsa")


# In[11]:


dis_gene = kg.kegg_link_graph("disease", "hsa", force_download = False)


# All `KEGGutils` graph nodes have a `nodetype` attribute which helps us differentiate between different objects in a graph:

# In[12]:


nx.get_node_attributes(dis_gene, "nodetype")['ds:H00773']


# In[13]:


dis_gene0.nodes['ds:H00773']


# the list of nodes that are linked to a particular given node can be obtained the same way as with every `networkx` graph

# In[14]:


dis_gene0['ds:H00773']


# In[15]:


dis_gene['ds:H00773']


# or using `kg.linked_nodes`, which returns the list of connected nodes, or optionally a dictionary with each node and its `nodetype`

# In[16]:


print(kg.linked_nodes.__doc__)


# In[17]:


ds_genes0= dis_gene0.linked_nodes("ds:H00773")
ds_genes0[:10]


# In[18]:


ds_genes = kg.linked_nodes(dis_gene, 'ds:H00773', return_dict = True)
ds_genes_list = kg.linked_nodes(dis_gene, 'ds:H00773', return_dict = False)


# in this case we obviously expect all the linked node to be human genes marked with `hsa`

# In[19]:


ds_genes_list[:10]


# we can find every enzyme associated with each of these genes using the gene-enzyme database from **KEGG**, which we download the same way as the previous one

# In[20]:


gene_enzyme0 = kg.KEGGlinkgraph(source_db = "hsa", target_db = "enzyme")


# In[21]:


gene_enzyme = kg.kegg_link_graph("hsa", "enzyme")


# it's possible to use the gene list we obtained before to narrow down a search on the complete KEGG *gene-enzyme* database with `kg.neighbor_graph()`: 

# In[22]:


print(kg.neighbor_graph.__doc__)


# In[23]:


ds_gene_enzyme0 = gene_enzyme0.neighbor_graph(ds_genes)


# In[24]:


ds_gene_enzyme = kg.neighbor_graph(gene_enzyme,
                                   ds_genes,
                                   keep_isolated_nodes = True,
                                   name = "ds_gene_enzyme")


# In[25]:


ds_gene_enzyme0.nodes


# In[26]:


ds_gene_enzyme.nodes


# the function returns a subgraph of `gene_enzyme` selecting both `ds_genes` nodes and their neighbors in the original graph.
# 
# Not all `ds_genes` are actually present in `gene_enzyme`, we can represent them as isolated nodes in the new graph using the `keep_isolated_nodes` option

# Graphic functionality is provided by `kg.draw()`

# In[27]:


print(kg.draw.__doc__)


# In[28]:


ds_gene_enzyme0.draw(layout = "random_layout")


# In[29]:


kg.draw(ds_gene_enzyme, layout = "random_layout")


# we see that the actual common enzymes between those obtained from `dis_gen` and those we found on `gene_enzyme` are very few

# In[30]:


set(ds_genes0) & set(gene_enzyme0.nodes)


# In[31]:


set(ds_genes.keys()) & set(gene_enzyme.nodes)


# Le'ts download the *enzyme-reaction* graph from **KEGG** for a last step

# In[32]:


enzyme_reaction0 = kg.KEGGlinkgraph(source_db = "enzyme", target_db = "reaction")


# In[33]:


enzyme_reaction = kg.kegg_link_graph("enzyme", "reaction")


# using the same mechanism as before we get the *enzyme-reaction* graph

# we can get every node of a particular `nodetype` in a graph using `kg.get_nodes_by_nodetype`:

# In[34]:


print(kg.get_nodes_by_nodetype.__doc__)


# In[35]:


ds_enzymes0 = ds_gene_enzyme0.list_by_nodetype("enzyme")
ds_enzymes0


# In[36]:


ds_enzymes = kg.get_nodes_by_nodetype(ds_gene_enzyme, "enzyme", return_dict = True)


# In[37]:


list(ds_enzymes.keys())


# the same `neighbor_graph()` function as before is used to calculated a subgraph of `enzyme_reaction`

# In[38]:


ds_enzyme_reaction0 = enzyme_reaction0.neighbor_graph(ds_enzymes0)


# In[39]:


ds_enzyme_reaction0.draw(layout = "kamada_kawai_layout")


# In[40]:


ds_enzyme_reaction = kg.neighbor_graph(enzyme_reaction,
                                       ds_enzymes,
                                       keep_isolated_nodes = True,
                                       name = "ds_enzyme_reaction" )


# In[41]:


kg.draw(ds_enzyme_reaction, layout = "kamada_kawai_layout")


# we can count all the reactions in the graph

# In[42]:


ds_reactions0 = ds_enzyme_reaction0.list_by_nodetype("reaction")
ds_reactions0


# In[43]:


ds_reactions = kg.get_nodes_by_nodetype(ds_enzyme_reaction, "reaction", return_dict = True)
list(ds_reactions.keys())


# we want now to create a projected graph in which we link two enzymes if they appear in the same reaction: we can use the `projected_graph()`function that does exactly that

# In[44]:


print(kg.projected_graph.__doc__)


# here we want to project the graph on the enzymes set that we obtain from `ds_gene_enzyme` using `get_nodetype_nodes()`

# In[45]:


ds_enzymes0 = ds_gene_enzyme0.list_by_nodetype("enzyme")


# In[46]:


ds_enzymes = kg.get_nodes_by_nodetype(ds_gene_enzyme, "enzyme", return_dict = True)


# In[47]:


list(ds_enzymes.keys())


# now we can finally project `ds_enzyme_reaction` onto the `ds_enzymes` set to get our enzyme correlation graph

# In[55]:


ds_enzyme_enzyme = kg.projected_graph(ds_enzyme_reaction, ds_enzymes, name = "ds_enzyme_enzyme")


# In[56]:


ds_enzyme_enzyme1 = ds_enzyme_reaction0.projected_graph(ds_enzymes)


# In[49]:


print(ds_enzyme_reaction0.source_db)


# In[57]:


ds_enzyme_enzyme0 = ds_enzyme_reaction0.projected_graph()


# In[ ]:


ds_enzyme_enzyme.edges


# In[ ]:


ds_enzyme_enzyme0.edges


# In[ ]:


kg.draw(ds_enzyme_enzyme)


# In[ ]:




