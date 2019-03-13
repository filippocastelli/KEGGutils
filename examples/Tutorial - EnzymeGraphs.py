#!/usr/bin/env python
# coding: utf-8

# # KEGGutils Tutorial: Enzymatic correlation graphs 
# ***

# ## Getting Started
# 
# ### Installing KEGGutils
# `KEGGutils` is available as a **PyPi** package  and can be easily installed using `pip` via `pip install KEGGutils`,
# this should be enough to get you going. 
# Alternatively you can visit the project's **Github** https://github.com/filippocastelli/KEGGutils and the **PyPi** page https://pypi.org/project/KEGGutils/ for a manual install.

# ### First things first
# Let's import `KEGGutils` and `networkx`, on which `KEGGutils` is based, with aliases

# In[1]:


import KEGGutils as kg
import networkx as nx


# check that the latest version of `KEGGutils` is installed

# In[2]:


kg.__version__


# we should be ready to go!
# 
# we can set a directory for our downloaded files, or we can use the default one:

# In[3]:


kg.download_dir


# please notice that the default directory is a runtime variable, remember to **change it every time `KEGGutils` is reloaded** (maybe we'll change that in the future).  

# to ensure that we use the latest available data we can remove all previously downloaded files, if it's your first time using `KEGGutils` you should have no problems.

# In[4]:


kg.delete_cached_files()


# ### The interesting stuff

# let's download the disease list from KEGG:

# In[5]:


kegg_diseases = kg.get_list("disease")


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


print(kg.kegg_graph.__doc__)


# In[10]:


dis_gene = kg.kegg_graph("disease", "hsa", force_download = False)


# All `KEGGutils` graph nodes have a `nodetype` attribute which helps us differentiate between different objects in a graph:

# In[11]:


nx.get_node_attributes(dis_gene, "nodetype")['ds:H00773']


# the list of nodes that are linked to a particular given node can be obtained the same way as with every `networkx` graph

# In[12]:


dis_gene['ds:H00773']


# or using `kg.linked_nodes`, which returns a dictionary with each node and its `nodetype`

# In[13]:


print(kg.linked_nodes.__doc__)


# In[14]:


ds_genes = kg.linked_nodes(dis_gene, 'ds:H00773')


# in this case we obviously expect all the linked node to be human genes marked with `hsa`

# In[15]:


list(ds_genes.items())[:10]


# we can find every enzyme associated with each of these genes using the gene-enzyme database from **KEGG**, which we download the same way as the previous one

# In[16]:


gene_enzyme = kg.kegg_graph("hsa", "enzyme")


# it's possible to use the gene list we obtained before to narrow down a search on the complete KEGG *gene-enzyme* database with `kg.neighbor_graph()`: 

# In[17]:


print(kg.neighbor_graph.__doc__)


# In[18]:


ds_gene_enzyme = kg.neighbor_graph(gene_enzyme,
                                   ds_genes,
                                   keep_isolated_nodes = True,
                                   name = "ds_gene_enzyme")


# the function returns a subgraph of `gene_enzyme` selecting both `ds_genes` nodes and their neighbors in the original graph.
# 
# Not all `ds_genes` are actually present in `gene_enzyme`, we can represent them as isolated nodes in the new graph using the `keep_isolated_nodes` option

# Graphic functionality is provided by `kg.draw()`

# In[19]:


print(kg.draw.__doc__)


# In[20]:


kg.draw(ds_gene_enzyme, layout = "random_layout")


# we see that the actual common enzymes between those obtained from `dis_gen` and those we found on `gene_enzyme` are very few

# In[21]:


set(ds_genes.keys()) & set(gene_enzyme.nodes)


# Le'ts download the *enzyme-reaction* graph from **KEGG** for a last step

# In[22]:


enzyme_reaction = kg.kegg_graph("enzyme", "reaction")


# using the same mechanism as before we get the *enzyme-reaction* graph

# we can get every node of a particular `nodetype` in a graph using `kg.get_nodetype_nodes`:

# In[23]:


print(kg.get_nodetype_nodes.__doc__)


# In[24]:


ds_enzymes = kg.get_nodetype_nodes(ds_gene_enzyme, "enzyme")


# In[25]:


list(ds_enzymes.keys())


# the same `neighbor_graph()` function as before is used to calculated a subgraph of `enzyme_reaction`

# In[26]:


ds_enzyme_reaction = kg.neighbor_graph(enzyme_reaction,
                                       ds_enzymes,
                                       keep_isolated_nodes = True,
                                       name = "ds_enzyme_reaction" )


# In[27]:


kg.draw(ds_enzyme_reaction, layout = "kamada_kawai_layout")


# we can count all the reactions in the graph

# In[28]:


ds_reactions = kg.get_nodetype_nodes(ds_enzyme_reaction, "reaction")
list(ds_reactions.keys())


# we want now to create a projected graph in which we link two enzymes if they appear in the same reaction: we can use the `projected_graph()`function that does exactly that

# In[29]:


print(kg.projected_graph.__doc__)


# here we want to project the graph on the enzymes set that we obtain from `ds_gene_enzyme` using `get_nodetype_nodes()`

# In[30]:


ds_enzymes = kg.get_nodetype_nodes(ds_gene_enzyme, "enzyme")


# In[31]:


list(ds_enzymes.keys())


# now we can finally project `ds_enzyme_reaction` onto the `ds_enzymes` set to get our enzyme correlation graph

# In[32]:


ds_enzyme_enzyme = kg.projected_graph(ds_enzyme_reaction, ds_enzymes, name = "ds_enzyme_enzyme")


# In[33]:


kg.draw(ds_enzyme_enzyme)


# In[3]:


get_ipython().system('jupyter nbconvert --to=script "Tutorial - EnzymeGraphs.ipnyb"')


# In[ ]:




