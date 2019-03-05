import networkx as nx
import requests

#%%
def populate_graph(graph, nodes_1, nodes_2, nodetype1, nodetype2):
    """Populates a pre-existing Graph given two list of nodes and two node labels
    
    Parameters:
        :graph (Graph): input graph
        :nodes_1(list): first list of nodes
        :nodes_2(list): second list of nodes
        :nodetype1(str): first nodetype
        :nodetype2(str): second nodetype
        """
    
    for i, nodo in enumerate(nodes_1):
        graph.add_node(nodo, nodetype = nodetype1)
        graph.add_node(nodes_2[i], nodetype = nodetype2)
        graph.add_edge(nodo, nodes_2[i])

def kegg_graph(target_db, source_db):
    """Returns a NetworkX Graph for a KEGG target-source database
    
    Parameters:
        :target_db (str): target category
        :source_db (str): source category
        
        both categories must be valid KEGG categories, see KEGG API Docs
        
    Returns:
        :Graph (NetworkX Graph): Graph with nodes of type target_db and source_db
    Example:

        >>> KEGGgraph = gen_graph("hsa", "disease")

        .. note:: gene category is represented with the corresponding KEGG organism code
        """
    
    url = kegg_url(target_db, source_db)
    print("Downloading database...")
    text = requests.get(url).text
    print("done")
    
    nodes1 = []
    nodes2 = []
    for line in text.splitlines():
        n1, n2 = line.strip().split('\t')
        nodes1.append(n1)
        nodes2.append(n2)

    graph = nx.Graph()
    populate_graph(graph, nodes1, nodes2, source_db, target_db)
    
    return graph
    
def get_organism_codes():
    """Returns all KEGG Organism name codes """
    
    org_url = "http://rest.kegg.jp/list/organism"
    
    org_codes = []
    print("Downloading kegg organism list...")
    organism_fulltext = requests.get(org_url).text
    print("done")
    
    for line in organism_fulltext.splitlines():
        T_identifier, kegg_code, description, hier= line.strip().split('\t')
        org_codes.append(kegg_code)
        
    return org_codes

def kegg_url(target_db, source_db):
    """Returns a KEGG database URL given two categories
    
    Parameters:
        :target_db (str): target category
        :source_db (str): source category
        
        both categories must be valid KEGG categories, see KEGG API Docs
        
    Returns:
        :url (str): url for the corresponding KEGG database
    Example:

        >>> kegg_url("hsa", "disease")
        'http://rest.kegg.jp/link/hsa/disease'

    .. warning:: 
        - gene category is represented with the corresponding KEGG organism code
    
        - target_db and source_db must be valid KEGG <database> names, or valid <org> names, see KEGG API Docs
        """
    db_categories = ["pathway", "brite", "module",
                     "ko", "genome", "vg",
                     "ag", "compound", "glycan", "reaction",
                     "rclass", "enzyme", "network", "variant",
                     "disease", "drug", "dgroup", "environ", "atc",
                     "jtc", "ndc", "yj", "pubmed"]
    
    #don't wanna download list it every time, too heavy
    #it's time for ugly global variables
    #just declaring doesnt actually create it
    global organism_names
    
    try:
        organism_names
    except NameError:
        organism_names = get_organism_codes()
    
    assert target_db != source_db, "Same key for target and source"
    assert all(key in db_categories+organism_names for key in [target_db, source_db]), "Invalid target or source KEGG database key"
            
    url = "http://rest.kegg.jp/link/"+target_db+"/"+source_db
    
    return url

def get_nodetype_nodes(kegg_graph, nodetype):
    """Given a KEGG graph returns list of nodes for a given nodetype
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via gen_graph()
        :nodetype (str): nodetype, is generally a <database> KEGG name
        
    Returns:
        :nodelist (list): list of nodes
    Example:
        >>> KEGG_graph = gen_graph("hsa", "disease")
        >>> nlist = nodelist(KEGG_graph, "hsa")
        >>> nlist[:10]
        ['hsa:7428',
         'hsa:4233',
         'hsa:2271',
         'hsa:201163',
         'hsa:7030',
         'hsa:894',
         'hsa:411',
         'hsa:1075',
         'hsa:2720',
         'hsa:2588']

    .. seealso:: gen_graph()
        """
        
    return [n for n in kegg_graph if kegg_graph.node[n]['nodetype'] == nodetype]

def connected_components(graph):
    """ Returns a list of connected components for a given graph, ordered by size"""
    subgraphs = list(nx.connected_component_subgraphs(graph))
    subgraphs = sorted(subgraphs, key=len, reverse=True)
    
    return subgraphs

def kegg_infos(item, verbose = False):
    """ Prints KEGG infos for a given database item 
    Parameters:
        :item (str): KEGG item you want infos about
        :verbose (Bool), False: if True get full KEGG description, if False get only first 4 lines
            Example:
        """
    
    url = "http://rest.kegg.jp/get/"+item
    infos = requests.get(url).text
    if verbose == False:
        infos = "\n".join(infos.splitlines()[1:4])
        
    print(infos)
        

