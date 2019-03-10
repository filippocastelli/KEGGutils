import networkx as nx
import matplotlib.pylab as plt
import requests
import os, glob

    #%%
# =============================================================================
# ERRORS AND EXCEPTIONS
# =============================================================================


class KeggUtilsGraphException(Exception):
    def __init__(self, graph, msg=None):
        self.graph = graph
        if msg is None:
            # Set some default useful error message
            msg = "There's a problem with graph: {}".format(graph.name)
            
        super(KeggUtilsGraphException, self).__init__(msg)

        
class NotAKeggGraphError(KeggUtilsGraphException):
    pass
        
class MissingNodetypeError(KeggUtilsGraphException):
    pass

class NoProjectedError(KeggUtilsGraphException):
    def __init__(self, graph, msg = None):
        self.graph = graph
        
        if msg is None:
            msg = "Graph {} has no projected graph".format(graph.name)
            
        super(NoProjectedError, self).__init__(graph, msg)

class KEGGOnlineError(Exception):
    def __init__(self, request, msg = None):
        self.request = request
        self.url = request.url
        self.status_code = request.status_code
        self.reason = request.reason
        
        if msg is None:
            msg = "Nework Request Error > url: {}, stat. code: {}, description: {}".format(self.url, self.status_code, self.reason)
        super(KEGGOnlineError, self).__init__(msg)
        
        
# =============================================================================
# DOWNLOADING FROM KEGG
# =============================================================================
download_dir = "./kegg_downloads/"    


def file_exists(filename):
    return os.path.isfile(download_dir + filename)

def mkdir(directory):
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass
    
def delete_cached_files():
    files_to_remove = glob.glob(download_dir + "*")
    
    print("> deleting the following files from {}".format(download_dir))
    print(*files_to_remove, sep='\n')
    
    for file in files_to_remove:
        os.remove(file)
    
   
def download_textfile(url, filename,
                      force_download = False,
                      verbose = True):
    
    start_message = "> Downloading {} from KEGG at {}".format(filename, url)
    end_message = "succesfully downloaded {}".format(filename)
    filepath = download_dir + filename
    if (not file_exists(filename)) or force_download:
        if verbose:
            print(start_message)
        request = get_online_request(url)
        text = request.text
        mkdir(download_dir)
        with open(filepath, 'w+') as text_file:
            text_file.write(text)
        if verbose == True:
            print(end_message)
#            print("\n")
    else:
        if verbose:
            print("> File {} already present in {}".format(filename, url))
            print("reading from cached file...")
#            print("\n")
        
        with open(filepath, 'r') as read_file:
            text =  read_file.read()
    return text
            
def get_online_request(url):
    """ Just an errored proxy for requests.get()"""
    request = requests.get(url) 
    if request.ok == False:
        raise KEGGOnlineError(request)   
    return request

def get_list(item):
    """Returns KEGG list of codes """
    
    url = "http://rest.kegg.jp/list/{}".format(item)
    filename = item+"_list"
    
    itemlist = []
    
    list_fulltext = download_textfile(url, filename)
    for line in list_fulltext.splitlines():
        entry, description = line.strip().split('\t')
        itemlist.append(entry)
        
    return itemlist
    
def get_organism_codes():
    """Returns all KEGG Organism name codes """
    
    org_url = "http://rest.kegg.jp/list/organism"
    org_filename = "organism_code_list"
    
    org_codes = []
    
    organism_fulltext = download_textfile(org_url, org_filename, verbose = False)
    
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
    


    organism_names = get_organism_codes()
    
    assert target_db != source_db, "Same key for target and source"
    assert all(key in db_categories+organism_names for key in [target_db, source_db]), "Invalid target or source KEGG database key"
            
    url = "http://rest.kegg.jp/link/"+target_db+"/"+source_db
    
    return url   

def get_infos(item, verbose = False):
    """ Prints KEGG infos for a given database item 
    Parameters:
        :item (str): KEGG item you want infos about
        :verbose (Bool), False: if True get full KEGG description, if False get only first 4 lines
            Example:
        """
    
    url = "http://rest.kegg.jp/get/"+item
    filename = item+"_description"
    
    infos = download_textfile(url, filename, verbose = False)
    if verbose == False:
        infos = "\n".join(infos.splitlines()[1:4])
        
    print("Infos on {} from KEGG:\n".format(item))
    print(infos)
    
def kegg_graph(source_db, target_db):
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
    graphname = "{}_to_{}".format(source_db, target_db)
    
    url = kegg_url(target_db, source_db)
    text = download_textfile(url, graphname)
    
    nodes1 = []
    nodes2 = []
    for line in text.splitlines():
        n1, n2 = line.strip().split('\t')
        nodes1.append(n1)
        nodes2.append(n2)

    graph = nx.Graph()
    graph.name = graphname
    populate_graph(graph, nodes1, nodes2, source_db, target_db)
    
    
    
    return graph

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
        

# =============================================================================
# GRAPH OPERATIONS
# =============================================================================

def has_nodetypes(graph):
    """Populates a pre-existing Graph given two list of nodes and two node labels
    
    Parameters:
        :graph (Graph): input graph
        
    Returns:
        :has_nodetype (bool): that's self explanatory bro
        
    """
    attributes = nx.get_node_attributes(graph, 'nodetype')
    if attributes == {}:
        return False
    else:
        return True

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
        
    if nodetype not in get_unique_nodetypes(kegg_graph):
        raise MissingNodetypeError(kegg_graph, "Requested nodetype is missing in graph nodes")
        
    node_list = [n for n in kegg_graph if kegg_graph.node[n]['nodetype'] == nodetype]
    
    return dict.fromkeys(node_list, nodetype)

def connected_components(graph):
    """ Returns a list of connected components for a given graph, ordered by size"""
    subgraphs = list(nx.connected_component_subgraphs(graph))
    subgraphs = sorted(subgraphs, key=len, reverse=True)
    
    return subgraphs


    
def get_unique_nodetypes(graph):
    """Given a KEGG graph returns list of its unique nodetypes
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via gen_graph()
        
    Returns:
        :nodetypes (list): list of unique nodetypes
    Example:
        >>> KEGG_graph = gen_graph("hsa", "disease")
        >>> nlist = get_unique_nodetypes(KEGG_graph)
        ['disease','hsa']

    .. seealso:: gen_graph()
        """
    if has_nodetypes(graph) == False:
        raise NotAKeggGraphError(graph, "Graph nodes are missing nodetype attribute")
        
    attributes = nx.get_node_attributes(graph, 'nodetype')
    unique_nodetypes = sorted(set(attributes.values()))
    
    return unique_nodetypes

def linked_nodes(graph, node):
    
    linked_nodes = list(graph[node])
    attributes = nx.get_node_attributes(graph, 'nodetype')

    linked_nodes_dict = dict(((k, attributes[k]) for k in linked_nodes))
    
    return linked_nodes_dict

def neighbor_graph(graph, node_dict, name = None, keep_isolated_nodes = False):
    """Neighbor Subgraph
    
    Given a Graph and a node list returns the subgraph generated with the nodes
    in the node list, the first neighbors of those nodes, and the edges between
    them
    
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via gen_graph()
        :nodelist (list): list of nodes for the nighbor graph
        :name (str): optional, name of the graph
        
    Returns:
        :neighbor_graph (Graph): graph of nodelist, first neighbors of those nodes\
        and edges between them
    .. seealso:: gen_graph()
    """
    
    
    nodeset = set()
    
    input_nodes_set = set(node_dict.keys())
    
    
    try:
        nx.get_node_attributes(graph, 'nodetype')
    except:
        NotAKeggGraphError(graph)
    
    nodeset.update(input_nodes_set)
    
    for node in input_nodes_set:
        try:
            nodeset.update(list(graph[node]))
        except KeyError:
            pass
        
    #we already have attributes for nodes in graph
    neighbor_graph = nx.Graph.copy(nx.subgraph(graph, nodeset))
    
    # ok now we want to add input nodes that are left behind
    if keep_isolated_nodes == True:
        difference_set = input_nodes_set - set(neighbor_graph.nodes)
        
        for node in difference_set:
            neighbor_graph.add_node(node, nodetype = node_dict[node])

    if name is None:
        name = "neighbor_graph_of_{}".format(graph.name)
        
    neighbor_graph.name = name
    
    return neighbor_graph

def projected_graph(graph, nodelist, multigraph = False, name = None):
    """Calculates the projected graph respect to a node list     
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via gen_graph()
        :nodelist (list): list of nodes
        :multigraph (bool): if True 
        
    Returns:
        :descendant_graph (Graph): graph of descendant nodes
    .. seealso:: gen_graph()
    """
    
    graphnodes_set = set(graph.nodes)
    nodelist_set = set(nodelist.keys())
    
    common_nodes = graphnodes_set & nodelist_set
    
    try:
        nodetype = graph.nodes[list(common_nodes)[0]]['nodetype']
    except IndexError:
        raise NoProjectedError(graph)
    
    disjoint_nodes = nodelist_set - set(get_nodetype_nodes(graph, nodetype))
    
    projected_graph = nx.Graph.copy(nx.projected_graph(graph, common_nodes, multigraph))
    
    for dis_node in disjoint_nodes:
        
        projected_graph.add_node(dis_node, nodetype = nodetype)
        
    if name == None:
        name = "{}_projected".format(graph.name)
    projected_graph.name = name
    
    return projected_graph

def graph_measures(graph):
    
    max_connected_component = connected_components(graph)[0]
    
    measures = {'totnodes': len(max_connected_component.nodes),
                'totedges': len(max_connected_component.edges),
                'ncliques': len(list(nx.enumerate_all_cliques(max_connected_component))),
                'radius':   nx.radius(max_connected_component),
                'diameter': nx.diameter(max_connected_component)}
    
    return measures

# =============================================================================
# PLOTTING
# =============================================================================
def draw(graph,
         title = None,
         layout = None,
         filename = None,
         return_ax = False):
    """Graph drawing made a bit easier
    
    Parameters:
        :graph (Graph): input graph, has to be generated via gen_graph()
        :layout (str): layout type, choose from 'bipartite_layout',\
        'circular_layout','kamada_kawai_layout','random_layout',\ 'shell_layout',\
        'spring_layout','spectral_layout'
        :filename (str): if a filename is selected saves the plot as filename.png
        :title (str): title for the graph
        :return_ax: if True returns ax for plot
        
    Returns:
        :ax (list): optional ax for the plot


        """
    default_layout = 'spring_layout'
    if layout is None: layout = default_layout
    graph_nodetypes = get_unique_nodetypes(graph)
    
    if len(graph_nodetypes) == 1:
        graph_nodetypes = graph_nodetypes*2
    
    if title is None: title = "{} > {} graph".format(graph_nodetypes[0], graph_nodetypes[1])
    
    layouts = {'circular_layout':         nx.circular_layout,
             'kamada_kawai_layout':     nx.kamada_kawai_layout,
             'random_layout':           nx.random_layout,
             'shell_layout':            nx.shell_layout,
             'spring_layout':           nx.spring_layout,
             'spectral_layout':         nx.spectral_layout
             }

    if layout not in layouts:
        print("layout {} not valid: using {} layout".format(layout, default_layout))
        layout = default_layout

    plt.figure()
    
    pos = layouts[layout](graph)
    nx.draw_networkx_nodes(graph,pos,node_size=700)
    nx.draw_networkx_edges(graph,pos)
    nx.draw_networkx_labels(graph,pos)
    
    if title is not None:
        plt.title(title)
        
    plt.axis('off')
    
    if filename is not None:
        plt.savefig('output.png')
           
    plt.show()
    
    if return_ax:
        ax = plt.gca()
        
        return ax
        
        
        
        

