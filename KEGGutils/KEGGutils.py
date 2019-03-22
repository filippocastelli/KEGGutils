import networkx as nx
import matplotlib.pylab as plt
from matplotlib import colors as mplcolors
import logging 
logging.basicConfig(level=logging.WARNING)

from KEGGutils.KEGGerrors import MissingNodetypeError,NotAKeggGraphError, NoProjectedError
from KEGGutils.KEGGhelpers import replace_dict_value, shift_pos, shorten_labels
from KEGGutils.KEGGapi import keggapi_link

# =============================================================================
# GRAPH OPERATIONS
# =============================================================================


def kegg_link_graph(source, target, force_download = False):
    """Returns a NetworkX bipartite link graph with nodes from source and target KEGG databases 
    
    Parameters
    ----------
    source : str
        source database
    target : str
        target database
    force_download : bool, optional
        if set to True overwrites pre-existing database file with the same name (the default is False, which [default_description])
    
    Returns
    -------
    graph
        bipartite link graph
    """


    graphname = "{}_to_{}".format(source, target)
    
    nodes1, nodes2 = keggapi_link(source, target, verbose = True)
    
    graph = nx.Graph()
    graph.name = graphname
    populate_graph(graph, nodes1, nodes2, source, target)
    
    return graph


def populate_graph(graph, nodes_1, nodes_2, nodetype1, nodetype2):
    """Populates a pre-existing Graph given two list of nodes and nodetypes
    
    Parameters
    ----------
    graph : graph
        graph to populate
    nodes_1 : list
        list of nodes to update
    nodes_2 : list
        list of nodes to update
    nodetype1 : str
        first nodetype
    nodetype2 : str
        second nodetype
    
    """
    for i, nodo in enumerate(nodes_1):
        graph.add_node(nodo, nodetype=nodetype1, label = nodo)
        graph.add_node(nodes_2[i], nodetype=nodetype2, label = nodes_2[i])
        graph.add_edge(nodo, nodes_2[i])

def has_nodetypes(graph):
    """Determines if graph nodes have the "nodetype" attribute
    
    Parameters:
        :graph (Graph): input graph
        
    Returns:
        :has_nodetype (bool): that's self explanatory bro
        
    """
    attributes = nx.get_node_attributes(graph, "nodetype")
    if attributes == {}:
        return False
    else:
        return True


def get_nodes_by_nodetype(kegg_graph, nodetype):
    """Given a KEGG graph returns all the nodes for a given nodetype
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via kegg_link_graph()
        :nodetype (str): nodetype, is generally a <database> KEGG name
        
    Returns:
        :nodedict (dict): dict of nodes and corresponding nodetype
    Example:
        >>> KEGG_graph = kegg_link_graph("hsa", "disease")
        >>> nodedict = get_nodes_by_nodetype(KEGG_graph, "hsa")
        >>> list(nodedict.items())[:5]
        [('hsa:7428', 'hsa'),
         ('hsa:4233', 'hsa'),
         ('hsa:2271', 'hsa'),
         ('hsa:201163', 'hsa'),
         ('hsa:7030', 'hsa')]

    .. seealso:: kegg_link_graph()
        """

    if nodetype not in get_unique_nodetypes(kegg_graph):
        raise MissingNodetypeError(
            kegg_graph, "Requested nodetype is missing in graph nodes"
        )

    node_list = [n for n in kegg_graph if kegg_graph.node[n]["nodetype"] == nodetype]

    return dict.fromkeys(node_list, nodetype)


def connected_components(graph):
    """ Returns a list of connected components for a given graph, ordered by size"""
    subgraphs = list(nx.connected_component_subgraphs(graph))
    subgraphs = sorted(subgraphs, key=len, reverse=True)

    return subgraphs


def get_unique_nodetypes(graph):
    """Given a KEGG graph returns list of its unique nodetypes
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via kegg_link_graph()
        
    Returns:
        :nodetypes (list): list of unique nodetypes
    Example:
        >>> KEGG_graph = kegg_link_graph("hsa", "disease")
        >>> nlist = get_unique_nodetypes(KEGG_graph)
        ['disease','hsa']

    .. seealso:: kegg_link_graph()
        """
    if has_nodetypes(graph) == False:
        raise NotAKeggGraphError(graph, "Graph nodes are missing nodetype attribute")

    attributes = nx.get_node_attributes(graph, "nodetype")
    unique_nodetypes = sorted(set(attributes.values()))

    return unique_nodetypes


def linked_nodes(graph, node):
    """Linked Nodes:
        Returns all nodes in graph linked to node
    
    Parameters:
        :graph (Graph): input graph, has to be generated via kegg_link_graph()
        :node (str): name of a node in graph
        
    Returns:
        :linked_nodes (dict): dict of linked nodes { node: nodetype}

    .. seealso:: kegg_link_graph()
        """

    linked_nodes = list(graph[node])
    attributes = nx.get_node_attributes(graph, "nodetype")

    linked_nodes_dict = dict(((k, attributes[k]) for k in linked_nodes))

    return linked_nodes_dict


def neighbor_graph(graph, node_dict, name=None, keep_isolated_nodes=False):
    """Neighbor Subgraph
    
    Given a Graph and a node list returns the subgraph generated with the nodes
    in the node dict, the first neighbors of those nodes, and the edges between
    them
    
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via kegg_link_graph()
        :node_dict (dict): dict of input nodes
        :name (str): optional, name of the graph
        
    Returns:
        :neighbor_graph (Graph): graph of nodelist, first neighbors of those nodes\
        and edges between them
    .. seealso:: kegg_link_graph()
    """

    nodeset = set()

    input_nodes_set = set(node_dict.keys())

    try:
        nx.get_node_attributes(graph, "nodetype")
    except:
        NotAKeggGraphError(graph)

    nodeset.update(input_nodes_set)

    for node in input_nodes_set:
        try:
            nodeset.update(list(graph[node]))
        except KeyError:
            pass

    # we already have attributes for nodes in graph
    neighbor_graph = nx.Graph.copy(nx.subgraph(graph, nodeset))

    # ok now we want to add input nodes that are left behind
    if keep_isolated_nodes == True:
        difference_set = input_nodes_set - set(neighbor_graph.nodes)

        for node in difference_set:
            neighbor_graph.add_node(node, nodetype=node_dict[node], label = node)

    if name is None:
        name = "neighbor_graph_of_{}".format(graph.name)

    neighbor_graph.name = name

    return neighbor_graph


def projected_graph(graph, nodelist, multigraph=False, name=None):
    """Calculates the projected graph respect to a node list     
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via kegg_link_graph()
        :nodelist (list): list of nodes
        :multigraph (bool): if True 
        :name (str): optional name of the graph
        
    Returns:
        :descendant_graph (Graph): graph of descendant nodes
    .. seealso:: kegg_link_graph()
    """

    graphnodes_set = set(graph.nodes)
    nodelist_set = set(nodelist.keys())

    common_nodes = graphnodes_set & nodelist_set

    try:
        nodetype = graph.nodes[list(common_nodes)[0]]["nodetype"]
    except IndexError:
        raise NoProjectedError(graph)

    disjoint_nodes = nodelist_set - set(get_nodes_by_nodetype(graph, nodetype))

    projected_graph = nx.Graph.copy(nx.projected_graph(graph, common_nodes, multigraph))

    for dis_node in disjoint_nodes:

        projected_graph.add_node(dis_node, nodetype=nodetype, label = dis_node)

    if name == None:
        name = "{}_projected".format(graph.name)
    projected_graph.name = name

    return projected_graph


def graph_measures(graph):

    max_connected_component = connected_components(graph)[0]

    measures = {
        "totnodes": len(max_connected_component.nodes),
        "totedges": len(max_connected_component.edges),
        "ncliques": len(list(nx.enumerate_all_cliques(max_connected_component))),
        "radius": nx.radius(max_connected_component),
        "diameter": nx.diameter(max_connected_component),
    }

    return measures


# =============================================================================
# PLOTTING
# =============================================================================
def draw(graph, title=None, layout=None, filename=None, return_ax=False, pos = None,
         font_size = 9,  alpha = 1.0, label_shift = (0,0), truncate_labels = 10):
    """Graph drawing made a bit easier
    
    Parameters:
        :graph (Graph): input graph, has to be generated via kegg_link_graph()
        :layout (str): layout type, choose from 'bipartite_layout',\
        'circular_layout','kamada_kawai_layout','random_layout',\ 'shell_layout',\
        'spring_layout','spectral_layout'
        :filename (str): if a filename is selected saves the plot as filename.png
        :title (str): title for the graph
        :return_ax: if True returns ax for plot
        
    Returns:
        :ax (list): optional ax for the plot


        """
    default_layout = "spring_layout"
    if layout is None:
        layout = default_layout
        
    node_groups = {}
    
    graph_nodetypes = get_unique_nodetypes(graph)
    
    base_colors = list(mplcolors.BASE_COLORS.keys())
    
    for i, nodetype in enumerate(graph_nodetypes):
        node_group = (get_nodes_by_nodetype(graph, nodetype).keys())
        node_groups.update({nodetype: (node_group, base_colors[i])})
        
    if title is None:
        if len(graph_nodetypes) == 1:
            title = "{} graph".format(graph_nodetypes[0])
        if len(graph_nodetypes) == 2:
            title = "{} > {} graph".format(graph_nodetypes[1], graph_nodetypes[0])
        else:
            title = "Graph plot"

    layouts = {
        "circular_layout": nx.circular_layout,
        "kamada_kawai_layout": nx.kamada_kawai_layout,
        "random_layout": nx.random_layout,
        "shell_layout": nx.shell_layout,
        "spring_layout": nx.spring_layout,
        "spectral_layout": nx.spectral_layout,
    }

    if layout not in layouts:
        logging.warning("layout {} not valid: using {} layout\nusing default layout".format(layout, default_layout))
        layout = default_layout

    plt.figure()    
    
    if pos is None:
        output_layout = layouts[layout](graph)
        pos = {}
        for key, value in output_layout.items():
            pos[key] = tuple(value)
            

    for nodetype, node_group in node_groups.items():
        nx.draw_networkx(graph, nodelist = node_group[0], pos = pos, node_color = node_group[1], with_labels = False, label = nodetype)
        
    nx.draw_networkx_edges(graph, pos)
    pos_labels = shift_pos(pos, label_shift)
    
    candidate_labels = nx.get_node_attributes(graph, "label")
    
    if candidate_labels != {}:
        if truncate_labels != False:
            labels = shorten_labels(candidate_labels, truncate_labels)
        else:
            labels = candidate_labels
    else:
#        labels = None
        nodelist = list(graph.nodes)
        labels = dict(zip(nodelist, nodelist))
    
    nx.draw_networkx_labels(graph, pos_labels,
                            labels = labels,
                            font_size = font_size,
                            alpha = alpha)
    
    plt.legend()
    if title is not None:
        plt.title(title)

    plt.axis("off")

    if filename is not None:
        plt.savefig("output.png")

    plt.show()

    if return_ax:
        ax = plt.gca()

        return ax

