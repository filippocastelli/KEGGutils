import networkx as nx
import matplotlib.pylab as plt
import matplotlib.image as mpimg
import requests
import os, glob
import json
import imghdr
import shutil
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
    def __init__(self, graph, msg=None):
        self.graph = graph

        if msg is None:
            msg = "Graph {} has no projected graph".format(graph.name)

        super(NoProjectedError, self).__init__(graph, msg)


class KEGGOnlineError(Exception):
    def __init__(self, request, msg=None):
        self.request = request
        self.url = request.url
        self.status_code = request.status_code
        self.reason = request.reason

        if msg is None:
            msg = "Nework Request Error > url: {}, stat. code: {}, description: {}".format(
                self.url, self.status_code, self.reason
            )
        super(KEGGOnlineError, self).__init__(msg)


class KEGGKeyError(Exception):
    def __init__(self, key, msg=None):
        self.key = key

        if msg is None:
            msg = "Invalid KEGG key: {}".format(self.key)
        super(KEGGKeyError, self).__init__(msg)


# =============================================================================
# KEGG API
# =============================================================================
download_dir = "./kegg_downloads/"

db_categories = [
    "pathway",
    "brite",
    "module",
    "ko",
    "genome",
    "vg",
    "ag",
    "compound",
    "glycan",
    "reaction",
    "rclass",
    "enzyme",
    "network",
    "variant",
    "disease",
    "drug",
    "dgroup",
    "environ",
    "atc",
    "jtc",
    "ndc",
    "yj",
    "pubmed",
    "hsa",
]


medicus = [
    "disease_ja",
    "drug_ja",
    "dgroup_ja",
    " environ_ja",
    "compound_ja",
    "brite_ja",
    "atc",
    "jtc",
    "ndc",
    "yj",
]

outside = ["pubmed", "ncbi-geneid", "ncbi-proteinid", "uniprot", "pubchem", "chebi"]

def msg_start_download(filename, url, verbose = True):
    if verbose == True:
        print("> Downloading {} from KEGG at {}".format(filename, url))
        
def msg_end_download(filename, verbose = True):
    if verbose == True:
        print("succesfully downloaded {}".format(filename))
        
def msg_file_already_exists(filename, download_dir, verbose):
    if verbose == True:
        print("> File {} already present in {}".format(filename, download_dir))
        print("reading from cached file...")

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
    print(*files_to_remove, sep="\n")

    for file in files_to_remove:
        os.remove(file)


def download_textfile(url, filename, force_download=False, verbose=True):

    filepath = download_dir + filename
    if (not file_exists(filename)) or force_download:
        msg_start_download(filename, url, verbose)
        
        request = get_online_request(url)
        text = request.text
        mkdir(download_dir)
        
        with open(filepath, "w+") as text_file:
            text_file.write(text)
            
        msg_end_download(filename, verbose)
        
    else:
        msg_file_already_exists(filename, download_dir, verbose)
        with open(filepath, "r") as read_file:
            text = read_file.read()
    return text

def download_json(url, filename, force_download=False, verbose=True):

    filepath = download_dir + filename
    if (not file_exists(filename)) or force_download:
        msg_start_download(filename, url, verbose)
        
        request = get_online_request(url)
        json_file = request.json()
        mkdir(download_dir)
        
        with open(filepath, "w+") as outfile:
            json.dump(json_file, outfile)
            
        msg_end_download(filename, verbose)
        
    else:            
        with open(filepath, "r") as read_file:
            data = json.load(read_file)
            
    return data

def download_pic(url, filename, force_download = False, verbose = False):
    
    filepath = download_dir + filename
    
    possible_filenames = {'gif': filename + '.gif', 
                           'png': filename + '.png'}
    
    possible_filepaths = {'gif': filepath + '.gif',
                           'png': filepath + '.png'}
    
    if all(not file_exists(filenames) for filenames in possible_filenames) or force_download:
        
        msg_start_download(filename, url, verbose)
        response = requests.get(url, stream=True)
        
        filetype = imghdr.what(response.raw)
        
        assert filetype in ["gif", "png"], "Filetype should always be png or gif"
        
        
        with open(possible_filepaths[filetype], 'wb') as out_file:
            shutil.copyfileobj(response.raw, out_file)
            
        img = mpgimg.imread(possible_filepaths[filetype])
        
        msg_end_download(filename, verbose)
        
        return img
        
    else:
        with open(filepath)


# =============================================================================
#  CODE TO FIX NO TIME NOW
# =============================================================================
def get_online_request(url):
    """ Just an errored proxy for requests.get()"""
    request = requests.get(url)
    if request.ok == False:
        raise KEGGOnlineError(request)
    return request

def process_request_text(fulltext, want_descr = False):
    itemlist = []
    descriptionlist = []
    
    for line in fulltext.splitlines():
        entry, description = line.strip().split("\t")
        itemlist.append(entry)
        descriptionlist.append(description)
    if want_descr == True:
        return itemlist, descriptionlist
    else:
        return itemlist
    
    
def push_backslash(stuff):
    stuff_url = ""
    
    if stuff is None:
        stuff = ""
    else:
        stuff_url = "/" + stuff
        
    return stuff, stuff_url
    

def keggapi_list(database, option = None, want_descriptions = False):
    """Returns KEGG list of codes """
    
    org_codes = get_organism_codes()
    if database not in db_categories:
        raise KEGGKeyError(database)
        
    if (option == "xl")&(database != "brite"):
        raise KEGGKeyError(database, msg = "option xl can only be used with brite argument")
        
    if (option in org_codes) & (database != "pathway") & (database != "module"):
        raise KEGGKeyError(database, msg = "only pathway and module list request are available for {}".format(option))
            
    option, optionurl = push_backslash(option)
    
    url = "http://rest.kegg.jp/list/{}{}".format(database, optionurl)
    filename = database + "_" + option + "_list"
    
    list_fulltext = download_textfile(url, filename)
    
    if want_descriptions == False:
        itemlist = process_request_text(list_fulltext, want_descr=want_descriptions)
    if want_descriptions == True:
        itemlist, descriptionlist = process_request_text(list_fulltext, want_descr=want_descriptions)

    if want_descriptions == False:
        return itemlist
    else:
        assert len(itemlist) == len(descriptionlist), "different item and description lengths, something's not working"
        return itemlist, descriptionlist


def keggapi_find(database, query, option = None, want_descriptions = False, verbose = False):
    
    options = ["formula", "exact_mass", "mol_weight"]
    
    if database not in db_categories:
        raise KEGGKeyError(database)
        
    if (database == "compound" or database == "drug") & (option is not None) & (option not in options):
        raise KEGGKeyError(database, msg = "only {} opts are available for {} database".format(options, database))
        
    query, queryurl = push_backslash(query)
    option, optionurl = push_backslash(option)
    
    url = "http://rest.kegg.jp/find/{}{}{}".format(database, queryurl, optionurl)
    
    filename = database + "_" + query + "_" + option
    fulltext = download_textfile(url, filename, verbose = verbose)
    if want_descriptions == False:
        itemlist = process_request_text(fulltext, want_descr=want_descriptions)
    if want_descriptions == True:
        itemlist, descriptionlist = process_request_text(fulltext, want_descr=want_descriptions)

    if want_descriptions == False:
        return itemlist
    else:
        assert len(itemlist) == len(descriptionlist), "different item and description lengths, something's not working"
        return itemlist, descriptionlist
    
    
def keggapi_get(dbentry, option = None, want_descriptions = False, verbose = True):
    
    options = ["aase","ntseq", "mol", "kcf","image","conf", "kgml","json"]
    
    if (option is not None) & (option not in options):
        raise KEGGKeyError(option, msg = "option {} invalid for GET".format(option))

    option, optionurl = push_backslash(option)
    
    url = "http://rest.kegg.jp/get/{}{}".format(dbentry, optionurl)
    
    if option == None:
        option = "description"

    filename = dbentry + "_" + option
    
    if option == "description":
        infos = download_textfile(url, filename, verbose = False)
        if verbose == False:
            infos = "\n".join(infos.splitlines()[1:4])
        print("Infos on {} from KEGG:\n".format(dbentry))
        print(infos)
        return infos
    elif option == "json":
        json_data = download_json(url, filename, verbose = verbose)
        return json_data
    elif option in ["mol", "kcf", "conf", "aaseq", "ntseq"]:
        text = download_textfile(url, filename, verbose = verbose)
        print(text)
        return text
    elif option == "image":
        response = requests.get(url, stream=True)
        
        with open('img.png', 'wb') as out_file:
            shutil.copyfileobj(response.raw, out_file)
            del response
        

    
    if want_descriptions == False:
        itemlist = process_request_text(fulltext, want_descr=want_descriptions)
    if want_descriptions == True:
        itemlist, descriptionlist = process_request_text(fulltext, want_descr=want_descriptions)

    if want_descriptions == False:
        return itemlist
    else:
        assert len(itemlist) == len(descriptionlist), "different item and description lengths, something's not working"
        return itemlist, descriptionlist


def get_organism_codes():
    """Returns all KEGG Organism name codes """

    org_url = "http://rest.kegg.jp/list/organism"
    org_filename = "organism_code_list"

    org_codes = []

    organism_fulltext = download_textfile(org_url, org_filename, verbose=False)

    for line in organism_fulltext.splitlines():
        T_identifier, kegg_code, description, hier = line.strip().split("\t")
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

    #    organism_names = get_organism_codes()

    if not target_db in db_categories:
        raise KEGGKeyError(target_db)
    if not source_db in db_categories:
        raise KEGGKeyError(source_db)

    if target_db == source_db:
        raise KEGGKeyError(
            source_db, "Same key for target and source: {}".format(source_db)
        )

    #    assert all(key in db_categories+organism_names for key in [target_db, source_db]), "Invalid target or source KEGG database key"

    url = "http://rest.kegg.jp/link/" + target_db + "/" + source_db

    return url


def get_infos(item, verbose=False):
    """ Prints KEGG infos for a given database item 
    Parameters:
        :item (str): KEGG item you want infos about
        :verbose (Bool), False: if True get full KEGG description, if False get only first 4 lines
        """

    url = "http://rest.kegg.jp/get/" + item
    filename = item + "_description"

    infos = download_textfile(url, filename, verbose=False)
    if verbose == False:
        infos = "\n".join(infos.splitlines()[1:4])

    print("Infos on {} from KEGG:\n".format(item))
    print(infos)


def kegg_graph(source_db, target_db, force_download=False):
    """Returns a NetworkX Graph for a KEGG target-source database
    
    Parameters:
        :target_db (str): target category
        :source_db (str): source category
        :force_download (bool): if set to True redownloads the graph from KEGG\
        everytime, if set to False checks if the graph was previously downloaded\
        and loads an offline copy
        
        both categories must be valid KEGG categories, see KEGG API Docs
        
    Returns:
        :Graph (NetworkX Graph): Graph with nodes of type target_db and source_db
    Example:

        >>> KEGGgraph = kegg_graph("hsa", "disease")

        .. note:: gene category is represented with the corresponding KEGG organism code
        """
    graphname = "{}_to_{}".format(source_db, target_db)

    url = kegg_url(target_db, source_db)
    text = download_textfile(url, graphname)

    nodes1 = []
    nodes2 = []
    for line in text.splitlines():
        n1, n2 = line.strip().split("\t")
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
        graph.add_node(nodo, nodetype=nodetype1)
        graph.add_node(nodes_2[i], nodetype=nodetype2)
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
    attributes = nx.get_node_attributes(graph, "nodetype")
    if attributes == {}:
        return False
    else:
        return True


def get_nodetype_nodes(kegg_graph, nodetype):
    """Given a KEGG graph returns all the nodes for a given nodetype
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via kegg_graph()
        :nodetype (str): nodetype, is generally a <database> KEGG name
        
    Returns:
        :nodedict (dict): dict of nodes and corresponding nodetype
    Example:
        >>> KEGG_graph = kegg_graph("hsa", "disease")
        >>> nodedict = get_nodetype_nodes(KEGG_graph, "hsa")
        >>> list(nodedict.items())[:5]
        [('hsa:7428', 'hsa'),
         ('hsa:4233', 'hsa'),
         ('hsa:2271', 'hsa'),
         ('hsa:201163', 'hsa'),
         ('hsa:7030', 'hsa')]

    .. seealso:: kegg_graph()
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
        :kegg_graph (Graph): input graph, has to be generated via kegg_graph()
        
    Returns:
        :nodetypes (list): list of unique nodetypes
    Example:
        >>> KEGG_graph = kegg_graph("hsa", "disease")
        >>> nlist = get_unique_nodetypes(KEGG_graph)
        ['disease','hsa']

    .. seealso:: kegg_graph()
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
        :graph (Graph): input graph, has to be generated via kegg_graph()
        :node (str): name of a node in graph
        
    Returns:
        :linked_nodes (dict): dict of linked nodes { node: nodetype}

    .. seealso:: kegg_graph()
        """

    linked_nodes = list(graph[node])
    attributes = nx.get_node_attributes(graph, "nodetype")

    linked_nodes_dict = dict(((k, attributes[k]) for k in linked_nodes))

    return linked_nodes_dict


def neighbor_graph(graph, node_dict, name=None, keep_isolated_nodes=False):
    """Neighbor Subgraph
    
    Given a Graph and a node list returns the subgraph generated with the nodes
    in the node list, the first neighbors of those nodes, and the edges between
    them
    
    
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via kegg_graph()
        :nodelist (list): list of nodes for the nighbor graph
        :name (str): optional, name of the graph
        
    Returns:
        :neighbor_graph (Graph): graph of nodelist, first neighbors of those nodes\
        and edges between them
    .. seealso:: kegg_graph()
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
            neighbor_graph.add_node(node, nodetype=node_dict[node])

    if name is None:
        name = "neighbor_graph_of_{}".format(graph.name)

    neighbor_graph.name = name

    return neighbor_graph


def projected_graph(graph, nodelist, multigraph=False, name=None):
    """Calculates the projected graph respect to a node list     
    Parameters:
        :kegg_graph (Graph): input graph, has to be generated via kegg_graph()
        :nodelist (list): list of nodes
        :multigraph (bool): if True 
        :name (str): optional name of the graph
        
    Returns:
        :descendant_graph (Graph): graph of descendant nodes
    .. seealso:: kegg_graph()
    """

    graphnodes_set = set(graph.nodes)
    nodelist_set = set(nodelist.keys())

    common_nodes = graphnodes_set & nodelist_set

    try:
        nodetype = graph.nodes[list(common_nodes)[0]]["nodetype"]
    except IndexError:
        raise NoProjectedError(graph)

    disjoint_nodes = nodelist_set - set(get_nodetype_nodes(graph, nodetype))

    projected_graph = nx.Graph.copy(nx.projected_graph(graph, common_nodes, multigraph))

    for dis_node in disjoint_nodes:

        projected_graph.add_node(dis_node, nodetype=nodetype)

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
def draw(graph, title=None, layout=None, filename=None, return_ax=False):
    """Graph drawing made a bit easier
    
    Parameters:
        :graph (Graph): input graph, has to be generated via kegg_graph()
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

    graph_nodetypes = get_unique_nodetypes(graph)

    nodetypes_dict = nx.get_node_attributes(graph, "nodetype")

    if len(graph_nodetypes) == 1:
        graph_nodetypes = graph_nodetypes * 2

    graph_colors = replace_dict_value(nodetypes_dict, graph_nodetypes[0], "b")
    graph_colors = replace_dict_value(graph_colors, graph_nodetypes[1], "r")

    if title is None:
        title = "{} > {} graph".format(graph_nodetypes[1], graph_nodetypes[0])

    layouts = {
        "circular_layout": nx.circular_layout,
        "kamada_kawai_layout": nx.kamada_kawai_layout,
        "random_layout": nx.random_layout,
        "shell_layout": nx.shell_layout,
        "spring_layout": nx.spring_layout,
        "spectral_layout": nx.spectral_layout,
    }

    if layout not in layouts:
        print("layout {} not valid: using {} layout".format(layout, default_layout))
        layout = default_layout

    plt.figure()

    pos = layouts[layout](graph)

    nx.draw_networkx_nodes(graph, pos, node_color=list(graph_colors.values()))
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_labels(graph, pos)

    if title is not None:
        plt.title(title)

    plt.axis("off")

    if filename is not None:
        plt.savefig("output.png")

    plt.show()

    if return_ax:
        ax = plt.gca()

        return ax


# =============================================================================
# MISC
# =============================================================================


def replace_dict_value(dictionary, old_value, new_value):
    """ Selectively replaces values in a dictionary
    
    Parameters:
        :dictionary(dict): input dictionary
        :old_value: value to be replaced
        :new_value: value to replace
        
    Returns:
        :output_dictionary (dict): dictionary with replaced values"""

    for key, value in dictionary.items():
        if value == old_value:
            dictionary[key] = new_value
    return dictionary
