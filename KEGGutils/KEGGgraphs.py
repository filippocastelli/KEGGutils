import networkx as nx
import logging

from KEGGutils.KEGGutils import get_nodes_by_nodetype, populate_graph, connected_components,\
linked_nodes, graph_measures, projected_graph, get_unique_nodetypes, draw, neighbor_graph
from KEGGutils.KEGGapi import keggapi_link, keggapi_info
from KEGGutils.KEGGerrors import KEGGDataBaseError, KEGGKeyError, KEGGgraphError


#keggapi_info_ = keggapi_info


class KEGGgraph(nx.Graph):
    """Base class from KEGGutils NetworkX compatible Graphs:
    Directly inherits from networkx.Graph()

    Methods
    -------
    list_by_nodetype(nodetype)
        returns a dict of nodes for the given nodetype in the form {key: nodetype}
    
    connected_components()
        returns a list of connected components
    
    linked_nodes(node)
        returns a dict of nodes in the graph connected to a givne node in the form {key: nodetype}

     graph_measures()
        returns a list of graph characteristics
    
    shortest_path(source_node, target_node):
        computes the shortest path between two nodes
    
    compose(othergraph, inplace = False)
        returns graph composition between self and one or more other graphs
    
    neighbor_graph(nodedict, keep_isolated_nodes, inplace = False)
        returns graph of internal neighbors to a given dict of nodes

    get_unique_nodetypes()
        returns unique nodetypes in graph

    draw(layout = None)
        plots the graph
    
    prune_isolated_nodes(inplace = False)
        removes isolated nodes
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
    def _find_arg_and_kick(self, args, to_find):
        """finds an arg in args using to_find as a substring"""

        for arg in args:
            if type(arg) == str:
                if arg.find(to_find) != -1:
                    logging.debug("found {} in args".format(arg))

                    newargs = [argument for argument in args if argument != arg]
                    args = newargs

                    return arg, newargs
        
    def list_by_nodetype(self, nodetype, return_dict = False):
        """Returns a list of nodes given a nodetype"""

        
        if return_dict != True:
            nodelist = get_nodes_by_nodetype(self, nodetype)
            
            #it's ugly af but removes duplicates in a line
            return list(dict.fromkeys(nodelist))
        else:
            return get_nodes_by_nodetype(self, nodetype, return_dict = True)
    
    def connected_components(self):
        """ returns a list of connected components"""
        return connected_components(self)
    
    def linked_nodes(self,node, return_dict = False):
        """ returns a list of linked nodes to node"""
        
        if return_dict == True:
            return linked_nodes(self, node, return_dict = True)
        else:
            return linked_nodes(self, node)
    
    def graph_measures(self):
        """ see KEGGutils.graph_measures()"""
        return graph_measures(self)
    
    def shortest_path(self, source_node, target_node):
        """ Returns the shortest path between two nodes"""
        return nx.shortest_path(self, source_node, target_node)
    
    def compose(self, othergraph, inplace = False):
        """Compose
        Calculates composition between two or mre graphs
        
        Parameters
        ----------
        othergraph : KEGGgraph
            Can be a single KEGGgraph, a list of graphs or a tuple of graphs, contains the other graphs you want to calculate composition respect to
        inplace : bool, optional
            If True the original graph calling the method is substituted by the composition (the default is False)
        
        Returns
        -------
        KEGGgraph
            composed graph
        """

        if type(othergraph) == self.__class__:
            composed_graph = nx.compose(self, othergraph)
        elif type(othergraph) is list or type(othergraph) is tuple :
            if not all((type(x) == self.__class__) for x in othergraph):
                raise KEGGgraphError([graph for graph in othergraph if not (type(graph) == self.__class__)][0])
            composed_graph = nx.compose_all(othergraph)
        
        if inplace ==True: 
            self.add_nodes_from(composed_graph.nodes(data = True)) 
            self.add_edges_from(composed_graph.edges(data = True))
            return self
        else:
            return composed_graph
    
    def neighbor_graph(self, nodelist, keep_isolated_nodes = True, inplace = False):
        """Neighbor Graph
        
        Parameters
        ----------
        nodelist : list
            contains the nodes we wish to calculate the neighbor graph respect to
        keep_isolated_nodes : bool, optional
            Isolated nodes can be automatically removed selecting False (the default is True)
        inplace : bool, optional
            Substitutes the original graph with the neighbor graph (the default is False, which [default_description])
        
        Returns
        -------
        graph
            neighbor graph
        """

        if set(nodelist) & set(self.nodes) == set():
            raise ValueError
            
        #i'll take the first intersection nodetype for all
        intersection_list = [x for x in nodelist if x in self.nodes] 
        nodetype = self.nodes[intersection_list[0]]['nodetype']
        nodedict = dict.fromkeys(nodelist, nodetype)
        
        ngraph = neighbor_graph(self, nodedict, keep_isolated_nodes = keep_isolated_nodes)
        
        if inplace == True:
            name = self.name
            self.clear() #removes nodes, edge, name
            self.add_nodes_from(ngraph.nodes(data = True))
            self.add_edges_from(ngraph.edges(data = True))
            self.name = name
            
            return self
        else:
            return ngraph
        
        
    def get_unique_nodetypes(self):
        """returns a list of  graph's unique nodetypes"""
        return get_unique_nodetypes(self)
        
    def draw(self, layout = None):
        """ plots the graph"""
        draw(self, title = self.name, layout = layout)

    def prune_isolated_nodes(self, inplace = False):
        """ removes isolated nodes"""
        isolates = nx.isolates(self)
        
        if inplace == True:
            self.remove_nodes_from(isolates)
            return self
        
        else:
            gcopy = self.copy()
            gcopy.remove_nodes_from(isolates)
            
            return gcopy
        

class KEGGlinkgraph(KEGGgraph):
    """KEGGlinkgraph class inherits from KEGGgraph
    Can be used to build bipartite graphs from KEGG API LINK functionality 
    
    Parameters
    ----------
    source_db : str
        Source Database
    target_db : str
        Target Database

    Properties
    ----------
    source_nodes : dict
        dict of source nodes {node: nodetype}
    target_nodes : dict
        dict of target nodes {node: nodetype}
    source_liked_db : list
        list of databases connected to source in KEGG
    target_linked_db : list
        list of databases connected to target in KEGG
    """

    
    source_db = None
    target_db = None
    
    source_nodes = {}
    target_nodes = {}
    
    source_linked_db = []
    target_linked_db = []
    
    
    
    def __init__(self, *args, **kwargs):
        
        if "source_db" in kwargs:
            self.source_db = kwargs.pop("source_db")
        if "target_db" in kwargs:
            self.target_db = kwargs.pop("target_db")
            
        super().__init__(*args, **kwargs)
        
        if (self.source_db is not None) and (self.target_db is not None):
            self.graph_init()
        
        
    def graph_init(self):   
        self.name = "{}_to_{}".format(self.source_db, self.target_db)
        
        nodes1, nodes2 = self._get_nodes_from_api(self.source_db, self.target_db)
        
        self.source_nodes = dict.fromkeys(nodes1, self.source_db)
        self.target_nodes = dict.fromkeys(nodes2, self.target_db)
        
        populate_graph(self, nodes1, nodes2, self.source_db, self.target_db)
        
        self.source_linked_db = self.source_infos(return_format = "dict", verbose = False)['linked db']
        self.target_linked_db = self.target_infos(return_format = "dict", verbose = False)['linked db']
        
        
    def source_infos(self,return_format = None, verbose = True):
        """ retrieves infos on the source database"""
        return keggapi_info(self.source_db, return_format = return_format, verbose = verbose)
    def target_infos(self,return_format = None, verbose = True):
        """ retrieves info on the target database"""
        return keggapi_info(self.target_db, return_format = return_format, verbose = verbose)
    
    def neighbor_graph(self, nodelist, keep_isolated_nodes = True, inplace = False):
        
        ngraph =  super(type(self) ,self).neighbor_graph(nodelist, keep_isolated_nodes = keep_isolated_nodes, inplace = False)
        
        ngraph.source_db = self.source_db
        ngraph.target_db = self.target_db
        
        ngraph.source_nodes = dict.fromkeys(ngraph.list_by_nodetype(self.source_db), self.source_db)
        ngraph.target_nodes = dict.fromkeys(ngraph.list_by_nodetype(self.target_db), self.target_db)
        
        if inplace == True:
            name = self.name
            self.clear() #removes nodes, edge, name, shouldnt touch other properties
            self.add_nodes_from(ngraph.nodes(data = True))
            self.add_edges_from(ngraph.edges(data = True))
            self.name = name
            return self
        else:
            return ngraph
        
        
    def projected_graph(self, nodelist = None, name = None):
        """Projects the link graph on a subset of nodes
        
        Parameters
        ----------
        nodelist : list, optional
            list of nodes you wish to project onto, if None automatically selects source nodes
        name : str, optional
            name of the graph
        
        
        Returns
        -------
        graph
            projection graph
        """
 
        assert nx.is_bipartite(self), "Graph is not bipartite, something's wrong with graph construction"
        if nodelist == None:
            nodedict = self.source_nodes
        else:
            if set(nodelist) & set(self.source_nodes.keys()) is set():
                raise KEGGKeyError(key = self.source_db,
                                   msg = "nodelist and target_nodes intersection ( {} nodetype) must not be null".format(self.source_db)) 
                
            intersection_nodes = [x for x in nodelist if x in self.nodes]
            nodetype = self.nodes[intersection_nodes[0]]['nodetype']
            nodedict = dict.fromkeys(nodelist, nodetype)
        

            
        proj_graph = projected_graph(self, nodedict = nodedict,
                                     multigraph = False,
                                     name = "projection of {} onto {}")
        if name == None:
            name = "projection of {} onto {} nodes".format(self.name, self.source_db)
            
        
        projected_kegggraph = KEGGgraph()

        projected_kegggraph.add_nodes_from(proj_graph.nodes(data = True)) 
        projected_kegggraph.add_edges_from(proj_graph.edges(data = True))
        projected_kegggraph.name = name
        
        return projected_kegggraph
        
    def _get_nodes_from_api(self, source, target):
        
        nodes1, nodes2 = keggapi_link(source, target, verbose = True)
        
        return nodes1, nodes2
    
    
    

class KEGGchain(KEGGgraph):
    """ KEGGchain
    Can be used to build chains of KEGGlinkgraphs, builds a composed graph of sequentially linked databases

    Init Arguments
    --------------
    chain : list
        list of databases, they have to be sequentially linked in KEGG

    Properties
    ----------

    chain_dbs : list
        list of databases
    chain : list
        list of chained KEGGlinkgraphs
    """
    chain_dbs = []
    chain = []


    def __init__(self, *args, **kwargs):
        if "chain" in kwargs:
            self.chain_dbs = kwargs.pop("chain")
            
        super().__init__(*args, **kwargs)
        
        if self.chain_dbs != []:    
            self.initchain()
            
            composed_graph = nx.compose_all(self.chain)
            self.add_nodes_from(composed_graph.nodes(data = True)) 
            self.add_edges_from(composed_graph.edges(data = True))
            
            self.name = ">".join(self.chain_dbs)+" chain"
        
        
    def initchain(self):
        self.chain = []
        for i, db in enumerate(self.chain_dbs):
            if i == len(self.chain_dbs) -1 :
                pass
            else:
                link = [self.chain_dbs[i], self.chain_dbs[i+1]]
                source_linkeddb = keggapi_info(link[0], return_format = "dict")['linked db']
                if link[1] not in source_linkeddb:
                    raise KEGGDataBaseError(db = link[1], msg = "KEGG database {} doesn't have a direct link to {}:\nplease compose your chain only using sequentially linkable databases".format(link[0], link[1]))
                self.chain.append(KEGGlinkgraph(source_db = link[0], target_db = link[1]))                
                
        

                
            

    
    
            
        
    
        