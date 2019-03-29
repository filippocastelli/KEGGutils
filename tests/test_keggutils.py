import unittest
from unittest.mock import patch
import networkx as nx

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import KEGGutils.KEGGutils as kgu
import KEGGutils.KEGGerrors as errors


hsa_nodes = ['hsa:9344', 'hsa:5894', 'hsa:673', 'hsa:5607']
enzyme_nodes = ['ec:2.7.8.2', 'ec:3.4.23.3', 'ec:2.3.3.10', 'ec:6.4.1.2']

testgraph = nx.Graph()
testgraph.add_nodes_from(hsa_nodes, nodetype = "hsa")
testgraph.add_nodes_from(enzyme_nodes, nodetype = "enzyme")
testgraph.add_edge(hsa_nodes[0], enzyme_nodes[0])


class KEGGutilsTest(unittest.TestCase):
    
    @patch('KEGGutils.KEGGutils.keggapi_link')
    def test_kegg_link_graph_all_nodes_are_correctly_copied(self, mocker):

        mocker.return_value = [hsa_nodes, enzyme_nodes]
        
        graph = kgu.kegg_link_graph("hsa", "enzyme")
        
        self.assertEqual(len(hsa_nodes) + len(enzyme_nodes), len(graph.nodes), "Not all nodes are being copied")
    
    @patch('KEGGutils.KEGGutils.keggapi_link')
    def test_kegg_link_graph_returns_bipartite_graphs(self, mocker):

        mocker.return_value = [hsa_nodes, enzyme_nodes]
        
        graph = kgu.kegg_link_graph("hsa", "enzyme")
        
        self.assertTrue(nx.is_bipartite(graph), "Graph is not bipartite")
        
    @patch('KEGGutils.KEGGutils.keggapi_link')
    def test_kegg_link_graph_adds_nodes(self, mock_keggapi_link):
        
        mock_keggapi_link.return_value = [hsa_nodes, enzyme_nodes]
        
        graph = kgu.kegg_link_graph("hsa", "enzyme")
        
        self.assertTrue(set(graph.nodes) == set(hsa_nodes + enzyme_nodes), "Different number of nodes")
    
    @patch('KEGGutils.KEGGutils.keggapi_link')
    def test_has_nodetypes_actuallyhas(self, mock_keggapi_link):
        
        mock_keggapi_link.return_value = [hsa_nodes, enzyme_nodes]
        
        graph = kgu.kegg_link_graph("hsa", "enzyme")
        
        self.assertTrue(kgu.has_nodetypes(graph))

    def test_has_nodetype_doesnthave(self):
        
        graph = nx.Graph()
        
        self.assertFalse(kgu.has_nodetypes(graph))
        
    @patch('KEGGutils.KEGGutils.keggapi_link')
    def test_get_nodes_by_nodetype(self, mock_keggapi_link):
        
        mock_keggapi_link.return_value = [hsa_nodes, enzyme_nodes]
        
        graph = kgu.kegg_link_graph("hsa", "enzyme")
        
        self.assertEqual(kgu.get_nodes_by_nodetype(graph, "hsa", return_dict = True), dict.fromkeys(hsa_nodes, "hsa"))
        
    @patch('KEGGutils.KEGGutils.keggapi_link')
    def test_get_nodes_by_nodetype_invalidgraph(self, mock_keggapi_link):
        
        mock_keggapi_link.return_value = [hsa_nodes, enzyme_nodes]
        
        graph = nx.Graph()
        
        with self.assertRaises(errors.NotAKeggGraphError):
            kgu.get_nodes_by_nodetype(graph, "hsa")
            
    @patch('KEGGutils.KEGGutils.keggapi_link')
    def test_get_nodes_by_nodetype_missingnodetype(self, mock_keggapi_link):
        
        mock_keggapi_link.return_value = [hsa_nodes, enzyme_nodes]
        
        graph = kgu.kegg_link_graph("hsa", "enzyme")
        
        with self.assertRaises(errors.MissingNodetypeError):
            kgu.get_nodes_by_nodetype(graph, "map")
            
    def test_get_nodes_by_nodetype_correct_nodetype_dict(self):
        hsa_nodetypes = dict.fromkeys(hsa_nodes, "hsa")
        self.assertEqual(kgu.get_nodes_by_nodetype(testgraph, "hsa", return_dict = True), hsa_nodetypes)
        
    def test_unique_nodetypes(self):
        self.assertEqual(kgu.get_unique_nodetypes(testgraph), ['enzyme', 'hsa'])
        
    def test_linked_nodes(self):
        self.assertEqual(kgu.linked_nodes(testgraph, hsa_nodes[0], return_dict = True), {enzyme_nodes[0]: "enzyme"})
        
    def test_neighbor_graph_nodes(self):
        self.assertEqual(set(kgu.neighbor_graph(testgraph, {hsa_nodes[0]: "hsa"}).nodes), set([hsa_nodes[0], enzyme_nodes[0]]))
    
    def test_neighbor_graph_edges(self):
        self.assertEqual(set(list(kgu.neighbor_graph(testgraph, {hsa_nodes[0]: 'hsa'}).edges)[0]), set([hsa_nodes[0], enzyme_nodes[0]]))
        
    