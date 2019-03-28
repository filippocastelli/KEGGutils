import unittest
from unittest.mock import patch
import networkx as nx

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import KEGGutils as kg
from KEGGutils import KEGGgraph, KEGGlinkgraph

from tests.http_mocker import mocked_requests_get


hsa_nodes = ['hsa:9344', 'hsa:5894', 'hsa:673', 'hsa:5607']
enzyme_nodes = ['ec:2.7.8.2', 'ec:3.4.23.3', 'ec:2.3.3.10', 'ec:6.4.1.2']

testgraph = nx.Graph()
testgraph.add_nodes_from(hsa_nodes, nodetype = "hsa")
testgraph.add_nodes_from(enzyme_nodes, nodetype = "enzyme")
testgraph.add_edge(hsa_nodes[0], enzyme_nodes[0])


class KEGGLinkGraphTests(unittest.TestCase):
    
    def setUp(self):
        kg.delete_cached_files(verbose = False)
        
    @patch('KEGGutils.KEGGgraphs.keggapi_link', return_value = [["s1", "s2", "s3"],["t1","t1","t2"] ])
    @patch('KEGGutils.KEGGgraphs.keggapi_info', return_value = {'linked db': ["enzyme"]})
    def test_kegglinkgraph_projection_has_edges(self, link_mocker, info_mocker):
        
        linkgraph = KEGGlinkgraph(source_db = "hsa", target_db = "enzyme")
        projected_graph = linkgraph.projected_graph()
        
        self.assertEqual(list(projected_graph.edges), [("s1", "s2")])
        
    @patch('KEGGutils.KEGGgraphs.keggapi_link', return_value = [hsa_nodes, enzyme_nodes])
    @patch('KEGGutils.KEGGutils.keggapi_link', return_value = [hsa_nodes, enzyme_nodes])
    @patch('KEGGutils.KEGGgraphs.keggapi_info', return_value = {'linked db': ["enzyme"]})
    def test_kegglinkgraph_same_nodes_as_original_function(self, link_mocker,linkmocker2, info_mocker):
        
        linkgraph = KEGGlinkgraph(source_db = "hsa", target_db = "enzyme")
        graph = kg.kegg_link_graph("hsa", "enzyme")
        
        print(linkgraph.nodes())
        print(graph.nodes())
        self.assertCountEqual(list(linkgraph.nodes()), list(graph.nodes()), "KEGGlinkgraph produces different nodes than kegg_link_graph method alone")
        
    @patch('KEGGutils.KEGGgraphs.keggapi_link', return_value = [hsa_nodes, enzyme_nodes])
    @patch('KEGGutils.KEGGgraphs.keggapi_info', return_value = {'linked db': ["enzyme"]})
    def test_kegglinkgraph_correct_nodetypes(self, link_mocker, info_mocker):

        linkgraph = KEGGlinkgraph(source_db = "hsa", target_db = "enzyme")

        self.assertEqual(linkgraph.source_nodes[hsa_nodes[0]], "hsa", "Source node dict contains wrong nodetype")
        self.assertEqual(linkgraph.target_nodes[enzyme_nodes[0]], "enzyme", "Target node dict contains wrong nodetype")
        
    @patch('KEGGutils.KEGGgraphs.keggapi_link', return_value = [hsa_nodes, enzyme_nodes])
    @patch('KEGGutils.KEGGgraphs.keggapi_info', return_value = {'linked db': ["enzyme"]})
    def test_kegglinkgraph_number_of_nodes(self, link_mocker, info_mocker):
        
        linkgraph = KEGGlinkgraph(source_db = "hsa", target_db = "enzyme")
        
        self.assertEqual((len(linkgraph.nodes)),(len(hsa_nodes) + len(enzyme_nodes)))
    
    @patch('KEGGutils.KEGGgraphs.keggapi_link', return_value = [hsa_nodes, enzyme_nodes])
    @patch('KEGGutils.KEGGgraphs.keggapi_info', return_value = {'linked db': ["enzyme"]})
    def test_1(self, linkmocker, infomocker):

        linkgraph = KEGGlinkgraph(source_db = "hsa", target_db = "enzyme")
        
        self.assertEqual(True, True, "ciao")
        
    