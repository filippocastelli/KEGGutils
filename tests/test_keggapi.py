import unittest
from unittest.mock import patch

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import KEGGutils.KEGGapi as kgapi
import KEGGutils.KEGGerrors as kgerrors


from tests.http_mocker import mocked_requests_get

class KEGGapiTest(unittest.TestCase):
    
    def setUp(self):
        kgapi.delete_cached_files(verbose = False)
    
    
    @patch('requests.get', side_effect = mocked_requests_get)
    def test_keggapi_info_returns_str(self, mockrequest):
        
        response_str = kgapi.keggapi_info("hsa", return_format = "str")

        self.assertEqual(response_str,  "hsa             descr_hsa")
        
    @patch('requests.get', side_effect = mocked_requests_get)
    def test_keggapi_info_returns_dict(self, mockrequest):
        
        response_dict = kgapi.keggapi_info("hsa", return_format = "dict")
        
        expected_dict = {"hsa": ["descr_hsa"]}
        self.assertEqual(response_dict, expected_dict)
    
    @patch('requests.get', side_effect = mocked_requests_get)
    def test_download_json_file_doesnt_exist(self, mockrequestsget):
        
        data = kgapi.download_json("http://rest.kegg.jp/get/br:br08301/json", "json_testing")
        self.assertEqual(data, {"key1": "value1"})
        
    @patch('KEGGutils.KEGGapi.get_organism_codes', return_value = ['hsa'])
    @patch('requests.get', side_effect = mocked_requests_get,)
    def test_keggapi_list_wantdescrFalse(self, mocked_orgcodes, mockrequestget):
        itemlist = kgapi.keggapi_list("hsa", want_descriptions = False)
        test_genes = set(["testgene1", "testgene2"])
        self.assertEqual(set(itemlist), set(test_genes))
        
        
    @patch('KEGGutils.KEGGapi.get_organism_codes', return_value = ['hsa'])
    @patch('requests.get', side_effect = mocked_requests_get,)
    def test_keggapi_list_wantdescrTrue(self, mocked_orgcodes, mockrequestget):
        itemlist, descriptionlist, _, _ = kgapi.keggapi_list("hsa", want_descriptions = True)
        test_genes = set(["testgene1", "testgene2"])
        test_descriptions = set(["testdescription1", "testdescription2"])
        
        self.assertEqual(set(itemlist), set(test_genes))
        self.assertEqual(set(descriptionlist), set(test_descriptions))
        
    @patch('requests.get', side_effect = mocked_requests_get)
    def test_download_textfile_file_doesnt_exist(self, mockedget):
        text = "testgene1\ttesttranscripttype1\ttesttranscriptpos1\ttestdescription1\ntestgene2\ttesttranscripttype2\ttesttranscriptpos2\ttestdescription2"
        self.assertEqual(kgapi.download_textfile("http://rest.kegg.jp/list/hsa", "textfile_testing", force_download = True), text)
        
    @patch('requests.get', side_effect = mocked_requests_get)
    def test_download_textfile_file_exists(self, mocker):
        filepath = kgapi.DOWNLOAD_DIR.joinpath("textfile-testing")
        text = "surprise_mofo"
        try:
            os.remove(filepath)
        except FileNotFoundError:
            pass
        
        filepath.write_text(text)
        
        self.assertEqual(kgapi.download_textfile("http://rest.kegg.jp/list/hsa", "textfile-testing"),text)
        
    @patch('requests.get', side_effect = mocked_requests_get)
    def test_download_textfile_raise_error_if_newline(self, mockrequest):
        
        with self.assertRaises(kgerrors.KEGGInvalidFileContent):
            kgapi.download_textfile(url =  "http://rest.kegg.jp/returninvalidtext", filename = "textfile_testing", force_download = True)
            

            
        
        
        
        
        
        
        
        