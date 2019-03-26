import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 


from tests.download_fresh_test_responses import unpickle_resp


# This method will be used by the mock to replace requests.get
def mocked_requests_get(*args, **kwargs):
    class MockResponse:
        def __init__(self, url = None, json_data = None, text = None, status_code = 404, reason = 'Not Found'):
            self.url = url
            self.json_data = json_data
            self.status_code = status_code
            self.text = text
            self.reason = reason
            
            self.ok = status_code == 200
            
            
            if (self.ok):
                self.reason = 'OK'
            
        def json(self):
            return self.json_data
        
    
    if args[0] == "http://rest.kegg.jp/get/br:br08301/json":
        return MockResponse(url = args[0],
                            json_data = {"key1": "value1"},
                            status_code = 200)
    elif args[0] == 'http://rest.kegg.jp/get/hsa05130/image':
        return unpickle_resp("testpng")
    elif args[0] == "http://rest.kegg.jp/get/C00002/image":
        return unpickle_resp("testgif")
    elif args[0] ==  "http://rest.kegg.jp/list/hsa":
        return MockResponse(url = args[0],
                            text = "testgene1\ttestdescription1\ntestgene2\ttestdescription2",
                            status_code = 200)
    elif args[0] == "http://rest.kegg.jp/returninvalidtext":
        return MockResponse(url = args[0],
                            text = '\n',
                            status_code = 200)
    elif args[0] ==  "http://rest.kegg.jp/info/hsa":
        return MockResponse(url = args[0],
                            text = "hsa             descr_hsa",
                            status_code = 200)
        
    elif args[0] ==  "http://rest.kegg.jp/link/enzyme/hsa":
        return MockResponse(url = args[0],
                            text = "node1\tnode2",
                            status_code = 200)
    return MockResponse()