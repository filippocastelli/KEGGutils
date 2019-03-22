import requests
import pickle

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

test_files_dir = "test_files_dir"

test_files_path = os.path.join(currentdir, test_files_dir)

def pickle_resp(url, filename, stream = False):
    print("> downloading: {}".format(filename))
    fpath = os.path.join(test_files_path, filename)
    
    response = requests.get(url, stream = stream)
    
    outfile = open(fpath, "wb")
    pickle.dump(response, outfile)
    outfile.close()
    
    return response


def unpickle_resp(filename):
    
    fpath = os.path.join(test_files_path, filename)
    
    infile = open(fpath, "rb")
    resp = pickle.load(infile)
    infile.close()
    
    return resp

def mkdir(directory):
    """ Creates a directory if it doesn't already exsist"""
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass

mkdir(test_files_dir)

def download_testfiles():
    
# =============================================================================
# GIF IMAGE
# =============================================================================
    
    url_gif = "http://rest.kegg.jp/get/C00002/image"
    pickle_resp(url_gif, "testgif")
    
# =============================================================================
# PNG IMAGE 
# =============================================================================
    
    url_png = "http://rest.kegg.jp/get/hsa05130/image"
    pickle_resp(url_png, "testpng")
    

# =============================================================================
# JSON
# =============================================================================
    
    url_json = "http://rest.kegg.jp/get/br:br08301/json"
    pickle_resp(url_json, "testjson")
    
# =============================================================================
# LIST OF GENES
# =============================================================================
    url_list = "http://rest.kegg.jp/list/hsa"
    pickle_resp(url_list, "testlist")