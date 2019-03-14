import os, glob, json, requests, shutil
import networkx as nx
import imghdr
import matplotlib.image as mpimg
import matplotlib.pylab as plt

from KEGGerrors import KEGGOnlineError, KEGGKeyError, KEGGInvalidFileContent
from KEGGhelpers import push_backslash

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


# =============================================================================
# INTERFACE MESSAGES
# =============================================================================

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



# =============================================================================
# FILE MANAGEMENT
# =============================================================================
def file_exists(filename):
    return os.path.isfile(get_fname_path(filename))

def get_fname_path(filename):
    return download_dir+ filename

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


# =============================================================================
# DOWNLOADS
# =============================================================================

# > URL REQUESTS
def request_image(url, fpath):
    response = requests.get(url, stream = True)
    
    with open(fpath, 'wb') as out_file:
        shutil.copyfileobj(response.raw, out_file)
        
    return response.raw


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
    
    if text == '\n':
        raise KEGGInvalidFileContent(filename, text)
        
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
    
    possible_filenames = {'gif': filename + '.gif', 
                       'png': filename + '.png'}

    if all(not file_exists(filenames) for filenames in possible_filenames.values()) or force_download:    
        
        msg_start_download(filename, url, verbose)
        
        tmp_path = get_fname_path("tmp_img")
        
        img = request_image(url, tmp_path)
        
        img_format = imghdr.what(tmp_path)
        
        assert img_format in ["gif", "png"], "Image format is wrong, something funny is happening in decoding probably"
        
        os.rename(tmp_path, get_fname_path(filename)+"."+img_format)
        
    else:
        for imgformat in ["gif", "png"]:
            try:
                img = mpimg.imread(get_fname_path(filename)+"."+imgformat)
            except:
                FileNotFoundError
    return img

# =============================================================================
# KEGG API COMMANDS
# =============================================================================
    

def keggapi_list(database, option = None, want_descriptions = False, force_download = False):
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
    
    list_fulltext = download_textfile(url, filename,force_download = force_download)
    
    if want_descriptions == False:
        itemlist = process_request_text(list_fulltext, want_descr=want_descriptions)
        return itemlist
    
    elif want_descriptions == True:
        itemlist, descriptionlist = process_request_text(list_fulltext, want_descr=want_descriptions)
        assert len(itemlist) == len(descriptionlist), "different length, funny innit"
        
        return itemlist, descriptionlist


def keggapi_find(database, query, option = None, want_descriptions = False, verbose = False, force_download = False):
    
    options = ["formula", "exact_mass", "mol_weight"]
    
    if database not in db_categories:
        raise KEGGKeyError(database)
        
    if (database == "compound" or database == "drug") & (option is not None) & (option not in options):
        raise KEGGKeyError(database, msg = "only {} opts are available for {} database".format(options, database))
        
    query, queryurl = push_backslash(query)
    option, optionurl = push_backslash(option)
    
    url = "http://rest.kegg.jp/find/{}{}{}".format(database, queryurl, optionurl)
    
    filename = database + "_" + query + "_" + option
    fulltext = download_textfile(url, filename, verbose = verbose, force_download = force_download)
    if want_descriptions == False:
        itemlist = process_request_text(fulltext, want_descr=want_descriptions)
        return itemlist
    elif want_descriptions == True:
        itemlist, descriptionlist = process_request_text(fulltext, want_descr=want_descriptions)
        assert len(itemlist) == len(descriptionlist), "different item and description lengths, something's not working"
        
        return itemlist, descriptionlist
    
    
def keggapi_get(dbentry, option = None, want_descriptions = False, verbose = True, force_download = False):
    
    options = ["aaseq","ntseq", "mol", "kcf","image","conf", "kgml","json"]
    
    if (option is not None) & (option not in options):
        raise KEGGKeyError(option, msg = "option {} invalid for GET".format(option))

    option, optionurl = push_backslash(option)
    
    url = "http://rest.kegg.jp/get/{}{}".format(dbentry, optionurl)
    
    if option == None:
        option = "description"

    filename = dbentry + "_" + option
    
    if option == "description":
        infos = download_textfile(url, filename, verbose = False, force_download = force_download)
        if verbose == False:
            infos = "\n".join(infos.splitlines()[1:4])
        print("Infos on {} from KEGG:\n".format(dbentry))
        print(infos)
        return infos
    elif option == "json":
        json_data = download_json(url, filename, verbose = verbose)
        return json_data
    elif option in ["mol", "kcf", "conf", "aaseq", "ntseq"]:
        text = download_textfile(url, filename, verbose = verbose, force_download = force_download)
        print(text)
        return text
    elif option == "image":
        img = download_pic(url, filename, verbose = True)
        plt.imshow(img)
        
        return img
    
    if want_descriptions == False:
        itemlist = process_request_text(text, want_descr=want_descriptions)
        return itemlist
    elif want_descriptions == True:
        itemlist, descriptionlist = process_request_text(text, want_descr=want_descriptions)
        assert len(itemlist) == len(descriptionlist), "different item and description lengths, something's not working"
        return itemlist, descriptionlist
    
def keggapi_link(source, target, verbose = True, force_download = False):
    
    if target not in db_categories:
        raise KEGGKeyError(target, msg = "source database {} is not a valid database".format(target))
    
    url = "http://rest.kegg.jp/link/{}/{}".format(target, source)
    
    filename = target + "_" + source + "_link"
    
    text = download_textfile(url, filename, force_download = force_download)
    
    link1, link2 = process_request_text(text, want_descr = True)
    
    
    return link1, link2
    

# =============================================================================
# KEGGAPI DERIVATE FUNCTIONS
# =============================================================================
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


def kegg_url(target_db, source_db): #deprecated
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

        
        
        
# =============================================================================
# MISC
# =============================================================================
        
