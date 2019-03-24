import os, glob, json, requests, shutil
import imghdr
import matplotlib.image as mpimg
import matplotlib.pylab as plt
import logging
import xml.etree.ElementTree as et
import pkg_resources
import pathlib


from slugify import slugify

RES_PATH = pkg_resources.resource_filename('KEGGutils', 'res/')
ORG_CODES = pkg_resources.resource_filename('KEGGutils', 'res/org_codes.txt')

CURRENT_DIR = pathlib.Path.cwd()
DOWNLOAD_DIR = CURRENT_DIR.joinpath("kegg_downloads")
DOWNLOAD_DIR.mkdir(exist_ok = True)

from KEGGutils.KEGGerrors import KEGGOnlineError, KEGGKeyError, KEGGInvalidFileContent
from KEGGutils.KEGGhelpers import push_backslash

#from KEGGutils.KEGGenums import KEGGoutside, KEGGorganisms, KEGGdatabases, KEGGmedicus


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
        
def msg_file_already_exists(filename, download_dir, verbose= True):
    if verbose == True:
        logging.warning("> File {} already present in {}\nreading from chached file...".format(filename, download_dir))
#        print("> File {} already present in {}".format(filename, download_dir))
#        print("reading from cached file...")



# =============================================================================
# FILE MANAGEMENT
# =============================================================================

def delete_cached_files():
    """Deletes all files in download_dir"""
    files_to_remove = glob.glob(download_dir + "*")

    print("> deleting the following files from {}".format(download_dir))
    print(*files_to_remove, sep="\n")

    for file in files_to_remove:
        os.remove(file)

def change_download_dir(newpath):
    """ Changes download directory to the one specified in newpath"""
    global DOWNLOAD_DIR
    DOWNLOAD_DIR = pathlib.Path(newpath).absolute()
    DOWNLOAD_DIR.mkdir(exist_ok = True)
    
def get_download_dir():
    """ Returns download directory path"""
    global DOWNLOAD_DIR
    return DOWNLOAD_DIR    

# =============================================================================
# DOWNLOADS
# =============================================================================

def is_kegg_up():
    """Sends a simple HTTP requests to see if KEGG is currently reachable"""
    
    resp = requests.head('http://rest.kegg.jp/info/kegg')
    
    return resp.status_code == 200

# > URL REQUESTS
def request_image(url, fpath):
    """Creates a requets for url, saves the response in fpath, returns response.raw"""
    response = requests.get(url, stream = True)
    
    if response.status_code == 200:
#        with open(fpath, 'wb') as out_file:
        with fpath.open(mode = "wb") as out_file:
            shutil.copyfileobj(response.raw, out_file)
    else:
        raise KEGGOnlineError
        
    return response.status_code == 200


def get_online_request(url):
    """ Just an errored proxy for requests.get()"""
    request = requests.get(url)
    if request.ok == False:
        raise KEGGOnlineError(request)
    return request

def process_request_text(fulltext, want_descr = False):
    """Preprocessing of item-description type text
    
    Separates couple of words in item \t description format
    
    Parameters:
        :fulltext (str): raw text
        :wannt_descr (bool): if False ignores descriptions and returns items only
        
    
    Returns:
        itemlist, *descriptionlist (list): list of items and descriptions (optional)"""
        
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
    """downloads a text file given an url and a filename
    
    Arguments:
        url {str} -- url for the request
        filename {[type]} -- desired file name
    
    Keyword Arguments:
        force_download {bool} -- If set to True replaces previous versions of the file (default: {False})
        verbose {bool} -- Display additional messages during download (default: {True})
    
    
    Returns:
        text [str] -- downloaded text
    """

    filename = slugify(filename)
    filepath = DOWNLOAD_DIR.joinpath(filename)


    if (not filepath.exists()) or force_download:
        msg_start_download(filename, url, verbose)
        
        request = get_online_request(url)
        text = request.text
        filepath.write_text(text)
            
        msg_end_download(filename, verbose)
        
    else:
        msg_file_already_exists(filename, download_dir, verbose)
        
        text = filepath.read_text()
    
    if text == '\n':
        raise KEGGInvalidFileContent(filename, text)
        
    return text

def download_json(url, filename, force_download=False, verbose=True):
    """ Downloads a json file
    
    Parameters
    ----------
    url : str
        url for the requests
    filename : str
        desired filename
    force_download : bool, optional
            if True replaces any previously downloaded file with the same name (the default is False)
    verbose : bool, optional
        displays additional messages (the default is True)
    
    Returns
    -------
    json
        json data
    """
    filename = slugify(filename)
    
    filepath = DOWNLOAD_DIR.joinpath(filename)
    
    if (filepath.exists()) or force_download:
        msg_start_download(filename, url, verbose)
        
        request = get_online_request(url)
        data = request.json()

        with filepath.open(mode = "w+") as outfile:
            json.dump(data, outfile)
            
        msg_end_download(filename, verbose)
        
    else:
        with filepath.open() as read_file:
            data = json.load(read_file)
            
    return data


def download_pic(url, filename, force_download = False, verbose = False):
    """Downloads a .gif or .png pic fom an url
    
    Parameters
    ----------
    url : str
        url for the request
    filename : str
        desired filename
    force_download : bool, optional
        if set to True replaces any previously downloaded file under the same filename (the default is False)
    verbose : bool, optional
        if set to True displays additional messages (the default is False)
    
    Returns
    -------
    image
        img data stream
    """

    filename = slugify(filename)
#    possible_filenames = {'gif': filename + '.gif', 
#                       'png': filename + '.png'}
    possible_paths = {'gif': DOWNLOAD_DIR.joinpath(filename+".gif"),
                      'png': DOWNLOAD_DIR.joinpath(filename+".png")}
    
    if all(not filepath.exists() for filepath in possible_paths.values()) or force_download:    
        
        msg_start_download(filename, url, verbose)
        
        tmp_path = DOWNLOAD_DIR.joinpath("tmp_img")
        
        imgformat = request_image(url, tmp_path)
        
        img_format = imghdr.what(tmp_path)
        
        assert img_format in ["gif", "png"], "Image format is wrong, something funny is happening in decoding probably"
        
        path = DOWNLOAD_DIR.joinpath("{}.{}".format(filename, img_format))
        tmp_path.rename(path)
        
        img = mpimg.imread(path)
        
    else:
        for imgformat in ["gif", "png"]:
            try:
                path = DOWNLOAD_DIR.joinpath("{}.{}".format(filename, img_format))
                img = mpimg.imread(path)
            except:
                FileNotFoundError
    return img

def download_xml(url, filename, force_download=False, verbose=True):
    """ Downloads a KGML xml file
    
    Parameters
    ----------
    url : str
        url for the requests
    filename : str
        desired filename
    force_download : bool, optional
            if True replaces any previously downloaded file with the same name (the default is False)
    verbose : bool, optional
        displays additional messages (the default is True)
    
    Returns
    -------
    tree
        XML tree
    """
    filename = slugify(filename)
    
    filepath = DOWNLOAD_DIR.joinpath(filename)
    
    if (not filepath.exists() ) or force_download:
        msg_start_download(filename, url, verbose)
        
        response = get_online_request(url)
#        treebytes = response.content
        treestr = response.text
        tree =  et.ElementTree(et.fromstring(treestr))
        
#        mkdir(download_dir)
        
        with filepath.open(mode = "w+") as outfile:
            outfile.write(treestr)
            
        msg_end_download(filename, verbose)
        
    else:
        msg_file_already_exists(filename, download_dir, verbose)
        tree = et.parse(filepath)
            
    return tree
# =============================================================================
# KEGG API COMMANDS
# =============================================================================
    

def keggapi_list(database, option = None, want_descriptions = False, force_download = False):
    """Interface for the KEGG API LIST command

    See https://www.kegg.jp/kegg/rest/keggapi.html for usage

    Parameters
    ----------
    database : str
        Database you wish to obtain a list about, one from pathway | brite | module | ko | genome | <org> | vg | ag | compound | glycan | reaction | rclass | enzyme | network | variant | disease | drug | dgroup | environ | organism | <medicus>
    option : str, optional
        'xl' option is applicable only to the 'brite' database for listing binary relation files (the default is None)
    want_descriptions : bool, optional
        If True returns descriptions for each item (the default is False)
    force_download : bool, optional
        If true replaces any pre-existing downloaded file (the default is False)
    
    Returns
    -------
    item, descritpions
        lists for items and their descriptions
    """
    
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
    """Interface for the KEGG API FIND command
    
    See https://www.kegg.jp/kegg/rest/keggapi.html for further info
    
    Parameters
    ----------
    database : str
        KEGG database you wish to query on, one from pathway | brite | module | ko | genome | genes | <org> | vg | ag | ligand | compound | glycan | reaction | rclass | enzyme | network | variant | disease | drug | dgroup | environ | <medicus>  
    query : str
        Desired query
    option : str, optional
        if database is "compound" or "drug", possible options are "formula", "exact mass" and "mol_weight" (the default is None)
    want_descriptions : bool, optional
        if set to True returns a list of descriptions for the found items (the default is False,)
    verbose : bool, optional
        if set to True displays additional messages (the default is False)
    force_download : bool, optional
        if set to True replaces any previously downloaded file under the same name (the default is False)
    
    Returns
    -------
    list, list
        list of items and descriptions
    """


    
    options = ["formula", "exact_mass", "mol_weight"]
    
    if database not in db_categories + ["genes"]:
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
    
    
def keggapi_get(dbentry, option = None, want_descriptions = False, verbose = False, force_download = False, show_result_image = True):
    """Interface for the KEGG API GET command

    for further info read https://www.kegg.jp/kegg/rest/keggapi.html
    
    Parameters
    ----------
    dbentry : str
        KEGG database entry you wish to GET, one from pathway | brite | module | ko | genome | <org> | vg | ag | compound | glycan | reaction | rclass | enzyme | network | variant | disease | drug | dgroup | environ | disease_ja | drug_ja | dgroup_ja | environ_ja | compound_ja
    option : str, optional
        one from aaseq | ntseq | mol | kcf | image | conf | kgml | json (the default is None])
    want_descriptions : bool, optional
        if True returns a list of descriptions for the requested items (the default is False)
    verbose : bool, optional
        is True displays additional messages  (the default is True)
    force_download : bool, optional
        if set to True replaces any file under the same filename (the default is False)
    show_result_image: boo, optional
        if set to True shows the downloaded image (the default is True)
    """

    
    options = ["aaseq","ntseq", "mol", "kcf","image","conf", "kgml","json"]
    
    if (option is not None) & (option not in options):
        raise KEGGKeyError(option, msg = "option {} invalid for GET".format(option))

    option, optionurl = push_backslash(option)
    
    url = "http://rest.kegg.jp/get/{}{}".format(dbentry, optionurl)
    
    if option == "":
        option = "description"

    filename = dbentry + "_" + option
    
    if option == "description":
        infos = download_textfile(url, filename, verbose = False, force_download = force_download)
        if verbose == False:
            infos = "\n".join(infos.splitlines()[1:4])
        print("Infos on {} from KEGG:\n".format(dbentry))
        print(infos)
        return
    elif option == "kgml":
        tree = download_xml(url, filename, force_download=force_download, verbose = verbose)
        return tree
    elif option == "json":
        json_data = download_json(url, filename, verbose = verbose)
        return json_data
    elif option in ["mol", "kcf", "conf", "ntseq"]:
        text = download_textfile(url, filename, verbose = verbose, force_download = force_download)
        print(text)
        return text
    elif option == "aaseq":
        text = download_textfile(url, filename, verbose = True, force_download = force_download)
        description = text.splitlines()[0]
        sequence = "".join(text.splitlines()[1:])
        if want_descriptions == True:
            return description, sequence
        else:
            return sequence
    elif option == "image":
        img = download_pic(url, filename, verbose = True)
        if show_result_image:
            plt.imshow(img)
        return img
    elif option == "kgml":
        raise NotImplementedError
    else:
        raise KEGGKeyError(key = dbentry, msg = "KEGG GET API request not recognized")

    
def keggapi_link(source, target, verbose = True, force_download = False):
    """Interface for the KEGG REST API LINK command 
    Given two different database names returns the linked relations between them
    
    for further info read https://www.kegg.jp/kegg/rest/keggapi.html
    
    Parameters
    ----------
    source : str
        source database name
    target : str
        target database name
    verbose : bool, optional
        displays additional infos during the download (the default is True)
    force_download : bool, optional
        forces overwriting over cached files (the default is False)
    
    Raises
    ------
    KEGGKeyError
        if a database key is invalid
    
    Returns
    -------
    link1 : list
        list of source nodes
    link2 : list
        list of target nodes
    """

    
    if target not in db_categories:
        raise KEGGKeyError(target, msg = "source database {} is not a valid database".format(target))
    
    url = "http://rest.kegg.jp/link/{}/{}".format(target, source)
    
    filename = target + "_" + source + "_link"
    
    text = download_textfile(url, filename, force_download = force_download)
    
    link1, link2 = process_request_text(text, want_descr = True)
    
    
    return link1, link2

def keggapi_conv(source, target, verbose = True, force_download = False):
    """ KEGG REST API interface to CONV command
    Converts KEGG codes to and from NCBI ProteinID, NCBI GeneID, Uniprot, CHEBI\
    and PubChem name standards
    
    for further info read https://www.kegg.jp/kegg/rest/keggapi.html
    
    Parameters
    ----------
    
    source str :
        source database or dbentry
    target str :
        target database or dbentry
    verbose bool : 
        if set to True displays additional messages (defaultis True)
    force_download bool: 
        forces overwriting cached files (default is False)
        
    Returns
    -------
    
    source_codes : list
        list of codes in the original database format
    target_codes : list
        list of codes in the target database format
    
    """
    
    org = get_organism_codes() + ["genes"]
    kegg_db = None
    outside_db = None
    keggdblist = ["compound", "glycan", "drug"]
    outsidedb1 = ["ncbi-geneid", "ncbi-proteinid", "uniprot"]
    outsidedb2 = ["pubchem", "chebi"]
    outsidedblist = outsidedb1 + outsidedb2
    
    
    #CHECK IF ARGUMENTS ARE VALID
    for entry in [source, target]:
        if entry in org + keggdblist:
            kegg_db = entry
        elif entry in outsidedblist:
            outside_db = entry
    
    if (kegg_db is None) ^ (outside_db is None):
        #second mode, source is a dbentry
        if target not in org + keggdblist + outsidedblist:
            raise KEGGKeyError(key = target)
        #further checks to implement
    else:
        if kegg_db in org:
            if outside_db not in outsidedb1:
                raise KEGGKeyError(key = outside_db)
        elif kegg_db in keggdblist:
            if outside_db not in outsidedb2:
                raise KEGGKeyError(key = outside_db)
        else:
            raise KEGGKeyError(key = kegg_db)


    url = "http://rest.kegg.jp/conv/{}/{}".format(target, source)
    
    filename = target + "_" + source + "_conv"
    
    text = download_textfile(url, filename, force_download = force_download)
    
    codes1, codes2 = process_request_text(text, want_descr = True)
    
    return codes1, codes2
        


def keggapi_info(database, verbose = True, force_download = False):
    """KEGG REST API interface for INFO command
    Displays information on a given database
    
    for further info read https://www.kegg.jp/kegg/rest/keggapi.html
    
    Parameters
    ----------
    database : str
        database of which you want to obtain infos on
    verbose : bool
        if set to False displays only the first 4 lines of text (default is True)
    force_download :  bool
        forces overwriting on previous cached files (default is False)
    """
    org = get_organism_codes()
    
    if database not in db_categories + org:
        raise KEGGKeyError(database, msg = "source database {} is not a valid database".format(database))
        
    url = "http://rest.kegg.jp/info/{}".format(database)
    
    filename = database+"_info"
    
    infos = download_textfile(url, filename, verbose=False)
    if verbose == False:
        infos = "\n".join(infos.splitlines()[1:4])

    print("Infos on {} from KEGG:\n".format(database))
    print(infos)
    
    
def keggapi_ddi(dbentry, force_download = False):
    """KEGG REST API interface for the DDI command
    lists drug-drug interactions for a given compound name
    
    Parameters
    ----------
    
    dbentry : str
        drug KEGG database entry
    force_download : bools
        forces overwriting over cached files (default is False)
        
    for further info read https://www.kegg.jp/kegg/rest/keggapi.html"""
    
    
        
    url = "http://rest.kegg.jp/ddi/{}".format(dbentry)
    
    filename = dbentry+"_ddi"
    
    text = download_textfile(url, filename, verbose=False)

    ddi_list = []
    
    for line in text.splitlines():
        drug1, drug2, ddi_code, interaction = line.strip().split("\t")
        ddi_list.append((drug1, drug2, ddi_code, interaction))
        
    return ddi_list
    
    

    
# =============================================================================
# KEGGAPI DERIVATE FUNCTIONS
# =============================================================================
def get_organism_codes(force_download = False):
    """Returns all KEGG Organism name codes """

    org_url = "http://rest.kegg.jp/list/organism"
    org_filename = "organism_code_list"

    org_codes = []
    
    
    if force_download:
        organism_fulltext = download_textfile(org_url, org_filename, verbose=False, force_download = force_download)
        
        for line in organism_fulltext.splitlines():
            T_identifier, kegg_code, description, hier = line.strip().split("\t")
            org_codes.append(kegg_code)
            
    else:
        with open("./res/org_codes.txt", "r+") as org_file:
            organism_fulltext = org_file.read()
            
        for line in organism_fulltext.splitlines():
            org_codes.append(line)
        


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

