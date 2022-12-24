import os, glob, json, requests, shutil
import imghdr
import matplotlib.image as mpimg
import matplotlib.pylab as plt
import logging
import xml.etree.ElementTree as et
import pkg_resources
import pathlib


from slugify import slugify

RES_PATH = pathlib.Path(pkg_resources.resource_filename("KEGGutils", "res/"))
ORG_CODES = pathlib.Path(pkg_resources.resource_filename("KEGGutils", "res/org_codes.txt"))

CURRENT_DIR = pathlib.Path.cwd()
DOWNLOAD_DIR = CURRENT_DIR.joinpath("kegg_downloads")

DOWNLOAD_DIR.mkdir(exist_ok=True)

from KEGGutils.KEGGerrors import KEGGOnlineError, KEGGKeyError, KEGGInvalidFileContent, KEGGInvalidContent
from KEGGutils.KEGGhelpers import push_backslash

# from KEGGutils.KEGGenums import KEGGoutside, KEGGorganisms, KEGGdatabases, KEGGmedicus


#download_dir = "./kegg_downloads/"

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


def msg_start_download(filename, url, verbose=True):
    if verbose == True:
        logging.info("> Downloading %s from KEGG at %s",filename, url)


def msg_end_download(filename, verbose=True):
    if verbose == True:
        logging.info("succesfully downloaded %s",filename)


def msg_file_already_exists(filename, verbose=True):
    if verbose == True:
        logging.warning(
            "> File %s already present in %s\nreading from chached file...",
                filename, DOWNLOAD_DIR
        )



# =============================================================================
# FILE MANAGEMENT
# =============================================================================


def delete_cached_files(verbose = True):
    """Deletes all files in download_dir"""
    
    files_to_remove = list(DOWNLOAD_DIR.glob("*"))
    files_to_remove = [file for file in files_to_remove if file.is_file()]
    
    if verbose == True:
        print("> deleting the following files from {}".format(str(DOWNLOAD_DIR)))
        print(*files_to_remove, sep="\n")

    for file in files_to_remove:
        os.remove(file)


def change_download_dir(newpath):
    """ Changes download directory to the one specified in newpath"""
    global DOWNLOAD_DIR
    DOWNLOAD_DIR = pathlib.Path(newpath).absolute()
    DOWNLOAD_DIR.mkdir(exist_ok=True)


def get_download_dir():
    """ Returns download directory path"""
    global DOWNLOAD_DIR
    return DOWNLOAD_DIR


# =============================================================================
# DOWNLOADS
# =============================================================================


def is_kegg_up():
    """Sends a simple HTTP requests to see if KEGG is currently reachable"""

    resp = requests.head("http://rest.kegg.jp/info/kegg")

    return resp.status_code == 200


# > URL REQUESTS
def request_image(url, fpath):
    """Creates a requets for url, saves the response in fpath, returns response.raw"""
    response = requests.get(url, stream=True)

    if response.status_code == 200:
        #        with open(fpath, 'wb') as out_file:
        with fpath.open(mode="wb") as out_file:
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


def process_request_text(fulltext, want_descr=False, mode = "bipartite_list"):
    """Preprocessing of item-description type text
    
    Separates couple of words in item \t description format
    
    Parameters
    ----------
    fulltext : str
        raw text
    mode : str
        parsing mode, valid ones are bipartite_list ! columns ! nested
    want_desc : bool
        valid only for "bipartite" mode, if False returns only first of the two lists

    Returns
    -------
    itemlist, descriptionlist : list
        list of items and their descriptions, only for bipartite_list mode
    parsed_dict : dict
        dictionary from parsed text, only for columns mode
        """
        
    validmodes = ["bipartite_list", "columns", "nested", "four_way_list"]
    
    if mode not in validmodes:
        raise ValueError
    
    # temporary fix for Oct22 API change involving /list/<org> requests
    # will be removed in KEGGutils 1.x release
    
    if mode == "four_way_list":
        itemlist = []
        descriptionlist = []
        transcripttypelist = []
        transcriptpositionlist = []
        
        for line in fulltext.splitlines():
            entry, transcript_type, transcript_position, description = line.strip().split("\t")
            itemlist.append(entry)
            descriptionlist.append(description)
            transcripttypelist.append(transcript_type)
            transcriptpositionlist.append(transcript_position)
        
        if want_descr == True:
            return itemlist, descriptionlist, transcripttypelist, transcriptpositionlist
        else:
            return itemlist
            
    if mode == "bipartite_list":
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
    if mode == "columns":
        parsed_dict = {}
        for line in fulltext.splitlines():
            split_line = line.split("  ")
            
            split_line = [entry.lstrip() for entry in split_line if entry != ""]
            
            if len(split_line) == 0:
                pass
            elif len(split_line) == 2:
                current_key, item = split_line
                
                parsed_dict[current_key] = []
                parsed_dict[current_key].append(item)
                #remove all empty entries
            elif len(split_line) == 1:
                
                try:
                    parsed_dict[current_key]
                except KeyError:
                    #shouldn't happen tho
                    raise KEGGInvalidContent(fulltext, msg = "text response cannot be parsed")
                    
                item = split_line[0]
                parsed_dict[current_key].append(item)
            else:
                raise NotImplementedError
                
        return parsed_dict
    
    
    if mode == "nested":
        
        lines = fulltext.splitlines()
        
        lastkey = None
        key_id = 0
        keys = []
        subdict = {}
        parsed_dict = {}
        for line in lines:
            
            entry_str = line[:12]
            content_str = line[12:]
            
            
            split_line = line.split("  ")
            
            split_line = [entry.lstrip() for entry in split_line if entry != ""]
            
            if entry_str.lstrip().rstrip() != "":
                entry_str_split = entry_str.split("  ")
                if entry_str_split[0] == "":
                    subkey = entry_str_split[1].lstrip().rstrip().lower()
#                    print("entering subkey {}".format(subkey))
                    flag = "s"
                else:
                    keytxt = entry_str_split[0].lstrip().rstrip()
                    if keytxt == "///":
                        pass
                    else:
                        key = keytxt.lower()
                        if key != lastkey or flag != "k":
                            parsed_dict.update({lastkey: subdict})
                            subdict = {}
                        if key in keys:
                            key = key + str(key_id)
                            key_id = key_id +1
#                        print("entering key {}".format(key))
                        keys.append(key)
                        flag = "k"
                        subkey = None
            
            if keytxt != "///":
                content = content_str.lstrip().rstrip()
                if flag == "k":
                    if "reference" in key:
                        subk = "reference_hook"
                    else:
                        subk = "descr"
                    subdict.update({subk : content})
                elif flag == "s":
                    subdict.update({subkey : content})
                
                lastkey = key
        
        parsed_dict.pop(None)
        return parsed_dict
                
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
    DOWNLOAD_DIR.mkdir(exist_ok=True)
    
    filename = slugify(filename)
    filepath = DOWNLOAD_DIR.joinpath(filename)

    if (not filepath.exists()) or force_download:
        msg_start_download(filename, url, verbose)

        request = get_online_request(url)
        text = request.text
        filepath.write_text(text)

        msg_end_download(filename, verbose)

    else:
        msg_file_already_exists(filename, verbose)

        text = filepath.read_text()

    if text == "\n":
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
    DOWNLOAD_DIR.mkdir(exist_ok=True)
    
    filename = slugify(filename)

    filepath = DOWNLOAD_DIR.joinpath(filename)

    if (not filepath.exists()) or force_download:
        msg_start_download(filename, url, verbose)

        request = get_online_request(url)
        data = request.json()

        with filepath.open(mode="w+") as outfile:
            json.dump(data, outfile)

        msg_end_download(filename, verbose)

    else:
        with filepath.open() as read_file:
            data = json.load(read_file)

    return data


def download_pic(url, filename, force_download=False, verbose=False):
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
    DOWNLOAD_DIR.mkdir(exist_ok=True)
    
    filename = slugify(filename)
    #    possible_filenames = {'gif': filename + '.gif',
    #                       'png': filename + '.png'}
    possible_paths = {
        "gif": DOWNLOAD_DIR.joinpath(filename + ".gif"),
        "png": DOWNLOAD_DIR.joinpath(filename + ".png"),
    }

    if (
        all(not filepath.exists() for filepath in possible_paths.values())
        or force_download
    ):

        msg_start_download(filename, url, verbose)

        tmp_path = DOWNLOAD_DIR.joinpath("tmp_img")
        
        request_image(url, tmp_path)

        img_format = imghdr.what(tmp_path)

        assert img_format in [
            "gif",
            "png",
        ], "Image format is wrong, something funny is happening in decoding probably"

        path = DOWNLOAD_DIR.joinpath("{}.{}".format(filename, img_format))
        tmp_path.rename(path)

        img = mpimg.imread(str(path))

    else:
        for imgformat in ["gif", "png"]:
            try:
                path = DOWNLOAD_DIR.joinpath("{}.{}".format(filename, imgformat))
                img = mpimg.imread(str(path))
            except FileNotFoundError:
                pass
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
    DOWNLOAD_DIR.mkdir(exist_ok=True)
    
    filename = slugify(filename)

    filepath = DOWNLOAD_DIR.joinpath(filename)

    if (not filepath.exists()) or force_download:
        msg_start_download(filename, url, verbose)

        response = get_online_request(url)
        #        treebytes = response.content
        treestr = response.text
        tree = et.ElementTree(et.fromstring(treestr))


        with filepath.open(mode="w+") as outfile:
            outfile.write(treestr)

        msg_end_download(filename, verbose)

    else:
        msg_file_already_exists(filename, verbose)
        tree = et.parse(filepath)

    return tree


# =============================================================================
# KEGG API COMMANDS
# =============================================================================


def keggapi_list(database, option=None, want_descriptions=False, force_download=False, return_url = False):
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
    return_url : bool, optional
        If True returns the interrogated URL
    
    Returns
    -------
    item, descritpions
        lists for items and their descriptions
    """

    org_codes = get_organism_codes()
    if database not in db_categories:
        raise KEGGKeyError(database)

    if (option == "xl") & (database != "brite"):
        raise KEGGKeyError(
            database, msg="option xl can only be used with brite argument"
        )

    if (option in org_codes) & (database != "pathway") & (database != "module"):
        raise KEGGKeyError(
            database,
            msg="only pathway and module list request are available for {}".format(
                option
            ),
        )
    
    option, optionurl = push_backslash(option)

    url = "http://rest.kegg.jp/list/{}{}".format(database, optionurl)
    if return_url == True:
        return url
    filename = database + "_" + option + "_list"

    list_fulltext = download_textfile(url, filename, force_download=force_download)

    if database in org_codes:
        # since October 2022 /list/<org> includes chromosomal positions
        # reponses are now in the form of <org>:<gene> <transcript_type> <position> <description>
        # this is not compatible with the current implementation of process_request_text
        itemlist, descriptionlist, transcripttypelist, transcriptpositionlist = process_request_text(list_fulltext, want_descr=True, mode="four_way_list")
        if want_descriptions:
            return itemlist, descriptionlist, transcripttypelist, transcriptpositionlist
        else:
            return itemlist
    else:
        if want_descriptions == False:
            itemlist = process_request_text(list_fulltext, want_descr=want_descriptions)
            return itemlist

        elif want_descriptions == True:
            itemlist, descriptionlist = process_request_text(
                list_fulltext, want_descr=want_descriptions
            )
            assert len(itemlist) == len(descriptionlist), "different length, funny innit"

            return itemlist, descriptionlist


def keggapi_find(
    database,
    query,
    option=None,
    want_descriptions=False,
    verbose=False,
    force_download=False,
    return_url = False,
):
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
    return_url : bool, optional
        If True returns the interrogated URL
    Returns
    -------
    list, list
        list of items and descriptions
    """

    options = ["formula", "exact_mass", "mol_weight"]

    if database not in db_categories + ["genes"]:
        raise KEGGKeyError(database)

    if (
        (database in ("compound", "drug"))
        & (option is not None)
        & (option not in options)
    ):
        raise KEGGKeyError(
            database,
            msg="only {} opts are available for {} database".format(options, database),
        )

    query, queryurl = push_backslash(query)
    option, optionurl = push_backslash(option)

    url = "http://rest.kegg.jp/find/{}{}{}".format(database, queryurl, optionurl)
    
    if return_url == True:
        return url
    
    filename = database + "_" + query + "_" + option
    fulltext = download_textfile(
        url, filename, verbose=verbose, force_download=force_download
    )
    if want_descriptions == False:
        itemlist = process_request_text(fulltext, want_descr=want_descriptions)
        return itemlist
    elif want_descriptions == True:
        itemlist, descriptionlist = process_request_text(
            fulltext, want_descr=want_descriptions
        )
        assert len(itemlist) == len(
            descriptionlist
        ), "different item and description lengths, something's not working"

        return itemlist, descriptionlist


def keggapi_get(
    dbentry,
    option=None,
    want_descriptions=False,
    verbose=False,
    force_download=False,
    show_result_image=True,
    return_dict = False,
    return_text = False,
    return_url = False,
):
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
    return_dic : bool, optional
        if set to True, returns description in dict format
    return_url : bool, optional
        If True returns the interrogated URL
        """

    options = ["aaseq", "ntseq", "mol", "kcf", "image", "conf", "kgml", "json"]

    if (option is not None) & (option not in options):
        raise KEGGKeyError(option, msg="option {} invalid for GET".format(option))

    option, optionurl = push_backslash(option)

    url = "http://rest.kegg.jp/get/{}{}".format(dbentry, optionurl)
    if return_url == True:
        return url
    
    if option == "":
        option = "description"

    filename = dbentry + "_" + option

    if option == "description":
        
        infos = download_textfile(
            url, filename, verbose=False, force_download=force_download
        )
        if return_dict == True:
            parsed_dict = process_request_text(infos, mode = "nested")
            return parsed_dict
        elif return_text == True:
            return infos
        if verbose == False:
            infos = "\n".join(infos.splitlines()[1:4])
        print("Infos on {} from KEGG:\n".format(dbentry))
        print(infos)

    elif option == "kgml":
        tree = download_xml(
            url, filename, force_download=force_download, verbose=verbose
        )
        return tree
    elif option == "json":
        json_data = download_json(url, filename, verbose=verbose)
        return json_data
    elif option in ["mol", "kcf", "conf", "ntseq"]:
        text = download_textfile(
            url, filename, verbose=verbose, force_download=force_download
        )
        print(text)
        return text
    elif option == "aaseq":
        text = download_textfile(
            url, filename, verbose=True, force_download=force_download
        )
        description = text.splitlines()[0]
        sequence = "".join(text.splitlines()[1:])
        if want_descriptions == True:
            return description, sequence
        else:
            return sequence
    elif option == "image":
        img = download_pic(url, filename, verbose=True)
        if show_result_image:
            plt.imshow(img)
        return img
    elif option == "kgml":
        raise NotImplementedError
    else:
        raise KEGGKeyError(key=dbentry, msg="KEGG GET API request not recognized")


def keggapi_link(source, target, verbose=True, force_download=False, return_url = False):
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
    return_url : bool, optional
        If True returns the interrogated URL
        
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
        raise KEGGKeyError(
            target, msg="source database {} is not a valid database".format(target)
        )

    url = "http://rest.kegg.jp/link/{}/{}".format(target, source)
    if return_url == True:
        return url
    
    
    filename = target + "_" + source + "_link"

    text = download_textfile(url, filename, force_download=force_download, verbose = verbose)

    link1, link2 = process_request_text(text, want_descr=True)

    return link1, link2


def keggapi_conv(source, target, verbose=True, force_download=False, return_url = False):
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
    return_url : bool, optional
        If True returns the interrogated URL        
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

    # CHECK IF ARGUMENTS ARE VALID
    for entry in [source, target]:
        if entry in org + keggdblist:
            kegg_db = entry
        elif entry in outsidedblist:
            outside_db = entry

    if (kegg_db is None) ^ (outside_db is None):
        # second mode, source is a dbentry
        if target not in org + keggdblist + outsidedblist:
            raise KEGGKeyError(key=target)
        # further checks to implement
    else:
        if kegg_db in org:
            if outside_db not in outsidedb1:
                raise KEGGKeyError(key=outside_db)
        elif kegg_db in keggdblist:
            if outside_db not in outsidedb2:
                raise KEGGKeyError(key=outside_db)
        else:
            raise KEGGKeyError(key=kegg_db)

    url = "http://rest.kegg.jp/conv/{}/{}".format(target, source)
    
    if return_url == True:
        return url
    filename = target + "_" + source + "_conv"

    text = download_textfile(url, filename, force_download=force_download, verbose = verbose)

    codes1, codes2 = process_request_text(text, want_descr=True)

    return codes1, codes2


def keggapi_info(database, verbose=True, force_download=False, return_format = None, return_url = False):
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
    retutn_format : str
        optional, specify a return format to return, str | dict (default is None)
        
    Returns
    -------
    info_str : str
        optional, plain text response of API INFO command
    info_dict : dict
        optional, parsed response of API INFO as a dictionary
    """
    
    valid_return_formats = (None, "str", "dict")
    
    if return_format not in valid_return_formats:
        raise ValueError("invalid {} format for keggapi_info return".format(return_format))
        
    org = get_organism_codes()

    if database not in db_categories + org:
        raise KEGGKeyError(
            database, msg="source database {} is not a valid database".format(database)
        )

    url = "http://rest.kegg.jp/info/{}".format(database)
    
    if return_url == True:
        return url
    
    filename = database + "_info"

    infos = download_textfile(url, filename, verbose=False, force_download = force_download)
    
    if verbose == True:
        logging.info("Infos on %s from KEGG:\n",database)
    if return_format == None:
        if verbose == False:
             print("\n".join(infos.splitlines()[1:4]))
        else:
            print(infos)
    elif return_format == "str":
        return infos
    elif return_format == "dict":
        processed_dict = process_request_text(infos, mode = "columns")
        return processed_dict
        


def keggapi_ddi(dbentry, force_download=False, return_url = False):
    """KEGG REST API interface for the DDI command
    lists drug-drug interactions for a given compound name
    
    Parameters
    ----------
    
    dbentry : str
        drug KEGG database entry
    force_download : bools
        forces overwriting over cached files (default is False)
    return_url : bool, optional
        If True returns the interrogated URL        
    for further info read https://www.kegg.jp/kegg/rest/keggapi.html"""

    url = "http://rest.kegg.jp/ddi/{}".format(dbentry)
    
    if return_url == True:
        return url
    
    
    
    filename = dbentry + "_ddi"

    text = download_textfile(url, filename, verbose=False, force_download = force_download)

    ddi_list = []

    for line in text.splitlines():
        drug1, drug2, ddi_code, interaction = line.split("\t")
        ddi_list.append((drug1, drug2, ddi_code, interaction))

    return ddi_list


# =============================================================================
# KEGGAPI DERIVATE FUNCTIONS
# =============================================================================
def get_organism_codes(force_download=False):
    """Returns all KEGG Organism name codes """

    org_url = "http://rest.kegg.jp/list/organism"
    org_filename = "organism_code_list"

    org_codes = []

    if force_download:
        organism_fulltext = download_textfile(
            org_url, org_filename, verbose=False, force_download=force_download
        )

        for line in organism_fulltext.splitlines():
            _, kegg_code, _, _ = line.strip().split("\t")
            org_codes.append(kegg_code)

    else:
        org_codespath = RES_PATH.joinpath("org_codes.txt")
#        with open("./res/org_codes.txt", "r+") as org_file:
#            organism_fulltext = org_file.read()
        organism_fulltext = org_codespath.read_text()

        for line in organism_fulltext.splitlines():
            org_codes.append(line)

    return org_codes


def kegg_url(target_db, source_db):  # deprecated
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

def get_references(item):
    """ Get References
    For a given item tries to find academic references
    
    Parameters:
    -----------
    input : str
        requested object identifier
    
    Returns:
    -------
    referencelist : list
        list of references
    """
    
    get_dict = keggapi_get(dbentry = item, return_dict = True)
    
    dict_keys = [key for key in get_dict.keys() if "reference" in key]
    
    return [get_dict[key] for key in dict_keys]