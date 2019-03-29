# =============================================================================
# ERRORS AND EXCEPTIONS
# =============================================================================

class KeggUtilsGraphException(Exception):
    def __init__(self, graph, msg=None):
        self.graph = graph
        if msg is None:
            msg = "There's a problem with graph: {}".format(graph.name)

        super(KeggUtilsGraphException, self).__init__(msg)


class NotAKeggGraphError(KeggUtilsGraphException):
    pass


class MissingNodetypeError(Exception):
    def __init__(self, nodetype, graph, msg=None):
        self.graph = graph
        self.nodetype = nodetype
        
        if msg is None:
            msg = "Nodetype {} is missing in graph {}".format(self.nodetype, self.graph.name)
            
        super().__init__(graph, msg)



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
        
class KEGGInvalidFileContent(Exception):
    def __init__(self, file, content, msg=None):
        
        self.file = file
        self.content = content

        if msg is None:
            msg = "Invalid content {} in file {}".format(self.content, self.file)
            
        super(KEGGInvalidFileContent, self).__init__(msg)

class KEGGDataBaseError(Exception):
    def __init__(self, db, msg=None):
        
        self.db = db


        if msg is None:
            msg = "Invalid KEGG database {}".format(self.db)
            
        super(KEGGDataBaseError, self).__init__(msg)
        
class KEGGInvalidContent(Exception):
    def __init__(self, content, msg=None):
        
        self.content = content

        if msg is None:
            msg = "Invalid content in {}".format(self.content)
            
        super(KEGGInvalidFileContent, self).__init__(msg)
        
class KGMLerror(Exception):
    def __init__(self, xml_file = None, tree = None, msg=None):
            
        self.xml_file = xml_file
        self.tree = tree

        if msg is None:
            msg = "Invalid tree {} in xml file {}".format(self.tree, self.xml_file)
            
        super(KGMLerror, self).__init__(msg)
        
class KEGGgraphError(Exception):
    def __init__(self, graph = None, msg = None):
        self.graph = graph
        
        if msg is None:
            msg = "Invalid KEGGgraph {}".format(self.graph)
            
        super(KEGGgraphError, self).__init__(msg)
        
    
class KEGGChainError(Exception):
    def __init__(self, chain, msg=None):
        self.chain = chain

        if msg is None:
            msg = "Invalid KEGG chain: {}".format(self.chain)
        super(KEGGChainError, self).__init__(msg)