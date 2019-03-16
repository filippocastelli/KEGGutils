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


class MissingNodetypeError(KeggUtilsGraphException):
    pass


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


class KEGGKeyError(KEGGOnlineError):
    def __init__(self, key, msg=None):
        self.key = key

        if msg is None:
            msg = "Invalid KEGG key: {}".format(self.key)
        super(KEGGKeyError, self).__init__(msg)
        
class KEGGInvalidFileContent(KEGGOnlineError):
    def __init__(self, file, content, msg=None):
        
        self.file = file
        self.content = content

        if msg is None:
            msg = "Invalid content {} in file {}".format(self.content, self.file)
            
        super(KEGGInvalidFileContent, self).__init__(msg)
        
