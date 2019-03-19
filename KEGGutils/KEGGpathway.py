import networkx as nx
import xml.etree.ElementTree as et
import logging
import matplotlib.pylab as plt

from KEGGutils import draw
from KEGGutils.KEGGapi import keggapi_get


class KEGGpathway(nx.DiGraph):
    """KEGG Pathway:
    Extends the NetworkX.DiGraph(), includes a parser to read KGML .xml files 
    
    Parameters
    ----------
    kgml_file : str
        location of a KGML .xml file to parse

    Properties
    ----------
    title : str
        title of the pathway
    labels: dict
        dictionary of node labels
    idcode: str
        id code of the pathway
    tree: XML ElementTree
        element tree for the KGML file
    pos: dict
        dictionary of graph node positions for plotting
    nodedict: dict
        dictionary of ENTRY occurrences, ordered by corresponding graph node
    kgml_file: str
        path to KGML .xml file to parse
    kegg_image: img
        pathway image downloaded from KEGG
    reactions: dict
        dictionary of reactions in the pathway, ordered by corresponding graph node
    genes: dict
        dictionary of genes in the pathway, ordered by corresponding graph node
    maps: dict
        dictionary of maps in the pathway, ordered by corresponding graph node
    orthologs: dict
        dictionary of orthologs in the pathway, ordered by corresponding graph node
    enzymes: dict
        dictionary of enzymes in the pathway, ordered by corresponding graph node
    groups: dict
        dictionary of groups in the pathway, ordered by corresponding graph node
    brite_nodes: dict
        dictionary of brite-style nodes in the pathway, ordered by corresponding graph node
    other_nodes: dict
        dictionary of unclassified nodes in the pathway, ordered by corresponding graph node
    relations: dict
        dictionary of relations in the pathway, ordered by corresponding graph edge

    elements: dict
        dictionary of the previous dictionaries ordered by nodetype
    element_keys: list
        keys for elements
    
    Methods:
    --------

    KEGGpathway.parse_kgml(kgml_file):
        parses a kgml file
    KEGGpathway.draw():
        Plots the pathway with networkx.draw_networkx() (note: KGML file format is missing lot of additional information, for visualization purposes is reccomended to use downloaded kegg images)
    KEGGpathway.download_img():
        Downloads the pathway visualization from KEGG and returns it
    KEGGpathway.list_by_nodetype(nodetype):
        Returns a list of corresponding nodetypes in the pathway.
    
    """

    title = ""
    labels = {}
    reactions = {}
    relations = {}
    genes = {}
    orthologs = {}
    maps = {}
    enzymes = {}
    groups = {}
    brite_nodes = {}
    other_nodes = {}

    elements = {
        "gene": genes,
        "reaction": reactions,
        "relation": relations,
        "ortholog": orthologs,
        "map": maps,
        "enzyme": enzymes,
        "group": groups,
        "brite": brite_nodes,
    }

    element_keys = list(elements.keys())

    idcode = ""
    pos = dict()
    nodedict = {}
    tree = None
    kgml_file = None
    kegg_image = []

    def __init__(self, *args, **kwargs):

        #        self.kgml_file, args = self._find_arg_and_kick(".xml")

        if "kgml_file" in kwargs:
            self.kgml_file = kwargs.pop("kgml_file")

        super().__init__(self, *args, **kwargs)

        if self.kgml_file is not None:
            self.parse_kgml(self.kgml_file)

    # =============================================================================
    # PUBLIC METHODS
    # =============================================================================
    def calc_pos(self):
        """Infers the position layout of the nodes
        Returns
        -------
        self.pos: dict
            dictionary of node positions ordered by graph node
        """

        for n in self.nodes():
            self.pos.update({n: self.node[n]["xy"]})

        return self.pos

    def parse_kgml(self, kgml_file):
        """parse_kgml parses a kgml .xml file
        
        Parameters
        ----------
        kgml_file : str
            path to .xml file
        
        """

        tree = et.parse(kgml_file)
        self.tree = tree

        self.title = tree.getroot().get("title")
        self.name = tree.getroot().get("name")
        self.idcode = tree.getroot().get("id")

        for entry in tree.getiterator("entry"):
            self._parse_entry(entry)
        for relation in tree.getiterator("relation"):
            self._parse_relation(relation)
        for reaction in tree.getiterator("reaction"):
            self._parse_reaction(reaction)

        _ = self.calc_pos()

    def draw(self):
        """draws the pathway graph with networkx.networkx_draw()
        """

        draw(graph=self, title=self.title, pos=self.pos)

    def download_img(self):
        """downloads the KEGG picture of the pathway
        Returns
        -------
        img: img
            KEGG picture of the pathway
        """

        self.kegg_image = keggapi_get(
            dbentry=self.name, option="image", show_result_image=False
        )
        plt.figure()
        plt.title(self.title)
        plt.imshow(self.kegg_image)

    def list_by_nodetype(self, nodetype):
        """Gets a list of nodetypes in the graph
        
        Parameters
        ----------
        nodetype : str
            desired nodetype
        
        Returns
        -------
        element_list
            list of elements in the graph by given nodetype
        """

        elementlist = []
        for item in list(self.elements[nodetype].values()):
            elementlist.append(item[nodetype])

        # remove duplicates
        elementlist = list(dict.fromkeys(elementlist))

        return elementlist

    # =============================================================================
    # PRIVATE METHODS
    # =============================================================================
    def _find_arg_and_kick(self, args, to_find):
        """finds an arg in args using to_find as a substring
        
        Parameters
        ----------
        args : list
            list of args
        to_find : str
            substring to identify the argument
        
        Returns
        -------
        arg: str
            found argument
        newargs: list
            list of arguments without arg
        """

        for arg in args:
            if type(arg) == str:
                if arg.find(to_find) != -1:
                    logging.debug("found {} in args".format(arg))

                    newargs = [argument for argument in args if argument != arg]
                    args = newargs

                    return arg, newargs

    def _parse_graphics(self, graphics):
        """parses graphics
        
        Parameters
        ----------
        graphics : xml entry
            xml entry for graphics, has to be an "entry" child
        
        Returns
        -------
        name : str
            name of the graphic entry
        x : int
            x coordinate
        y : int
            y coordinate
        
        """

        x = int(graphics.get("x"))
        y = -int(graphics.get("y"))
        name = graphics.get("name")

        return name, x, y

    def _parse_substrate_product(self, substrate):
        """parses substrate and product entries
        
        Parameters
        ----------
        substrate ( or product ) : xml entry
            has to be a child of "reaction"
        
        Returns
        -------
       id_ : int
            ID of substrate or product
        name : str
            name of substrate or product entry
        """

        id_ = int(substrate.get("id"))
        name = substrate.get("name")

        return id_, name

    def _parse_entry(self, entry):
        """parses ENTRY xml Elements
        
        Parameters
        ----------
        entry : xml element
            xml ENTRY type element, see KEGG KGML documentation
        """

        g_name = None
        g_x = None
        g_y = None
        component_id = None

        node_id = entry.get("id")
        node_name = entry.get("name")
        node_type = entry.get("type")
        node_kegglink = entry.get("link")  # implied
        #    node_reaction = entry.get('reaction') #implied

        for child in entry.getchildren():
            if child.tag == "graphics":
                # this should happen only once
                g_name, g_x, g_y = self._parse_graphics(child)
            if child.tag == "component":
                component_id = child.get("component")

        node_title = g_name if g_name is not None else node_name
        self.labels[node_id] = node_title

        xy = (g_x, g_y) if ((g_x is not None) and (g_y is not None)) else None

        nodes_to_add = node_name.split()
        for i, node in enumerate(nodes_to_add):
            node_index = node_id
            if i > 1:
                node_index = node_id + "_" + str(i - 1)

            xy = (
                (g_x - 10 * i, 10 * i + g_y)
                if ((g_x is not None) and (g_y is not None))
                else None
            )
            self.add_node(
                node_index,
                name=node,
                label=node_title,
                nodetype=node_type,
                xy=xy,
                kegglink=node_kegglink,
                component_id=component_id,
            )

            self.nodedict[node_index] = (node, node_title, node_type)

            for element in self.elements.keys():
                if node_type == element:
                    self.elements[element][node_index] = {
                        element: node,
                        "description": node_title,
                        "nodetype": node_type,
                        "kegglink": node_kegglink,
                    }

    def _parse_reaction(self, reaction):
        """Parses REACTION xml elements
        
        Parameters
        ----------
        reaction : xml element
            see KEGG KGML documentation
        
        """

        reaction_id = reaction.get("id")
        reaction_name = reaction.get("name")
        reaction_type = reaction.get("type")

        reaction_alt = None
        substrates = []
        products = []
        alts = []

        for child in reaction.getchildren():
            if child.tag == "substrate":
                substrate_id, substrate_name = self._parse_substrate_product(child)
                substrates.append((substrate_id, substrate_name))
            if child.tag == "product":
                product_id, product_name = self._parse_substrate_product(child)
                products.append((product_id, product_name))

            if (reaction_alt is None) and (child.find("alt") is not None):
                reaction_alt = child.find("alt").get("name")
                alts.append(reaction_alt)

        self.reactions[reaction_id] = {
            "id": reaction_id,
            "name": reaction_name,
            "type": reaction_type,
            "substrates": substrates,
            "products": products,
        }

    def _parse_relation(self, relation):
        """Parses RELATION xml elements
        
        Parameters
        ----------
        relation : xml element
            see KEGG KGML documentation
        
        """
        relation_type = relation.get("type")
        relation_entry1 = relation.get("entry1")
        relation_entry2 = relation.get("entry2")
        subtypes = []

        for subtype in relation.getchildren():
            subtype_name = subtype.get("name")

            if subtype_name in ("compound", "hidden compound"):
                subtype_value = int(subtype.get("value"))
            else:
                subtype_value = subtype.get("value")
            subtypes.append((subtype_name, subtype_value))

        all_entry1 = [
            i for i in list(self.nodes) if i.startswith(relation_entry1 + "_")
        ]
        all_entry1.insert(0, relation_entry1)
        all_entry2 = [
            i for i in list(self.nodes) if i.startswith(relation_entry2 + "_")
        ]
        all_entry2.insert(0, relation_entry2)

        for e1 in all_entry1:
            for e2 in all_entry2:
                self.add_edge(e1, e2, relation_type=relation_type, subtypes=subtypes)

                self.relations[e1 + "to" + e2] = {
                    "nodes": (self.node[e1]["name"], self.node[e2]["name"]),
                    "relation_type": relation_type,
                    "subtypes": subtypes,
                }
