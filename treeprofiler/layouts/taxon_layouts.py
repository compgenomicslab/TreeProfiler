from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, LegendFace
from collections import  OrderedDict
from treeprofiler.src import utils 
paried_color = ["red", "darkblue", "darkgreen", "darkyellow", "violet", "mediumturquoise", "sienna", "lightCoral", "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal"]

#collapse in layout
#kingdom, phylum, class, order, family, genus, species, subspecies
def get_level(node, level=0):
    if node.is_root:
        return level
    else:
        return get_level(node.up, level + 1)

def summary(nodes):
    "Return a list of names summarizing the given list of nodes"
    return list(OrderedDict((first_name(node), None) for node in nodes).keys())

def first_name(tree):
    "Return the name of the first node that has a name"
    
    sci_names = []
    for node in tree.traverse('preorder'):
        if node.is_leaf:
            sci_name = node.props.get('sci_name')
            sci_names.append(sci_name)

    return next(iter(sci_names))

class TaxaClade(TreeLayout):
    def __init__(self, name, level, rank, color_dict, legend=True):
        super().__init__(name, aligned_faces=True)

        self.activate = False
        self.name = name
        self.column = level
        self.rank = rank
        self.color_dict = color_dict
        self.legend = legend

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.rank,
                                    variable='discrete',
                                    colormap=self.color_dict,
                                    )

    # def set_node_style(self, node):
    #     if not node.is_root and node.props.get('rank') == self.rank:
    #         if node.props.get('sci_name'):
    #             node.sm_style["bgcolor"] = self.color_dict[node.props.get('sci_name')] # highligh clade

    def set_node_style(self, node):
        named_lineage = node.props.get('named_lineage', None)
        if named_lineage:
            for clade, color in self.color_dict.items():
                if clade in named_lineage:
                    node.sm_style["hz_line_color"] = color
                    node.sm_style["hz_line_width"] = 2
                    node.sm_style["vt_line_color"] = color
                    node.sm_style["vt_line_width"] = 2
                    #node.sm_style["draw_descendants"] = False
                    node.sm_style["outline_color"] = color
                    break

        if not node.is_leaf:
            # Collapsed face
            if node.props.get('rank') == self.rank:
                if node.props.get('sci_name') is not None:
                    sci_name = node.props.get('sci_name')
                    color = self.color_dict.get(sci_name, 'gray')
                    node.add_face(TextFace(sci_name, padding_x=2, color = color),
                        position="branch_right", column=1, collapsed_only=True)
            

        # if not node.is_root and node.props.get('rank') == self.rank: 
        #     if node.props.get('sci_name'):
        #         color = self.color_dict[node.props.get('sci_name')]
        #         node.sm_style["hz_line_color"] = color
        #         node.sm_style["hz_line_width"] = 2
        #         node.sm_style["vt_line_color"] = color
        #         node.sm_style["vt_line_width"] = 2
        #         #node.sm_style["draw_descendants"] = False
        #         node.sm_style["outline_color"] = color
                

            

class LayoutSciName(TreeLayout):
    def __init__(self, name="Scientific name", color_dict={}):
        super().__init__(name, aligned_faces=True)
        self.color_dict = color_dict

    def set_node_style(self, node):
        if node.is_leaf:
            sci_name = node.props.get('sci_name')
            prot_id = node.name

            rank_colordict = self.color_dict.get(node.props.get('rank'),'')
            if rank_colordict:
                color = rank_colordict.get(sci_name, 'gray')
            else:
                color = 'gray'
            node.add_face(TextFace(sci_name, color = color, padding_x=2, min_fsize=4, max_fsize=25),
                column=0, position="branch_right")

            if len(prot_id) > 40:
                prot_id = prot_id[0:37] + " ..."
           
            #node.add_face(TextFace(prot_id, color = 'Gray', padding_x=2), column = 2, position = "aligned")
        else:
            # Collapsed face
            names = summary(node.children)
            texts = names if len(names) < 6 else (names[:3] + ['...'] + names[-2:])
            for i, text in enumerate(texts):
                sci_name = node.props.get('sci_name')
                rank_colordict = self.color_dict.get(node.props.get('rank'),'')
                if rank_colordict:
                    color = rank_colordict.get(sci_name, 'gray')
                else:
                    color = 'gray'
                node.add_face(TextFace(text, padding_x=2, color = color, min_fsize=4, max_fsize=25),
                        position="branch_right", column=1, collapsed_only=True)
                        
class TaxaRectangular(TreeLayout):
    def __init__(self, name="Last common ancestor", rank=None, color_dict={}, rect_width=20, column=0, legend=True ):
        super().__init__(name, aligned_faces=True)

        self.active = True
        self.rank = rank
        self.color_dict = color_dict
        self.rect_width = rect_width
        self.column = column

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title='TaxaRectangular_'+self.rank,
                                    variable='discrete',
                                    colormap=self.color_dict,
                                    )

    def set_node_style(self, node):
        node_rank = node.props.get('rank')
        node_sciname = node.props.get('sci_name')
        if node_sciname and (node_rank == self.rank):
            lca = node_sciname
            color = self.color_dict.get(lca, 'lightgray')
            level = get_level(node, level=self.column)
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node_sciname}</b><br>'
            if lca:
                tooltip += f'rank: {node_rank}<br>'
                tooltip += f'sci_name: {lca}<br>'

            lca_face = RectFace(self.rect_width, None, text = lca, color=color, padding_x=1, padding_y=1, tooltip=tooltip)
            lca_face.rotate_text = True
            node.add_face(lca_face, position='aligned', column=level)
            node.add_face(lca_face, position='aligned', column=level,
                collapsed_only=True)
        # else:
        #     named_lineage = node.props.get('named_lineage', None)
        #     lca = next((elem for elem in named_lineage if elem in self.taxa_list), None)
        #     if lca:
        #         color = self.color_dict.get(lca, 'lightgray')
        #         #level = get_level(node, level=self.column)
        #         level = self.column
        #         tooltip = ""
        #         if node.name:
        #             tooltip += f'<b>{node.name}</b><br>'
        #         if lca:
        #             tooltip += f'rank: {node_rank}<br>'
        #             tooltip += f'sci_name: {lca}<br>'

        #         lca_face = RectFace(self.rect_width, None, text = lca, color=color, padding_x=1, padding_y=1, tooltip=tooltip)
        #         lca_face.rotate_text = True
        #         node.add_face(lca_face, position='aligned', column=level)
        #         node.add_face(lca_face, position='aligned', column=level,
        #             collapsed_only=True)

class TaxaCollapse(TreeLayout):
    def __init__(self, name="Last common ancestor", rank=None, color_dict={}, rect_width=20, column=0, padding_x=1, padding_y=0, legend=True ):
        super().__init__(name, aligned_faces=True)

        self.active = True
        self.rank = rank
        self.color_dict=color_dict
        self.rect_width = rect_width
        self.column = column
        self.padding_x = padding_x
        self.padding_y = padding_y

        self.taxa_list =  list(color_dict.keys())

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        super().set_tree_style(tree, tree_style)
        text = TextFace(" ", min_fsize=10, max_fsize=15, padding_x=self.padding_x, width=self.rect_width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title='TaxaRectangular_'+self.rank,
                                    variable='discrete',
                                    colormap=self.color_dict,
                                    )
    
    def set_node_style(self, node):
        node_rank = node.props.get('rank')
        node_sciname = node.props.get('sci_name')
        named_lineage = node.props.get('named_lineage', None)
        #lca = next((elem for elem in named_lineage if elem in self.taxa_list), None)
        lca_value = node.props.get('lca')
        if lca_value:
            lca_dict = utils.string_to_dict(lca_value)
            lca = lca_dict.get(self.rank, None)
            if lca:
                color = self.color_dict.get(lca, 'lightgray')
                #level = get_level(node, level=self.column)
                level = self.column
                tooltip = ""
                if node.name:
                    tooltip += f'<b>{lca}</b><br>'
                if lca:
                    tooltip += f'rank: {self.rank}<br>'
                    tooltip += f'sci_name: {lca}<br>'

                lca_face = RectFace(self.rect_width, None, text = lca, color=color, padding_x=1, padding_y=1, tooltip=tooltip)
                lca_face.rotate_text = True
                node.sm_style["draw_descendants"] = False
                
                node.add_face(lca_face, position='aligned', column=level)
                node.add_face(lca_face, position='aligned', column=level,
                    collapsed_only=True)
                node.add_face(TextFace(lca, color = color, padding_x=2),
                column=1, position="branch_right", collapsed_only=True)

            # elif node_sciname and (node_rank == self.rank):
            #     lca = node_sciname
            #     color = self.color_dict.get(lca, 'lightgray')
            #     level = self.column
            #     tooltip = ""
            #     if node.name:
            #         tooltip += f'<b>{node.name}</b><br>'
            #     if lca:
            #         tooltip += f'rank: {node_rank}<br>'
            #         tooltip += f'sci_name: {lca}<br>'

            #     lca_face = RectFace(self.rect_width, None, text = lca, color=color, padding_x=1, padding_y=1, tooltip=tooltip)
            #     lca_face.rotate_text = True
            #     node.sm_style["draw_descendants"] = False
            #     node.add_face(lca_face, position='aligned', column=level)
            #     node.add_face(lca_face, position='aligned', column=level,
            #         collapsed_only=True)

class LayoutEvolEvents(TreeLayout):
    def __init__(self, name="Evolutionary events", 
            prop="evoltype",
            speciation_color="blue", 
            duplication_color="red", node_size = 2,
            legend=True):
        super().__init__(name)
        
        self.prop = prop
        self.speciation_color = speciation_color
        self.duplication_color = duplication_color
        self.node_size = node_size
        self.legend = legend

        self.active = True

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            colormap = { "Speciation event": self.speciation_color,
                         "Duplication event": self.duplication_color }
            tree_style.add_legend(title=self.name, 
                    variable="discrete",
                    colormap=colormap)

    def set_node_style(self, node):
        if not node.is_leaf:
            if node.props.get(self.prop, "") == "S":
                node.sm_style["fgcolor"] = self.speciation_color
                node.sm_style["size"] = self.node_size

            elif node.props.get(self.prop, "") == "D":
                node.sm_style["fgcolor"] = self.duplication_color
                node.sm_style["size"] = self.node_size