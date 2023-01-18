from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, LegendFace

paried_color = ["red", "darkblue", "darkgreen", "darkyellow", "violet", "mediumturquoise", "sienna", "lightCoral", "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal"]

#collapse in layout
#kingdom, phylum, class, order, family, genus, species, subspecies
def get_level(node, level=0):
    if node.is_root():
        return level
    else:
        return get_level(node.up, level + 1)

class TaxaClade(TreeLayout):
    def __init__(self, name, level, rank, color_dict):
        super().__init__(name, aligned_faces=True)

        self.activate = False
        self.name = name
        self.column = level
        self.rank = rank
        self.color_dict = color_dict
    
    def set_node_style(self, node):
        if not node.is_root() and node.props.get('rank') == self.rank:
            if node.props.get('sci_name'):
                #text_face = TextFace(node.props.get('sci_name'), color='black')
                #face_name = OutlineFace(node.props.get('sci_name'), collapsing_height= float("inf"))
                node.sm_style["bgcolor"] = self.color_dict[node.props.get('sci_name')] # highligh clade
                #node.sm_style["draw_descendants"] = False
                #node.add_face(text_face, column = self.column, position = "aligned")
                #node.add_face(text_face, column = self.column, position = "aligned", collapsed_only=True)
                #node.add_face(face_name, column = 5, position = 'branch_right', collapsed_only=True)

class TaxaRectangular(TreeLayout):
    def __init__(self, name="Last common ancestor",
            rect_width=15, column=0):
        super().__init__(name, aligned_faces=True)

        self.active = True

        self.rect_width = rect_width
        self.column = column

    def set_node_style(self, node):
        if node.props.get('sci_name'):
            lca = node.props.get('sci_name')
            color = node.props.get('sci_name_color', 'lightgray')
            
            level = get_level(node, level=self.column)
            # lca_face = RectFace(self.rect_width, float('inf'), 
            #         color = color, 
            #         text = lca,
            #         fgcolor = "white",
            #         padding_x = 1, padding_y = 1)
            lca_face = RectFace(self.rect_width, float('inf'), text = lca, color=color, padding_x=1, padding_y=1)
            lca_face.rotate_text = True
            node.add_face(lca_face, position='aligned', column=level)
            node.add_face(lca_face, position='aligned', column=level,
                collapsed_only=True)

# def taxa_layout(rank, color_dict=None):
#     def layout_fn(node):
#         if not node.is_root() and node.props.get('rank') == rank:
#             node.sm_style["bgcolor"] = color_dict[node.props.get('sci_name')] # highligh clade

#             #print(node.props.get('sci_name'))
#             #node.sm_style["hz_line_color"] = paried_color[count]
#             # node.sm_style["hz_line_width"] = 2
            
#             # children = node.children
#             # if children:
#             #     for child in children:
#             #         if child:
#             #             child.sm_style["hz_line_color"] = paried_color[count]
#             # count += 1
#     return layout_fn
#     return

def collapse_kingdom():
    def layout_fn(node):
        color = "red"
        if not node.is_root() and  node.props.get('rank') == 'superkingdom':
            face_name = TextFace(node.props.get('sci_name'), color=color)
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = color
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "level1_kingdom"
    return layout_fn
    return

def collapse_phylum():
    def layout_fn(node):
        color="orange"
        if not node.is_root() and  node.props.get('rank') == 'phylum':
            face_name = TextFace(node.props.get('sci_name'), color=color)
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = color
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "level2_phylum"
    return layout_fn
    return

def collapse_class():
    def layout_fn(node):
        color="yellow"
        if not node.is_root() and  node.props.get('rank') == 'class':
            face_name = TextFace(node.props.get('sci_name'), color=color)
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = color
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "level3_class"
    return layout_fn
    return

def collapse_order():
    def layout_fn(node):
        color="green"
        if not node.is_root() and  node.props.get('rank') == 'order':
            face_name = TextFace(node.props.get('sci_name'), color=color)
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = color
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "level4_order"
    return layout_fn
    return 

def collapse_family():
    def layout_fn(node):
        color="blue"
        if not node.is_root() and  node.props.get('rank') == 'family':
            face_name = TextFace(node.props.get('sci_name'), color=color)
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = color
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "level5_family"
    return layout_fn
    return 

def collapse_genus():
    def layout_fn(node):
        color="indigo"
        if not node.is_root() and  node.props.get('rank') == 'genus':
            face_name = TextFace(node.props.get('sci_name'), color="indigo")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "indigo"
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "level6_genus"
    return layout_fn
    return

def collapse_species():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'species':
            face_name = TextFace(node.props.get('sci_name'), color="violet")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "violet"
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "level7_species"
    return layout_fn
    return