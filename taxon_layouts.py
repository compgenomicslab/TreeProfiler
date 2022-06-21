from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace

#collapse in layout
#kingdom, phylum, class, order, family, genus, species, subspecies

def collapse_kingdom():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'superkingdom':
            face_name = TextFace(node.props.get('sci_name'), color="red")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "red"
            node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_kingdom"
    return layout_fn
    return

def collapse_phylum():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'phylum':
            face_name = TextFace(node.props.get('sci_name'), color="orange")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "orange"
            node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_phylum"
    return layout_fn
    return

def collapse_class():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'class':
            face_name = TextFace(node.props.get('sci_name'), color="yellow")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "yellow"
            node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_class"
    return layout_fn
    return

def collapse_order():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'order':
            face_name = TextFace(node.props.get('sci_name'), color="green")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "green"
            node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_order"
    return layout_fn
    return 

def collapse_family():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'family':
            face_name = TextFace(node.props.get('sci_name'), color="blue")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "blue"
            node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_genus"
    return layout_fn
    return 

def collapse_genus():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'genus':
            face_name = TextFace(node.props.get('sci_name'), color="indigo")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "indigo"
            node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_genus"
    return layout_fn
    return

def collapse_species():
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'species':
            face_name = TextFace(node.props.get('sci_name'), color="violet")
            node.sm_style["draw_descendants"] = False
            node.sm_style["outline_color"] = "violet"
            node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_species"
    return layout_fn
    return