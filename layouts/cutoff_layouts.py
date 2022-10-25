from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace

#collapse in layout
#kingdom, phylum, class, order, family, genus, species, subspecies

def collapse_cutoff(prop, cutoff):
    def layout_fn(node):
        if node.props.get(prop):
            if not node.is_root() and float(node.props.get(prop)) < float(cutoff):
                #node = node.up
                #face_name = TextFace(node.props.get('sci_name'), color="red")
                node.sm_style["draw_descendants"] = False
                node.sm_style["outline_color"] = "red"
                #node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_" + prop
    return layout_fn
    return