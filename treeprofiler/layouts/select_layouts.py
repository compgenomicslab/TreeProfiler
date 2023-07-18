from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import (RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, \
                            SelectedFace, SelectedCircleFace, SelectedRectFace)

# branch thicken, background highlighted to purple
def select_layout(prop, level, color='purple'):
    def layout_fn(node):
        if node.is_leaf and node.props.get(prop):
            
            prop_text = node.props.get(prop)
            #prop_face = SelectedCircleFace(name='prop')
            node.add_face(prop_face, column=level, position = "branch_right")

    return layout_fn
    return     

# hightlighted as rectangular

# hightlighted as circle