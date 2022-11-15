from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace


def text_layout(prop, level, colour_dict=None):
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            
            prop_text = node.props.get(prop)
            if prop_text:
                if colour_dict:
                    prop_face = TextFace(prop_text, color=colour_dict[prop_text])
                else:
                    prop_face = TextFace(prop_text, color='blue')
            node.add_face(prop_face, column = level, position = "branch_right")
            node.sm_style["bgcolor"] = 'black' # highligh clade
            # while (node):
            #         node = node.up
            #         if node:
            #             node.sm_style["hz_line_width"] = 5
    return layout_fn
    return

def label_layout(prop, level, colour_dict=None):
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            prop_text = node.props.get(prop)
            if prop_text:
                if colour_dict:
                    node.sm_style["hz_line_color"] = colour_dict[prop_text]
                    node.sm_style["hz_line_width"] = 2
                    # while (node):
                    #     node = node.up
                    #     if node:
                    #         node.sm_style["hz_line_color"] = colour_dict[prop_text]
                    #         node.sm_style["hz_line_width"] = 2
            #node.sm_style["bgcolor"] = 'black' # highligh clade
    return layout_fn
    return

def rectangular_layout(prop, level, colour_dict=None):
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            prop_text = node.props.get(prop)
            if prop_text:
                if colour_dict:
                    label_rect = RectFace(width=50,height=50, color=colour_dict[prop_text], padding_x=1, padding_y=1)
                    node.add_face(label_rect, column = level,  position = 'branch_right')
    return layout_fn
    return