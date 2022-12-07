from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace

#paried_color = ["red", "darkblue", "darkgreen", "darkyellow", "violet", "mediumturquoise", "sienna", "lightCoral", "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal"]

def text_layout(prop, level, colour_dict=None, internal_rep='counter'):
    internal_prop = prop+'_'+internal_rep
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            
            prop_text = node.props.get(prop)
            if prop_text:
                if colour_dict:
                    prop_face = TextFace(prop_text, color=colour_dict[prop_text])
                else:
                    prop_face = TextFace(prop_text, color='blue')
            node.add_face(prop_face, column = level, position = "branch_right")

        elif node.is_leaf() and node.props.get(internal_prop):
            piechart_face = get_piechartface(node, internal_prop, colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
            node.add_face(piechart_face, column = level+5, position = "branch_right", collapsed_only=True)

        elif node.props.get(internal_prop):
            piechart_face = get_piechartface(node, internal_prop, colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
            node.add_face(piechart_face, column = level+5, position = "branch_right", collapsed_only=True)
        
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
        elif node.is_leaf() and node.props.get(internal_prop):
            piechart_face = get_piechartface(node, internal_prop, colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
            node.add_face(piechart_face, column = level+5, position = "branch_right", collapsed_only=True)

        elif node.props.get(internal_prop):
            piechart_face = get_piechartface(node, internal_prop, colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
    return layout_fn
    return

def rectangular_layout(prop, level, colour_dict=None, internal_rep='counter'):
    internal_prop = prop+'_'+internal_rep
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            prop_text = node.props.get(prop)
            if prop_text:
                if colour_dict:
                    label_rect = RectFace(width=50,height=50, color=colour_dict[prop_text], padding_x=1, padding_y=1)
                    node.add_face(label_rect, column = level,  position = 'aligned')
        elif node.is_leaf() and node.props.get(internal_prop):
            piechart_face = get_piechartface(node, internal_prop, colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
            node.add_face(piechart_face, column = level+5, position = "branch_right", collapsed_only=True)
            
        elif node.props.get(internal_prop):
            piechart_face = get_piechartface(node, internal_prop, colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
    return layout_fn
    return

def get_piechartface(node, prop, colour_dict=None):
    piechart_data = []
    counter_props = node.props.get(prop).split('||')
    for counter_prop in counter_props:
        k, v = counter_prop.split('--')
        piechart_data.append([k,float(v),colour_dict[k],None])

    if piechart_data:
        piechart_face = PieChartFace(radius=50, data=piechart_data)
        return piechart_face
    else:
        return None