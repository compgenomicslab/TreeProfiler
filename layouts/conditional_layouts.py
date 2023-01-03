from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import (RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, \
                            SelectedFace, SelectedCircleFace, SelectedRectFace)
from utils import to_code, call, counter_call

# branch thicken, background highlighted to purple
def highlight_layout(conditions, level, prop2type={}, color='purple'):
    conditional_output = to_code(conditions)
    def layout_fn(node):
        final_call = False
        for condition in conditional_output:
            #normal
            
            op = condition[1]
            if op == 'in':
                value = condition[0]
                prop = condition[2]
                datatype = prop2type[prop]
                final_call = call(node, prop, datatype, op, value)

            elif ':' in condition[0] :
                internal_prop, leaf_prop = condition[0].split(':')
                value = condition[2]
                datatype = prop2type[internal_prop]
                final_call = counter_call(node, internal_prop, leaf_prop, datatype, op, value)
            else:
                prop = condition[0]
                value = condition[2]
                datatype = prop2type[prop]
                final_call = call(node, prop, datatype, op, value)
            
            if final_call == False:
                break
            else:
                continue
        
        if final_call:
            
            prop_face = SelectedRectFace(name='prop')
            node.sm_style["bgcolor"] = color # highligh clade
            #node.sm_style["hz_line_width"] = 5
            #node.add_face(prop_face, column=level, position = "branch_right")
            while (node):
                node = node.up
                if node:
                    node.sm_style["hz_line_width"] = 5
                    #node.sm_style["hz_line_color"] = color
                    #node.add_face(OutlineFace(color=color), column=level, collapsed_only=True)
    return layout_fn
    return     

# conditional collapse layouts
def collapsed_by_layout(conditions, level, prop2type={}, color='red'):
    conditional_output = to_code(conditions)
    def layout_fn(node):
        final_call = False
        for condition in conditional_output:
            #normal
            
            op = condition[1]
            if op == 'in':
                value = condition[0]
                prop = condition[2]
                datatype = prop2type[prop]
                final_call = call(node, prop, datatype, op, value)

            elif ':' in condition[0] :
                internal_prop, leaf_prop = condition[0].split(':')
                value = condition[2]
                datatype = prop2type[internal_prop]
                final_call = counter_call(node, internal_prop, leaf_prop, datatype, op, value)
            else:
                prop = condition[0]
                value = condition[2]
                datatype = prop2type[prop]
                final_call = call(node, prop, datatype, op, value)
            
            if final_call == False:
                break
            else:
                continue
        if final_call:
            if not node.is_root():
                node.sm_style["draw_descendants"] = False
                node.sm_style["outline_color"] = color
    return layout_fn
    return

# for boolean layouts
from distutils.util import strtobool
from utils import check_nan

class LayoutBinary(TreeLayout):
    def __init__(self, name, level, color, prop_colour_dict, bool_prop, reverse=False):
        super().__init__(name)
        self.aligned_faces = True
        self.bool_prop = bool_prop
        self.column = level
        self.color = color
        self.prop_colour_dict = prop_colour_dict
        self.internal_prop = bool_prop+'_counter'
        self.reverse = reverse

    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.bool_prop):
            prop_bool = node.props.get(self.bool_prop)
            if not check_nan(prop_bool):
                str2bool = strtobool(prop_bool)
                if self.reverse:
                    if not bool(str2bool):
                        prop_face = CircleFace(radius=200, color=self.color, padding_x=1, padding_y=1)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                    else:
                        prop_face = CircleFace(radius=200, color='white', padding_x=1, padding_y=1)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                else:
                    if bool(str2bool):
                        prop_face = CircleFace(radius=200, color=self.color, padding_x=1, padding_y=1)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                    else:
                        prop_face = CircleFace(radius=200, color='white', padding_x=1, padding_y=1)
                        node.add_face(prop_face, column=self.column, position = "aligned")
            else:
                prop_face = TextFace('NaN', color='black')
                node.add_face(prop_face, column=self.column, position = "aligned")
        elif node.is_leaf() and node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.prop_colour_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column+2, position = "aligned", collapsed_only=True)

        elif node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.prop_colour_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column+2, position = "aligned", collapsed_only=True)

def boolean_layout(prop, level, color, prop_colour_dict, internal_rep='counter', reverse=False):
    internal_prop = prop+'_'+internal_rep
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            prop_text = node.props.get(prop)
            
            if reverse:
                if not bool(strtobool(prop_text)):
                    prop_face = CircleFace(radius=100, color=color, padding_x=1, padding_y=1)
                    node.add_face(prop_face, column=level, position = "branch_right")
                else:
                    prop_face = CircleFace(radius=100, color='white', padding_x=1, padding_y=1)
                    node.add_face(prop_face, column=level, position = "branch_right")
            else:
                if bool(strtobool(prop_text)):
                    # lca_face = RectFace(15, float('inf'), 
                    # color = color, 
                    # #text = lca,
                    # fgcolor = "white",
                    # padding_x = 1, padding_y = 1)
                    prop_face = CircleFace(radius=100, color=color)
                    node.add_face(prop_face, column=level, position = "aligned")
                else:
                    prop_face = CircleFace(radius=100, color='white')
                    node.add_face(prop_face, column=level, position = "aligned")
        elif node.is_leaf() and node.props.get(internal_prop):
            piechart_face = get_piechartface(node, internal_prop, prop_colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
            node.add_face(piechart_face, column = level+5, position = "branch_right", collapsed_only=True)

        elif node.props.get(internal_prop):
            
            piechart_face = get_piechartface(node, internal_prop, prop_colour_dict)
            node.add_face(piechart_face, column = level, position = "branch_top")
            node.add_face(piechart_face, column = level+5, position = "branch_right", collapsed_only=True)
    
    layout_fn.aligned_faces = True
    return layout_fn


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

# hightlighted as rectangular

# hightlighted as circle