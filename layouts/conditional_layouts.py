from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import (RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, \
                            SelectedFace, SelectedCircleFace, SelectedRectFace)
import re
import operator

operator_dict = {
                '<':operator.lt,
                '<=':operator.le,
                '=':operator.eq,
                '!=':operator.ne,
                '>':operator.gt,
                '>=':operator.ge,
                }

def call(node, prop, datatype, operator_string, right_value):

    if datatype == 'str':
        if operator_string != '=' and operator_string != '!=':
            return False
        else:
            left_value = node.props.get(prop)
            if left_value:
                return operator_dict[operator_string](left_value, right_value)
    
    elif datatype == 'num':
        left_value = node.props.get(prop)
        if left_value:
            return operator_dict[operator_string](left_value, right_value)
        else:
            return False

    #left_value = node.props.get(prop)
    
    # if left_value:
    #     return operator_dict[operator_string](left_value, right_value)
    # else:
    #     return False

def to_code(string):
    conditional_output = []
    operators = list(operator_dict.keys())
    r = re.compile( '|'.join( '(?:{})'.format(re.escape(o)) for o in sorted(operators, reverse=True, key=len)) )

    #codes = string.split(',')
    # code = code.replace(",", " and ") # ',' means and
    # code = code.replace(";", " or ") # ';' means or
    # code = code.replace("=", " = ")
    # code = code.replace('>', " > ")
    # code = code.replace('>=', ' >= ')

    condition_strings = string.split(',')
    for condition_string in condition_strings:
        ops = r.findall(condition_string)
        for op in ops:
            condition_string = re.sub(op, ' '+op+' ', condition_string)
            left_value, op, right_value = condition_string.split(None,2)
            conditional_output.append([left_value, op, right_value])

    return conditional_output

# branch thicken, background highlighted to purple
def highlight_layout(conditions, level, prop2type={}, color='purple'):
    conditional_output = to_code(conditions)
    def layout_fn(node):
        final_call = False
        for condition in conditional_output:
            prop = condition[0]
            op = condition[1]
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
            node.add_face(prop_face, column=level, position = "branch_right")
            while (node):
                node = node.up
                if node:
                    node.sm_style["hz_line_width"] = 5
                    #node.sm_style["hz_line_color"] = color
                    #node.add_face(OutlineFace(color=color), column=level, collapsed_only=True)
    return layout_fn
    return     


# for boolean layouts
# hightlighted as rectangular

# hightlighted as circle