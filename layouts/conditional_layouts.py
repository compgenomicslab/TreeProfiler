from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import (RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, \
                            SelectedFace, SelectedCircleFace, SelectedRectFace, LegendFace)
from layouts.general_layouts import get_piechartface, get_heatmapface
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
    def __init__(self, name=None, level=1, color='red', prop_colour_dict=None, \
        bool_prop=None, reverse=False, radius=25, padding_x=2, padding_y=2, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.bool_prop = bool_prop
        self.column = level
        self.color = color
        self.prop_colour_dict = prop_colour_dict
        self.internal_prop = bool_prop+'_counter'
        self.reverse = reverse
        self.radius = radius
        self.padding_x = padding_x
        self.padding_y = padding_y
        self.legend = legend
        self.width = 70
        self.height = 50

    # def set_tree_style(self, tree, tree_style):
    #     super().set_tree_style(tree, tree_style)
    #     text = TextFace(self.name, max_fsize=11, padding_x=1)
    #     tree_style.aligned_panel_header.add_face(text, column=self.column)
    def update_header_width(self):
        return

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.bool_prop, min_fsize=5, max_fsize=10, padding_x=self.padding_x, width=70, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)
        if self.legend:
            if self.prop_colour_dict:
                if self.reverse:
                    title = 'ReverseBinary_' + self.bool_prop
                else:
                    title = 'Binary_' + self.bool_prop
                tree_style.add_legend(title=title,
                                    variable='discrete',
                                    colormap={self.bool_prop:self.color}
                                    )
                # tree_style.add_legend(title=self.internal_prop,
                #                     variable='discrete',
                #                     colormap=self.prop_colour_dict
                #                     )
                

    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.bool_prop):
            prop_bool = node.props.get(self.bool_prop)
            
            if not check_nan(prop_bool):
                str2bool = strtobool(prop_bool)
                tooltip = ""
                if node.name:
                    tooltip += f'<b>{node.name}</b><br>'
                if self.bool_prop:
                    tooltip += f'<br>{self.bool_prop}: {node.props.get(self.bool_prop)}<br>'

                if self.reverse:
                    if not bool(str2bool):
                        
                        #prop_face = CircleFace(radius=self.radius, color=self.color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        prop_face = RectFace(width=self.width, height=self.height, color=self.color,  padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                    else:
                        #prop_face = CircleFace(radius=self.radius, color='white', padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        prop_face = RectFace(width=self.width, height=self.height, color='white',  padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                else:
                    if bool(str2bool):
                        #prop_face = CircleFace(radius=self.radius, color=self.color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        prop_face = RectFace(width=self.width, height=self.height, color=self.color,  padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                    else:
                        #prop_face = CircleFace(radius=self.radius, color='white', padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        prop_face = RectFace(width=self.width, height=self.height, color='white',  padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        node.add_face(prop_face, column=self.column, position = "aligned")
            else:
                #prop_face = CircleFace(radius=self.radius, color='grey', padding_x=self.padding_x, padding_y=self.padding_y)
                prop_face = TextFace('NaN', min_fsize=5, max_fsize=10, padding_x=self.padding_x, width=self.width)
                node.add_face(prop_face, column=self.column, position = "aligned")
        
        elif node.is_leaf() and node.props.get(self.internal_prop):
            # piechart_face = get_piechartface(node, self.internal_prop, self.prop_colour_dict, self.radius)
            # node.add_face(piechart_face, column = self.column, position = "branch_top")
            heatmapFace = get_heatmapface(node, self.internal_prop, self.color, width=self.width, height=self.height)
            node.add_face(heatmapFace, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            # if node.is_root():
            #     piechart_face = get_piechartface(node, self.internal_prop, self.prop_colour_dict, self.radius*0.15)
            # else:
            #     piechart_face = get_piechartface(node, self.internal_prop, self.prop_colour_dict, self.radius)
            #node.add_face(piechart_face, column = self.column, position = "branch_top")
            heatmapFace = get_heatmapface(node, self.internal_prop, self.color, width=self.width, height=self.height)
            node.add_face(heatmapFace, column = self.column, position = "aligned", collapsed_only=True)


# hightlighted as rectangular

# hightlighted as circle