from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import (RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, \
                            SelectedFace, SelectedCircleFace, SelectedRectFace, LegendFace)
from treeprofiler.layouts.general_layouts import get_heatmapface, get_aggregated_heatmapface
from treeprofiler.src.utils import to_code, call, counter_call, check_nan
# for boolean layouts
try:
    from distutils.util import strtobool
except ImportError:
    from treeprofiler.src.utils import strtobool

# branch thicken, background highlighted to purple
class LayoutHighlight(TreeLayout):
    def __init__(self, name, color2conditions, column, prop2type=None, legend=True, width=70, padding_x=1, padding_y=0):
        super().__init__(name)
        self.name = name
        self.aligned_faces = True
        self.prop2type = prop2type
        self.color2conditions = color2conditions

        self.min_fsize = 5 
        self.max_fsize = 15
        self.width = 70
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            colormap = {','.join(v) if isinstance(v, list) else v: k for k, v in self.color2conditions.items()}
            tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=colormap
                                    )

        for color, conditions in self.color2conditions.items():
            conditional_output = to_code(conditions)
            for node in tree.traverse():
                final_call = False
                for condition in conditional_output:
                    op = condition[1]
                    if op == 'in':
                        value = condition[0]
                        prop = condition[2]
                        datatype = self.prop2type.get(prop)
                        final_call = call(node, prop, datatype, op, value)

                    elif ':' in condition[0] :
                        internal_prop, leaf_prop = condition[0].split(':')
                        value = condition[2]
                        datatype = self.prop2type[internal_prop]
                        final_call = counter_call(node, internal_prop, leaf_prop, datatype, op, value)
                    else:
                        prop = condition[0]
                        value = condition[2]
                        datatype = self.prop2type.get(prop)
                        final_call = call(node, prop, datatype, op, value)
                    if final_call == False:
                        break
                    else:
                        continue
                
                if final_call:
                    #prop_face = SelectedRectFace(name='prop')
                    node.add_prop(f'hl_{conditions}', color)  # highligh clade
                    node.add_prop(f'hl_{conditions}_endnode', True)
                    while (node):
                        node = node.up
                        if node:
                            node.add_prop(f'hl_{conditions}', True)
                            #node.sm_style["hz_line_width"] = 5
        return

    def set_node_style(self, node):
        for color, conditions in self.color2conditions.items():
            if node.props.get(f'hl_{conditions}'):
                node.sm_style["hz_line_width"] = 5
                node.sm_style["hz_line_color"] = color
                node.sm_style["outline_color"] = color
                if node.props.get(f'hl_{conditions}_endnode'):
                    node.sm_style["bgcolor"] = color

class LayoutCollapse(TreeLayout):
    def __init__(self, name, color2conditions, column, prop2type=None, legend=True, width=70, padding_x=1, padding_y=0):
        super().__init__(name)
        self.name = name
        self.aligned_faces = True
        self.prop2type = prop2type
        self.color2conditions = color2conditions

        self.min_fsize = 5 
        self.max_fsize = 15
        self.width = 70
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            colormap = {','.join(v) if isinstance(v, list) else v: k for k, v in self.color2conditions.items()}
            tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=colormap
                                    )

        for color, conditions in self.color2conditions.items():
            conditional_output = to_code(conditions)
            for node in tree.traverse():
                final_call = False
                for condition in conditional_output:
                    op = condition[1]
                    if op == 'in':
                        value = condition[0]
                        prop = condition[2]
                        datatype = self.prop2type.get(prop)
                        final_call = call(node, prop, datatype, op, value)

                    elif ':' in condition[0] :
                        internal_prop, leaf_prop = condition[0].split(':')
                        value = condition[2]
                        datatype = self.prop2type[internal_prop]
                        final_call = counter_call(node, internal_prop, leaf_prop, datatype, op, value)
                    else:
                        prop = condition[0]
                        value = condition[2]
                        datatype = self.prop2type.get(prop)
                        final_call = call(node, prop, datatype, op, value)
                        
                    if final_call == False:
                        break
                    else:
                        continue
                
                if final_call:
                    #prop_face = SelectedRectFace(name='prop')
                    node.add_prop(f'cl_{conditions}', color)  # highligh clade
                    node.add_prop(f'cl_{conditions}_endnode', True)
                    while (node):
                        node = node.up
                        if node:
                            node.add_prop(f'hl_{conditions}', True)
                            #node.sm_style["hz_line_width"] = 5
        return

    def set_node_style(self, node):
        for color, conditions in self.color2conditions.items():
            if not node.is_root:
                if node.props.get(f'cl_{conditions}'):
                    node.sm_style["draw_descendants"] = False
                    node.sm_style["outline_color"] = color
                    if node.props.get(f'cl_{conditions}_endnode'):
                        node.sm_style["draw_descendants"] = False
                        node.sm_style["outline_color"] = color

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
                datatype = prop2type.get(prop)
                final_call = call(node, prop, datatype, op, value)

            elif ':' in condition[0] :
                internal_prop, leaf_prop = condition[0].split(':')
                value = condition[2]
                datatype = prop2type[internal_prop]
                final_call = counter_call(node, internal_prop, leaf_prop, datatype, op, value)
            else:
                prop = condition[0]
                value = condition[2]
                datatype = prop2type.get(prop)
                final_call = call(node, prop, datatype, op, value)
            
            if final_call == False:
                break
            else:
                continue
        if final_call:
            if not node.is_root:
                node.sm_style["draw_descendants"] = False
                node.sm_style["outline_color"] = color
    return layout_fn
    return

class LayoutBinary(TreeLayout):
    def __init__(self, name=None, level=1, color='#E60A0A', \
            bool_prop=None, reverse=False, aggregate=False, \
            max_count=0, \
            radius=25, padding_x=1, padding_y=0, width=70, \
            legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.bool_prop = bool_prop
        self.column = level
        self.color = color
        self.negative_color = '#EBEBEB'
        self.internal_prop = bool_prop+'_counter'
        self.reverse = reverse
        self.aggregate = aggregate
        self.max_count = max_count
        self.radius = radius
        self.padding_x = padding_x
        self.padding_y = padding_y
        self.legend = legend
        self.width = width
        self.height = None
        self.min_fsize = 5
        self.max_fsize = 10
        
        
    # def set_tree_style(self, tree, tree_style):
    #     super().set_tree_style(tree, tree_style)
    #     text = TextFace(self.name, max_fsize=11, padding_x=1)
    #     tree_style.aligned_panel_header.add_face(text, column=self.column)
    def update_header_width(self):
        return

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.bool_prop, min_fsize=10, max_fsize=15, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)
        if self.legend:
            
            if self.reverse:
                title = 'ReverseBinary_' + self.bool_prop
                colormap = {
                    "False": self.color,
                    "True" : self.negative_color,
                    "NA": 'white'
                }
                tree_style.add_legend(title=title,
                                    variable='discrete',
                                    colormap=colormap,
                                    )
            else:
                title = 'Binary_' + self.bool_prop
                colormap = {
                    "True": self.color,
                    "False" : self.negative_color,
                    "NA": 'white'
                }
                tree_style.add_legend(title=title,
                                    variable='discrete',
                                    colormap=colormap,
                                    )

                
    def set_node_style(self, node):
        # need to correct
        if node.is_leaf and node.props.get(self.bool_prop):
            #if node.props.get(self.bool_prop):
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
                        #prop_face = CircleFace(radius=self.radius, color=self.negative_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        prop_face = RectFace(width=self.width, height=self.height, color=self.negative_color,  padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                else:
                    if bool(str2bool):
                        #prop_face = CircleFace(radius=self.radius, color=self.color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        prop_face = RectFace(width=self.width, height=self.height, color=self.color,  padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        node.add_face(prop_face, column=self.column, position = "aligned")
                    else:
                        #prop_face = CircleFace(radius=self.radius, color=self.negative_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        prop_face = RectFace(width=self.width, height=self.height, color=self.negative_color,  padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                        node.add_face(prop_face, column=self.column, position = "aligned")
            else: #mising
                prop_face = RectFace(width=self.width, height=self.height, text="NA", color=self.negative_color,  padding_x=self.padding_x, padding_y=self.padding_y, stroke_color=self.negative_color, tooltip=None)
                node.add_face(prop_face, column=self.column, position = "aligned")
            # else:
            #     prop_face = RectFace(width=self.width, height=self.height, text="NA", color=self.negative_color,  padding_x=self.padding_x, padding_y=self.padding_y, stroke_color=self.negative_color, tooltip=None)
            #     node.add_face(prop_face, column=self.column, position = "aligned")
        
        elif node.is_leaf and node.props.get(self.internal_prop):
            if self.aggregate:
                heatmapFace = get_aggregated_heatmapface(node, self.internal_prop, max_color=self.color, width=self.width, height=self.height, padding_x=self.padding_x, max_count=self.max_count)
            else:
                heatmapFace = get_heatmapface(node, self.internal_prop, max_color=self.color, width=self.width, height=self.height, padding_x=self.padding_x,)
            node.add_face(heatmapFace, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            if self.aggregate:
                heatmapFace = get_aggregated_heatmapface(node, self.internal_prop, max_color=self.color, width=self.width, height=self.height, padding_x=self.padding_x, max_count=self.max_count)
            else:
                heatmapFace = get_heatmapface(node, self.internal_prop, max_color=self.color, width=self.width, height=self.height, padding_x=self.padding_x,reverse=self.reverse)
            node.add_face(heatmapFace, column = self.column, position = "aligned", collapsed_only=True)
        # else:
        #     prop_face = RectFace(width=self.width, height=self.height, text="NA", color=self.negative_color,  padding_x=self.padding_x, padding_y=self.padding_y, stroke_color=self.negative_color, tooltip=None)
        #     node.add_face(prop_face, column=self.column, position = "aligned")
