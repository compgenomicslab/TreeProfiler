import matplotlib as mpl
import numpy as np

from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import TextFace, Face, ScaleFace, LegendFace, RectFace
from ete4.smartview.renderer.draw_helpers import *
from treeprofiler.src.utils import random_color, add_suffix

import colorsys



__all__ = [ "LayoutBarplot" ]

def heatmap_gradient(hue, intensity, granularity):
    min_lightness = 0.35 
    max_lightness = 0.9
    base_value = intensity

    # each gradient must contain 100 lightly descendant colors
    colors = []   
    rgb2hex = lambda rgb: '#%02x%02x%02x' % rgb
    l_factor = (max_lightness-min_lightness) / float(granularity)
    l = min_lightness
    while l <= max_lightness:
        l += l_factor
        rgb =  rgb2hex(tuple(map(lambda x: int(x*255), colorsys.hls_to_rgb(hue, l, base_value))))
        colors.append(rgb)
        
    colors.append("#ffffff")
    return list(reversed(colors))

def color_gradient(c1, c2, mix=0):
    """ Fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1) """
    # https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def swap_pos(pos, angle):
    if abs(angle) >= pi / 2:
        if pos == 'branch_top':
            pos = 'branch_bottom'
        elif pos == 'branch_bottom':
            pos = 'branch_top'
    return pos

class LayoutPlot(TreeLayout):
    def __init__(self, name=None, prop=None, width=200, size_prop=None, 
            color_prop=None, color_gradient=None, color="red", colors=None,
            position="aligned", column=0, padding_x=10, padding_y=0, size_range=[], 
            internal_rep='avg', scale=True, legend=True, active=True):
        super().__init__(name, 
                aligned_faces=True if position == "aligned" else False,
                legend=legend, active=active)
    
        self.width = width
        self.position = position
        self.column = column
        
        self.scale = scale
        self.padding_x = padding_x
        self.padding_y = padding_y

        self.internal_rep = internal_rep
        self.prop = prop
        # if not (size_prop or color_prop):
            # raise InvalidUsage("Either size_prop or color_prop required")

        self.size_prop = size_prop
        self.color_prop = color_prop
        self.size_range = size_range
        self.size_range = size_range

        self.color = color
        self.colors = colors
        self.color_gradient = color_gradient
        if self.color_prop and not self.color_gradient:
            self.color_gradient = ("#FFF", self.color)
    
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        def update_vals(metric, node):
            p, minval, maxval, uniqvals = vals[metric]
            prop = node.props.get(p)
            try:
                prop = float(prop) # prop by default is string 
                if type(prop) in [int, float]:
                    vals[metric][1] = min(minval, prop)
                    vals[metric][2] = max(maxval, prop)
                elif prop is None or prop == "":
                    return
                else:
                    uniqvals.add(prop)
            except:
                pass
        
        # only when size_range is not provided, calculate min and max values
        if self.size_range == []:
            vals = { 
                    "size": [ self.size_prop, 0, 0, set() ],    # min, max, unique
                    "color":  [ self.color_prop,  0, 0, set() ] # min, max, unique
                    }
                
            for node in tree.traverse():
                if self.size_prop:
                    update_vals("size", node)

                if self.color_prop:
                    update_vals("color", node)
            
            if self.size_prop:
                self.size_range = vals["size"][1:3]
                            # if self.color_prop:
            #     unique = vals["color"][3]
            #     if len(unique):
            #         colors = self.colors or random_color(num=len(unique))
            #         if type(colors) == dict:
            #             self.colors = colors.copy()
            #         else:
            #             colors = list(colors)
            #             self.colors = {}
            #             for idx, value in enumerate(unique):
            #                 self.colors[value] = colors[idx % len(colors)]
            #         if self.legend:
            #             tree_style.add_legend(title=self.name,
            #                     variable="discrete",
            #                     colormap=self.colors)
            #     else:
            #         self.color_range = vals["color"][1:3]
            #         if self.legend:
            #             tree_style.add_legend(title=self.name, 
            #                     variable="continuous",
            #                     value_range=self.color_range,
            #                     color_range=self.color_gradient)

    # def get_size(self, node, prop):
    #     if self.size_range != [0, 0]:
    #         minval, maxval = self.size_range
    #     else:
    #         minval, maxval = 0,1
    #     return float(node.props.get(prop, 0)) / float(maxval) * self.width
    
    def get_size(self, node, prop):
        if not self.size_prop:
            return self.width
        minval, maxval = self.size_range
        return float(node.props.get(prop, 0)) / maxval * self.width

    def get_color(self, node):
        if not self.color_prop:
            return self.color
        prop = node.props.get(self.color_prop)
        if prop is None:
            return None
        if self.color_range:
            minval, maxval = self.color_range
            return color_gradient(*self.color_gradient, (prop - minval) / maxval)
        else:
            return self.colors.get(prop)

    def get_legend(self):
        return self.legend

class LayoutBarplot(LayoutPlot):
    def __init__(self, name=None, prop=None, width=200, size_prop=None,
            color_prop=None, position="aligned", column=0, 
            color_gradient=None, color=None, colors=None, 
            padding_x=10, padding_y=0, scale=True, legend=True, active=True, 
            internal_rep='avg', scale_size=None, size_range=None, scale_range=None):

        name = name or f'Barplot_{size_prop}_{color_prop}'
        super().__init__(name=name, prop=prop, width=width, size_prop=size_prop,
                color_prop=color_prop, position=position, column=column,
                color_gradient=color_gradient, colors=colors, color=color,
                padding_x=padding_x, padding_y=padding_y,  scale=scale, size_range=size_range, 
                legend=legend, active=active,
                internal_rep=internal_rep)

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.scale and self.size_range:
            self.scale_width = self.width
            self.scale_range = self.size_range
            scale = ScaleFace(width=self.width, scale_range=self.size_range, 
                    formatter='%.2f',
                    padding_x=self.padding_x, padding_y=2)
            text = TextFace(self.prop, max_fsize=15, padding_x=self.padding_x, rotation=315)
            tree_style.aligned_panel_header.add_face(scale, column=self.column)
            tree_style.aligned_panel_header.add_face(text, column=self.column)
        
        if self.legend:
            if self.color:
                colormap = {self.prop: self.color
                            }
            else:
                colormap = self.colors

            tree_style.add_legend(title=self.prop,
                                    variable='discrete',
                                    colormap=colormap
                                    )

    def get_color(self, node, color_prop, color_dict):
        if color_dict and color_prop:
            prop = node.props.get(color_prop)
            if prop is None:
                return None
            else:
                return color_dict.get(prop, None)
        else:
            return self.color

    def set_node_style(self, node):
        internal_prop = self.prop + '_' + self.internal_rep  
        if node.props.get(self.prop) is not None:
            if node.is_leaf:
                width = self.get_size(node, self.prop)
                color = self.get_color(node, self.color_prop, self.colors)
                tooltip = ""
                if node.name:
                    tooltip += f'<b>{node.name}</b><br>'
                if self.size_prop:
                    tooltip += f'<br>{self.prop}: {node.props.get(self.prop)}<br>'
                if self.color_prop:
                    tooltip += f'<br>{self.color_prop}: {color}<br>'
                face = RectFace(width, None, color=color, 
                    tooltip=tooltip, padding_x=self.padding_x, padding_y=self.padding_y)
                node.add_face(face, position=self.position, column=self.column,
                        collapsed_only=False)
        
        elif node.is_leaf and node.props.get(internal_prop):
            width = self.get_size(node, internal_prop)
            color = self.get_color(node, self.color_prop, self.colors)
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.size_prop:
                tooltip += f'<br>{self.prop}: {node.props.get(internal_prop)}<br>'
            if self.color_prop:
                tooltip += f'<br>{self.color_prop}: {color}<br>'
            face = RectFace(width, None, color=color, 
                    tooltip=tooltip, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(face, position=self.position, column=self.column,
                    collapsed_only=False)

        elif node.props.get(internal_prop):
            width = self.get_size(node, internal_prop)
            color = self.get_color(node, self.color_prop, self.colors)
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.size_prop:
                tooltip += f'<br>{self.size_prop}: {node.props.get(internal_prop)}<br>'
            if self.color_prop:
                tooltip += f'<br>{self.color_prop}: {color}<br>'
            face = RectFace(width, None, color=color, 
                    tooltip=tooltip, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(face, position=self.position, column=self.column,
                    collapsed_only=True)

class LayoutHeatmap(TreeLayout):
    def __init__(self, name=None, column=0, width=70, height=None, 
            padding_x=1, padding_y=0, heatmap_prop=None, internal_rep=None,
            value_color=None, value_range=[], color_range=None, minval=0, maxval=None, 
            absence_color="#EBEBEB",
            legend=True):

        super().__init__(name)
        self.aligned_faces = True

        self.heatmap_prop = heatmap_prop
        self.internal_prop = add_suffix(heatmap_prop, internal_rep)
        self.column = column
        self.value_color = value_color
        self.value_range = value_range
        self.color_range = color_range
        self.absence_color = absence_color
        self.maxval = maxval
        self.minval = minval

        self.width = width
        self.height = height
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)

        text = TextFace(self.heatmap_prop,  padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            tree_style.add_legend(title=self.heatmap_prop,
                                    variable='continuous',
                                    value_range=self.value_range ,
                                    color_range=[
                                        self.color_range.get(20),
                                        self.color_range.get(10),
                                        self.color_range.get(1),
                                    ]
                                    )
    def set_node_style(self, node):
        heatmap_num = node.props.get(self.heatmap_prop)
        if heatmap_num is not None and heatmap_num != 'NaN':
            heatmap_num = float(heatmap_num)
            if node.is_leaf:
                # heatmap
                tooltip = ""
                if node.name:
                    tooltip += f'<b>{node.name}</b><br>'
                if self.heatmap_prop:
                    tooltip += f'<br>{self.heatmap_prop}: {heatmap_num}<br>'
                
                gradient_color = self.value_color.get(heatmap_num)
                
                if gradient_color:
                    identF = RectFace(width=self.width, height=self.height, 
                    color=gradient_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                    node.add_face(identF, column = self.column,  position = 'aligned')
            
        elif node.is_leaf and node.props.get(self.internal_prop):
            heatmap_num = node.props.get(self.internal_prop)
            heatmap_num = float(heatmap_num)
            # heatmap
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.heatmap_prop:
                tooltip += f'<br>{self.heatmap_prop}: {heatmap_num}<br>'

            gradient_color = self.value_color.get(heatmap_num)

            if gradient_color:
                identF = RectFace(width=self.width, height=self.height, 
                color=gradient_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)
        
        elif node.props.get(self.internal_prop):
            heatmap_num = node.props.get(self.internal_prop)
            heatmap_num = float(heatmap_num)
            # heatmap
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.heatmap_prop:
                tooltip += f'<br>{self.heatmap_prop}: {heatmap_num}<br>'
            
            gradient_color = self.value_color.get(heatmap_num)

            if gradient_color:
                identF = RectFace(width=self.width, height=self.height,
                color=gradient_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)

        else:
            # heatmap
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.heatmap_prop:
                tooltip += f'<br>{self.heatmap_prop}: {heatmap_num}<br>'

            identF = RectFace(width=self.width, height=self.height, text=heatmap_num, color=self.absence_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=None)
            node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=False)

class LayoutHeatmapOld(TreeLayout):
    def __init__(self, name=None, column=0, width=70, height=None, padding_x=1, padding_y=0, \
            internal_rep=None, prop=None, maxval=100, minval=0, mean_val=0, std_val=0, \
            norm_method='min-max', color_dict=None, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.num_prop = prop
        self.column = column
        self.color_dict = color_dict
        self.maxval = maxval
        self.minval = minval
        self.mean_val = mean_val
        self.std_val = std_val
        self.norm_method = norm_method

        self.internal_prop = add_suffix(prop, internal_rep)
        self.width = width
        self.height = height
        self.padding_x = padding_x
        self.padding_y = padding_y
        self.min_fsize = 5
        self.max_fsize = 15
        
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.num_prop, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            tree_style.add_legend(title=self.num_prop,
                                    variable='continuous',
                                    value_range=[self.minval, self.maxval],
                                    color_range=[
                                        self.color_dict[20], 
                                        self.color_dict[10], 
                                        self.color_dict[1]
                                    ]
                                    )

    def min_max_normalize(self, value):
        return (value - self.minval) / (self.maxval - self.minval)

    def mean_normalize(self, value):
        return (value - self.mean_val) / (self.maxval - self.minval)

    def z_score_normalize(self, value):
        return (value - self.mean_val) / self.std_val

    def _get_color(self, search_value, norm_method='min-max'):
        num = len(self.color_dict)
        search_value = float(search_value)
        if norm_method == "min-max":
            normalized_value = self.min_max_normalize(search_value)
            index_values = np.linspace(0, 1, num)
        elif norm_method == "mean":
            normalized_value = self.mean_normalize(search_value)
            index_values = np.linspace(-1, 1, num)
        elif norm_method == "zscore":
            normalized_value = self.z_score_normalize(search_value)
            index_values = np.linspace(-3, 3, num)
        else:
            raise ValueError("Unsupported normalization method.")
        index = np.abs(index_values - normalized_value).argmin() + 1
        #index = np.abs(index_values - search_value).argmin() + 1
        return self.color_dict.get(index, "")

    def set_node_style(self, node):
        if node.props.get(self.num_prop) is not None:
            if node.is_leaf:
                # heatmap
                tooltip = ""
                if node.name:
                    tooltip += f'<b>{node.name}</b><br>'
                if self.num_prop:
                    tooltip += f'<br>{self.num_prop}: {node.props.get(self.num_prop)}<br>'

                gradient_color = self._get_color(node.props.get(self.num_prop), norm_method=self.norm_method)
                if gradient_color:
                    identF = RectFace(width=self.width, height=self.height, text="%.2f" % (float(node.props.get(self.num_prop))), \
                    color=gradient_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                else:   # for miss data
                    identF = RectFace(width=self.width, height=self.height, text="NA", 
                    color="", 
                    padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)

                node.add_face(identF, column = self.column,  position = 'aligned')
        
        elif node.is_leaf and node.props.get(self.internal_prop):
            # heatmap
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.num_prop:
                tooltip += f'<br>{self.internal_prop}: {node.props.get(self.internal_prop)}<br>'
            
            gradient_color = self._get_color(node.props.get(self.internal_prop), norm_method=self.norm_method)
            if gradient_color:
                identF = RectFace(width=self.width, height=self.height, text="%.2f" % (float(node.props.get(self.internal_prop))), \
                color=gradient_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            else:   # for miss data
                identF = RectFace(width=self.width, height=self.height, text="NA", 
                color="", 
                padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)

        elif node.props.get(self.internal_prop):
            # heatmap
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.num_prop:
                tooltip += f'<br>{self.internal_prop}: {node.props.get(self.internal_prop)}<br>'
            
            gradient_color = self._get_color(node.props.get(self.internal_prop), norm_method=self.norm_method)
            if gradient_color:
                identF = RectFace(width=self.width, height=self.height, text="%.2f" % (float(node.props.get(self.internal_prop))), \
                color=gradient_color, padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            else:   # for miss data
                identF = RectFace(width=self.width, height=self.height, text="NA", 
                color="", 
                padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)

class LayoutBranchScore(TreeLayout):
    def __init__(self, name, color_dict, score_prop, internal_rep=None, \
    value_range=None, color_range=None, show_score=False, legend=True, active=True):
        super().__init__(name)
        self.aligned_faces = True
        self.score_prop = score_prop
        if internal_rep:
            self.internal_prop = add_suffix(score_prop, internal_rep)
        else:
            self.internal_prop = None
        self.color_dict = color_dict
        self.legend = legend
        self.absence_color = "black"
        self.value_range = value_range
        self.color_range = color_range
        self.show_score = show_score
        self.line_width = 3
        self.active = active

    def set_tree_style(self, tree, tree_style):
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.name,
                                    variable='continuous',
                                    value_range=self.value_range,
                                    color_range=self.color_range,
                                    )

    # def _get_color_for_score(self, search_value):
    #     num = len(self.color_dict)
    #     index_values = np.linspace(self.value_range[0], self.value_range[1], num)
    #     index = np.abs(index_values - search_value).argmin() + 1
    #     return self.color_dict.get(index, self.absence_color)

    def set_node_style(self, node):
        prop_score = node.props.get(self.score_prop)
        if prop_score is not None:
            prop_score = float(prop_score)
            node.sm_style["hz_line_color"] = self.color_dict.get(prop_score)
            node.sm_style["hz_line_width"] = self.line_width
            node.sm_style["vt_line_color"] = self.color_dict.get(prop_score)
            node.sm_style["vt_line_width"] = self.line_width
            node.sm_style["outline_color"] = self.color_dict.get(prop_score)
            
            if self.show_score:                
                node.add_face(
                    TextFace("%.2f" % (float(prop_score)), color=self.color_dict.get(prop_score)),
                    position="branch_bottom")

        elif node.is_leaf and node.props.get(self.internal_prop):
            prop_score = node.props.get(self.internal_prop)
            prop_score = float(prop_score)
            node.sm_style["hz_line_color"] = self.color_dict.get(prop_score)
            node.sm_style["hz_line_width"] = self.line_width
            node.sm_style["vt_line_color"] = self.color_dict.get(prop_score)
            node.sm_style["vt_line_width"] = self.line_width
            node.sm_style["outline_color"] = self.color_dict.get(prop_score)
            
            if self.show_score:                
                node.add_face(
                    TextFace("%.2f" % (float(prop_score)), color=self.color_dict.get(prop_score)),
                    position="branch_bottom")

        elif node.props.get(self.internal_prop):
            prop_score = node.props.get(self.internal_prop)
            prop_score = float(prop_score)
            node.sm_style["hz_line_color"] = self.color_dict.get(prop_score)
            node.sm_style["hz_line_width"] = self.line_width
            node.sm_style["vt_line_color"] = self.color_dict.get(prop_score)
            node.sm_style["vt_line_width"] = self.line_width
            node.sm_style["outline_color"] = self.color_dict.get(prop_score)

            if self.show_score:                
                node.add_face(
                    TextFace("%.2f" % (float(prop_score)), color=self.color_dict.get(prop_score)),
                    position="branch_bottom")

class LayoutBubbleNumerical(TreeLayout):
    def __init__(self, name=None, prop=None, position="aligned", 
            column=0, color=None, max_radius=10, abs_maxval=None,
            padding_x=2, padding_y=0, 
            scale=True, legend=True, active=True, 
            internal_rep='avg'):

        name = name or f'Barplot_{size_prop}_{color_prop}'
        super().__init__(name)

        self.aligned_faces = True
        self.num_prop = prop
        self.internal_prop = add_suffix(prop, internal_rep)
        
        self.column = column
        self.position = position
        self.color = color
        self.positive_color = "#ff0000"
        self.negative_color = "#0000ff"
        self.internal_rep = internal_rep
        self.max_radius = float(max_radius)
        self.abs_maxval = float(abs_maxval)
        self.fgopacity = 0.7
        self.padding_x = padding_x
        self.padding_y = padding_y
        
        self.legend = legend
        self.active = active

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.num_prop, min_fsize=5, max_fsize=15, padding_x=self.padding_x, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)
        colormap = {
            "positive": self.positive_color,
            "negative": self.negative_color
            }
        if self.legend:
            tree_style.add_legend(title=self.num_prop,
                                    variable='discrete',
                                    colormap=colormap
                                    )

    def _get_bubble_size(self, search_value):
        search_value = abs(float(search_value))
        bubble_size = search_value / self.abs_maxval * self.max_radius
        return bubble_size


    def set_node_style(self, node):
        number = node.props.get(self.num_prop)
        if number is not None:
            # Ensure number is converted to float
            number = float(number)
            
            # Set bubble size and color based on the number's value
            bubble_size = self._get_bubble_size(number)
            bubble_color = self.positive_color if number > 0 else self.negative_color
            
            # Apply styles to the node
            node.sm_style["size"] = bubble_size
            node.sm_style["fgcolor"] = bubble_color
            node.sm_style["fgopacity"] = self.fgopacity

        elif node.is_leaf and node.props.get(self.internal_prop):
            # Since it's a leaf node with internal properties, use internal_prop as number
            number = float(node.props.get(self.internal_prop))
            
            # Set bubble size and color based on the number's value
            bubble_size = self._get_bubble_size(number)
            bubble_color = self.positive_color if number > 0 else self.negative_color

            # Apply styles to the node
            node.sm_style["size"] = bubble_size
            node.sm_style["fgcolor"] = bubble_color
            node.sm_style["fgopacity"] = self.fgopacity

        elif node.props.get(self.internal_prop):
            # Handle non-leaf nodes with internal properties
            number = float(node.props.get(self.internal_prop))
            
            # Set bubble size and color based on the number's value
            bubble_size = self._get_bubble_size(number)
            bubble_color = self.positive_color if number > 0 else self.negative_color

            # Apply styles to the node
            node.sm_style["size"] = bubble_size
            node.sm_style["fgcolor"] = bubble_color
            node.sm_style["fgopacity"] = self.fgopacity
