import matplotlib as mpl
import numpy as np

from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import TextFace, Face, ScaleFace, LegendFace, RectFace
from ete4.smartview.renderer.draw_helpers import *
from ete4.treeview.svg_colors import random_color

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
            position="aligned", column=0, padding_x=10, padding_y=0,
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
            

        if self.color_prop:
            unique = vals["color"][3]
            if len(unique):
                colors = self.colors or random_color(num=len(unique))
                if type(colors) == dict:
                    self.colors = colors.copy()
                else:
                    colors = list(colors)
                    self.colors = {}
                    for idx, value in enumerate(unique):
                        self.colors[value] = colors[idx % len(colors)]
                if self.legend:
                    tree_style.add_legend(title=self.name,
                            variable="discrete",
                            colormap=self.colors)
            else:
                self.color_range = vals["color"][1:3]
                if self.legend:
                    tree_style.add_legend(title=self.name, 
                            variable="continuous",
                            value_range=self.color_range,
                            color_range=self.color_gradient)

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
            color_gradient=None, color="red", colors=None, 
            padding_x=10, padding_y=0, scale=True, legend=True, active=True, 
            internal_rep='avg', scale_size=None, scale_range=None):

        name = name or f'Barplot_{size_prop}_{color_prop}'
        super().__init__(name=name, prop=prop, width=width, size_prop=size_prop,
                color_prop=color_prop, position=position, column=column,
                color_gradient=color_gradient, colors=colors, color=color,
                padding_x=padding_x, padding_y=padding_y,scale=scale, legend=legend, active=active,
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
            colormap = { self.prop: self.color
                        }
            tree_style.add_legend(title=self.prop,
                                    variable='discrete',
                                    colormap=colormap
                                    )
        
    def set_node_style(self, node):
        internal_prop = self.prop + '_' + self.internal_rep  
        if node.is_leaf and node.props.get(self.prop):
            width = self.get_size(node, self.prop)
            color = self.color
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
            color = self.color
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
            color = self.color
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
    def __init__(self, name=None, column=0, width=70, height=None, padding_x=1, padding_y=0, \
        internal_rep=None, prop=None, maxval=100, minval=0, min_color="#ffffff", \
        max_color="#971919", legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.num_prop = prop
        self.column = column
        #self.colour_dict = colour_dict
        self.min_color = min_color
        self.max_color = max_color
        self.maxval = maxval
        self.minval = minval
        self.internal_prop = prop+'_'+internal_rep
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
            colormap = { self.num_prop: {self.max_color},
                        }
            tree_style.add_legend(title=self.num_prop,
                                    variable='continuous',
                                    colormap=colormap,
                                    value_range=[self.minval, self.maxval],
                                    color_range=[self.max_color, self.min_color ]
                                    )

    def set_node_style(self, node):
        c1 = self.min_color
        c2 = self.max_color #red
        if node.is_leaf and node.props.get(self.num_prop):
            # heatmap
            
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.num_prop:
                tooltip += f'<br>{self.num_prop}: {node.props.get(self.num_prop)}<br>'

            try:
                ratio = float(node.props.get(self.num_prop)) / self.maxval
                gradient_color = color_gradient(c1, c2, mix=ratio)
                identF = RectFace(width=self.width, height=self.height, text="%.2f" % (float(node.props.get(self.num_prop))), color=gradient_color, 
                    padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
               
            except ValueError: # for miss data
                identF = RectFace(width=self.width, height=self.height, text="NA", color=c1, 
                padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)

            node.add_face(identF, column = self.column,  position = 'aligned')
        
        elif node.is_leaf and node.props.get(self.internal_prop):
            # heatmap
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.num_prop:
                tooltip += f'<br>{self.internal_prop}: {node.props.get(self.internal_prop)}<br>'
            try:
                ratio = float(node.props.get(self.internal_prop)) / self.maxval
                gradient_color = color_gradient(c1, c2, mix=ratio)
                identF = RectFace(width=self.width, height=self.height, text="%.2f" % (float(node.props.get(self.internal_prop))), color=gradient_color, 
                    padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            except ValueError: # for miss data
                identF = RectFace(width=self.width, height=self.height, text="NA", color=c1, 
                padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            node.add_face(identF, column = self.column,  position = 'aligned')

        elif node.props.get(self.internal_prop):
            # heatmap
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.num_prop:
                tooltip += f'<br>{self.internal_prop}: {node.props.get(self.internal_prop)}<br>'
            try:
                ratio = float(node.props.get(self.internal_prop)) / self.maxval
                gradient_color = color_gradient(c1, c2, mix=ratio)
                identF = RectFace(width=self.width, height=self.height, text="%.2f" % (float(node.props.get(self.internal_prop))), color=gradient_color, 
                    padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            except ValueError: # for miss data
                identF = RectFace(width=self.width, height=self.height, text="NA", color=c1, 
                padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)
            # identF = RectFace(width=50,height=50,text="NA", color=c1, 
            #     padding_x=1, padding_y=1)
            # node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)
