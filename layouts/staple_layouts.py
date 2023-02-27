import matplotlib as mpl
import numpy as np

from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import TextFace, Face, ScaleFace, LegendFace
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
            position="aligned", column=0, padding_x=10, 
            internal_rep='avg', scale=True, legend=True, active=True):
        super().__init__(name, 
                aligned_faces=True if position == "aligned" else False,
                legend=legend, active=active)
    
        self.width = width
        self.position = position
        self.column = column
        
        self.scale = scale
        self.padding_x = padding_x

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
            padding_x=10, scale=True, legend=True, active=True, 
            internal_rep='avg', scale_size=None, scale_range=None):

        name = name or f'Barplot_{size_prop}_{color_prop}'
        super().__init__(name=name, prop=prop, width=width, size_prop=size_prop,
                color_prop=color_prop, position=position, column=column,
                color_gradient=color_gradient, colors=colors, color=color,
                padding_x=padding_x, scale=scale, legend=legend, active=active,
                internal_rep=internal_rep)

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        
        if self.scale and self.size_range:
            self.scale_width = self.width
            self.scale_range = self.size_range
            scale = ScaleFace(width=self.width, scale_range=self.size_range, 
                    formatter='%.2f',
                    padding_x=self.padding_x, padding_y=2)
            text = TextFace(self.prop, max_fsize=11, padding_x=self.padding_x, rotation=315)
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
        # width = self.get_size(node)
        # color = self.get_color(node)
        
        if node.is_leaf() and node.props.get(self.prop):
            width = self.get_size(node, self.prop)
            
            
            #print(self.scale_width, self.scale_range)
            #color = self.get_color(node)
            color = self.color
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.size_prop:
                tooltip += f'<br>{self.prop}: {node.props.get(self.prop)}<br>'
            if self.color_prop:
                tooltip += f'<br>{self.color_prop}: {color}<br>'
            face = RectFace(width, None, color=color, 
                   tooltip=tooltip, padding_x=self.padding_x)
            node.add_face(face, position=self.position, column=self.column,
                    collapsed_only=False)
        
        elif node.is_leaf() and node.props.get(internal_prop):
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
                    tooltip=tooltip, padding_x=self.padding_x)
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
                    tooltip=tooltip, padding_x=self.padding_x)
            node.add_face(face, position=self.position, column=self.column,
                    collapsed_only=True)

class LayoutHeatmap(TreeLayout):
    def __init__(self, name=None, column=0, width=70, height=50, internal_rep=None, \
        prop=None, maxval=100, minval=0, min_color="#ffffff", max_color="#ff0000",\
        legend=True):
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
        self.padding_x = 2
        self.padding_y = 2
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.num_prop, min_fsize=5, max_fsize=10, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            colormap = { self.num_prop: {self.max_color},
                        }
            tree_style.add_legend(title=self.num_prop,
                                    variable='continuous',
                                    colormap=colormap,
                                    value_range=[self.minval, self.maxval],
                                    color_range=[self.min_color, self.max_color]
                                    )

    def set_node_style(self, node):
        c1 = self.min_color
        c2 = self.max_color #red
        if node.is_leaf() and node.props.get(self.num_prop):
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
                identF = RectFace(width=self.width, height=self.height, text="NaN", color=c1, 
                padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)

            node.add_face(identF, column = self.column,  position = 'aligned')
        
        elif node.is_leaf() and node.props.get(self.internal_prop):
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
                identF = RectFace(width=self.width, height=self.height, text="NaN", color=c1, 
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
                identF = RectFace(width=self.width, height=self.height, text="NaN", color=c1, 
                padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
            node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)
            # identF = RectFace(width=50,height=50,text="NaN", color=c1, 
            #     padding_x=1, padding_y=1)
            # node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)

        ## old way
        # if node.is_leaf() and node.props.get(self.num_prop):
            # redgradient = heatmap_gradient(0.95, 0.6, 10)
            # relative_abundance = float(node.props.get(self.num_prop))
            # tooltip = ""
            # if node.name:
            #     tooltip += f'<b>{node.name}</b><br>'
            # if self.num_prop:
            #     tooltip += f'<br>{self.num_prop}: {node.props.get(self.num_prop)}<br>'
            
            # try:
            #     color_idx = int(relative_abundance*10)
            #     color = redgradient[color_idx]
                
            #     identF = RectFace(width=50,height=50,text="%.2f" % (relative_abundance), color=color, 
            #     padding_x=1, padding_y=1, tooltip=tooltip)
            # except ValueError: # for miss data
            #     color = redgradient[0]
                
            #     identF = RectFace(width=50,height=50,text="NaN", color=color, 
            #     padding_x=1, padding_y=1, tooltip=tooltip)
            # node.add_face(identF, column = self.column,  position = 'aligned')
            
        # elif node.is_leaf() and node.props.get(self.internal_prop):
        #     # heatmap
        #     redgradient = heatmap_gradient(0.95, 0.6, 10)
        #     relative_abundance = float(node.props.get(self.internal_prop))
        #     color_idx = int(relative_abundance*10)
        #     color = redgradient[color_idx]
        #     tooltip = ""
        #     if node.name:
        #         tooltip += f'<b>{node.name}</b><br>'
        #     if self.internal_prop:
        #         tooltip += f'<br>{self.internal_prop}: {node.props.get(self.internal_prop)}<br>'
            
        #     identF = RectFace(width=50,height=50,text="%.1f" % (relative_abundance), color=color, 
        #     padding_x=1, padding_y=1)
        #     node.add_face(identF, column = self.column,  position = 'aligned')

        # elif node.props.get(self.internal_prop):
        #     # heatmap
        #     redgradient = heatmap_gradient(0.95, 0.6, 10)
        #     relative_abundance = float(node.props.get(self.internal_prop))
        #     color_idx = int(relative_abundance*10)
            
        #     color = redgradient[color_idx]
        #     tooltip = ""
        #     if node.name:
        #         tooltip += f'<b>{node.name}</b><br>'
        #     if self.internal_prop:
        #         tooltip += f'<br>{self.internal_prop}: {node.props.get(self.internal_prop)}<br>'
            
        #     identF = RectFace(width=50,height=50,text="%.1f" % (relative_abundance), color=color, 
        #     padding_x=1, padding_y=1, tooltip=tooltip)
        #     #face_name = TextFace(node.props.get('name'), color="red")
        #     #face_name = TextFace("%.1f" % (relative_abundance), color=color)
        #     node.add_face(identF, column = self.column,  position = 'aligned', collapsed_only=True)
class RectFace(Face):
    def __init__(self, width, height, color='gray',
            opacity=0.7,
            text=None, fgcolor='black', # text color
            min_fsize=6, max_fsize=15,
            ftype='sans-serif',
            tooltip=None,
            name="",
            padding_x=0, padding_y=0):

        Face.__init__(self, name=name, padding_x=padding_x, padding_y=padding_y)

        self.width = width
        self.height = height
        self.stretch = True
        self.color = color
        self.opacity = opacity
        # Text related
        self.text = str(text) if text is not None else None
        self.rotate_text = False
        self.fgcolor = fgcolor
        self.ftype = ftype
        self.min_fsize = min_fsize
        self.max_fsize = max_fsize

        self.tooltip = tooltip

    def __name__(self):
        return "RectFace"

    def compute_bounding_box(self, 
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before):

        if drawer.TYPE == 'circ':
            pos = swap_pos(pos, point[1])

        box = super().compute_bounding_box( 
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before)

        x, y, dx, dy = box
        zx, zy = self.zoom
        zx = 1 if self.stretch\
                and pos.startswith('aligned')\
                and drawer.TYPE != 'circ'\
                else zx

        r = (x or 1e-10) if drawer.TYPE == 'circ' else 1

        def get_dimensions(max_width, max_height):
            if not (max_width or max_height):
                return 0, 0
            if (type(max_width) in (int, float) and max_width <= 0) or\
               (type(max_height) in (int, float) and max_height <= 0):
                return 0, 0

            width = self.width / zx if self.width is not None else None
            height = self.height / zy if self.height is not None else None

            if width is None:
                return max_width or 0, min(height or float('inf'), max_height)
            if height is None:
                return min(width, max_width or float('inf')), max_height

            hw_ratio = height / width

            if max_width and width > max_width:
                width = max_width
                height = width * hw_ratio
            if max_height and height > max_height:
                height = max_height
                if not self.stretch or drawer.TYPE == 'circ':
                    width = height / hw_ratio

            height /= r  # in circular drawer
            return width, height

        max_dy = dy * r  # take into account circular mode

        if pos == 'branch_top':
            width, height = get_dimensions(dx, max_dy)
            box = (x, y + dy - height, width, height) # container bottom

        elif pos == 'branch_bottom':
            width, height = get_dimensions(dx, max_dy)
            box = (x, y, width, height) # container top

        elif pos == 'branch_right':
            width, height = get_dimensions(dx, max_dy)
            box = (x, y + (dy - height) / 2, width, height)

        elif pos.startswith('aligned'):
            width, height = get_dimensions(None, dy)
            # height = min(dy, (self.height - 2 * self.padding_y) / zy)
            # width = min(self.width - 2 * self.padding_x) / zx

            if pos == 'aligned_bottom':
                y = y + dy - height
            elif pos == 'aligned_top':
                y = y
            else:
                y = y + (dy - height) / 2

            box = (x, y, width, height)

        self._box = Box(*box)
        return self._box

    def draw(self, drawer):
        self._check_own_variables()

        circ_drawer = drawer.TYPE == 'circ'
        style = {'fill': self.color, 'opacity': self.opacity}
        if self.text and circ_drawer:
            rect_id = get_random_string(10)
            style['id'] = rect_id

        yield draw_rect(self._box,
                self.name,
                style=style,
                tooltip=self.tooltip)

        if self.text:
            x, y, dx, dy = self._box
            zx, zy = self.zoom

            r = (x or 1e-10) if circ_drawer else 1
            if self.rotate_text:
                rotation = 90
                self.compute_fsize(dy * zy / (len(self.text) * zx) * r,
                                   dx * zx / zy, zx, zy)

                text_box = Box(x + (dx - self._fsize / (2 * zx)) / 2,
                        y + dy / 2,
                        dx, dy)
            else:
                rotation = 0
                self.compute_fsize(dx / len(self.text), dy, zx, zy)
                text_box = Box(x + dx / 2,
                        y + (dy - self._fsize / (zy * r)) / 2,
                        dx, dy)
            text_style = {
                'max_fsize': self._fsize,
                'text_anchor': 'middle',
                'ftype': f'{self.ftype}, sans-serif', # default sans-serif
                }

            if circ_drawer:
                offset = dx * zx + dy * zy * r / 2
                # Turn text upside down on bottom
                if y + dy / 2 > 0:
                    offset += dx * zx + dy * zy * r
                text_style['offset'] = offset

            yield draw_text(text_box,
                    self.text,
                    rotation=rotation,
                    anchor=('#' + str(rect_id)) if circ_drawer else None,
                    style=text_style)
    