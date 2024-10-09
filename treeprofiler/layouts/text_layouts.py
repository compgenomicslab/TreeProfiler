from collections import OrderedDict, namedtuple
from math import pi
import random, string

from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import Face, RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from ete4.smartview.renderer.draw_helpers import draw_text, draw_line, draw_array, draw_rect
from treeprofiler.layouts.general_layouts import get_piechartface, get_stackedbarface
from treeprofiler.src.utils import random_color, add_suffix
"""
label_layout, colorbranch_layout, rectangular_layout   
"""
Box = namedtuple('Box', 'x y dx dy')  # corner and size of a 2D shape

class AlignLinkFace(Face):
    def __init__(self, width=70, height=None,
            stroke_color='gray', stroke_width=0.5,
            line_type=1, opacity=0.8):
        """Line types: 0 solid, 1 dotted, 2 dashed"""

        Face.__init__(self, padding_x=0, padding_y=0)
        self.line = None
        self.width = width
        self.height = height
        self.stroke_color = stroke_color
        self.stroke_width = stroke_width
        self.type = line_type;
        self.opacity = opacity

        self.always_drawn = True

    def __name__(self):
        return "AlignLinkFace"

    def _compute_bounding_box(self,
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before):

        if drawer.NPANELS > 1 and drawer.viewport and pos == 'branch_right':
            x, y = point
            dx, dy = size
            p1 = (x + bdx + dx_before, y + dy/2)
            if drawer.TYPE == 'rect':
                p2 = (drawer.viewport.x + drawer.viewport.dx, y + dy/2)
            else:
                aligned = sorted(drawer.tree_style.aligned_grid_dxs.items())
                # p2 = (drawer.node_size(drawer.tree)[0], y + dy/2)
                if not len(aligned):
                    return Box(0, 0, 0, 0)
                p2 = (aligned[0][1] - bdx, y + dy/2)
                if p1[0] > p2[0]:
                    return Box(0, 0, 0, 0)
                p1, p2 = cartesian(p1), cartesian(p2)

            self.line = (p1, p2)

        return Box(0, 0, 0, 0) # Should not take space

    def get_box(self):
        return Box(0, 0, 0, 0) # Should not take space

    def fits(self):
        return True

    def _draw(self, drawer):

        if drawer.NPANELS < 2:
            return None

        style = {
                'type': self.type,
                'stroke': self.stroke_color,
                'stroke-width': self.stroke_width,
                'opacity': self.opacity,
                }
        if drawer.panel == 0 and drawer.viewport and\
          (self.node.is_leaf or self.node.is_collapsed)\
          and self.line:
            p1, p2 = self.line
            yield draw_line(p1, p2, 'align-link', style=style)

    def get_random_string(self, length):
        """ Generates random string to nameless trees """
        letters = string.ascii_lowercase
        result_str = ''.join(random.choice(letters) for i in range(length))
        return result_str


    def compute_bounding_box(self,
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before):

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

        if pos == 'branch_right':
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
        style = {
            'fill': self.stroke_color,
            'opacity': self.opacity,
            'stroke': self.stroke_color,
            'stroke-width': self.stroke_width
            }
        if circ_drawer:
            rect_id = self.get_random_string(10)
            style['id'] = rect_id

        yield draw_rect(self._box,
                self.name,
                style=style,
            )


class LayoutText(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, width=70, min_fsize=5, max_fsize=15, padding_x=1, padding_y=0, legend=True, aligned_faces=True):
        super().__init__(name, aligned_faces=aligned_faces)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'
        self.legend = legend
        self.width = width
        self.height = None
        self.min_fsize = min_fsize 
        self.max_fsize = max_fsize
        self.absence_color = "#EBEBEB"
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        if node.is_leaf and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                if self.color_dict:
                    prop_face = TextFace(prop_text, color=self.color_dict.get(prop_text, 'black'),min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width )
                else:
                    prop_face = TextFace(prop_text, color='black', min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width )
            node.add_face(prop_face, column=self.column, position="aligned")
            
        elif node.is_leaf and node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=True)
            # piechart_face = get_piechartface(node, self.internal_prop, self.color_dict, radius=25)
            # node.add_face(piechart_face, column = self.column, position = "branch_right", collapsed_only=False)
            # node.add_face(piechart_face, column = self.column, position = "branch_right", collapsed_only=True)
        else:
            #prop_face = CircleFace(radius=self.radius, color='grey', padding_x=self.padding_x, padding_y=self.padding_y)
            prop_face = RectFace(width=self.width, height=self.height, color=self.absence_color, \
                    padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
            node.add_face(prop_face, column=self.column, position="aligned", collapsed_only=True)

class LayoutColorbranch(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, legend=True, width=70, padding_x=1, padding_y=0):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'
        self.legend = legend
        self.height = None
        self.absence_color = "#EBEBEB"
        self.width = width
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=5, max_fsize=15, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        prop_text = node.props.get(self.text_prop)
        if prop_text is not None:
            if type(prop_text) == list:
                prop_text = ",".join(prop_text)
            else:
                pass

            if self.color_dict:
                if node.is_leaf:
                    node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                    padding_x=self.padding_x),column=0, position="branch_right")
                node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                padding_x=self.padding_x),column=self.column, position="branch_right", collapsed_only=True)

                node.sm_style["hz_line_color"] = self.color_dict.get(prop_text,"")
                node.sm_style["hz_line_width"] = 3
                node.sm_style["vt_line_color"] = self.color_dict.get(prop_text,"")
                node.sm_style["vt_line_width"] = 3
                node.sm_style['outline_color'] = self.color_dict.get(prop_text,"")
                node.add_face(RectFace(width=self.width, height=None, color=self.absence_color, \
                    padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None),column=self.column, position="aligned")
        
        elif node.is_leaf and node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=False)


        if node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=True)
        
        else:
            prop_face = RectFace(width=self.width, height=self.height, color=self.absence_color, \
                    padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
            node.add_face(prop_face, column=self.column, position="aligned", collapsed_only=True)

class LayoutRect(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, width=70, height=None, padding_x=1, padding_y=0, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.absence_color = "#EBEBEB"
        self.internal_prop = add_suffix(text_prop, 'counter')
        self.legend = legend
        self.width = width
        self.height = height
        self.min_fsize = 5
        self.max_fsize = 15
        self.padding_x = padding_x
        self.padding_y = padding_y
    
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            if self.color_dict:
                self.color_dict['NA'] = self.absence_color
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )
            else:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap={'NA':self.absence_color}
                                    )
    def set_node_style(self, node):
        if node.is_leaf:
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                tooltip = ""
                if node.name:
                    tooltip += f'<b>{node.name}</b><br>'
                if self.text_prop:
                    tooltip += f'<br>{self.text_prop}: {prop_text}<br>'
                
                if self.color_dict:
                    color = self.color_dict.get(str(prop_text),"")
                    prop_face = RectFace(width=self.width, height=self.height, color=color, \
                        padding_x=self.padding_x , padding_y=self.padding_y, tooltip=tooltip)
                    node.add_face(prop_face, column=self.column, position="aligned")
            else:
                #prop_face = CircleFace(radius=self.radius, color='grey', padding_x=self.padding_x, padding_y=self.padding_y)
                prop_face = RectFace(width=self.width, height=self.height, text="NA", color=self.absence_color, \
                        padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
                node.add_face(prop_face, column=self.column, position="aligned")
        
        elif node.is_leaf and node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=True)

        else:
            
            prop_face = RectFace(width=self.width, height=self.height, color=self.absence_color, \
                    padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
            node.add_face(prop_face, column=self.column, position="aligned", collapsed_only=True)

class LayoutPiechart(TreeLayout):
    def __init__(self, name, color_dict, text_prop, radius=20, padding_x=1, padding_y=0, legend=True, aligned_faces=True):
        super().__init__(name, aligned_faces=aligned_faces)
        self.aligned_faces = True

        self.text_prop = text_prop+"_counter"
        self.internal_prop = text_prop+'_counter'
        self.color_dict = color_dict
        self.radius = radius
        self.padding_x = padding_x
        self.padding_y = padding_y
        self.legend = legend
    
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        if not node.is_leaf:
            if node.props.get(self.internal_prop):
                piechart_face = get_piechartface(node, self.internal_prop, self.color_dict, radius=self.radius)
                node.add_face(piechart_face, column = 1, position = "branch_right", collapsed_only=False)
                node.add_face(piechart_face, column = 1, position = "branch_right", collapsed_only=True)

class LayoutBackground(TreeLayout):
    def __init__(self, name, color_dict, text_prop, width=70, column=0, 
        padding_x=1, padding_y=0, legend=True, aligned_faces=True):
        super().__init__(name, aligned_faces=aligned_faces)

        self.aligned_faces = True

        self.text_prop = text_prop
        self.internal_prop = text_prop+'_counter'
        self.color_dict = color_dict
        self.width = width
        self.column = column
        self.padding_x = padding_x
        self.padding_y = padding_y
        self.absence_color = "#EBEBEB"
        self.legend = legend
    
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            if self.color_dict:
                self.color_dict['NA'] = self.absence_color
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        prop_text = node.props.get(self.text_prop)
        if prop_text:
            
            if type(prop_text) == list:
                prop_text = ",".join(prop_text)
            else:
                pass
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.text_prop:
                tooltip += f'<br>{self.text_prop}: {prop_text}<br>'
            
            if self.color_dict:
                color = self.color_dict.get(str(prop_text), self.absence_color)
                align_link_face = AlignLinkFace(width=self.width*2, height=None,
                    stroke_color=color, stroke_width=2, line_type=0, opacity=0.7)
                node.sm_style["bgcolor"] = color
                node.sm_style["fgopacity"] = 0.7
                node.sm_style['outline_color'] = color
                if node.is_leaf:
                    node.add_face(align_link_face,
                        position='branch_right',
                        column=0,
                        )
        elif node.is_leaf and node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=False)

        if node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=True)
        # else:
        #     prop_face = RectFace(width=self.width, height=None, color=self.absence_color, \
        #             padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
        #     node.add_face(prop_face, column=self.column, position="aligned", collapsed_only=True)

class LayoutBubbleCategorical(TreeLayout):
    def __init__(self, name=None, prop=None, position="branch_right", 
            column=0, color_dict=None, 
            max_radius=1, padding_x=2, padding_y=0, 
            scale=True, legend=True, active=True):

        name = name or f'Barplot_{size_prop}_{color_prop}'
        super().__init__(name)

        self.aligned_faces = True
        self.prop = prop
        self.internal_prop = add_suffix(prop, 'counter')
        
        self.column = column
        self.position = position
        self.color_dict = color_dict
        self.absence_color = "#EBEBEB"
        self.max_radius = float(max_radius)
        self.fgopacity = 0.8

        self.padding_x = padding_x
        self.padding_y = padding_y
   
        self.legend = legend
        self.active = active

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            if self.color_dict:
                self.color_dict['NA'] = self.absence_color
                tree_style.add_legend(title=self.prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )
            else:
                tree_style.add_legend(title=self.prop,
                                    variable='discrete',
                                    colormap={'NA':self.absence_color}
                                    )

    def set_node_style(self, node):
        prop_text = node.props.get(self.prop)
        if prop_text is not None:
            if type(prop_text) == list:
                prop_text = ",".join(prop_text)
            else:
                pass
            if self.color_dict:
                bubble_color = self.color_dict.get(prop_text, self.absence_color)
                bubble_size = self.max_radius
                # node.sm_style["fgcolor"] = bubble_color
                # node.sm_style["size"] = bubble_size
                # node.sm_style["fgopacity"] = self.fgopacity
                prop_face = CircleFace(radius=bubble_size, color=bubble_color, 
                padding_x=self.padding_x, padding_y=self.padding_y)
                node.add_face(prop_face, column=0, 
                position="branch_right", collapsed_only=False)

        
