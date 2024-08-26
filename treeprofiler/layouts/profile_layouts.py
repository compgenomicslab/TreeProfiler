from io import StringIO 
from collections import OrderedDict, namedtuple
import numpy as np
import math
import re

from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import (RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, \
                            SelectedFace, SelectedCircleFace, SelectedRectFace, LegendFace,
                            SeqFace, Face, AlignmentFace)
from ete4.smartview.renderer.draw_helpers import draw_text, draw_line, draw_array
from ete4 import SeqGroup
from treeprofiler.layouts.general_layouts import get_piechartface, get_heatmapface
from treeprofiler.src.utils import get_consensus_seq


Box = namedtuple('Box', 'x y dx dy')  # corner and size of a 2D shape

profilecolors = {
    'A':"#301515" ,
    'R':"#145AFF" ,
    'N':"#00DCDC" ,
    'D':"#E60A0A" ,
    'C':"#E6E600" ,
    'Q':"#00DCDC" ,
    'E':"#E60A0A" ,
    'G':"#EBEBEB" ,
    'H':"#8282D2" ,
    'I':"#0F820F" ,
    'S':"#FA9600" ,
    'K':"#145AFF" ,
    'M':"#E6E600" ,
    'F':"#3232AA" ,
    'P':"#DC9682" ,
    'L':"#0F820F" ,
    'T':"#FA9600" ,
    'W':"#B45AB4" ,
    'Z':"#FF69B4" ,
    'V':"#0F820F" ,
    'B':"#FF69B4" ,
    'Y':"#3232AA" ,
    'X':"#BEA06E",
    # '.':"#FFFFFF",
    '-':"#cccce8",
    # '-': "#EBEBEB",
    }
gradientscolor = {
    'z': '#D3D3D3', # absence lightgray 
    'x': '#000000', # black
    '-': '#ffffff', # white
    'a': '#ffede5', # lightest -> darkest reds (a->t) 
    'b': '#fee5d9', 'c': '#fedbcc', 'd': '#fdcdb9',
    'e': '#fcbfa7', 'f': '#fcaf93', 'g': '#fca082',
    'h': '#fc9070', 'i': '#fc8161', 'j': '#fb7252',
    'k': '#f96044', 'l': '#f44f39', 'm': '#f03d2d',
    'n': '#e32f27', 'o': '#d52221', 'p': '#c7171c',
    'q': '#b81419', 'r': '#aa1016', 's': '#960b13',
    't': '#7e0610'}

# Draw categorical/numerical matrix as MSA using ProfileAlignmentFace
class LayoutPropsMatrix(TreeLayout):
    def __init__(self, name="Profile", matrix_type='categorical', alignment=None, \
            matrix_props=None, width=None, profiles=None, 
            poswidth=20, height=20, column=0, range=None, \
            summarize_inner_nodes=False, value_range=[], value_color={}, \
            legend=True, active=True):

        super().__init__(name, active=active)
        self.alignment = SeqGroup(alignment) if alignment else None
        self.matrix_type = matrix_type
        self.matrix_props = matrix_props
        self.profiles = profiles

        if width:
            self.width = width
        else:
            self.width = poswidth * len(matrix_props)

        self.height = height
        self.column = column
        self.aligned_faces = True

        self.length = len(next(self.alignment.iter_entries())[1]) if self.alignment else None
        self.scale_range = range or (0, self.length)
        self.value_range = value_range
        self.value_color = value_color

        self.summarize_inner_nodes = summarize_inner_nodes
        self.legend = legend

    def set_tree_style(self, tree, tree_style):
        if self.length:
            #face = TextScaleFace(width=self.width, scale_range=self.scale_range, 
            #                    headers=self.matrix_props, padding_y=0, rotation=270)
            face = MatrixScaleFace(width=self.width, scale_range=(0, self.length), padding_y=0)
            header = self.matrix_props[0]
            title = TextFace(header, min_fsize=5, max_fsize=12, 
                    padding_x=self.width/2, padding_y=2, width=self.width/2)
            tree_style.aligned_panel_header.add_face(face, column=self.column)
            tree_style.aligned_panel_header.add_face(title, column=self.column)
            if self.matrix_type == 'categorical':
                colormap = {value: profilecolors[letter] for value, letter in self.value_color.items()}
                tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=colormap,
                                    )

            if self.matrix_type == 'numerical':
                max_color = gradientscolor['t']
                min_color = gradientscolor['a']
                color_range = [max_color, min_color]
                tree_style.add_legend(title=self.name,
                                    variable='continuous',
                                    value_range=self.value_range,
                                    color_range=color_range,
                                    )
    def _get_seq(self, node):
        if self.alignment:
            return self.alignment.get_seq(node.name)
        return node.props.get("seq", None)

    def get_seq(self, node):
        if node.is_leaf:
            return self._get_seq(node)
        if self.summarize_inner_nodes:
            # TODO: summarize inner node's seq
            matrix = ''
            for leaf in node.leaves():
                matrix += ">"+leaf.name+"\n"
                matrix += self._get_seq(leaf)+"\n"
            
            try:
                if self.mode == "numerical":
                    consensus_seq = get_consensus_seq(matrix, 0.1)
                elif self.mode == 'profiles':
                    consensus_seq = get_consensus_seq(matrix, 0.7)
                else:
                    consensus_seq = get_consensus_seq(matrix, 0.7)
                return str(consensus_seq)
            except ValueError:
                return None
        else:
            first_leaf = next(node.leaves())
            return self._get_seq(first_leaf)
    
    def set_node_style(self, node):
        seq = self.get_seq(node)
        if len(self.profiles) > 1:
            poswidth = self.width / (len(self.profiles)-1 )
        else:
            poswidth = self.width

        if seq:
            seqFace = ProfileAlignmentFace(seq, gap_format=None, seqtype='aa', 
            seq_format=self.matrix_type, width=self.width, height=self.height, 
            poswidth=poswidth,
            fgcolor='black', bgcolor='#bcc3d0', gapcolor='gray',
            gap_linewidth=0.2,
            max_fsize=12, ftype='sans-serif', 
            padding_x=0, padding_y=0)
            node.add_face(seqFace, column=self.column, position='aligned', 
                    collapsed_only=(not node.is_leaf))

# Draw presence/absence matrix as MSA using ProfileAlignmentFace
class LayoutProfile(TreeLayout):
    def __init__(self, name="Profile", mode='profiles',
            alignment=None, seq_format='profiles', profiles=None, 
            width=None, poswidth=20, height=20,
            column=0, range=None, summarize_inner_nodes=False, 
            value_range=[], value_color={}, legend=True,
            active=True):
        super().__init__(name, active=active)
        self.alignment = SeqGroup(alignment) if alignment else None
        self.mode = mode
        
        if width:
            self.width = width
        else:
            self.width = poswidth * len(profiles)
        #total_width = self.seqlength * self.poswidth
        # if self.width:
        #     self.w_scale = self.width / total_width
        # else:
        #     self.width = total_width
        
        self.height = height
        self.column = column
        self.aligned_faces = True
        self.seq_format = seq_format
        self.profiles = profiles

        self.length = len(next(self.alignment.iter_entries())[1]) if self.alignment else None
        self.scale_range = range or (0, self.length)
        self.value_range = value_range
        self.value_color = value_color

        self.summarize_inner_nodes = summarize_inner_nodes
        self.legend = legend

    def set_tree_style(self, tree, tree_style):
        if self.length:
            face = TextScaleFace(width=self.width, scale_range=self.scale_range, 
                                headers=self.profiles, padding_y=0, rotation=270)
            tree_style.aligned_panel_header.add_face(face, column=self.column)

        if self.legend:
            if self.mode == 'profiles':
                color_dict = {}
                for i in range(len(self.profiles)):
                    profile_val = self.profiles[i]
                    #profile_color = profilecolors[list(profilecolors.keys())[i % len(profilecolors)]]
                    color_dict[profile_val] = ''

                tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=color_dict,
                                    )

    def _get_seq(self, node):
        if self.alignment:
            return self.alignment.get_seq(node.name)
        return node.props.get("seq", None)

    def get_seq(self, node):
        if node.is_leaf:
            return self._get_seq(node)
        if self.summarize_inner_nodes:
            # TODO: summarize inner node's seq
            matrix = ''
            for leaf in node.leaves():
                matrix += ">"+leaf.name+"\n"
                matrix += self._get_seq(leaf)+"\n"
            
            try:
                if self.mode == "numerical":
                    consensus_seq = get_consensus_seq(matrix, 0.1)
                elif self.mode == 'profiles':
                    consensus_seq = get_consensus_seq(matrix, 0.7)
                else:
                    consensus_seq = get_consensus_seq(matrix, 0.7)
                return str(consensus_seq)
            except ValueError:
                return None
        else:
            first_leaf = next(node.leaves())
            return self._get_seq(first_leaf)
    
    def set_node_style(self, node):
        
        seq = self.get_seq(node)
        if len(self.profiles) > 1:
            poswidth = self.width / (len(self.profiles)-1 )
        else:
            poswidth = self.width

        if seq:
            seqFace = ProfileAlignmentFace(seq, gap_format=None, seqtype='aa', 
            seq_format=self.seq_format, width=self.width, height=self.height, 
            poswidth=poswidth,
            fgcolor='black', bgcolor='#bcc3d0', gapcolor='gray',
            gap_linewidth=0.2,
            max_fsize=12, ftype='sans-serif', 
            padding_x=0, padding_y=0)
            node.add_face(seqFace, column=self.column, position='aligned', 
                    collapsed_only=(not node.is_leaf)) 

# Draw presence/absence, categorical/numerical matrix as drawing array using ProfileFace
class LayoutPropsMatrixOld(TreeLayout):
    def __init__(self, name="Profile", matrix=None, matrix_type='categorical', \
            matrix_props=None, is_list=False, width=None, poswidth=20, height=20,
            column=0, range=None, summarize_inner_nodes=False, value_range=[], \
            value_color={}, legend=True, active=True):
        super().__init__(name, active=active)
        self.matrix = matrix
        self.matrix_type = matrix_type
        self.matrix_props = matrix_props
        self.is_list = is_list

        if width:
            self.width = width
        else:
            self.width = poswidth * len(matrix_props)

        self.height = height
        self.column = column
        self.aligned_faces = True

        self.length = len(next((value for value in self.matrix.values() if value != [None]), None)) if any(value != [None] for value in self.matrix.values()) else 0
        self.scale_range = range or (0, self.length)
        self.value_range = value_range
        self.value_color = value_color

        self.summarize_inner_nodes = summarize_inner_nodes
        self.legend = legend

    def set_tree_style(self, tree, tree_style):
        if self.length:
            if self.is_list:
                # first not None list to set the column
                ncols = len(next((value for value in self.matrix.values() if value != [None]), None)) if any(value != [None] for value in self.matrix.values()) else 0
                if ncols > 1:
                    total_width = self.width * (ncols-1)
                else:
                    total_width = self.width
                face = MatrixScaleFace(width=total_width, scale_range=(0, ncols), padding_y=0)
                header = self.matrix_props
                title = TextFace(header, min_fsize=5, max_fsize=12, 
                    padding_x=0, padding_y=2, width=self.width)
                tree_style.aligned_panel_header.add_face(face, column=self.column)
                tree_style.aligned_panel_header.add_face(title, column=self.column)
                
            else:
                face = TextScaleFace(width=self.width, scale_range=self.scale_range, 
                                    headers=self.matrix_props, padding_y=0, rotation=270)
                tree_style.aligned_panel_header.add_face(face, column=self.column)

        if self.legend:
            if self.matrix_type == 'numerical':
                keys_list = list(self.value_color.keys())  
                middle_index = len(keys_list) // 2  
                middle_key = keys_list[middle_index]  
                middle_value = self.value_color[middle_key]  

                if self.value_range:
                    color_gradient = [
                        self.value_color[self.value_range[1]], 
                        middle_value,
                        self.value_color[self.value_range[0]]
                        ]
                    tree_style.add_legend(title=self.name,
                                    variable="continuous",
                                    value_range=self.value_range,
                                    color_range=color_gradient,
                                    )
            if self.matrix_type == 'categorical':
                tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=self.value_color,
                                    )
    def _get_array(self, node):
        if self.matrix:
            return self.matrix.get(node.name)

    # def get_array(self, node):
    #     if node.is_leaf:
    #         return self._get_array(node)
    #     else:
    #         first_leaf = next(node.leaves())
    #         return self._get_array(first_leaf)

    def get_array(self, node):
        if self.matrix.get(node.name):
            return self.matrix.get(node.name)
        else:
            first_leaf = next(node.leaves())
            return self._get_array(first_leaf)

    def set_node_style(self, node):
        array = self.get_array(node)
        #array = self.get_array(node)
        if array:
            if not self.is_list:
                if len(self.matrix_props) > 1:
                    poswidth = self.width / (len(self.matrix_props) - 1)
                else:
                    poswidth = self.width
                if array:
                    profileFace = ProfileFace(array, self.value_color, gap_format=None, \
                    seq_format=self.matrix_type, width=self.width, height=self.height, \
                    poswidth=poswidth, tooltip=True)
                    node.add_face(profileFace, column=self.column, position='aligned', \
                        collapsed_only=(not node.is_leaf))
            else:
                poswidth = self.width * len(array)

                profileFace = ProfileFace(array, self.value_color, gap_format=None, \
                    seq_format=self.matrix_type, width=poswidth, height=self.height, \
                    poswidth=poswidth, tooltip=True)
                node.add_face(profileFace, column=self.column, position='aligned', \
                    collapsed_only=(not node.is_leaf))

class LayoutPropsMatrixBinary(TreeLayout):
    def __init__(self, name="Binary_profiling", matrix=None,  
            matrix_props=None, is_list=False, width=None, 
            poswidth=20, height=20,
            column=0, range=None, summarize_inner_nodes=False, 
            value_range=[], 
            value_color={}, legend=True, active=True):
        super().__init__(name, active=active)
        self.matrix = matrix
        self.matrix_props = matrix_props
        self.is_list = is_list

        if width:
            self.width = width
        else:
            self.width = poswidth * len(matrix_props)

        self.height = height
        self.column = column
        self.aligned_faces = True

        self.length = len(next((value for value in self.matrix.values() if value != [None]), None)) if any(value != [None] for value in self.matrix.values()) else 0
        self.scale_range = range or (0, self.length)
        self.value_range = value_range
        self.value_color = value_color

        self.summarize_inner_nodes = summarize_inner_nodes
        self.legend = legend

    def set_tree_style(self, tree, tree_style):
        if self.length:
            if self.is_list:
                # first not None list to set the column
                ncols = len(next((value for value in self.matrix.values() if value != [None]), None)) if any(value != [None] for value in self.matrix.values()) else 0
                face = MatrixScaleFace(width=self.width, scale_range=(0, ncols), padding_y=0)
                header = self.matrix_props[0]
                title = TextFace(header, min_fsize=5, max_fsize=12, 
                    padding_x=0, padding_y=2, width=self.width)
                tree_style.aligned_panel_header.add_face(face, column=self.column)
                tree_style.aligned_panel_header.add_face(title, column=self.column)
                
            else:
                face = TextScaleFace(width=self.width, scale_range=self.scale_range, 
                                    headers=self.matrix_props, padding_y=0, rotation=270)
                tree_style.aligned_panel_header.add_face(face, column=self.column)

        if self.legend:
            keys_list = list(self.value_color.keys())  
            middle_index = len(keys_list) // 2  
            middle_key = keys_list[middle_index]  
            middle_value = self.value_color[middle_key]  

            if self.value_range:
                color_gradient = [
                    self.value_color[self.value_range[1]], 
                    middle_value,
                    self.value_color[self.value_range[0]]
                    ]
                    
            tree_style.add_legend(title=self.name,
                                variable='continuous',
                                value_range=self.value_range,
                                color_range=color_gradient
                                )

    def _get_array(self, node):
        if self.matrix:
            return self.matrix.get(node.name)

    # def get_array(self, node):
    #     if node.is_leaf:
    #         return self._get_array(node)
    #     else:
    #         first_leaf = next(node.leaves())
    #         return self._get_array(first_leaf)

    def get_array(self, node):
        if self.matrix.get(node.name):
            return self.matrix.get(node.name)
        else:
            first_leaf = next(node.leaves())
            return self._get_array(first_leaf)

    def set_node_style(self, node):
        array = self.get_array(node)

        if len(self.matrix_props) > 1:
            poswidth = self.width / (len(self.matrix_props)-1 )
        else:
            poswidth = self.width
        
        if array:
            profileFace = ProfileFace(array, self.value_color, gap_format=None, \
            seq_format='numerical', width=self.width, height=self.height, \
            poswidth=poswidth, tooltip=True)
            node.add_face(profileFace, column=self.column, position='aligned', \
                collapsed_only=(not node.is_leaf))



#Faces
class TextScaleFace(Face):
    def __init__(self, name='', width=None, color='black',
            scale_range=(0, 0), headers=None, tick_width=100, line_width=1,
            formatter='%.0f', 
            min_fsize=10, max_fsize=15, ftype='sans-serif',
            padding_x=0, padding_y=0, rotation=0):

        Face.__init__(self, name=name,
                padding_x=padding_x, padding_y=padding_y)

        self.width = width
        self.height = None
        self.range = scale_range
        self.headers = headers

        self.color = color
        self.min_fsize = min_fsize
        self.max_fsize = max_fsize
        self._fsize = max_fsize
        self.ftype = ftype
        self.formatter = formatter
        self.rotation=rotation
        
        self.tick_width = tick_width
        self.line_width = line_width

        self.vt_line_height = 10

    def __name__(self):
        return "ScaleFace"

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

        x, y, _, dy = box
        zx, zy = self.zoom

        self.viewport = (drawer.viewport.x, drawer.viewport.x + drawer.viewport.dx)

        self.height = (self.line_width + 10 + self.max_fsize) / zy
        height = min(dy, self.height)
        if pos == "aligned_bottom":
            y = y + dy - height

        self._box = Box(x, y, self.width / zx, height)
        return self._box

    def draw(self, drawer):
        x0, y, _, dy = self._box
        zx, zy = self.zoom
        p1 = (x0, y + dy - 5 / zy)
        p2 = (x0 + self.width, y + dy - self.vt_line_height / (2 * zy))
        if drawer.TYPE == 'circ':
            p1 = cartesian(p1)
            p2 = cartesian(p2)
        # yield draw_line(p1, p2, style={'stroke-width': self.line_width,
        #                                'stroke': self.color})


        #nticks = round((self.width * zx) / self.tick_width)
        
        if len(self.headers) > 1:
            nticks = len(self.headers)
        else:
            nticks = 1
        dx = self.width / nticks
        range_factor = (self.range[1] - self.range[0]) / self.width

        if self.viewport:
            sm_start = round(max(self.viewport[0] - self.viewport_margin - x0, 0) / dx)
            sm_end = nticks - round(max(x0 + self.width - (self.viewport[1] +
                self.viewport_margin), 0) / dx)
        else:
            sm_start, sm_end = 0, nticks
            
        for i in range(sm_start, sm_end + 1):
            x = x0 + i * dx + dx/2 
            
            number = range_factor * i * dx
            if number == 0:
                text = "0"
            else:
                text = self.formatter % number if self.formatter else str(number)

            #text = text.rstrip('0').rstrip('.') if '.' in text else text
            try:
                text = self.headers[i]
                self.compute_fsize(self.tick_width / len(text), dy, zx, zy)
                text_style = {
                    'max_fsize': self._fsize,
                    'text_anchor': 'left', # left, middle or right
                    'ftype': f'{self.ftype}, sans-serif', # default sans-serif
                    }
                    
                
                text_box = Box(x,
                        y,
                        # y + (dy - self._fsize / (zy * r)) / 2,
                        dx, dy)
                yield draw_text(text_box, text, style=text_style,rotation=self.rotation)

                p1 = (x, y + dy - self.vt_line_height / zy)
                p2 = (x, y + dy)

                yield draw_line(p1, p2, style={'stroke-width': self.line_width,
                                               'stroke': self.color})
            except IndexError:
                break
                
class ProfileAlignmentFace(Face):
    def __init__(self, seq, bg=None,
            gap_format='line', seqtype='aa', seq_format='profiles',
            width=None, height=None, # max height
            fgcolor='black', bgcolor='#bcc3d0', gapcolor='gray',
            gap_linewidth=0.2,
            max_fsize=12, ftype='sans-serif', poswidth=5,
            padding_x=0, padding_y=0):

        Face.__init__(self, padding_x=padding_x, padding_y=padding_y)

        self.seq = seq
        self.seqlength = len(self.seq)
        self.seqtype = seqtype

        self.autoformat = True  # block if 1px contains > 1 tile

        self.seq_format = seq_format
        self.gap_format = gap_format
        self.gap_linewidth = gap_linewidth
        self.compress_gaps = False

        self.poswidth = poswidth
        self.w_scale = 1
        self.width = width    # sum of all regions' width if not provided
        self.height = None  # dynamically computed if not provided

        total_width = self.seqlength * self.poswidth
        if self.width:
            self.w_scale = self.width / total_width
        else:
            self.width = total_width
        self.bg = profilecolors
        # self.fgcolor = fgcolor
        # self.bgcolor = bgcolor
        self.gapcolor = gapcolor

        # Text
        self.ftype = ftype
        self._min_fsize = 8
        self.max_fsize = max_fsize
        self._fsize = None

        self.blocks = []
        self.build_blocks()

    def __name__(self):
        return "AlignmentFace"

    def get_seq(self, start, end):
        """Retrieves sequence given start, end"""
        return self.seq[start:end]

    def build_blocks(self):
        pos = 0
        for reg in re.split('([^-]+)', self.seq):
            if reg:
                if not reg.startswith("-"):
                    self.blocks.append([pos, pos + len(reg) - 1])
                pos += len(reg)

        self.blocks.sort()

    def compute_bounding_box(self,
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before):

        if pos != 'branch_right' and not pos.startswith('aligned'):
            raise InvalidUsage(f'Position {pos} not allowed for SeqMotifFace')

        box = super().compute_bounding_box(
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before)

        x, y, _, dy = box

        zx, zy = self.zoom
        zx = 1 if drawer.TYPE != 'circ' else zx

            # zx = drawer.zoom[0]
            # self.zoom = (zx, zy)

        if drawer.TYPE == "circ":
            self.viewport = (0, drawer.viewport.dx)
        else:
            self.viewport = (drawer.viewport.x, drawer.viewport.x + drawer.viewport.dx)

        self._box = Box(x, y, self.width / zx, dy)
        return self._box

    # def get_fsize(self, dx, dy, zx, zy, max_fsize=None):
    #     return min([dx * zx * CHAR_HEIGHT, abs(dy * zy), max_fsize or 4])

    def draw(self, drawer):
        def get_height(x, y):
            r = (x or 1e-10) if drawer.TYPE == 'circ' else 1
            default_h = dy * zy * r
            h = min([self.height or default_h, default_h]) / zy
            # h /= r
            return y + (dy - h) / 2, h

        # Only leaf/collapsed branch_right or aligned
        x0, y, dx, dy = self._box
        zx, zy = self.zoom
        zx = drawer.zoom[0] if drawer.TYPE == 'circ' else zx


        if self.gap_format in ["line", "-"]:
            p1 = (x0, y + dy / 2)
            p2 = (x0 + self.width, y + dy / 2)
            if drawer.TYPE == 'circ':
                p1 = cartesian(p1)
                p2 = cartesian(p2)
            yield draw_line(p1, p2, style={'stroke-width': self.gap_linewidth,
                                           'stroke': self.gapcolor})
        vx0, vx1 = self.viewport
        too_small = (self.width * zx) / (self.seqlength) < 1

        posw = self.poswidth * self.w_scale
        viewport_start = vx0 - self.viewport_margin / zx
        viewport_end = vx1 + self.viewport_margin / zx
        sm_x = max(viewport_start - x0, 0)
        sm_start = round(sm_x / posw)
        w = self.seqlength * posw
        sm_x0 = x0 if drawer.TYPE == "rect" else 0
        sm_end = self.seqlength - round(max(sm_x0 + w - viewport_end, 0) / posw)
        
        # if too_small:
        #     # for start, end in self.blocks:
        #     #     if end >= sm_start and start <= sm_end:
        #     #         bstart = max(sm_start, start)
        #     #         bend = min(sm_end, end)
        #     #         bx = x0 + bstart * posw
        #     #         by, bh = get_height(bx, y)
        #     #         box = Box(bx, by, (bend + 1 - bstart) * posw, bh)
                    
        #     #         yield [ "pixi-block", box ]
            
        #     # total position of columns
        #     seq = self.get_seq(sm_start, sm_end)

        #     # starting point
        #     sm_x = sm_x if drawer.TYPE == 'rect' else x0
        #     # get height and how high the box should be
        #     y, h = get_height(sm_x, y)
        #     # create a box
        #     sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)

        #     yield draw_array(sm_box,[gradientscolor[x] for x in seq])
        if self.seq_format == "numerical":
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)
            
            # fsize = self.get_fsize(dx / len(seq), dy, zx, zy, 20)
            # style = {
            #     'fill': "black",
            #     'max_fsize': fsize,
            #     'ftype': 'sans-serif', # default sans-serif
            #    }

            yield [ f'pixi-gradients', sm_box, seq]
            # yield draw_array(sm_box,[gradientscolor[x] for x in seq])
            #yield draw_text(sm_box, for i in seq, "jjj", style=style)
            
        elif self.seq_format == "categorical": 
            
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)
            # aa_type = "notext"
            # yield [ f'pixi-aa_{aa_type}', sm_box, seq ]
            yield draw_array(sm_box, [profilecolors[x] for x in seq])

        else: # when is "profiles":
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)
            
            if self.seq_format == 'profiles' or posw * zx < self._min_fsize:
                aa_type = "notext"
                tooltip = f'<p>{seq}</p>'
                style = {
                    'fill': "black",
                    'max_fsize': 14,
                    'ftype': 'sans-serif', # default sans-serif
                   }
                # yield draw_array(sm_box, [gradientscolor[x] for x in seq], tooltip=tooltip)            
                # yield [ f'pixi-aa_{aa_type}', sm_box, seq ]
                yield [ f'pixi-gradients', sm_box, seq]
            # else:
            #     aa_type = "text"
            #     yield [ f'pixi-aa_{aa_type}', sm_box, seq ]
            # sm_x0 = sm_x0 + posw/2 - zx*2 # centering text in the middle of the box
            # for i in range(len(seq)):
            #     sm_box = Box(sm_x+sm_x0+(posw * i), y, posw, h)
            #     yield draw_text(sm_box, seq[i], "jjj", style=style)

class ProfileFace(Face):
    def __init__(self, seq, value2color=None,
            gap_format='line', seq_format='categorical', # profiles, numerical, categorical
            width=None, height=None, # max height
            gap_linewidth=0.2,
            max_fsize=12, ftype='sans-serif', poswidth=5,
            padding_x=0, padding_y=0, tooltip=True):

        Face.__init__(self, padding_x=padding_x, padding_y=padding_y)

        self.seq = seq
        self.seqlength = len(self.seq)
        self.value2color = value2color
        self.absence_color = '#EBEBEB'

        self.autoformat = True  # block if 1px contains > 1 tile

        self.seq_format = seq_format
        self.gap_format = gap_format
        self.gap_linewidth = gap_linewidth
        self.compress_gaps = False

        self.poswidth = poswidth
        self.w_scale = 1
        self.width = width    # sum of all regions' width if not provided
        self.height = None  # dynamically computed if not provided
        self.tooltip = tooltip

        total_width = self.seqlength * self.poswidth
        if self.width:
            self.w_scale = self.width / total_width
        else:
            self.width = total_width


        # Text
        self.ftype = ftype
        self._min_fsize = 8
        self.max_fsize = max_fsize
        self._fsize = None

        self.blocks = []
        self.build_blocks()
    
    def __name__(self):
        return "ProfileFace"

    def get_seq(self, start, end):
        """Retrieves sequence given start, end"""
        return self.seq[start:end]

    def build_blocks(self):
        pos = 0
        for reg in self.seq:
            reg = str(reg)
            if reg:
                if not reg.startswith("-"):
                    self.blocks.append([pos, pos + len(reg) - 1])
                pos += len(reg)

        self.blocks.sort()

    def compute_bounding_box(self,
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before):

        if pos != 'branch_right' and not pos.startswith('aligned'):
            raise InvalidUsage(f'Position {pos} not allowed for Profile')

        box = super().compute_bounding_box(
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before)

        x, y, _, dy = box

        zx, zy = self.zoom
        zx = 1 if drawer.TYPE != 'circ' else zx

            # zx = drawer.zoom[0]
            # self.zoom = (zx, zy)

        if drawer.TYPE == "circ":
            self.viewport = (0, drawer.viewport.dx)
        else:
            self.viewport = (drawer.viewport.x, drawer.viewport.x + drawer.viewport.dx)

        self._box = Box(x, y, self.width / zx, dy)
        return self._box
    
    def get_rep_numbers(self, num_array, rep_num):
        def find_closest(numbers, target):
            # Find the number in 'numbers' that is closest to calculated target
            return min(numbers, key=lambda x: abs(x - target))

        # get the representative numbers from given array
        seg_size = math.ceil(len(num_array) / rep_num)
        rep_elements = []

        for i in range(rep_num):
            start_index = int(i * seg_size)
            end_index = int((i + 1) * seg_size)
            if i == rep_num - 1:
                end_index = len(num_array)

            segment = num_array[start_index:end_index]
            if segment:
                if len(segment) != 0:
                    segment_average = sum(segment) / len(segment)
                else:
                    segment_average = 0
                
                # Find the cloest number to the average
                closest_number = find_closest(segment, segment_average)
                rep_elements.append(closest_number)
        return rep_elements
            
    def draw(self, drawer):
        def get_height(x, y):
            r = (x or 1e-10) if drawer.TYPE == 'circ' else 1
            default_h = dy * zy * r
            h = min([self.height or default_h, default_h]) / zy
            # h /= r
            return y + (dy - h) / 2, h

        # Only leaf/collapsed branch_right or aligned
        x0, y, dx, dy = self._box
        zx, zy = self.zoom
        zx = drawer.zoom[0] if drawer.TYPE == 'circ' else zx


        if self.gap_format in ["line", "-"]:
            p1 = (x0, y + dy / 2)
            p2 = (x0 + self.width, y + dy / 2)
            if drawer.TYPE == 'circ':
                p1 = cartesian(p1)
                p2 = cartesian(p2)
            yield draw_line(p1, p2, style={'stroke-width': self.gap_linewidth,
                                           'stroke': self.gapcolor})
        vx0, vx1 = self.viewport
        too_small = (self.width * zx) / (self.seqlength) < 1

        posw = self.poswidth * self.w_scale
        viewport_start = vx0 - self.viewport_margin / zx
        viewport_end = vx1 + self.viewport_margin / zx
        sm_x = max(viewport_start - x0, 0)
        sm_start = round(sm_x / posw)
        w = self.seqlength * posw
        sm_x0 = x0 if drawer.TYPE == "rect" else 0
        sm_end = self.seqlength - round(max(sm_x0 + w - viewport_end, 0) / posw)
        
        # total width of the matrix: self.width
        # total number of column: self.seqlength
        # at least 1px per column

        if self.seq_format == "numerical":
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)

            if too_small: # only happens in data-matrix visualization            
                #ncols_per_px = math.ceil(self.seqlength / (zx * sm_box.dx)) #jordi's idea
                ncols_per_px = math.ceil(self.width * zx)
                rep_elements = self.get_rep_numbers(seq, ncols_per_px)
                if self.tooltip:
                    tooltip = f'<p>{seq}</p>'
                else:
                    tooltip = ''
                yield draw_array(sm_box, [self.value2color[x] if x is not None else self.absence_color for x in rep_elements], tooltip=tooltip)
            else:
                # fsize = self.get_fsize(dx / len(seq), dy, zx, zy, 20)
                # style = {
                #     'fill': "black",
                #     'max_fsize': fsize,
                #     'ftype': 'sans-serif', # default sans-serif
                #    }
                if self.tooltip:
                    tooltip = f'<p>{seq}</p>'
                else:
                    tooltip = ''
                # yield draw_text(sm_box, for i in seq, "jjj", style=style)
                yield draw_array(sm_box, [self.value2color[x] if x is not None else self.absence_color for x in seq], tooltip=tooltip)
            
            
        if self.seq_format == "categorical":
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)
            #tooltip = f'<p>{seq}</p>'
            yield draw_array(sm_box, [self.value2color[x] if x is not None else self.absence_color for x in seq])

class MatrixScaleFace(Face):
    def __init__(self, name='', width=None, color='black',
            scale_range=(0, 0), tick_width=80, line_width=1,
            formatter='%.0f',
            min_fsize=6, max_fsize=12, ftype='sans-serif',
            padding_x=0, padding_y=0):

        Face.__init__(self, name=name,
                padding_x=padding_x, padding_y=padding_y)

        self.width = width
        self.height = None
        self.range = scale_range
        self.columns = scale_range[1]

        self.color = color
        self.min_fsize = min_fsize
        self.max_fsize = max_fsize
        self._fsize = max_fsize
        self.ftype = ftype
        self.formatter = formatter

        self.tick_width = tick_width
        self.line_width = line_width

        self.vt_line_height = 10

    def __name__(self):
        return "ScaleFace"

    def compute_bounding_box(self,
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before):

        if drawer.TYPE == 'circ' and abs(point[1]) >= pi/2:
            pos = swap_pos(pos)

        box = super().compute_bounding_box(
            drawer,
            point, size,
            dx_to_closest_child,
            bdx, bdy,
            bdy0, bdy1,
            pos, row,
            n_row, n_col,
            dx_before, dy_before)

        x, y, _, dy = box
        zx, zy = self.zoom

        self.viewport = (drawer.viewport.x, drawer.viewport.x + drawer.viewport.dx)

        self.height = (self.line_width + 10 + self.max_fsize) / zy

        height = min(dy, self.height)

        if pos == "aligned_bottom":
            y = y + dy - height

        self._box = Box(x, y, self.width / zx, height)
        return self._box

    def draw(self, drawer):
        x0, y, _, dy = self._box
        zx, zy = self.zoom

        p1 = (x0, y + dy - 5 / zy)
        p2 = (x0 + self.width, y + dy - self.vt_line_height / (2 * zy))

        # count the middle point of each column
        if self.columns > 1:
            half_width_col = self.width / (self.columns-1) / 2
            p1 = (x0 + half_width_col, y + dy - 5 / zy)
            p2 = (x0 + (self.width + half_width_col), y + dy - self.vt_line_height / (2 * zy))
            
            if drawer.TYPE == 'circ':
                p1 = cartesian(p1)
                p2 = cartesian(p2)
            yield draw_line(p1, p2, style={'stroke-width': self.line_width,
                                        'stroke': self.color})
        else:
            half_width_col = self.width / 2

        


        #nticks = round((self.width * zx) / self.tick_width)
        if self.columns > 1:
            nticks = self.columns - 1
        else:
            nticks = 1
        
        
        dx = self.width / nticks
        range_factor = (self.range[1] - self.range[0]) / self.width

        if self.viewport:
            sm_start = round(max(self.viewport[0] - self.viewport_margin - x0, 0) / dx)
            sm_end = nticks - round(max(x0 + self.width - (self.viewport[1] +
                self.viewport_margin), 0) / dx)
        else:
            sm_start, sm_end = 0, nticks

        if self.columns > 1:
            for i in range(sm_start, sm_end + 1):

                x = x0 + i * dx
                
                # number = range_factor * i * dx

                # if number == 0:
                #     text = "0"
                # else:
                #     #actual_number = number + 1
                #     text = self.formatter % number if self.formatter else str(number)

                # text = text.rstrip('0').rstrip('.') if '.' in text else text

                text = str(i+1)
                self.compute_fsize(self.tick_width / len(text), dy, zx, zy)
                text_style = {
                    'max_fsize': self._fsize,
                    'text_anchor': 'middle',
                    'ftype': f'{self.ftype}, sans-serif', # default sans-serif
                    }
                text_box = Box(x+ half_width_col,
                        y,
                        # y + (dy - self._fsize / (zy * r)) / 2,
                        dx, dy)

                # column index starts from 1
                yield draw_text(text_box, text, style=text_style)

                # vertical line tick
                p1 = (x + half_width_col, y + dy - self.vt_line_height / zy)
                p2 = (x + half_width_col, y + dy)

                yield draw_line(p1, p2, style={'stroke-width': self.line_width,
                                            'stroke': self.color})
        else:
            x = x0
            text = str(1)
            self.compute_fsize(self.tick_width / len(text), dy, zx, zy)
            text_style = {
                'max_fsize': self._fsize,
                'text_anchor': 'middle',
                'ftype': f'{self.ftype}, sans-serif', # default sans-serif
                }
            text_box = Box(x+ half_width_col,
                    y,
                    # y + (dy - self._fsize / (zy * r)) / 2,
                    dx, dy)

            # column index starts from 1
            yield draw_text(text_box, text, style=text_style)

            # vertical line tick
            p1 = (x + half_width_col, y + dy - self.vt_line_height / zy)
            p2 = (x + half_width_col, y + dy)

            yield draw_line(p1, p2, style={'stroke-width': self.line_width,
                                        'stroke': self.color})