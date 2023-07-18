from io import StringIO 
from collections import OrderedDict, namedtuple
import numpy as np
import re

from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import (RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace, \
                            SelectedFace, SelectedCircleFace, SelectedRectFace, LegendFace,
                            SeqFace, Face, ScaleFace, AlignmentFace)
from ete4.smartview.renderer.draw_helpers import draw_text, draw_line, draw_array
from ete4 import SeqGroup
from treeprofiler.layouts.general_layouts import get_piechartface, get_heatmapface, color_gradient
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
    }
gradientscolor = {
    'a': '#F62D2D', 'b': '#F74444', 'c': '#F85B5B', 'd': '#F97373', 'e': '#FA8A8A', 'f': '#FBA1A1',
    'g': '#FCB9B9', 'h': '#FDD0D0', 'i': '#FEE7E7', 'j': '#e6e6f3', 'k': '#FFFFFF', 'l': '#E4E8F5',
    'm': '#C9D1EB', 'n': '#AFBBE1', 'o': '#94A4D7', 'p': '#7A8ECD', 'q': '#5F77C3', 'r': '#4561B9',
    's': '#2A4AAF', 't': '#1034A6', '-': "#EBEBEB"
}

class LayoutProfile(TreeLayout):
    def __init__(self, name="Profile", mode='single',
            alignment=None, seq_format='profiles', profiles=None, width=None, poswidth=20, height=20,
            column=0, range=None, summarize_inner_nodes=False, value_range=[], value_color={}, legend=True):
        super().__init__(name)
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
            #face = ScaleFace(width=self.width, scale_range=self.scale_range, padding_y=0)
            tree_style.aligned_panel_header.add_face(face, column=self.column)

        if self.legend:
            if self.mode == 'numerical':
                if self.value_range:
                    color_gradient = [gradientscolor['a'], gradientscolor['k'], gradientscolor['t']]
                    tree_style.add_legend(title=self.name,
                                    variable="continuous",
                                    value_range=self.value_range,
                                    color_range=color_gradient,
                                    )
            if self.mode == 'single':
                color_dict = {}
                for val, letter in self.value_color.items():
                    color_dict[val] = profilecolors[letter]
                tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=color_dict,
                                    )
            if self.mode == 'multi':
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
                    consensus_seq = get_consensus_seq(StringIO(matrix), 0.1)
                elif self.mode == 'multi':
                    consensus_seq = get_consensus_seq(StringIO(matrix), 0.7)
                else:
                    consensus_seq = get_consensus_seq(StringIO(matrix), 0.7)
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

class LayoutGOslim(TreeLayout):
    def __init__(self, name=None, column=1, min_color="#ffffff", max_color="red", go_propfile=[], goslim_prop=None, padding_x=2, padding_y=2, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.go_propfile = go_propfile
        self.goslim_prop = goslim_prop
        self.internal_prop = goslim_prop +'_counter' #GOslims_counter
        self.column = column
        self.max_color = max_color
        self.min_color = min_color
        self.padding_x = padding_x
        self.padding_y = padding_y
        self.legend = legend
        self.width = 70
        self.height = None
        self.min_fsize = 5
        self.max_fsize = 10

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        header = f'{self.go_propfile[1]}({self.go_propfile[0]})'
        text = TextFace(header, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=70, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)
    
    def set_node_style(self, node):
        entry, desc = self.go_propfile
        if node.is_leaf:
            tooltip = ""
            if node.name:
                tooltip += f'<b>{node.name}</b><br>'
            if self.go_propfile:
                tooltip += f'<br>{self.goslim_prop}: {entry}, {desc}<br>'


            if node.props.get(self.goslim_prop):
                goslims = node.props.get(self.goslim_prop)
                if entry in goslims:
                    profiling_face = RectFace(width=self.width, height=self.height, color=self.max_color,  
                    padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                    node.add_face(profiling_face, column=self.column, position = "aligned")
                else:
                    profiling_face = RectFace(width=self.width, height=self.height, color=self.min_color,  
                    padding_x=self.padding_x, padding_y=self.padding_y, tooltip=tooltip)
                    node.add_face(profiling_face, column=self.column, position = "aligned")
                # if entry in goslims[0]:
                #     index = goslims[0].index(entry)
                #     relative_count = goslims[2][index]
                #     c1 = 'white'
                #     c2 = self.color
                #     gradient_color = color_gradient(c1, c2, mix=relative_count)
                #     prop_face = RectFace(width=self.width, height=self.height, color=gradient_color,  padding_x=self.padding_x, padding_y=self.padding_y)
                #     node.add_face(prop_face, column=self.column, position = "aligned")
            elif node.props.get(self.internal_prop):
                profiling_face = self.get_profile_gradientface(node, entry, self.internal_prop, self.max_color, 
                                width=self.width, height=self.height, padding_x=self.padding_x, padding_y=self.padding_y)
                node.add_face(profiling_face, column=self.column, position = "aligned")
        else: 
            if node.props.get(self.internal_prop):
                profiling_face = self.get_profile_gradientface(node, entry, self.internal_prop, self.max_color, 
                                    width=self.width, height=self.height, padding_x=self.padding_x, padding_y=self.padding_y)
                node.add_face(profiling_face, column = self.column, position = "aligned", collapsed_only=True)

    
    def get_profile_gradientface(self, node, target_go, internal_prop, color, width, height, padding_x, padding_y, min_color="#ffffff", tooltip=None):
        counter_props = node.props.get(internal_prop).split('||')
    
        total = 0
        positive = 0

        for counter_prop in counter_props:
            k, v = counter_prop.split('--')
            if k == target_go:
                positive = float(v)
            total += float(v)

        total = int(total)
        ratio = positive / total
        if ratio < 0.05 and ratio != 0: # show minimum color for too low
            ratio = 0.05
        c1 = min_color
        c2 = color
        gradient_color = color_gradient(c1, c2, mix=ratio)
        text = f"{positive} / {total}"
        # gradientFace = RectFace(width=100,height=50,text="%.1f" % (ratio*100), color=gradient_color, 
        #         padding_x=1, padding_y=1)

        if not tooltip:
            if node.name:
                tooltip = f'<b>{node.name}</b><br>'
            else:
                tooltip = ''
            if target_go:
                tooltip += '<br>{}: {} / {} <br>'.format(target_go, positive, total)

        gradientFace = RectFace(width=width, height=height, 
                                #text=text, 
                                color=gradient_color, 
                                padding_x=padding_x, padding_y=padding_y, tooltip=tooltip)
        return gradientFace

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

        if too_small or self.seq_format == "[]":
            for start, end in self.blocks:
                if end >= sm_start and start <= sm_end:
                    bstart = max(sm_start, start)
                    bend = min(sm_end, end)
                    bx = x0 + bstart * posw
                    by, bh = get_height(bx, y)
                    box = Box(bx, by, (bend + 1 - bstart) * posw, bh)
                    
                    yield [ "pixi-block", box ]

        elif self.seq_format == "gradients":
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)
            #yield [ f'pixi-gradients', sm_box, seq ]
            yield draw_array(sm_box,[gradientscolor[x] for x in seq])

        elif self.seq_format == "categories":
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)
            #yield [ f'pixi-gradients', sm_box, seq ]
            yield draw_array(sm_box, [profilecolors[x] for x in seq])

        else:
            seq = self.get_seq(sm_start, sm_end)
            sm_x = sm_x if drawer.TYPE == 'rect' else x0
            y, h = get_height(sm_x, y)
            sm_box = Box(sm_x+sm_x0, y, posw * len(seq), h)
            if self.seq_format == 'profiles' or posw * zx < self._min_fsize:
                aa_type = "notext"
                tooltip = f'<p>{seq}</p>'
                yield draw_array(sm_box, [gradientscolor[x] for x in seq], tooltip=tooltip)
            else:
                aa_type = "text"
                yield [ f'pixi-aa_{aa_type}', sm_box, seq ]
            