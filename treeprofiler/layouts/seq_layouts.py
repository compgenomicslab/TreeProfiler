
from ete4 import SeqGroup
from ete4.smartview import TreeLayout
from ete4.smartview import AlignmentFace, SeqMotifFace, Face
from ete4.smartview.renderer import draw_helpers
from ete4.smartview.renderer.draw_helpers import draw_line, draw_text, Box

#from utils import get_consensus_seq
from pathlib import Path
from io import StringIO
import json


__all__ = [ "LayoutAlignment" ]
DOMAIN2COLOR = 'pfam2color.json' # smart2color.json

def get_colormap():
    with open(Path(__file__).parent / DOMAIN2COLOR) as handle:
        _pfam2color = json.load(handle)
    return _pfam2color

class LayoutAlignment(TreeLayout):
    def __init__(self, name="Alignment",
            alignment=None, alignment_prop=None, format='seq', width=700, height=15,
            column=0, scale_range=None, window=[], summarize_inner_nodes=True, 
            aligned_faces=True):
        super().__init__(name, aligned_faces=aligned_faces)
        #self.alignment = SeqGroup(alignment) if alignment else None
        self.alignment_prop = alignment_prop
        self.width = width
        self.height = height
        self.column = column
        self.aligned_faces = True
        self.format = format

        #self.length = len(next(self.alignment.iter_entries())[1]) if self.alignment else None
        self.scale_range = (0, scale_range) or (0, self.length)
        self.window = window
        self.summarize_inner_nodes = summarize_inner_nodes

    def set_tree_style(self, tree, tree_style):
        if self.scale_range:
            if self.window:
                face = ScaleFace(width=self.width, scale_range=self.window, padding_y=10)
                tree_style.aligned_panel_header.add_face(face, column=self.column)
            else:
                face = ScaleFace(width=self.width, scale_range=self.scale_range, padding_y=10)
                tree_style.aligned_panel_header.add_face(face, column=self.column)
    
    def get_seq(self, node, window=[]):
        seq = node.props.get(self.alignment_prop, None)
        return seq

    def set_node_style(self, node):
        seq = self.get_seq(node)
        if seq:
            seq = str(seq) # convert Bio.seq.seq to string seq
            if self.window:
                start, end = self.window
                seq = seq[start:end]
            seqFace = AlignmentFace(seq, seq_format=self.format, bgcolor='grey',
                    width=self.width, height=self.height)
            node.add_face(seqFace, column=self.column, position='aligned',
                    collapsed_only=(not node.is_leaf)) 


def get_alnface(seq_prop, level):
    def layout_fn(node):
        if node.is_leaf:
            seq = node.props.get(seq_prop)
            seq_face = AlignmentFace(seq, seqtype='aa',
            gap_format='line', seq_format='[]',
            width=None, height=None, # max height
            fgcolor='black', bgcolor='#bcc3d0', gapcolor='gray',
            gap_linewidth=0.2,
            max_fsize=12, ftype='sans-serif',
            padding_x=0, padding_y=0)
            node.add_face(seq_face, position="aligned", column=level)
    return layout_fn

class LayoutDomain(TreeLayout):
    def __init__(self, prop, name,
            column=10, colormap={},
            min_fsize=4, max_fsize=15,
            padding_x=5, padding_y=0):
        super().__init__(name or "Domains layout")
        self.prop = prop
        self.column = column
        self.aligned_faces = True
        self.min_fsize = min_fsize
        self.max_fsize = max_fsize
        self.padding = draw_helpers.Padding(padding_x, padding_y)
        if not colormap:
            self.colormap = get_colormap()
        else:
            self.colormap = colormap

    def get_doms(self, node):
        pair_delimiter = "@"
        item_seperator = "||"
        dom_list = []
        if node:
            dom_prop = node.props.get(self.prop, [])
            if dom_prop:
                
                dom_list = [dom.split(pair_delimiter) for dom in dom_prop.split(item_seperator)]
                #print(dom_list)
                return dom_list
            #return node.props.get(self.prop, [])
        # else:
        #     first_node = next(node.leaves())
        #     print("here", first_node.props.get(self.prop, []))
        #     return first_node.props.get(self.prop, [])

    def parse_doms(self, dom_list):
        doms = []
        for name, start, end in dom_list:
            color = self.colormap.get(name, "lightgray")
            dom = [int(start), int(end), "()", 
                   None, None, color, color,
                   "arial|30|black|%s" %(name)]
            doms.append(dom)
        return doms

    def set_node_style(self, node):
        dom_list = self.get_doms(node)
        if dom_list:
            doms = self.parse_doms(dom_list)
            fake_seq = '-' * int(node.props.get("len_alg", 0))
            
            if doms or fake_seq:
                seqFace = SeqMotifFace(seq=fake_seq, motifs=doms, width=400,
                        height=30)
                node.add_face(seqFace, column=self.column, 
                        position="aligned",
                        collapsed_only=(not node.is_leaf))
        else:
            print("no domain found for node %s" % node.name)

class ScaleFace(Face):
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
        if drawer.TYPE == 'circ':
            p1 = cartesian(p1)
            p2 = cartesian(p2)
        yield draw_line(p1, p2, style={'stroke-width': self.line_width,
                                       'stroke': self.color})


        nticks = round((self.width * zx) / self.tick_width)
        dx = self.width / nticks
        range_factor = (self.range[1] - self.range[0]) / self.width

        if self.viewport:
            sm_start = round(max(self.viewport[0] - self.viewport_margin - x0, 0) / dx)
            sm_end = nticks - round(max(x0 + self.width - (self.viewport[1] +
                self.viewport_margin), 0) / dx)
        else:
            sm_start, sm_end = 0, nticks

        for i in range(sm_start, sm_end + 1):
            x = x0 + i * dx
            number = range_factor * i * dx + self.range[0]
            
            if number == 0:
                text = "0"
            else:
                text = self.formatter % number if self.formatter else str(number)

            text = text.rstrip('0').rstrip('.') if '.' in text else text

            self.compute_fsize(self.tick_width / len(text), dy, zx, zy)
            text_style = {
                'max_fsize': self._fsize,
                'text_anchor': 'middle',
                'ftype': f'{self.ftype}, sans-serif', # default sans-serif
                }
            text_box = Box(x,
                    y,
                    # y + (dy - self._fsize / (zy * r)) / 2,
                    dx, dy)

            yield draw_text(text_box, text, style=text_style)

            p1 = (x, y + dy - self.vt_line_height / zy)
            p2 = (x, y + dy)

            yield draw_line(p1, p2, style={'stroke-width': self.line_width,
                                           'stroke': self.color})