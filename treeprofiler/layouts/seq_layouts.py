
from ete4 import SeqGroup
from ete4.smartview import TreeLayout
from ete4.smartview import AlignmentFace, SeqMotifFace, ScaleFace
from ete4.smartview.renderer import draw_helpers
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
            column=0, scale_range=None, summarize_inner_nodes=True, aligned_faces=True):
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
        self.summarize_inner_nodes = summarize_inner_nodes

    def set_tree_style(self, tree, tree_style):
        if self.scale_range:
            face = ScaleFace(width=self.width, scale_range=self.scale_range, padding_y=10)
            tree_style.aligned_panel_header.add_face(face, column=self.column)
    
    def _get_seq(self, node):
        return node.props.get(self.alignment_prop, None)
        # if self.alignment:
        #     return self.alignment.get_seq(node.name)
        # else:
        #     return node.props.get(alignment_prop, None)

    def get_seq(self, node):
        if node.is_leaf:
            return self._get_seq(node)

        if self.summarize_inner_nodes:
            return self._get_seq(node)
        else:
            first_leaf = next(node.leaves())
            return self._get_seq(first_leaf)
    
    def set_node_style(self, node):
        seq = self.get_seq(node)

        if seq:
            seq = str(seq) # convert Bio.seq.seq to string seq
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