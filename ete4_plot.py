#!/usr/bin/env python
from ete4 import Tree, PhyloTree
from ete4.parser.newick import NewickError
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
#from ete4.smartview.renderer.layouts.staple_layouts import LayoutBarplot
from layouts import staple_layouts, taxon_layouts
from collections import defaultdict
import csv
import sys

TREEFILE = sys.argv[1]
#METADATAFILE = sys.argv[2]
#OUTFILE = sys.argv[3]



t = Tree(TREEFILE, format=1)
t.explore(tree_name=None, layouts=[], show_leaf_name=True, 
            show_branch_length=True, show_branch_support=True, port=5000,
            custom_api={}, custom_route={})    