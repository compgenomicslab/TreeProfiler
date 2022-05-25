from ete4 import Tree, PhyloTree
from ete4.parser.newick import NewickError
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview.renderer.faces import RectFace, TextFace, AttrFace, CircleFace, SeqMotifFace, ScaleFace
from ete4.smartview.renderer.layouts.ncbi_taxonomy_layouts import LayoutLastCommonAncestor
from ete4 import GTDBTaxa, NCBITaxa
from ete4 import random_color

import colorsys
from collections import defaultdict
import csv
import sys

#NEWICK = '/home/deng/Projects/metatree_drawer/metatreedrawer/demo/tree_novel.nw' #phylotree.nw
METADATA = './demo/novelfam_itol_taxon.txt' #emapper_annotations.tsv
LAYOUTTEMPLATE = '/home/deng/Projects/metatree_drawer/metatreedrawer/ete4layout_template.txt'

NEWICK = sys.argv[1]

# Pick two colors. Values from 0 to 1. See "hue" at
# http://en.wikipedia.org/wiki/HSL_and_HSV
def color_gradient(hue, intensity, granularity):
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
        rgb =  rgb2hex(tuple(map(lambda x: int(x*255), 
                                 colorsys.hls_to_rgb(hue, l, base_value))))
        colors.append(rgb)
        
    colors.append("#ffffff")
    return colors

def ete4_parse(newick):
    try:
        tree = PhyloTree(newick)
    except NewickError:
        try:
            tree = PhyloTree(newick, format=1)            
        except NewickError:
            tree = PhyloTree(newick, format=1, quoted_node_names=True)

    # Correct 0-dist trees
    has_dist = False
    for n in tree.traverse(): 
        if n.dist > 0: 
            has_dist = True
            break
    if not has_dist: 
        for n in tree.iter_descendants(): 
            n.dist = 1

    return tree

def parse_metadata(metadata):
    metatable = []
    tsv_file = open(metadata)
    read_tsv = csv.DictReader(tsv_file, delimiter="\t")

    for row in read_tsv:
        metatable.append(row)
    tsv_file.close()

    columns = parse_csv(metadata)
    return metatable, columns

def parse_csv(metadata):
    columns = defaultdict(list) # each value in each column is appended to a list
    with open(metadata) as f:
        reader = csv.DictReader(f, delimiter="\t") # read rows into a dictionary format
        for row in reader: # read a row as {column1: value1, column2: value2,...}
            for (k,v) in row.items(): # go over each column name and value 
                columns[k].append(v) # append the value into the appropriate list
                                    # based on column name k
    return columns

def load_metadata_to_tree(tree, metadata):
    annotations, matrix = parse_metadata(metadata)

    columns = list(matrix.keys()) # load header of columns 
    for annotation in annotations:
        gene_name = next(iter(annotation.items()))[1] #gene name must be on first column
        try:
            target_node = tree.search_nodes(name=gene_name)[0]
            for _ in range(1, len(columns)):
                
                # if columns[_] == 'seed_ortholog': # only for emapper annotations
                #     taxid, gene = annotation[columns[_]].split('.', 1)
                #     target_node.add_prop('taxid', taxid)
                #     target_node.add_prop('gene', gene)

                target_node.add_prop(columns[_], annotation[columns[_]])
        except:
            pass

    return tree, matrix

# add layouts to leaf
def get_level(node, level=1):
    if node.is_root():
        return level
    else:
        return get_level(node.up, level + 1)

def get_layout_text(prop, color, column):
    def layout_new(node):
        nstyle = NodeStyle()
        if node.is_leaf():
            # Modify the aspect of the root node
            #nstyle["fgcolor"] = "green" # yellow
            #level = get_level(node)
            #nstyle["size"] = 5
            node.set_style(nstyle)
            node.add_face(TextFace(f'{node.props.get(prop)}',
                        color=color), 
                        column=column, position='aligned')
    return layout_new

def get_layout_lca_rects(column):
    def layout_fn(node):
       
        if node.props.get('sci_name'):
            lca = node.props.get('sci_name')
            color = node.props.get('sci_name_color', 'lightgray')
            
            level = get_level(node, level=column)
            lca_face = RectFace(self.rect_width, float('inf'), 
                    color = color, 
                    text = lca,
                    fgcolor = "white",
                    padding_x = 1, padding_y = 1)
            lca_face.rotate_text = True
            node.add_face(lca_face, position='aligned', column=level)
            node.add_face(lca_face, position='aligned', column=level,
                collapsed_only=True)

    layout_fn.__name__ = 'Last common ancestor'
    layout_fn.contains_aligned_face = True
    return layout_fn

def get_layout_numeric(prop, column):
    redgradient = color_gradient(0.95, 0.6, 10)
    
    def layout_new(node):
        nstyle = NodeStyle()
        if node.props.get(prop):
            # Modify the aspect of the root node
            #nstyle["fgcolor"] = "green" # yellow
            #level = get_level(node)
            normalize_value = float(node.props.get(prop))
            color_idx = int(normalize_value*10)
            color = redgradient[10-color_idx]
            nstyle["fgcolor"] = color 
            nstyle["size"] = color_idx
            
            node.set_style(nstyle)
            # node.add_face(TextFace(f'{node.props.get(prop)}',
            #             color=color), 
            #             column=column, position='aligned')
    return layout_new

def set_layouts(matrix, aligned_faces=True):
    props = list(matrix.keys())[1:3] # exclude the first column, which is name
    layouts = []
    column = 0 
    column_colors = random_color(num=len(props))
    for prop in props:
        layout = TreeLayout(name=prop, ns=get_layout_numeric(prop, column), aligned_faces=True)
        
        # try:
        #     #print(column_colors[column])
        #     #raw_array = list(map(float, matrix[prop]))
        #     layout = TreeLayout(name=prop, ns=get_layout_numeric(prop, column), aligned_faces=True)
        # except:
        #     print
        #     #layout = TreeLayout(name=prop, ns=get_layout_text(prop, column_colors[column], column), aligned_faces=True)
        #     pass
        #layout = TreeLayout(name=prop, ns=get_layout_text(prop, column_colors[column], column), aligned_faces=True)
        
        layouts.append(layout)
        column += 1
    return layouts

def annotate_tree(newick, sp_delimiter=None, sp_field=0):
    gtdb = GTDBTaxa()
    
    def return_spcode(leaf):
        try:
            return leaf.name.split(sp_delimiter)[sp_field]
        except IndexError:
            return leaf.name

    #tree = PhyloTree(newick)
    # extract sp codes from leaf names
    tree.set_species_naming_function(return_spcode)
    
    gtdb.annotate_tree(tree, taxid_attr="name")

    # Annotate tree for smartview visualization
    set_sci_names = set()
    for n in tree.traverse():
        set_sci_names.add(n.props.get('sci_name'))

    colors = random_color(num=len(set_sci_names), l=0.5, s=0.5)
    sci_name2colors = defaultdict()

    for ele in set_sci_names:
        if ele in sci_name2colors.keys():
            continue
        else:
            sci_name2colors[ele] = colors.pop()

    all_props = set()
    name2n = {}
    for n in tree.traverse():
        name2n[n.name] = n
        sci_name = n.props.get('sci_name')
        sci_name_color = sci_name2colors[sci_name]
        n.add_prop('sci_name_color', sci_name_color)
        all_props.update(set(n.props.keys()))


    all_props.discard("_speciesFunction")
    return tree

tree = ete4_parse(NEWICK)
# gtdb= GTDBTaxa()
# gtdb.annotate_tree(tree, taxid_attr="name")
#tree = annotate_tree(tree)
#tree = annotate_tree(tree.write(properties=[]))

if METADATA:
    tree, matrix = load_metadata_to_tree(tree, METADATA)
    #print(tree.write(format=1, properties=[]))
    layouts = set_layouts(matrix)

#layouts = []
#layouts.append(LayoutLastCommonAncestor())
tree.explore(tree_name='example',layouts=layouts)