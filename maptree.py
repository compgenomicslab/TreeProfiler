from ete4 import Tree, PhyloTree
from ete4.parser.newick import NewickError
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
#from ete4.smartview.renderer.layouts.staple_layouts import LayoutBarplot
from layouts import staple_layouts, taxon_layouts
from layouts.taxon_layouts import *
from collections import defaultdict
import csv
import sys

NEWICK = sys.argv[1]
METADATA = sys.argv[2]
#METADATA = ''

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
    count = 0
    for annotation in annotations:
        gene_name = next(iter(annotation.items()))[1] #gene name must be on first column
        #GB_GCA_011332645.1
        
        # has to be leaf
        #target_nodes = list(filter(lambda n: n.name==gene_name and n.is_leaf(), tree.traverse()))
        
        # all the nodes including internal node
        target_nodes = list(filter(lambda n: n.name==gene_name, tree.traverse()))


        if target_nodes:
            target_node = target_nodes[0]
            count += 1
            for _ in range(1, len(columns)):
            # if columns[_] == 'seed_ortholog': # only for emapper annotations
            #     taxid, gene = annotation[columns[_]].split('.', 1)
            #     target_node.add_prop('taxid', taxid)
            #     target_node.add_prop('gene', gene)
                target_node.add_prop(columns[_], annotation[columns[_]])
        else:
            pass
        #target_node = tree.search_nodes(name=gene_name)[0] # funny 
    print(count)
       

    return tree, matrix

from ete4 import GTDBTaxa, NCBITaxa
def annotate_taxa(tree, db="GTDB", taxid_attr="name", sp_delimiter='.', sp_field=0):
    def return_spcode(leaf):
        try:
            return leaf.name.split(sp_delimiter)[sp_field]
        except IndexError:
            return leaf.name

    if db == "GTDB":
        gtdb = GTDBTaxa()
        gtdb.annotate_tree(tree,  taxid_attr=taxid_attr)
    elif db == "NCBI":
        ncbi = NCBITaxa()
        # extract sp codes from leaf names
        tree.set_species_naming_function(return_spcode)
        ncbi.annotate_tree(tree, taxid_attr="species")

    # tree.annotate_gtdb_taxa(taxid_attr='name')
    # assign internal node as sci_name
    for n in tree.traverse():
        if not n.is_leaf():
            nleaves = str(len(n))
            n.add_prop("nleaves", nleaves)
        if n.name:
            pass
        else:
            n.name = n.props.get("sci_name", "")

    return tree

from collections import Counter,defaultdict
import numpy as np

def children_prop_array(node, prop):
    array = [n.props.get(prop) for n in node.iter_leaves() if n.props.get(prop)] 
    return array

def get_stats(array):
    if array:
        np_array = np.float_(array)
        return np_array.sum(), np_array.mean()
    else:
        return None
    #return np_array.sum(), np_array.mean(), np_array.max(), np_array.min()

#load prop
def get_calculation(tree, props):

    for node in tree.traverse():
        for prop in props:
            try:
                node_sum, node_mean = get_stats(children_prop_array(node, prop))
                if node.is_leaf():
                    pass
                else:
                    node.add_prop(prop+"_sum", node_sum)
                    node.add_prop(prop+"_mean", node_mean)
            except Exception as e: 
                pass
    return tree

# get config
def get_config(props, num_col=[]):
    config = {}
    #props = props[1:] # ignore first column which is leaf name
    for idx, prop in enumerate(props):
        if idx not in num_col and idx != 0:
            config[prop] = 0 #TextLayout
        elif idx in num_col:
            #config[prop] = 3 #barplot
            config[prop] = 4 #heatmap
    return config
#########################################run#############################################
# load and clean tree which will clean the original internal nodes

clean_newick = ete4_parse(NEWICK).write(properties=[])
tree = PhyloTree(clean_newick)
#tree = ete4_parse(NEWICK)

# Metadata annotation
if METADATA:
    tree, matrix = load_metadata_to_tree(tree, METADATA)
    props = list(matrix.keys())
    layouts_config = get_config(props, [1,2,3,4,5])

# Taxonomic annotation
if "Taxon" in layouts_config.keys():
    tree = annotate_taxa(tree, taxid_attr="Taxon")
else:
    tree = annotate_taxa(tree, taxid_attr="name")
#########################################calculation#############################################
# abundance calculation sum
# prop = props[1]
# tree = get_calculation(tree, prop)

# relative abundance for heatmap
heatmap_prop = props[1:]
tree = get_calculation(tree, heatmap_prop)

#################################################layout##########################################
# demo collapse
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
import colorsys

def new_collapse_class(prop, level):
    def layout_fn(node):
        if not node.is_root() and  node.props.get('rank') == 'species':
            
            node.sm_style["draw_descendants"] = False
            node.sm_style["bgcolor"] ="#FFFFFF"
            node.sm_style["outline_color"] = "deepyellow"
            #node.sm_style["shape"] = "square"

            #face_name = TextFace(node.props.get('sci_name')+"||"+str(node.props.get('random_abundance_sum')), color="red")
            if not node.is_leaf():
                nleaves = str(len(node))
                face_name = TextFace(node.props.get('sci_name')+"("+nleaves+")", color="red")            
            else:
                face_name = TextFace(node.props.get('sci_name'), color="red")
            node.add_face(face_name, column = 8,  position = 'aligned', collapsed_only=True)

            # heatmap
            relative_abundance = float(node.props.get(prop))
            color_idx = int(relative_abundance*10)
            color = redgradient[color_idx]

            identF = RectFace(width=50,height=50,text="%.1f" % (relative_abundance*100), color=color, 
            padding_x=1, padding_y=1)
            #face_name = TextFace(node.props.get('name'), color="red")
            #face_name = TextFace("%.1f" % (relative_abundance*100), color=color)
            node.add_face(identF, column = level,  position = 'branch_right')
            node.add_face(identF, column = level,  position = 'branch_right', collapsed_only=True)

    layout_fn.name = "level5_species"
    return layout_fn
    return

def new_collapse_heatmap(prop, level):
    def layout_fn(node):
        if not node.is_root() and node.props.get(prop):

            # heatmap
            relative_abundance = float(node.props.get(prop))
            color_idx = int(relative_abundance*10)
            color = redgradient[color_idx]

            identF = RectFace(width=50,height=50,text="%.1f" % (relative_abundance*100), color=color, 
            padding_x=1, padding_y=1)
            #face_name = TextFace(node.props.get('name'), color="red")
            #face_name = TextFace("%.1f" % (relative_abundance*100), color=color)
            node.add_face(identF, column = level,  position = 'branch_right')
            node.add_face(identF, column = level,  position = 'branch_right', collapsed_only=True)

    layout_fn.name = "level5_species"
    return layout_fn
    return

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
        rgb =  rgb2hex(tuple(map(lambda x: int(x*255), colorsys.hls_to_rgb(hue, l, base_value))))
        colors.append(rgb)
        
    colors.append("#ffffff")
    return list(reversed(colors))

def heatmap_layout(prop, level):
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            relative_abundance = float(node.props.get(prop))
            color_idx = int(relative_abundance*10)
            color = redgradient[color_idx]

            identF = RectFace(width=50,height=50,text="%.1f" % (relative_abundance*100), color=color, 
            padding_x=1, padding_y=1)
            #face_name = TextFace(node.props.get('name'), color="red")
            #face_name = TextFace("%.1f" % (relative_abundance*100), color=color)
            node.add_face(identF, column = level,  position = 'branch_right')
            node.add_face(identF, column = level,  position = 'branch_right', collapsed_only=True)
    return layout_fn
    return

def text_layout(prop, level, color='blue'):
    def layout_fn(node):
        if node.is_leaf() and node.props.get(prop):
            prop_text = node.props.get(prop)
            prop_face = TextFace(prop_text, color=color)
            node.add_face(prop_face, column = level, position = "branch_right")
            node.sm_style["bgcolor"] = 'black' # highligh clade
            # while (node):
            #         node = node.up
            #         if node:
            #             node.sm_style["hz_line_width"] = 5
    return layout_fn
    return

redgradient = color_gradient(0.95, 0.6, 10)



#print(tree.search_nodes(name="s__EX4484-205 sp002255045"))
# for n in tree.iter_leaves():
#     print(n.props.get("sample1"))

def get_layouts(config):
    layouts = []
    level = 10
    mapping = {
        0: "TextLayout",
        1: "NodeStyleLayout",
        2: "TaxonLayout",
        3: "NumericLayout",
        4: "HeatmapLayout",
    }

    for prop, code in config.items():
        layout_name = mapping[code]
        if layout_name == "TextLayout":
            layout = TreeLayout(name=prop, ns=text_layout(prop, level))
            level += 1
        # elif layout_name == "TaxonLayout":
        #     pass
        elif layout_name == "NumericLayout":
            layout = LayoutBarplot(name=prop, size_prop=prop)
            
        elif layout_name == "HeatmapLayout":
            layout = TreeLayout(name=prop, ns=heatmap_layout(prop, level))
            level += 1
        layouts.append(layout)

    taxon_layouts = [
        TreeLayout(name='level1_kingdom', ns=collapse_kingdom()),
        TreeLayout(name='level2_phylum', ns=collapse_phylum()),
        TreeLayout(name='level3_class', ns=collapse_class()),
        TreeLayout(name='level4_order', ns=collapse_order()),
        TreeLayout(name='level5_family', ns=collapse_family()),
        TreeLayout(name='level6_genus', ns=collapse_genus()),
        TreeLayout(name='level7_species', ns=collapse_species()),
    ]
    layouts = layouts + taxon_layouts
    return layouts

layouts = get_layouts(layouts_config)

layouts = [
#     # TreeLayout(name="collapse_cutoff", ns=collapse_cutoff('random_fraction', 0.70)),
#     #TreeLayout(name='level3_class', ns=new_collapse_class()),
#     #TreeLayout(name="taxa_face", ns=taxa_rect_layout()),
    
    
#     # collapse_heatmap
    # TreeLayout(name="collapse_heatmap_sample1", ns=new_collapse_heatmap("sample1_mean", 5)),
    # TreeLayout(name="collapse_heatmap_sample2", ns=new_collapse_heatmap("sample2_mean", 6)),
    
#     ### Default Text Layout
#     # TextAlign layout
#     TreeLayout(name='abundance', ns=text_layout('abundance', 10)),
#     # NodeStyle layout
    

#     ### Bar plot of count
#     # metadata_p__Thermoproteota.txt
#     # LayoutBarplot(name="abundance", size_prop="abundance"),
    
#     ### taxa default ####
#     # metadata_p__Thermoproteota*
#     # TreeLayout(name='level1_kingdom', ns=collapse_kingdom()),
#     # TreeLayout(name='level2_phylum', ns=collapse_phylum()),
#     # TreeLayout(name='level3_class', ns=collapse_class()),
#     # TreeLayout(name='level4_order', ns=collapse_order()),
#     # TreeLayout(name='level5_family', ns=collapse_family()),
#     # TreeLayout(name='level6_genus', ns=collapse_genus()),
#     # TreeLayout(name='level7_species', ns=collapse_species()),

#     ### relative abundance ####
#     # metadata_p__Thermoproteota_relative.txt
    # TreeLayout(name="sample1",ns=heatmap_layout("sample1", 5)), 
    # TreeLayout(name="sample2",ns=heatmap_layout("sample2", 6)), 
    # TreeLayout(name="sample3",ns=heatmap_layout("sample3", 7)),
    # TreeLayout(name="sample4",ns=heatmap_layout("sample4", 8)),
    # TreeLayout(name="sample5",ns=heatmap_layout("sample5", 9)), 
    

]


tree.explore(tree_name='example',layouts=layouts, port=5000)
