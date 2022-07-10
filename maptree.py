from ete4 import Tree, PhyloTree
from ete4.parser.newick import NewickError
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace

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
        if n.name:
            pass
        else:
            n.name = n.props.get("sci_name", "")

    return tree


# load and clean tree which will clean the original internal nodes
clean_newick = ete4_parse(NEWICK).write(properties=[])
tree = PhyloTree(clean_newick)
#tree = ete4_parse(NEWICK)

# Taxonomic annotation
tree = annotate_taxa(tree)

# Metadata annotation
if METADATA:
   tree, matrix = load_metadata_to_tree(tree, METADATA)

#collapse
def collapse_cutoff(prop, cutoff):
    def layout_fn(node):
        if node.props.get(prop):
            if not node.is_root() and float(node.props.get(prop)) < float(cutoff):
                #node = node.up
                #face_name = TextFace(node.props.get('sci_name'), color="red")
                node.sm_style["draw_descendants"] = False
                node.sm_style["outline_color"] = "red"
                #node.add_face(face_name, column = 5,  position = 'aligned', collapsed_only=True)
    layout_fn.name = "collapse_" + prop
    return layout_fn
    return


from taxon_layouts import *
layouts = [
    # TreeLayout(name="collapse_cutoff", ns=collapse_cutoff('random_fraction', 0.70)),
    # TreeLayout(name='collapse_kingdom', ns=collapse_kingdom()),
    # TreeLayout(name='collapse_phylum', ns=collapse_phylum()),
    # TreeLayout(name='collapse_class', ns=collapse_class()),
    # TreeLayout(name='collapse_order', ns=collapse_order()),
    # TreeLayout(name='collapse_family', ns=collapse_family()),
    # TreeLayout(name='collapse_genus', ns=collapse_genus()),
    # TreeLayout(name='collapse_species', ns=collapse_species()),
]


tree.explore(tree_name='example',layouts=layouts, port=5000)
