from ete4 import Tree, PhyloTree
from ete4.parser.newick import NewickError
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
import sys

NEWICK = sys.argv[1]
#METADATA = 
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

def collapse_tree(tree, rank="genus"):

    return

# load and clean tree
clean_newick = ete4_parse(NEWICK).write(properties=[])
tree = PhyloTree(clean_newick)

#tree = PhyloTree(NEWICK, format=1, quoted_node_names=True)
# def return_spcode(leaf):
#     try:
        
#         return leaf.name.split(":")[1]
#     except IndexError:
#         return leaf.name

#     # extract sp codes from leaf names
# tree.set_species_naming_function(return_spcode)

# Taxonomic annotation
from ete4 import GTDBTaxa
gtdb = GTDBTaxa()
gtdb.annotate_tree(tree,  taxid_attr="name")
#tree.annotate_gtdb_taxa(taxid_attr='name')


# collapse
def processable_node(node):
    if node.props.get("rank") == "subspecies":
        return True
    else:
        return False

#print(tree.write(is_leaf_fn=processable_node, properties=[]))

# for leaf in tree.iter_leaves(is_leaf_fn=processable_node):
#     print(leaf)

t2 = ete4_parse(tree.write(is_leaf_fn=processable_node, properties=[], format=1))
t2.explore(tree_name='example',layouts=[])