from ete4 import Tree, PhyloTree
from ete4.parser.newick import NewickError

NEWICK = '/home/deng/Projects/gtdb_mapper/novelfam/tree_novel.nw'
METADATA = '/home/deng/Projects/gtdb_mapper/novelfam/novelfam_itol_taxon.txt'
LAYOUTTEMPLATE = '/home/deng/Projects/gtdb_mapper/novelfam/ete4layout_template.txt'

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

def get_layout_metadata(template):
    return layout_info

tree = ete4_parse(NEWICK)