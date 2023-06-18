
import sys
from ete4 import Tree

# load gtdb trees and combine
# ar122_r202.tree & bac120_r202.tree
tree = Tree(sys.argv[1], format = 1,quoted_node_names = True)
tree_bact = Tree(sys.argv[2], format=1,quoted_node_names = True)

root_node_bac = next(n.get_tree_root() for n in tree_bact.traverse())
root_node_arch = next(n.get_tree_root() for n in tree.traverse())
root_node_arch.add_child(root_node_bac)
sys.stdout.write(root_node_arch.write())
#root_node_arch.write(outfile='merge_gtdb.tree',format_root_node=True)
