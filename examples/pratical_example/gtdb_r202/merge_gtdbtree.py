
import sys
from ete4 import Tree

# load gtdb trees and combine
# ar122_r202.tree & bac120_r202.tree
tree = Tree(open(sys.argv[1]))
tree_bact = Tree(open(sys.argv[2]))

root_node_bac = next(n.root for n in tree_bact.traverse())
root_node_arch = next(n.root for n in tree.traverse())
root_node_arch.add_child(root_node_bac)
sys.stdout.write(root_node_arch.write(parser=0))
#root_node_arch.write(outfile='merge_gtdb.tree',format_root_node=True)
