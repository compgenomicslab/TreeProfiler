import sys
from ete4 import Tree

filename = sys.argv[1]

t = Tree(open(filename), parser=1)

print("Target tree internal node if Taxa_0 and Taxa_1 contains the following properties:  ")
print(t.common_ancestor(['Taxa_0', 'Taxa_1']).props)
