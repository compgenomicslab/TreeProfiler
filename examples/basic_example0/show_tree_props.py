import sys
from ete4 import Tree

filename = sys.argv[1]

t = Tree(open(filename), parser=1)

print("Target tree internal node Root contains the following properties:  ")
print(t.props)

print("Target tree leaf node Taxa_0 contains the following propertiies:  ")
print(t['Taxa_0'].props)