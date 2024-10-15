import sys
from ete4 import Tree

filename = sys.argv[1]

t = Tree(open(filename), parser=1)

print("Target tree internal node Root contains the following properties:  ")
print(t.props)

print(f"Target tree leaf node {next(t.leaves()).name}contains the following propertiies:  ")
print(next(t.leaves()).props)