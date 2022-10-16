#!/usr/bin/env python

from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree

from collections import defaultdict
from collections import Counter
from scipy import stats
import itertools
import math
import numpy as np
import csv
import sys


TREEFILE = sys.argv[1]
METADATAFILE = sys.argv[2]

def parse_csv(input_file):
    metadata = {}
    with open(input_file, 'r') as f:
       
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames
        
        node_header, node_props = headers[0], headers[1:]
        
        for row in reader:
            nodename = row[node_header]
            del row[node_header]
            metadata[nodename] = dict(row)
        
    #print(metadata)
    return metadata, node_props

def ete4_parse(newick):
    try:
        tree = Tree(newick)
    except NewickError:
        try:
            tree = Tree(newick, format=1)            
        except NewickError:
            tree = Tree(newick, format=1, quoted_node_names=True)

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

def load_metadata_to_tree(tree, metadata_dict):
    for node, props in metadata_dict.items():
        hits = tree.search_nodes(name=node) # including internal nodes
        if hits:
            target_node = hits[0]
            for key,value in props.items():
                target_node.add_prop(key, value)
        else:
            pass
    return tree

def merge_annotations(nodes, target_props):
    internal_props = {}

    for target_prop in target_props:
        
        try:
            prop_array = np.array(children_prop_array(nodes, target_prop),dtype=np.float64)
            n, (smin, smax), sm, sv, ss, sk = stats.describe(prop_array)

            internal_props[add_suffix(target_prop, 'sum')] = np.sum(prop_array)
            internal_props[add_suffix(target_prop, 'min')] = smin
            internal_props[add_suffix(target_prop, 'max')] = smax
            internal_props[add_suffix(target_prop, 'mean')] = sm
            if math.isnan(sv) == False:
                internal_props[add_suffix(target_prop, 'variance')] = sv
            else:
                internal_props[add_suffix(target_prop, 'variance')] = 0

        except ValueError: # for strings
            prop_list = children_prop_array(nodes, target_prop)
            internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])

    return internal_props

def add_suffix(name, suffix, delimiter='_'):
    return str(name) + delimiter + str(suffix)

def children_prop_array(nodes, prop):
    array = [n.props.get(prop) for n in nodes if n.props.get(prop)] 
    return array

def main():
    metadata_dict, node_props = parse_csv(METADATAFILE)
    tree = ete4_parse(TREEFILE)
    annotated_tree = load_metadata_to_tree(tree, metadata_dict)
    
    node2leaves = annotated_tree.get_cached_content()
    for node in annotated_tree.traverse("postorder"):
        if node.is_leaf():
            pass
        else:
            internal_props = merge_annotations(node2leaves[node], node_props)
            node.add_prop('leaves', len(node2leaves[node]))
            for key,value in internal_props.items():
                node.add_prop(key, value)
            
    # print(annotated_tree.write(properties=[]))
    # annotated_tree
    return annotated_tree


output_tree = main()
print(output_tree.write(properties=[]))
output_tree.explore(tree_name='example',layouts=[], port=5000)