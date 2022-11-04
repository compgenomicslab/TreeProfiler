#!/usr/bin/env python


from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa

from argparse import ArgumentParser
import argparse
from collections import defaultdict
from collections import Counter
from scipy import stats
import b64pickle
import itertools
import math
import numpy as np
import csv
import sys


__author__ = 'Ziqi DENG'
__license__ = "GPL v2"
__email__ = 'dengziqi1234@gmail.com'
__version__ = '0.0.1'
__date__ = '01-11-2022'
__description__ = ('A program for profiling metadata on target '
                    'tree and conduct summary analysis')

def read_args():
    """
    Parse the input parameters
    Return the parsed arguments.
    """
    parser = ArgumentParser(description=
        "metatreeprofiler.py (ver. "+__version__+
        " of "+__date__+")." + __description__+ " Authors: "+
        __author__+" ("+__email__+")",
        formatter_class=argparse.RawTextHelpFormatter)

    # input parameters group
    group = parser.add_argument_group(title='input parameters',
        description="Input parameters")
    group.add_argument('-t', '--tree',
        type=str,
        required=False,
        help="Input tree, .nw file, customized tree input")
    group.add_argument('--taxatree',
        type=str,
        required=False,
        help="<kingdom|phylum|class|order|family|genus|species|subspecies> reference tree from taxonomic database")
    group.add_argument('--taxadb',
        type=str,
        #default='GTDB',
        required=False,
        help="<NCBI|GTDB> for taxonomic profiling or fetch taxatree default [GTDB]")    
    group.add_argument('-d', '--metadata',
        type=str,
        required=True,
        help="<metadata.csv> .csv, .tsv. mandatory input")
    group.add_argument('--delimiter',
        type=str,
        default='\t',
        required=False,
        help="delimiter of metadata columns. default [tab]")
    group.add_argument('--text_column',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as categorical data")
    group.add_argument('--num_column',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as numerical data")
    group.add_argument('--taxon_column',
        type=str,
        required=False,
        help="<col1> name of columns which need to be read as taxon data")

    group = parser.add_argument_group(title='Analysis arguments',
        description="Analysis parameters")
    group.add_argument('--rank_limit',
        type=str,
        required=False,
        help="TAXONOMIC_LEVEL prune annotate tree by rank limit")
    group.add_argument('--taxonomic_profile',
        default=False,
        action='store_true',
        required=False,
        help="TAXONOMIC_LEVEL prune annotate tree by rank limit")

    group = parser.add_argument_group(title='Plot arguments',
        description="Plot parameters")
    group.add_argument('--TextLayouts',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be plot as Textlayouts")
    group.add_argument('--RectangularLayouts',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be plot as Textlayouts")
    group.add_argument('--HeatmapLayouts',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as HeatmapLayouts")
    group.add_argument('--BarplotLayouts',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as BarplotLayouts")

    group = parser.add_argument_group(title='Output arguments',
        description="Output parameters")
    group.add_argument('--interactive',
        default=True,
        action='store_false',
        help="run interactive session")
    group.add_argument('--plot',
        type=str,
        required=False,
        help="output as pdf")
    group.add_argument('-o', '--outtree',
        type=str,
        required=False,
        help="output annotated tree")
    group.add_argument('--outtsv',
        type=str,
        required=False,
        help="output annotated tsv file")


    args = parser.parse_args()
    return args

def parse_csv(input_file, delimiter='\t'):
    metadata = {}
    with open(input_file, 'r') as f:
       
        reader = csv.DictReader(f, delimiter=delimiter)
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

def merge_annotations(nodes, target_props, dtype='str'):
    internal_props = {}

    for target_prop in target_props:
        if dtype == 'str':
            prop_list = children_prop_array(nodes, target_prop)
            internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])
            
        elif dtype == 'num':
            prop_array = np.array(children_prop_array(nodes, target_prop),dtype=np.float64)
            n, (smin, smax), sm, sv, ss, sk = stats.describe(prop_array)

            internal_props[add_suffix(target_prop, 'sum')] = np.sum(prop_array)
            internal_props[add_suffix(target_prop, 'min')] = smin
            internal_props[add_suffix(target_prop, 'max')] = smax
            internal_props[add_suffix(target_prop, 'mean')] = sm
            if math.isnan(sv) == False:
                internal_props[add_suffix(target_prop, 'std')] = sv
            else:
                internal_props[add_suffix(target_prop, 'std')] = 0
        # try:
        #     prop_array = np.array(children_prop_array(nodes, target_prop),dtype=np.float64)
        #     n, (smin, smax), sm, sv, ss, sk = stats.describe(prop_array)

        #     internal_props[add_suffix(target_prop, 'sum')] = np.sum(prop_array)
        #     internal_props[add_suffix(target_prop, 'min')] = smin
        #     internal_props[add_suffix(target_prop, 'max')] = smax
        #     internal_props[add_suffix(target_prop, 'mean')] = sm
        #     if math.isnan(sv) == False:
        #         internal_props[add_suffix(target_prop, 'std')] = sv
        #     else:
        #         internal_props[add_suffix(target_prop, 'std')] = 0

        # except ValueError: # for strings
        #     prop_list = children_prop_array(nodes, target_prop)
        #     internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])

    return internal_props

def add_suffix(name, suffix, delimiter='_'):
    return str(name) + delimiter + str(suffix)

def children_prop_array(nodes, prop):
    array = [n.props.get(prop) for n in nodes if n.props.get(prop)] 
    return array

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

def taxatree_prune(tree, rank_limit='subspecies'):
    rank_limit = rank_limit.lower()
    
    ex = False
    while not ex:
        ex = True
        for n in tree.iter_leaves():
            if n.props.get('rank') != rank_limit:
                n.detach()
                ex = False

    return tree

def tree2table(tree, internal_node=True, props=[], outfile='tree2table.csv'):
    node2leaves = {}
    leaf2annotations = {}
    with open(outfile, 'w', newline='') as csvfile:
        fieldnames = props
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for node in tree.traverse():
            if internal_node:
                output_row = dict(node.props)
                for k, prop in output_row.items():
                    if type(prop) == list:
                        output_row[k] = '|'.join(str(v) for v in prop)
                writer.writerow(output_row)
            else:
                if node.is_leaf():
                    output_row = dict(node.props)
                    for k, prop in output_row.items():
                        if type(prop) == list:
                            output_row[k] = '|'.join(str(v) for v in prop)
                    writer.writerow(output_row)
                else:
                    pass

    return 

def main():
    args = read_args()
    print(args)
    # parse csv to metadata table
    if args.metadata:
        metadata_dict, node_props = parse_csv(args.metadata, delimiter=args.delimiter)

    # parse tree
    if args.tree:
        tree = ete4_parse(args.tree)
    elif args.taxa and args.taxadb:
        tree = ''

    if args.text_column:
        text_column = args.text_column.split(',')
    else:
        text_column = []

    if args.num_column:
        num_column = args.num_column.split(',')
    else:
        num_column = []

    # if COLFORMAT:
    #     col_format = ['num','num','num','str','str','str']
    #     #col_format = COLFORMAT.split(',')
    #     pos_map = defaultdict(list)
    #     for pos, ele in enumerate(col_format):
    #         pos_map[ele].append(pos)

    #     if col_format:
    #         str_column = pos_map['str']
    #         numerical_column = pos_map['num']        

    # parse csv to metadata table
    
        # if str_column:
        #     str_props =  [node_props[i] for i in str_column]
        # if numerical_column:
        #     num_props =  [node_props[i] for i in numerical_column]

    rest_column = list(set(node_props) - set(text_column) - set(num_column))
    
    
    # load annotations to leaves
    annotated_tree = load_metadata_to_tree(tree, metadata_dict)
    
    # merge annotations
    node2leaves = annotated_tree.get_cached_content()
    for node in annotated_tree.traverse("postorder"):
        internal_props = {}
        if node.is_leaf():
            pass
        else:
            if text_column:
                internal_props_text = merge_annotations(node2leaves[node], text_column, dtype='str')
                internal_props.update(internal_props_text)
            
            if num_column:
                internal_props_num = merge_annotations(node2leaves[node], num_column, dtype='num')
                internal_props.update(internal_props_num)

            if rest_column:
                internal_props_rest = merge_annotations(node2leaves[node], rest_column, dtype='str')
                internal_props.update(internal_props_rest)
            
            #internal_props = {**internal_props_text, **internal_props_num, **internal_props_rest}
            for key,value in internal_props.items():
                node.add_prop(key, value)
    
    ### Anslysis ###
    # taxa annotations
    if args.taxonomic_profile:
        if not args.taxadb:
            print('Please specify which taxa db using --taxadb <GTDB|NCBI>')
            asdasds
        else:
            annotated_tree = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name")
        # if args.taxon_column:
        #     annotated_tree = annotate_taxa(annotated_tree, taxid_attr=taxon_column)
        # else:
        #     annotated_tree = annotate_taxa(annotated_tree, taxid_attr="name")
    
    # collapse tree by rank
    if args.rank_limit:
        annotated_tree = taxatree_prune(annotated_tree, rank_limit=rank_limit)

    annotated_tree.explore(tree_name='example',layouts=[], port=5000)
    return annotated_tree

if __name__ == '__main__':
    main()

#output_tree = main()


# # # write to pickle
# with open(OUTPUTTREE+'.ete', 'w') as f:
#     f.write(b64pickle.dumps(output_tree, encoder='pickle', pack=False))

# # read from pickle
# with open(OUTPUTTREE+'.ete', 'r') as f:
#     file_content = f.read()
#     print(b64pickle.loads(file_content, encoder='pickle', unpack=False))

# write to newick tree
#output_tree.write(outfile=OUTPUTTREE, properties=[], format=1)

# write to tsv file
# fieldnames = [
#     'name', 'support', 'dist',  'leaves', 'sample1_sum', 'sample1_min', 'sample1_max', 'sample1_mean', 'sample1_variance', 'sample2_sum', 'sample2_min', 'sample2_max', 'sample2_mean', 'sample2_variance', 'sample3_sum', 'sample3_min', 'sample3_max', 'sample3_mean', 'sample3_variance', 'sample4_sum', 'sample4_min', 'sample4_max', 'sample4_mean', 'sample4_variance', 'sample5_sum', 'sample5_min', 'sample5_max', 'sample5_mean', 'sample5_variance', 'random_type_counter', 'taxid', 'sci_name', 'common_name', 'lineage', 'rank', 'named_lineage', 'nleaves', 'sample1', 'sample2', 'sample3', 'sample4','sample5','random_type' 
# ]
# tree2table(output_tree, internal_node=True, props=fieldnames)


# interactive explore
#output_tree.explore(tree_name='example',layouts=[], port=5000)