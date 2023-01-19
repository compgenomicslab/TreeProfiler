#!/usr/bin/env python


from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from layouts import text_layouts, taxon_layouts, staple_layouts, conditional_layouts
from plot import plot

from argparse import ArgumentParser
import argparse
from collections import defaultdict
from collections import Counter
from scipy import stats
import colorsys
import random
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


#colors_50 = ["#E41A1C","#C72A35","#AB3A4E","#8F4A68","#735B81","#566B9B","#3A7BB4","#3A85A8","#3D8D96","#419584","#449D72","#48A460","#4CAD4E","#56A354","#629363","#6E8371","#7A7380","#87638F","#93539D","#A25392","#B35A77","#C4625D","#D46A42","#E57227","#F67A0D","#FF8904","#FF9E0C","#FFB314","#FFC81D","#FFDD25","#FFF12D","#F9F432","#EBD930","#DCBD2E","#CDA12C","#BF862B","#B06A29","#A9572E","#B65E46","#C3655F","#D06C78","#DE7390","#EB7AA9","#F581BE","#E585B8","#D689B1","#C78DAB","#B791A5","#A8959F","#999999"]
paried_color = ["red", "darkblue", "lightgreen", "darkyellow", "violet", "mediumturquoise", "sienna", "lightCoral", "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal"]

def read_args():
    """
    Parse the input parameters
    Return the parsed arguments.
    """
    parser = ArgumentParser(description=
        "treeprofiler.py (ver. "+__version__+
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
    group.add_argument('-d', '--metadata',
        type=str,
        required=False,
        help="<metadata.csv> .csv, .tsv. mandatory input")
    group.add_argument('--annotated_tree',
        default=False,
        action='store_true',
        required=False,
        help="inputtree already annotated by treeprofileer")
    group.add_argument('--tree_type',
        type=str,
        default='newick',
        required=False,
        help="statistic calculation to perform for numerical data in internal nodes, [newick, ete]")
    group.add_argument('--no_colnames',
        default=False,
        action='store_true',
        required=False,
        help="metadata table doesn't contain columns name")
    group.add_argument('--text_prop',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as categorical data")
    group.add_argument('--num_prop',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as numerical data")
    group.add_argument('--bool_prop',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as boolean data")
    group.add_argument('--text_prop_idx',
        type=str,
        required=False,
        help="1,2,3 or [1-5] index of columns which need to be read as categorical data")
    group.add_argument('--num_prop_idx',
        type=str,
        required=False,
        help="1,2,3 or [1-5] index columns which need to be read as numerical data")
    group.add_argument('--bool_prop_idx',
        type=str,
        required=False,
        help="1,2,3 or [1-5] index columns which need to be read as boolean data")
    group.add_argument('--taxatree',
        type=str,
        required=False,
        help="<kingdom|phylum|class|order|family|genus|species|subspecies> reference tree from taxonomic database")
    group.add_argument('--taxadb',
        type=str,
        default='GTDB',
        required=False,
        help="<NCBI|GTDB> for taxonomic profiling or fetch taxatree default [GTDB]")    
    group.add_argument('--taxon_column',
        type=str,
        required=False,
        help="<col1> name of columns which need to be read as taxon data")
    group.add_argument('--taxon_delimiter',
        type=str,
        default=';',
        required=False,
        help="delimiter of taxa columns. default [;]")
    group.add_argument('--taxa_field',
        type=int,
        default=0,
        required=False,
        help="field of taxa name after delimiter. default 0")
    group.add_argument('--taxonomic_profile',
        default=False,
        action='store_true',
        required=False,
        help="Determine if you need taxonomic profile on tree")

    group = parser.add_argument_group(title='Analysis arguments',
        description="Analysis parameters")
    group.add_argument('--num_stat',
        default='all',
        type=str,
        required=False,
        help="statistic calculation to perform for numerical data in internal nodes, [all, sum, avg, max, min, std] ")  
    group.add_argument('--internal_plot_measure',
        default='avg',
        type=str,
        required=False,
        help="statistic measures to be shown in numerical layout for internal nodes, [default: avg]")  

    group.add_argument('--counter_stat',
        default='raw',
        type=str,
        required=False,
        help="statistic calculation to perform for categorical data in internal nodes, raw count or in percentage [raw, relative] ")  
    
    group.add_argument('--rank_limit',
        type=str,
        required=False,
        help="TAXONOMIC_LEVEL prune annotate tree by rank limit")
    group.add_argument('--pruned_by', 
        type=str,
        required=False,
        action='append',
        help='target tree pruned by customized conditions')
    group.add_argument('--collapsed_by', 
        type=str,
        required=False,
        action='append',
        help='target tree collapsed by customized conditions')
    group.add_argument('--highlighted_by', 
        type=str,
        required=False,
        action='append',
        help='target tree highlighted by customized conditions')
    
    group = parser.add_argument_group(title='basic treelayout arguments',
        description="treelayout parameters")
    group.add_argument('--drawer',
        type=str,
        required=False,
        help="Circular or Rectangular")
    group.add_argument('--collapse_level',
        type=str,
        required=False,
        help="default collapse level, default is 10") 
    group.add_argument('--ultrametric',
        default=False,
        action='store_true',
        required=False,
        help="ultrametric tree")

    group = parser.add_argument_group(title='Plot arguments',
        description="Plot parameters")
    group.add_argument('--BinaryLayout',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as BinaryLayout")
    group.add_argument('--RevBinaryLayout',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as RevBinaryLayout")

    group.add_argument('--ColorbranchLayout',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as Textlayouts")
    group.add_argument('--LabelLayout',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as LabelLayout")
    group.add_argument('--RectangularLayout',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as RectangularLayout")
    
    
    group.add_argument('--HeatmapLayout',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as HeatmapLayout")
    group.add_argument('--BarplotLayout',
        type=str,
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as BarplotLayouts")
    
    group.add_argument('--TaxonLayout',
        default=False,
        action='store_true',
        help="activate TaxonLayout")

    group = parser.add_argument_group(title='Output arguments',
        description="Output parameters")
    group.add_argument('--interactive',
        default=False,
        action='store_true',
        help="run interactive session")
    group.add_argument('--port',
        type=str,
        default=5000,
        help="run interactive session on custom port")
    group.add_argument('--plot',
        type=str,
        required=False,
        help="output as pdf")
    group.add_argument('--ete4out',
        default=False,
        action='store_true',
        help="export intermediate tree in ete4")
    group.add_argument('-o', '--outtree',
        type=str,
        required=False,
        help="output annotated tree")
    group.add_argument('--outtsv',
        type=str,
        required=False,
        help="output annotated tsv file")
    group.add_argument('--out_color_dict',
        action="store_true", 
        required=False,
        help="print color dictionary of each property")

    args = parser.parse_args()
    return args

# 
def parse_csv(input_file, delimiter='\t', no_colnames=False):
    """
    Takes tsv table as input
    Return 
    metadata, as dictionary of dictionaries for each node's metadata
    node_props, a list of property names(column names of metadata table)
    columns, dictionary of property name and it's values
    """
    metadata = {}
    columns = defaultdict(list)
    with open(input_file, 'r') as f:
        if no_colnames:
            fields_len = len(next(f).split(delimiter))
            headers = ['col'+str(i) for i in range(fields_len)]
            reader = csv.DictReader(f, delimiter=delimiter, fieldnames=headers)
        else:
            reader = csv.DictReader(f, delimiter=delimiter)
            headers = reader.fieldnames
        node_header, node_props = headers[0], headers[1:]
        
        for row in reader:
            nodename = row[node_header]
            del row[node_header]
            row = {k: 'NaN' if not v else v for k, v in row.items() } ## replace empty to NaN
            metadata[nodename] = dict(row)
            for (k,v) in row.items(): # go over each column name and value 
                columns[k].append(v) # append the value into the appropriate list
                                    # based on column name k

    return metadata, node_props, columns


def ete4_parse(newick):
    #parse tree via ete4
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

def load_metadata_to_tree(tree, metadata_dict, prop2type={}, taxon_column=None, taxon_delimiter=';'):
    name2leaf = {}
    # preload all leaves to save time instead of search in tree
    for leaf in tree.iter_leaves():
        name2leaf[leaf.name] = leaf

    # load all metadata to leaf nodes
    for node, props in metadata_dict.items():

        #hits = tree.get_leaves_by_name(node)
        #hits = tree.search_nodes(name=node) # including internal nodes
        if node in name2leaf.keys():
            target_node = name2leaf[node]
            for key,value in props.items():
                if key == taxon_column:
                    taxon_prop = value.split(taxon_delimiter)[-1]
                    target_node.add_prop(key, taxon_prop)
                elif key in prop2type and prop2type[key]=='num':
                    if math.isnan(float(value)):
                        target_node.add_prop(key, value)
                    else:
                        target_node.add_prop(key, float(value))  
                else:
                    target_node.add_prop(key, value)
        else:
            pass
        
        # if hits:
            
        #     target_node = hits[0]
        #     for key,value in props.items():
        #         if key == taxon_column:
        #             taxon_prop = value.split(taxon_delimiter)[-1]
        #             target_node.add_prop(key, taxon_prop)
        #         else:
        #             target_node.add_prop(key, value)
        # elif len(hits)>1:
        #     print('repeat')
        #     break
        # else:
        #     pass
    return tree

# def merge_annotations(nodes, target_props, dtype='str'):
#     internal_props = {}

#     for target_prop in target_props:
#         if dtype == 'str':
#             prop_list = children_prop_array(nodes, target_prop)
#             internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])
            
#         elif dtype == 'num':
#             prop_array = np.array(children_prop_array(nodes, target_prop),dtype=np.float64)
#             n, (smin, smax), sm, sv, ss, sk = stats.describe(prop_array)

#             internal_props[add_suffix(target_prop, 'sum')] = np.sum(prop_array)
#             internal_props[add_suffix(target_prop, 'min')] = smin
#             internal_props[add_suffix(target_prop, 'max')] = smax
#             internal_props[add_suffix(target_prop, 'avg')] = sm
#             if math.isnan(sv) == False:
#                 internal_props[add_suffix(target_prop, 'std')] = sv
#             else:
#                 internal_props[add_suffix(target_prop, 'std')] = 0

#     return internal_props

def merge_text_annotations(nodes, target_props, counter_stat='raw'):
    internal_props = {}
    for target_prop in target_props:
        if counter_stat == 'raw':
            prop_list = children_prop_array(nodes, target_prop)
            internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(str(key), value, '--') for key, value in sorted(dict(Counter(prop_list)).items())])

        elif counter_stat == 'relative':
            prop_list = children_prop_array(nodes, target_prop)
            counter_line = []
            #print(dict(Counter(prop_list)))
            total = sum(dict(Counter(prop_list)).values())
            #print(total)
            for key, value in sorted(dict(Counter(prop_list)).items()):
                #print(key, value)
                rel_val = '{0:.2f}'.format(float(value)/total)
                counter_line.append(add_suffix(key, rel_val, '--'))
            internal_props[add_suffix(target_prop, 'counter')] = '||'.join(counter_line)
            #internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])

        else:
            print('Invalid stat method')
            break
    
    return internal_props

# def merge_bool_annotations(nodes, target_props, counter_stat='raw'):
#     internal_props = {}
#     for target_prop in target_props:
#         if counter_stat == 'raw':
#             prop_list = children_prop_array(nodes, target_prop)
#             counter_line = []
#             for key, value in dict(Counter(prop_list)).items():
#                 counter_line.append(add_suffix(key, value, '--'))
            
#             internal_props[add_suffix(target_prop, 'counter')] = '||'.join(counter_line)
#             # internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(str(key), value, '--') for key, value in dict(Counter(prop_list)).items()])

#         # elif counter_stat == 'relative':
#         #     prop_list = children_prop_array(nodes, target_prop)
#         #     internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])

#         else:
#             print('Invalid stat method')
#             break
    
#     return internal_props

def merge_num_annotations(nodes, target_props, num_stat='all'):
    internal_props = {}
    for target_prop in target_props:
        prop_array = np.array(children_prop_array(nodes, target_prop),dtype=np.float64)
        prop_array = prop_array[~np.isnan(prop_array)] # remove nan data
        if prop_array.any():
            n, (smin, smax), sm, sv, ss, sk = stats.describe(prop_array)

            if num_stat == 'all':
                internal_props[add_suffix(target_prop, 'avg')] = sm
                internal_props[add_suffix(target_prop, 'sum')] = np.sum(prop_array)
                internal_props[add_suffix(target_prop, 'max')] = smax
                internal_props[add_suffix(target_prop, 'min')] = smin
                if math.isnan(sv) == False:
                    internal_props[add_suffix(target_prop, 'std')] = sv
                else:
                    internal_props[add_suffix(target_prop, 'std')] = 0
            
            elif num_stat == 'avg':
                internal_props[add_suffix(target_prop, 'avg')] = sm
            elif num_stat == 'sum':
                internal_props[add_suffix(target_prop, 'sum')] = np.sum(prop_array)
            elif num_stat == 'max':
                internal_props[add_suffix(target_prop, 'max')] = smax
            elif num_stat == 'min':
                internal_props[add_suffix(target_prop, 'min')] = smin
            elif num_stat == 'std':
                if math.isnan(sv) == False:
                    internal_props[add_suffix(target_prop, 'std')] = sv
                else:
                    internal_props[add_suffix(target_prop, 'std')] = 0
            else:
                print('Invalid stat method')
                break
        else:
            break
    #print(internal_props)
    return internal_props

def add_suffix(name, suffix, delimiter='_'):
    return str(name) + delimiter + str(suffix)

def children_prop_array(nodes, prop):
    array = [n.props.get(prop) for n in nodes if n.props.get(prop)] 
    return array

def annotate_taxa(tree, db="GTDB", taxid_attr="name", sp_delimiter='.', sp_field=0):
    global rank2values
    def return_spcode(leaf):
        #print(leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field])
        try:
            return leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field]
        except IndexError:
            return leaf.props.get(taxid_attr)

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
    rank2values = defaultdict(list)
    for n in tree.traverse():
        if db == 'NCBI':
            n.del_prop('_speciesFunction')
        if n.props.get('rank') and n.props.get('rank') != 'Unknown':
            rank2values[n.props.get('rank')].append(n.props.get('sci_name',''))
        
        if n.name:
            pass
        else:
            n.name = n.props.get("sci_name", "")
    
    return tree, rank2values

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

from utils import to_code, call, counter_call
def conditional_prune(tree, conditions_input, prop2type):
    conditional_output = []
    for line in conditions_input:
        single_one = to_code(line)
        conditional_output.append(single_one)

    ex = False
    while not ex:
        ex = True
        for n in tree.traverse():
            if not n.is_root():
                final_call = False
                for or_condition in conditional_output:
                    for condition in or_condition:
                        op = condition[1]
                        if op == 'in':
                            value = condition[0]
                            prop = condition[2]
                            datatype = prop2type[prop]
                            final_call = call(n, prop, datatype, op, value)
                        elif ":" in condition[0]:
                            internal_prop, leaf_prop = condition[0].split(':')
                            value = condition[2]
                            datatype = prop2type[internal_prop]
                            final_call = counter_call(n, internal_prop, leaf_prop, datatype, op, value)
                        else:
                            prop = condition[0]
                            value = condition[2]
                            prop = condition[0]
                            value = condition[2]
                            datatype = prop2type[prop]
                            final_call = call(n, prop, datatype, op, value)
                        if final_call == False:
                            break
                        else:
                            continue
                    if final_call:
                        #print('cut', n.name)
                        n.detach()
                        ex = False
                    else:
                        pass
            else:
                if n.dist == 0: 
                    n.dist = 1
    return tree

def tree2table(tree, internal_node=True, props=[], outfile='tree2table.csv'):
    node2leaves = {}
    leaf2annotations = {}
    if not props:
        leaf = tree.get_farthest_leaf()[0]
        props = list(leaf.props)
        
    with open(outfile, 'w', newline='') as csvfile:
        fieldnames = props
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for node in tree.traverse():
            if node.is_leaf():
                output_row = dict(node.props)
                for k, prop in output_row.items():
                    if type(prop) == list:
                        output_row[k] = '|'.join(str(v) for v in prop)
                writer.writerow(output_row)
            else:
                if internal_node:
                    output_row = dict(node.props)
                    for k, prop in output_row.items():
                        if type(prop) == list:
                            output_row[k] = '|'.join(str(v) for v in prop)
                    writer.writerow(output_row)
                else:
                    pass
    return 

def wrtie_color(color_dict):
    with open('color_dict.txt','w') as f:
        for sub_dict in color_dict:
            for key,value in sub_dict.items():
                if type(value) != dict:
                    f.write('PROPERTY'+'\t'+ key+'\n')
                    f.write('COLOR'+'\t'+ value+'\n')
                    f.write('\n')
                else:
                    f.write('PROPERTY'+'\t'+ key+'\n')
                    for sub_k,sub_v in value.items():
                        f.write('VAR'+'\t'+ sub_k+'\n')
                        f.write('COLOR'+'\t'+ sub_v+'\n')
                    f.write('\n')
    return

def get_layouts(argv_input, layout_name, level, internal_rep): 
    props = []
    layouts = []
    prop_color_dict = {} # key = value, value = color id

    # identify range [1-5], index 1,2,3 and column names
    for i in argv_input.split(','):
        if i[0] == '[' and i[-1] == ']':
            column_start, column_end = get_range(i)
            for j in range(column_start, column_end+1):
                props.append(node_props[j-1])
        else:
            try:
                i = int(i)
                props.append(node_props[i-1])
            except ValueError:
                props.append(i)

    # load layout for each prop
    for idx, prop in enumerate(props):
        color_dict = {} # key = value, value = color id
        if layout_name in ['binary', 'revbinary']:
            
            if columns:
                prop_values = sorted(list(set(columns[prop])))
            else:
                prop_values = sorted(list(set(children_prop_array(annotated_tree, prop))))
            nvals = len(prop_values)

            for i in range(0, nvals): # only positive, negative, NaN, three options
                color_dict[prop_values[i]] = paried_color[i]
            
            color = random_color(h=None)
            if layout_name == 'binary':
                layout = conditional_layouts.LayoutBinary('Binary_'+prop, level, color, color_dict, prop, reverse=False)

            elif layout_name == 'revbinary':
                layout = conditional_layouts.LayoutBinary('ReverseBinary_'+prop, level, color, color_dict, prop, reverse=True)
            internal_prop = prop + '_' + internal_rep
            prop_color_dict[internal_prop] = color_dict
            prop_color_dict[prop] = color

        # numerical layouts
        elif layout_name == 'heatmap':
            layout =  staple_layouts.LayoutHeatmap('Heatmap_'+prop, level, internal_rep, prop)
        
        elif layout_name == 'barplot':
            if prop2type[prop] == 'num':
                size_prop = prop+'_'+internal_rep # using internal prop to set the range in case rank_limit cut all the leaves
            else:
                size_prop = prop

            layout =  staple_layouts.LayoutBarplot(name='Barplot_'+prop, prop=prop, \
                                        color=paried_color[level], size_prop=size_prop, 
                                        column=level, internal_rep=internal_rep
                                        )

            prop_color_dict[prop] = paried_color[level]

        # categorical layouts
        elif layout_name in ['label','rectangular', 'colorbranch']:
            
            if columns:
                prop_values =  sorted(list(set(columns[prop])))
            else:
                prop_values = sorted(list(set(children_prop_array(annotated_tree, prop))))
            nvals = len(prop_values)

            for i in range(0, nvals):
                if nvals <= 14:
                    color_dict[prop_values[i]] = paried_color[i]
                else:
                    color_dict[prop_values[i]] = random_color(h=None)
            
            if layout_name == 'label':
                layout = text_layouts.LayoutText('Label_'+prop, level, color_dict, text_prop = prop)
                #layout = TreeLayout(name=prop+'_'+layout_name, ns=text_layouts.text_layout(prop, level, color_dict, internal_rep))
            
            elif layout_name == 'rectangular':
                layout = text_layouts.LayoutRect('Rectangular_'+prop, level, color_dict, text_prop = prop)
            
            elif layout_name == 'colorbranch':
                layout = text_layouts.LayoutColorbranch('Colorbranch_'+prop, level, color_dict, text_prop = prop)
            
            prop_color_dict[prop] = color_dict
        layouts.append(layout)
        level += 1
        
    return layouts, level, prop_color_dict

def get_range(input_range):
    column_range = input_range[input_range.find("[")+1:input_range.find("]")]
    column_start, column_end = [int(i) for i in column_range.split('-')]
    #column_list_idx = [i for i in range(column_start, column_end+1)]
    return column_start, column_end

def random_color(h=None):
    """Generates a random color in RGB format."""
    if not h:
        h = random.random()
    s = 0.5
    l = 0.5
    return _hls2hex(h, l, s)
 
def _hls2hex(h, l, s):
    return '#%02x%02x%02x' %tuple(map(lambda x: int(x*255),
                                    colorsys.hls_to_rgb(h, l, s)))

def main():
    import time
    global prop2type
    global text_prop, num_prop, bool_prop
    global annotated_tree, node_props, columns

    # get params
    args = read_args()

    total_color_dict = []
    layouts = []
    level = 2 # level 1 is the leaf name

    # parse csv to metadata table
    start = time.time()
    
    if args.metadata:
        if args.no_colnames:
            # property key will be named col1, col2, col3, ... if without headers
            metadata_dict, node_props, columns = parse_csv(args.metadata, no_colnames=args.no_colnames)
        else:
            metadata_dict, node_props, columns = parse_csv(args.metadata)
    else: # annotated_tree
        node_props=[]
        columns = {}
        
        
    #code goes here
    end = time.time()
    print('Time for parse_csv to run: ', end - start)
    
    # parse input tree
    if args.tree:
        if args.tree_type == 'newick':
            tree = ete4_parse(args.tree)
        elif args.tree_type == 'ete':
            with open(args.tree, 'r') as f:
                file_content = f.read()
                tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
        
    # if refer tree from taxadb, input tree will be ignored
    elif args.taxa and args.taxadb:
        tree = ''

    if args.text_prop:
        text_prop = args.text_prop.split(',')
    else:
        text_prop = []

    if args.num_prop:
        num_prop = args.num_prop.split(',')
    else:
        num_prop = []
    
    if args.bool_prop:
        bool_prop = args.bool_prop.split(',')
    else:
        bool_prop = []

    if args.text_prop_idx:
        text_prop_idx = []
        for i in args.text_prop_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                text_prop_start, text_prop_end = get_range(i)
                for j in range(text_prop_start, text_prop_end+1):
                    text_prop_idx.append(j)
            else:
                text_prop_idx.append(int(i))

        text_prop = [node_props[index-1] for index in text_prop_idx]
    
    if args.num_prop_idx:
        num_prop_idx = []
        for i in args.num_prop_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                num_prop_start, num_prop_end = get_range(i)
                for j in range(num_prop_start, num_prop_end+1):
                    num_prop_idx.append(j)
            else:
                num_prop_idx.append(int(i))

        num_prop = [node_props[index-1] for index in num_prop_idx]

    if args.bool_prop_idx:
        bool_prop_idx = []
        for i in args.bool_prop_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                bool_prop_start, bool_prop_end = get_range(i)
                for j in range(bool_prop_start, bool_prop_end+1):
                    bool_prop_idx.append(j)
            else:
                bool_prop_idx.append(int(i))

        bool_prop = [node_props[index-1] for index in bool_prop_idx]

    #rest_prop = []
    rest_prop = list(set(node_props) - set(text_prop) - set(num_prop) - set(bool_prop))
    
    # output datatype of each property of each tree node including internal nodes
    prop2type = {
        # start with leaf name
        'name':'str',
        'dist':'num',
        'support':'num',
        # taxonomic features
        'rank':'str',
        'sci_name':'str',
        'taxid':'str',
        'lineage':'str',
        'named_lineage':'str'
        } 

    for prop in text_prop+bool_prop:
        prop2type[prop] = 'str'
        prop2type[prop+'_counter'] = 'str'
    for prop in num_prop:
        prop2type[prop] = 'num'
        prop2type[prop+'_avg'] = 'num'
        prop2type[prop+'_sum'] = 'num'
        prop2type[prop+'_max'] = 'num'
        prop2type[prop+'_min'] = 'num'
        prop2type[prop+'_std'] = 'num'
    for prop in rest_prop:
        prop2type[prop] = 'str'
        
    

    # load annotations to leaves
    start = time.time()
    
    taxon_column = []
    # load all metadata to leaf nodes
    if not args.annotated_tree:
        if args.taxon_column: # to identify taxon column as taxa property from metadata
            taxon_column.append(args.taxon_column)
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type, taxon_column=args.taxon_column, taxon_delimiter=args.taxon_delimiter)
        else:
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type)
    else:
        # if args.tree_type == 'newick':
            
        # with open(args.annotated_tree, 'r') as f:
        #     file_content = f.read()
        #     annotated_tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
        annotated_tree = tree

    end = time.time()
    print('Time for load_metadata_to_tree to run: ', end - start)
    
    # stat method
    if args.counter_stat:
        counter_stat = args.counter_stat

    if args.num_stat:
        num_stat = args.num_stat

    # merge annotations depends on the column datatype
    start = time.time()
    if not args.annotated_tree:
        #pre load node2leaves to save time
        node2leaves = annotated_tree.get_cached_content()
        count = 0
        for node in annotated_tree.traverse("postorder"):
            internal_props = {}
            if node.is_leaf():
                pass
            else:
                
                if text_prop:
                    internal_props_text = merge_text_annotations(node2leaves[node], text_prop, counter_stat=counter_stat)
                    internal_props.update(internal_props_text)

                if num_prop:
                    internal_props_num = merge_num_annotations(node2leaves[node], num_prop, num_stat=num_stat)
                    internal_props.update(internal_props_num)

                if bool_prop:
                    internal_props_bool = merge_text_annotations(node2leaves[node], bool_prop, counter_stat=counter_stat)
                    internal_props.update(internal_props_bool)

                # deprecated
                # if rest_column:
                #     internal_props_rest = merge_text_annotations(node2leaves[node], rest_column, counter_stat=counter_stat)
                #     internal_props.update(internal_props_rest)
                
                #internal_props = {**internal_props_text, **internal_props_num, **internal_props_rest}
                #print(internal_props.items())
                for key,value in internal_props.items():
                    node.add_prop(key, value)
    else:
        pass
    end = time.time()
    print('Time for merge annotations to run: ', end - start)
    
    
    # taxa annotations
    start = time.time()
    if args.taxonomic_profile:
        if not args.taxadb:
            print('Please specify which taxa db using --taxadb <GTDB|NCBI>')
        else:
            if args.taxadb == 'GTDB':
                if args.taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr=args.taxon_column)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name")
            elif args.taxadb == 'NCBI':
                if args.taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr=args.taxon_column, sp_delimiter=args.taxon_delimiter, sp_field=args.taxa_field)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name", sp_delimiter=args.taxon_delimiter, sp_field=args.taxa_field)

    else:
        rank2values = {}
    end = time.time()
    print('Time for annotate_taxa to run: ', end - start)
    ### Anslysis settings###

    # # prune tree by rank
    if args.rank_limit:
        annotated_tree= taxatree_prune(annotated_tree, rank_limit=args.rank_limit)

    # # prune tree by condition 
    # if args.pruned_by: # need to be wrap with quotes
    #     condition_strings = args.pruned_by
    #     annotated_tree= conditional_prune(annotated_tree, condition_strings, prop2type)

    # collapse tree by condition 
    if args.collapsed_by: # need to be wrap with quotes
        condition_strings = args.collapsed_by
        for condition in condition_strings:
            c_layout = TreeLayout(name='Collapsed_by_'+condition, \
                                    ns=conditional_layouts.collapsed_by_layout(condition, prop2type = prop2type, level=level))
            layouts.append(c_layout)

    # label node by condition
    if args.highlighted_by: # need to be wrap with quotes
        condition_strings = args.highlighted_by
        for condition in condition_strings:
            s_layout = TreeLayout(name='Highlighted_by_'+condition, \
                                    ns=conditional_layouts.highlight_layout(condition, prop2type = prop2type, level=level))
            layouts.append(s_layout)
    
    #### Layouts settings ####
    # Taxa layouts
    if args.TaxonLayout:
        taxon_color_dict = {}
        taxa_layouts = []
        
        # generate a rank2values dict for pre taxonomic annotated tree
        if not rank2values:
            rank2values = defaultdict(list)
            for n in tree.traverse():
                if n.props.get('rank') and n.props.get('rank') != 'Unknown':
                    rank2values[n.props.get('rank')].append(n.props.get('sci_name',''))
        else:
            pass
        
        # assign color for each value of each rank
        for rank, value in sorted(rank2values.items()):
            color_dict = {} 
            nvals = len(value)
            for i in range(0, nvals):
                if nvals <= 14:
                    color_dict[value[i]] = paried_color[i]
                else:
                    color_dict[value[i]] = random_color(h=None)

            layout = taxon_layouts.TaxaClade(name='TaxaClade_'+rank, level=level, rank = rank, color_dict=color_dict)
            taxa_layouts.append(layout)
            taxon_color_dict[rank] = color_dict
        layouts = layouts + taxa_layouts
        level += 1
        total_color_dict.append(taxon_color_dict)

    # numerical representative mearsure 
    if args.num_stat != 'all':
        internal_num_rep = args.num_stat
    else:
        internal_num_rep = args.internal_plot_measure

    # get layouts
    if args.HeatmapLayout:
        heatmap_layouts, level, _ = get_layouts(args.HeatmapLayout, 'heatmap', level, internal_num_rep)
        layouts.extend(heatmap_layouts)

    if args.BarplotLayout:
        barplot_layouts, level,color_dict = get_layouts(args.BarplotLayout, 'barplot', level, internal_num_rep)
        layouts.extend(barplot_layouts)
        total_color_dict.append(color_dict)

    # categorical and boolean 
    if args.ColorbranchLayout:
        colorbranch_layouts, level, color_dict = get_layouts(args.ColorbranchLayout, 'colorbranch', level, 'counter')
        layouts.extend(colorbranch_layouts)
        total_color_dict.append(color_dict)

    if args.RectangularLayout:
        rectangular_layouts, level, color_dict = get_layouts(args.RectangularLayout, 'rectangular', level, 'counter')
        layouts.extend(rectangular_layouts)
        total_color_dict.append(color_dict)
        
    if args.LabelLayout:
        label_layouts, level, color_dict = get_layouts(args.LabelLayout, 'label', level, 'counter')
        layouts.extend(label_layouts)
        total_color_dict.append(color_dict)

    if args.BinaryLayout:
        label_layouts, level, color_dict = get_layouts(args.BinaryLayout, 'binary', level, 'counter')
        layouts.extend(label_layouts)
        total_color_dict.append(color_dict)

    if args.RevBinaryLayout:
        label_layouts, level, color_dict = get_layouts(args.RevBinaryLayout, 'revbinary', level, 'counter')
        layouts.extend(label_layouts)
        total_color_dict.append(color_dict)
    #print(total_color_dict)
    #### prune at the last step in case of losing leaves information
    # prune tree by rank
    # if args.rank_limit:
    #     annotated_tree= taxatree_prune(annotated_tree, rank_limit=args.rank_limit)
    # prune tree by condition 
    if args.pruned_by: # need to be wrap with quotes
        condition_strings = args.pruned_by
        annotated_tree= conditional_prune(annotated_tree, condition_strings, prop2type)

    #### Output #####
    if args.out_color_dict:
        wrtie_color(total_color_dict)
    
    if args.ete4out:
        with open('./annotated_tree.ete', 'w') as f:
            f.write(b64pickle.dumps(annotated_tree, encoder='pickle', pack=False))
    
    if args.outtree:
        annotated_tree.write(outfile=args.outtree, properties = [], format=1)

    if args.interactive:
        if args.annotated_tree:
            if args.tree_type == 'ete':
                # exisiting props in internal node
                existing_internal_props = list(annotated_tree.props.keys())
                # exisiting props in leaf node
                existing_leaf_props = list(annotated_tree.get_farthest_leaf()[0].props.keys()) 
                popup_prop_keys = list(set(existing_internal_props+existing_leaf_props))
            elif args.tree_type == 'newick':
                # props which add in the arguments
                required_internal_props = list(prop2type.keys()) 
                # exisiting prop in leaf node
                existing_leaf_props = list(annotated_tree.get_farthest_leaf()[0].props.keys()) 
                popup_prop_keys = list(set(required_internal_props + existing_leaf_props))
        else:
            # all the metadata to the leaves, no internal
            popup_prop_keys = list(prop2type.keys())
        annotated_tree.explore(tree_name='example',layouts=layouts, port=args.port, popup_prop_keys=sorted(popup_prop_keys))
    elif args.plot:
        plot(annotated_tree, layouts, args.port, args.plot)
    if args.outtsv:
        tree2table(annotated_tree, internal_node=True, outfile=args.outtsv)

    
    return annotated_tree

if __name__ == '__main__':
    main()