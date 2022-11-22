#!/usr/bin/env python


from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from layouts import text_layouts, taxon_layouts, staple_layouts, heatmap_layouts, conditional_layouts

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


#colours_50 = ["#E41A1C","#C72A35","#AB3A4E","#8F4A68","#735B81","#566B9B","#3A7BB4","#3A85A8","#3D8D96","#419584","#449D72","#48A460","#4CAD4E","#56A354","#629363","#6E8371","#7A7380","#87638F","#93539D","#A25392","#B35A77","#C4625D","#D46A42","#E57227","#F67A0D","#FF8904","#FF9E0C","#FFB314","#FFC81D","#FFDD25","#FFF12D","#F9F432","#EBD930","#DCBD2E","#CDA12C","#BF862B","#B06A29","#A9572E","#B65E46","#C3655F","#D06C78","#DE7390","#EB7AA9","#F581BE","#E585B8","#D689B1","#C78DAB","#B791A5","#A8959F","#999999"]
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
        required=True,
        help="<metadata.csv> .csv, .tsv. mandatory input")
    group.add_argument('--no_colnames',
        default=False,
        action='store_true',
        required=False,
        help="metadata table doesn't contain columns name")
    group.add_argument('--text_column',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as categorical data")
    group.add_argument('--num_column',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as numerical data")
    group.add_argument('--bool_column',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as boolean data")
    group.add_argument('--text_column_idx',
        type=str,
        required=False,
        help="1,2,3 or 1-5 index of columns which need to be read as categorical data")
    group.add_argument('--num_column_idx',
        type=str,
        required=False,
        help="1,2,3 or 1-5 index columns which need to be read as numerical data")
    group.add_argument('--bool_column_idx',
        type=str,
        required=False,
        help="1,2,3 or 1-5 index columns which need to be read as boolean data")
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

    group.add_argument('--text_stat',
        default='raw',
        type=str,
        required=False,
        help="statistic calculation to perform for categorical data in internal nodes, raw count or in percentage [raw, relative] ")  
    
    group.add_argument('--rank_limit',
        type=str,
        required=False,
        help="TAXONOMIC_LEVEL prune annotate tree by rank limit")
    group.add_argument('--collapsed_by', 
        type=str,
        required=False,
        help='target tree collapsed by customized conditions')
    group.add_argument('--highlighted_by', 
        type=str,
        required=False,
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
        help="<col1,col2> names of columns which need to be plot as BinaryLayout")
    group.add_argument('--RevBinaryLayout',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be plot as RevBinaryLayout")

    group.add_argument('--ColorbranchLayout',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be plot as Textlayouts")
    group.add_argument('--LabelLayout',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be plot as LabelLayout")
    group.add_argument('--RectangularLayout',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be plot as RectangularLayout")
    
    
    group.add_argument('--HeatmapLayout',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as HeatmapLayout")
    group.add_argument('--BarplotLayout',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as BarplotLayouts")
    
    
    group.add_argument('--TaxonLayout',
        type=str,
        required=False,
        help="<col1,col2> names of columns which need to be read as TaxonLayouts")

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

def parse_csv(input_file, delimiter='\t', no_colnames=False):
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
            metadata[nodename] = dict(row)
            for (k,v) in row.items(): # go over each column name and value 
                columns[k].append(v) # append the value into the appropriate list
                                    # based on column name k

    return metadata, node_props, columns

def parse_csv_to_column(metadata):
    columns = defaultdict(list) # each value in each column is appended to a list
    with open(metadata) as f:
        reader = csv.DictReader(f, delimiter="\t") # read rows into a dictionary format
        for row in reader: # read a row as {column1: value1, column2: value2,...}
            for (k,v) in row.items(): # go over each column name and value 
                columns[k].append(v) # append the value into the appropriate list
                                    # based on column name k
    return columns

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

def load_metadata_to_tree(tree, metadata_dict, taxon_column=None, taxon_delimiter=';'):
    for node, props in metadata_dict.items():
        hits = tree.search_nodes(name=node) # including internal nodes
        if hits:
            target_node = hits[0]
            for key,value in props.items():
                if key == taxon_column:
                    taxon_prop = value.split(taxon_delimiter)[-1]
                    target_node.add_prop(key, taxon_prop)
                else:
                    target_node.add_prop(key, value)
        else:
            pass
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

def merge_text_annotations(nodes, target_props, text_stat='raw'):
    internal_props = {}
    for target_prop in target_props:
        if text_stat == 'raw':
            prop_list = children_prop_array(nodes, target_prop)
            internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])

        elif text_stat == 'relative':
            prop_list = children_prop_array(nodes, target_prop)
            internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])

        else:
            print('Invalid stat method')
            break
    return internal_props

def merge_num_annotations(nodes, target_props, num_stat='all'):
    internal_props = {}
    for target_prop in target_props:
        prop_array = np.array(children_prop_array(nodes, target_prop),dtype=np.float64)
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
    return internal_props

def add_suffix(name, suffix, delimiter='_'):
    return str(name) + delimiter + str(suffix)

def children_prop_array(nodes, prop):
    array = [n.props.get(prop) for n in nodes if n.props.get(prop)] 
    return array

def annotate_taxa(tree, db="GTDB", taxid_attr="name", sp_delimiter='.', sp_field=0):
    # def return_spcode(leaf):
    #     try:
    #         return leaf.name.split(sp_delimiter)[sp_field]
    #     except IndexError:
    #         return leaf.name

    if db == "GTDB":
        gtdb = GTDBTaxa()
        gtdb.annotate_tree(tree,  taxid_attr=taxid_attr)
    elif db == "NCBI":
        ncbi = NCBITaxa()
        # extract sp codes from leaf names
        #tree.set_species_naming_function(return_spcode)
        ncbi.annotate_tree(tree, taxid_attr=taxid_attr)

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

def get_layouts(argv_input, layout_name, level, internal_rep):
    props = []
    layouts = []

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
        # numerical layouts
        if layout_name == 'heatmap':
            layout =  TreeLayout(name=prop+'_'+layout_name, ns=heatmap_layouts.heatmap_layout(prop, level, internal_rep))
        
        elif layout_name == 'barplot':
            layout =  staple_layouts.LayoutBarplot(name=prop+'_'+layout_name, prop=prop, \
                                    color_prop=paried_color[level], size_prop=prop, 
                                    column=level, internal_rep=internal_rep
                                    )
        
        # categorical layouts
        elif layout_name in ['label','rectangular', 'colorbranch']:
            colour_dict = {} # key = value, value = colour id
            prop_values = list(set(columns[prop]))
            #prop_values = list(set(children_prop_array(annotated_tree, prop)))
            nvals = len(prop_values)

            for i in range(0, nvals):
                if nvals <= 14:
                    colour_dict[prop_values[i]] = paried_color[i]
                else:
                    colour_dict[prop_values[i]] = random_color(h=None)
            
            if layout_name == 'label':
                layout = TreeLayout(name=prop+'_'+layout_name, ns=text_layouts.text_layout(prop, level, colour_dict, internal_rep))
            
            elif layout_name == 'rectangular':
                layout = TreeLayout(name=prop+'_'+layout_name, ns=text_layouts.rectangular_layout(prop, level, colour_dict, internal_rep))
            
            elif layout_name == 'colorbranch':
                layout = TreeLayout(name=prop+'_'+layout_name, ns=text_layouts.label_layout(prop, level, colour_dict, internal_rep))

        elif layout_name in ['binary', 'revbinary']:
            prop_colour_dict = {} # key = value, value = colour id
            prop_values = list(set(columns[prop]))
            nvals = len(prop_values)
            for i in range(0, nvals):
                prop_colour_dict[prop_values[i]] = paried_color[i]

            color = random_color(h=None)
            #color = random_color(h=None)
            if layout_name == 'binary':
                layout = TreeLayout(name=prop+'_'+layout_name, ns=conditional_layouts.boolean_layout(prop, level, color, prop_colour_dict, internal_rep))

            elif layout_name == 'revbinary':
                layout = TreeLayout(name=prop+'_'+layout_name, ns=conditional_layouts.boolean_layout(prop, level, color, prop_colour_dict ,internal_rep, reverse=True))
        layouts.append(layout)
        level += 1
    return layouts, level

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
    global annotated_tree, node_props, columns
    args = read_args()

    layouts = []
    level = 2 # level 1 is the leaf name

    # parse csv to metadata table
    if args.metadata:
        if args.no_colnames:
            metadata_dict, node_props, columns = parse_csv(args.metadata, no_colnames=args.no_colnames)
        else:
            metadata_dict, node_props, columns = parse_csv(args.metadata)
            

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
    
    if args.bool_column:
        bool_column = args.bool_column.split(',')
    else:
        bool_column = []

    if args.text_column_idx:
        text_column_idx = []
        for i in args.text_column_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                text_column_start, text_column_end = get_range(i)
                for j in range(text_column_start, text_column_start+1):
                    text_column_idx.append(j)
            else:
                text_column_idx.append(int(i))

        text_column = [node_props[index-1] for index in text_column_idx]

    if args.num_column_idx:
        num_column_idx = []
        for i in args.num_column_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                num_column_start, num_column_end = get_range(i)
                for j in range(num_column_start, num_column_start+1):
                    num_column_idx.append(j)
            else:
                num_column_idx.append(int(i))

        num_column = [node_props[index-1] for index in num_column_idx]

    if args.bool_column_idx:
        bool_column_idx = []
        for i in args.bool_column_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                bool_column_start, bool_column_end = get_range(i)
                for j in range(bool_column_start, bool_column_end+1):
                    bool_column_idx.append(j)
            else:
                bool_column_idx.append(int(i))

        bool_column_idx = [node_props[index-1] for index in bool_column_idx]


    rest_column = list(set(node_props) - set(text_column) - set(num_column) - set(bool_column))
    
    # load annotations to leaves
    if args.taxon_column:
        annotated_tree = load_metadata_to_tree(tree, metadata_dict, args.taxon_column, args.taxon_delimiter)
    else:
        annotated_tree = load_metadata_to_tree(tree, metadata_dict)
    
    # stat method
    if args.text_stat:
        text_stat = args.text_stat

    if args.num_stat:
        num_stat = args.num_stat

    # merge annotations
    node2leaves = annotated_tree.get_cached_content()
    for node in annotated_tree.traverse("postorder"):
        internal_props = {}
        if node.is_leaf():
            pass
        else:
            if text_column:
                internal_props_text = merge_text_annotations(node2leaves[node], text_column, text_stat=text_stat)
                internal_props.update(internal_props_text)
            
            if num_column:
                internal_props_num = merge_num_annotations(node2leaves[node], num_column, num_stat=num_stat)
                internal_props.update(internal_props_num)

            if bool_column:
                internal_props_bool = merge_text_annotations(node2leaves[node], text_column, text_stat=text_stat)
                internal_props.update(internal_props_bool)

            if rest_column:
                internal_props_rest = merge_text_annotations(node2leaves[node], rest_column, text_stat=text_stat)
                internal_props.update(internal_props_rest)
            
            #internal_props = {**internal_props_text, **internal_props_num, **internal_props_rest}
            for key,value in internal_props.items():
                node.add_prop(key, value)
    
    # taxa annotations
    if args.taxonomic_profile:
        if not args.taxadb:
            print('Please specify which taxa db using --taxadb <GTDB|NCBI>')
        else:
            if args.taxon_column:
                annotated_tree = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr=args.taxon_column)
            else:
                annotated_tree = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name")
        # if args.taxon_column:
        #     annotated_tree = annotate_taxa(annotated_tree, taxid_attr=taxon_column)
        # else:
        #     annotated_tree = annotate_taxa(annotated_tree, taxid_attr="name")
    
    ### Anslysis settings###

    # prune tree by rank
    if args.rank_limit:
        annotated_tree= taxatree_prune(annotated_tree, rank_limit=args.rank_limit)

    # collapse tree by condition 
    if args.collapsed_by: # need to be wrap with quotes
        print(args.collapsed_by)

    # label node by condition
    if args.highlighted_by: # need to be wrap with quotes
        condition_strings = args.highlighted_by

        prop2type = {'name':'str'} # start with leaf name
        for prop in text_column+rest_column:
            prop2type[prop] = 'str'
        for prop in num_column:
            prop2type[prop] = 'num'

        s_layout = TreeLayout(name='highlighted', \
                                ns=conditional_layouts.highlight_layout(condition_strings, prop2type = prop2type, level=level))
        layouts.append(s_layout)
        
    #### Layouts settings ####
    # numerical
    if args.num_stat != 'all':
        internal_num_rep = args.num_stat
    else:
        internal_num_rep = args.internal_plot_measure

    if args.HeatmapLayout:
        heatmap_layouts, level = get_layouts(args.HeatmapLayout, 'heatmap', level, internal_num_rep)
        layouts.extend(heatmap_layouts)

    if args.BarplotLayout:
        barplot_layouts, level = get_layouts(args.BarplotLayout, 'barplot', level, internal_num_rep)
        layouts.extend(barplot_layouts)

    # boolean 
    
    # categorical
    if args.ColorbranchLayout:
        colorbranch_layouts, level = get_layouts(args.ColorbranchLayout, 'colorbranch', level, 'counter')
        layouts.extend(colorbranch_layouts)

    if args.RectangularLayout:
        rectangular_layouts, level = get_layouts(args.RectangularLayout, 'rectangular', level, 'counter')
        layouts.extend(rectangular_layouts)
        
    if args.LabelLayout:
        label_layouts, level = get_layouts(args.LabelLayout, 'label', level, 'counter')
        layouts.extend(label_layouts)

    if args.BinaryLayout:
        label_layouts, level = get_layouts(args.BinaryLayout, 'binary', level, 'counter')
        layouts.extend(label_layouts)

    if args.RevBinaryLayout:
        label_layouts, level = get_layouts(args.RevBinaryLayout, 'revbinary', level, 'counter')
        layouts.extend(label_layouts)

    # Taxa layouts
    if args.TaxonLayout:
        taxon_prop = args.TaxonLayout
        # taxa_layouts = [
        #     TreeLayout(name='level3_class', ns=taxon_layouts.class_layout()),
        # ]


        taxa_layouts = [
            TreeLayout(name='level1_kingdom', ns=taxon_layouts.collapse_kingdom()),
            TreeLayout(name='level2_phylum', ns=taxon_layouts.collapse_phylum()),
            TreeLayout(name='level3_class', ns=taxon_layouts.collapse_class()),
            TreeLayout(name='level4_order', ns=taxon_layouts.collapse_order()),
            TreeLayout(name='level5_family', ns=taxon_layouts.collapse_family()),
            TreeLayout(name='level6_genus', ns=taxon_layouts.collapse_genus()),
            TreeLayout(name='level7_species', ns=taxon_layouts.collapse_species()),
        ]


        layouts = layouts + taxa_layouts
    #### Output #####
    if not args.interactive:
        if args.outtree:
            annotated_tree.write(outfile=args.plot, format=1)
        elif args.plot:
            annotated_tree.explore(tree_name='example',layouts=[], port=5000)
    else:
        annotated_tree.explore(tree_name='example',layouts=layouts, port=5000)
    
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