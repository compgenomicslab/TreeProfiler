#!/usr/bin/env python
import numbers
import math
import colorsys
import random
import sys
import os

from collections import defaultdict
from itertools import islice
from io import StringIO
import numpy as np


from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from treeprofiler.tree_image import get_image
from treeprofiler.layouts import (
    text_layouts, taxon_layouts, staple_layouts, 
    conditional_layouts, seq_layouts, profile_layouts)
from treeprofiler.src import b64pickle
from treeprofiler.src.utils import (
    ete4_parse, taxatree_prune, conditional_prune,
    children_prop_array, children_prop_array_missing, 
    flatten, get_consensus_seq)

paried_color = ["red", "darkblue", "lightgreen", "sienna", "lightCoral", "violet", "mediumturquoise",   "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal", "darkyellow"]

DESC = "plot tree"

def poplulate_plot_args(plot_args_p):
    """
    Parse the input parameters
    Return the parsed arguments.
    """
    group = plot_args_p.add_argument_group(title='Conditional display arguments',
        description="Conditional display  parameters")
    
    group.add_argument('--internal-plot-measure',
        default='avg',
        type=str,
        required=False,
        help="statistic measures to be shown in numerical layout for internal nodes, [default: avg]")  

    group.add_argument('--collapsed-by', 
        type=str,
        required=False,
        action='append',
        help='target tree nodes collapsed by customized conditions')  
    group.add_argument('--highlighted-by', 
        type=str,
        required=False,
        action='append',
        help='target tree nodes highlighted by customized conditions')
        
    # group = plot_args_p.add_argument_group(title='Basic treelayout arguments',
    #     description="treelayout parameters")

    # group.add_argument('--drawer',
    #     type=str,
    #     required=False,
    #     help="Circular or Rectangular")
    # group.add_argument('--collapse_level',
    #     type=str,
    #     required=False,
    #     help="default collapse level, default is 10") 
    # group.add_argument('--ultrametric',
    #     default=False,
    #     action='store_true',
    #     required=False,
    #     help="ultrametric tree")

    group = plot_args_p.add_argument_group(title="Properties' layout arguments",
        description="Prop layout parameters")
    group.add_argument('--column-width',
        type=int,
        default=20,
        help="customize column width of each layout.[default: 20]"
    )
    group.add_argument('--barplot-width',
        type=int,
        default=200,
        help="customize barplot width of barplot layout.[default: 200]"
    )
    # group.add_argument('--profiling_width',
    #     type=int,
    #     default=None,
    #     help="customize profiling width of each profiling layout."
    # )
    group.add_argument('--padding-x',
        type=int,
        default=1,
        help="customize horizontal column padding distance of each layout.[default: 1]"
    )
    group.add_argument('--padding-y',
        type=int,
        default=0,
        help="customize vertical padding distance of each layout.[default: 0]"
    )
    group.add_argument('--binary-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of properties which need to be plot as binary-layout which highlights the postives")
    group.add_argument('--revbinary-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of properties which need to be plot as revbinary-layout which highlights the negatives")
    group.add_argument('--colorbranch-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of properties where branches will be colored based on different values.")
    group.add_argument('--label-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of properties where values will be displayed on the aligned panel.")
    group.add_argument('--rectangle-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names  of properties where values will be label as rectangular color block on the aligned panel.")
    group.add_argument('--heatmap-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of numerical properties which need to be read as heatmap-layout")
    group.add_argument('--barplot-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of numerical properties which need to be read as barplot_layouts")
    group.add_argument('--taxonclade-layout',
        default=False,
        action='store_true',
        help="Activate taxonclade_layout which clades will be colored based on taxonomy of each node.")
    group.add_argument('--taxonrectangle-layout',
        default=False,
        action='store_true',
        help="Activate taxonrectangle-layout which taxonomy of each node will be display as rectangular blocks in aligned panel.")
    group.add_argument('--emapper-layout',
        default=False,
        action='store_true',
        help="Activate emapper_layout which display all the annotation from EggNOG-mapper.") #domain_layout
    group.add_argument('--domain-layout',
        default=False,
        action='store_true',
        help="Activate domain_layout which display protein domain annotation in sequence.") #domain_layout
    group.add_argument('--alignment-layout',
        default=False,
        action='store_true',
        help="Display Multiple Sequence Alignment layout in aligned panel.")
    group.add_argument('--profiling-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of properties which need to be convert to presence-absence profiling matrix of each value")
    group.add_argument('--multi-profiling-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names of properties containing values as list which need to be convert to presence-absence profiling matrix")
    group.add_argument('--categorical-matrix-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names which need to be plot as categorical_matrix_layout for categorical values")
    group.add_argument('--numerical-matrix-layout',
        nargs='+',
        required=False,
        help="<prop1,prop2> names which need to be plot as numerical_matrix_layout for numerical values ")

    group = plot_args_p.add_argument_group(title='Visualizing output arguments',
        description="Visualizing output parameters")
    # group.add_argument('--interactive',
    #     default=False,
    #     action='store_true',
    #     help="run interactive session")
    group.add_argument('--verbose',
        action="store_false", 
        required=False,
        help="show detail on prompt when visualizing taget tree.")
    group.add_argument('--port',
        type=str,
        default=5000,
        help="run interactive session on custom port.[default: 5000]")
    group.add_argument('--plot',
        type=str,
        required=False,
        help="output as pdf")
    group.add_argument('--out-colordict',
        action="store_true", 
        required=False,
        help="print color dictionary of each property")


### visualize tree
def run(args):
    global prop2type, properties, tree
    node_props=[]
    properties = {}
    rank2values = {}

    total_color_dict = []
    layouts = []
    level = 1 # level 1 is the leaf name

    # checking file and output exists
    if not os.path.exists(args.tree):
       raise FileNotFoundError(f"Input tree {args.tree} does not exist.")

    #parse input tree
    if args.tree:
        if args.input_type == 'newick':
            try:
                tree = ete4_parse(open(args.tree), internal_parser=args.internal_parser)
            except Exception as e:
                print(e)
                sys.exit(1)
        # elif args.input_type == 'nexus':
        #    try:
        #         tree = ete4_parse(open(args.tree), internal_parser=args.internal_parser)
        #     except Exception as e:
        #         print(e)
        #         sys.exit(1)
        elif args.input_type == 'ete':
            try:
                with open(args.tree, 'r') as f:
                    file_content = f.read()
                    tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
            except ValueError as e:
                print(e)
                print("In valid ete format.")
                sys.exit(1)
    
    #rest_prop = []
    if args.prop2type:
        prop2type = {}
        with open(args.prop2type, 'r') as f:
            for line in f:
                line = line.rstrip()
                prop, value = line.split('\t')
                prop2type[prop] = eval(value)
        
        popup_prop_keys = list(prop2type.keys()) 

    else:
        prop2type = {# start with leaf name
                'name':str,
                'dist':float,
                'support':float,
                'rank': str,
                'sci_name': str,
                'taxid': str,
                'lineage':str,
                'named_lineage': str,
                'evoltype': str,
                'dup_sp': str,
                'dup_percent': float,
                }
        popup_prop_keys = list(prop2type.keys()) 

        if args.input_type == 'ete':
            leafa, _, leafb, _ = tree._get_farthest_and_closest_leaves()
            leaf_prop2type = get_prop2type(leafa)
            leaf_prop2type.update(get_prop2type(leafb))
            internal_node_prop2type = get_prop2type(tree)
            prop2type.update(leaf_prop2type)
            prop2type.update(internal_node_prop2type)
        
        # elif args.input_type == 'newick':
        #     popup_prop_keys = list(prop2type.keys()) 
    
    
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
    # numerical representative mearsure 
    internal_num_rep = args.internal_plot_measure

    # Get the input arguments in order
    input_order = []
    for arg in sys.argv[1:]:
        if arg.startswith('-') and arg.endswith('layout'):
            input_order.append(arg[2:])
        else:
            continue
    
    for layout in input_order:
        if layout == 'heatmap-layout':
            heatmap_layouts, level = get_heatmap_layouts(tree, args.heatmap_layout, level, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y, internal_rep=internal_num_rep)
            layouts.extend(heatmap_layouts)

        if layout == 'label-layout':
            label_layouts, level, color_dict = get_label_layouts(tree, args.label_layout, level, prop2type=prop2type, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)

        if layout == 'colorbranch-layout':
            colorbranch_layouts, level, color_dict = get_colorbranch_layouts(tree, args.colorbranch_layout, level, prop2type=prop2type, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y)
            layouts.extend(colorbranch_layouts)
            total_color_dict.append(color_dict)

        if layout == 'rectangle-layout':
            rectangle_layouts, level, color_dict = get_rectangle_layouts(tree, args.rectangle_layout, level, prop2type=prop2type, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y)
            layouts.extend(rectangle_layouts)
            total_color_dict.append(color_dict)

        if layout == 'binary-layout':
            label_layouts, level, color_dict = get_binary_layouts(tree, args.binary_layout, level, prop2type=prop2type, column_width=args.column_width, reverse=False, padding_x=args.padding_x, padding_y=args.padding_y)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)

        if layout == 'revbinary-layout':
            label_layouts, level, color_dict = get_binary_layouts(tree, args.revbinary_layout, level, prop2type=prop2type, column_width=args.column_width, reverse=True,  padding_x=args.padding_x, padding_y=args.padding_y)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)
        
        if layout == 'barplot-layout':
            barplot_layouts, level,color_dict = get_barplot_layouts(tree, args.barplot_layout, level, prop2type, column_width=args.barplot_width, padding_x=args.padding_x, padding_y=args.padding_y, internal_rep=internal_num_rep)
            layouts.extend(barplot_layouts)
            total_color_dict.append(color_dict)

        if layout == 'alignment-layout':
            #fasta_file = args.alignment_layout
            lengh = len(max(children_prop_array(tree, 'alignment'),key=len))
            aln_layout = seq_layouts.LayoutAlignment(name='Alignment_layout', 
                        alignment_prop='alignment', column=level, scale_range=lengh,
                        summarize_inner_nodes=True)
            layouts.append(aln_layout)

        if layout == 'domain-layout':
            domain_layout = seq_layouts.LayoutDomain(name="Domain_layout", prop='dom_arq')
            layouts.append(domain_layout)
        
        # presence-absence profiling based on categorical data
        if layout == 'profiling-layout':
            profiling_props = args.profiling_layout
            for profiling_prop in profiling_props:
                matrix, all_values = single2profile(tree, profiling_prop)
                profile_layout = profile_layouts.LayoutProfile(name=f'Profiling_{profiling_prop}', mode='multi',
                    alignment=matrix, seq_format='profiles', profiles=all_values, column=level, summarize_inner_nodes=True, poswidth=args.column_width)
                level += 1
                layouts.append(profile_layout)
        
        # presence-absence profiling based on list data
        if layout == 'multi-profiling-layout':
            profiling_props = args.multi_profiling_layout
            for profiling_prop in profiling_props:
                matrix, all_values = multiple2profile(tree, profiling_prop)
                profile_layout = profile_layouts.LayoutProfile(name=f'Profiling_{profiling_prop}', mode='multi',
                alignment=matrix, seq_format='profiles', profiles=all_values, column=level, summarize_inner_nodes=False, poswidth=args.column_width)
                level += 1
                layouts.append(profile_layout)
        
        # categorical matrix
        if layout == 'categorical-matrix-layout':
            profiling_props = args.categorical_matrix_layout
            matrix, value2color = props2matrix(tree, profiling_props, dtype=str)
            profile_layout = profile_layouts.LayoutProfile(name='categorical_matrix_layout', mode='single',
                alignment=matrix, seq_format='categories', profiles=profiling_props, value_color=value2color, column=level, poswidth=args.column_width)
            level += 1
            layouts.append(profile_layout)

        # numerical matrix
        if layout == 'numerical-matrix-layout':
            profiling_props = args.numerical_matrix_layout
            matrix, maxval, minval = props2matrix(tree, profiling_props)
            #profile_layout = TreeLayout(name='numerical_profiling_layout', ns=get_alnface(alignment, level), aligned_faces = True)
            profile_layout = profile_layouts.LayoutProfile(name='numerical_matrix_layout', mode='numerical', 
                alignment=matrix, seq_format='gradients', profiles=profiling_props, value_range=[minval, maxval], column=level, poswidth=args.column_width)
            level += 1
            layouts.append(profile_layout)

    # emapper layout 
    if args.emapper_layout:
        text_props = [
            'seed_ortholog',
            'max_annot_lvl',
            'COG_category',
            'Description',
            'Preferred_name',
            ]
        #label_layouts, level, _ = get_layouts(tree, text_props, 'rectangular', level, 'counter', prop2type=prop2type)
        label_layouts, level, _ = get_rectangle_layouts(tree, text_props, level, prop2type=prop2type, column_width=args.column_width)
            
        layouts.extend(label_layouts)
        
        num_props = [
            #'evalue',
            'score'
        ]
        #barplot_layouts, level, _ = get_layouts(tree, num_props, 'barplot', level, internal_num_rep, prop2type=prop2type)
        barplot_layouts, level, _ = get_barplot_layouts(tree, num_props, level, prop2type, column_width=args.barplot_width, internal_rep=internal_num_rep)   
        layouts.extend(barplot_layouts)
        
        multiple_text_props = [
            'eggNOG_OGs', #28PAR@1|root,2QVY3@2759|Eukaryota
            'GOs', #GO:0000002,GO:0000003
            'KEGG_ko', #ko:K04451,ko:K10148
            'KEGG_Pathway', #ko01522,ko01524
            'KEGG_Module', #M00118 
            'KEGG_Reaction', #R00497
            'KEGG_rclass', #RC00141
            'EC', #1.18.6.1,1.3.7.14,1.3.7.15
            'BRITE', #ko00000,ko00001,ko03000
            'KEGG_TC', #3.A.1.133.1 

            # Domains
            'CAZy',
            'BiGG_Reaction',
            'PFAMs'
        ]
        
        for multiple_text_prop in multiple_text_props:
            matrix, all_values = multiple2profile(tree, multiple_text_prop)
            multiple_text_prop_layout = profile_layouts.LayoutProfile(name="Profiling_"+multiple_text_prop, 
            mode='multi', alignment=matrix, profiles=all_values, column=level)
            level += 1
            layouts.append(multiple_text_prop_layout)
            
    # Taxa layouts
    if args.taxonclade_layout or args.taxonrectangle_layout:
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
            
            if args.taxonclade_layout:
                taxa_layout = taxon_layouts.TaxaClade(name='TaxaClade_'+rank, level=level, rank = rank, color_dict=color_dict)
                taxa_layouts.append(taxa_layout)

            if args.taxonrectangle_layout:
                taxa_layout = taxon_layouts.TaxaRectangular(name = "TaxaRect_"+rank, rank=rank ,color_dict=color_dict, column=level)
                taxa_layouts.append(taxa_layout)
                #level += 1
            taxon_color_dict[rank] = color_dict
            
        #taxa_layouts.append(taxon_layouts.TaxaRectangular(name = "Last Common Ancester", color_dict=taxon_color_dict, column=level))
        taxa_layouts.append(taxon_layouts.LayoutSciName(name = 'Taxa Scientific name', color_dict=taxon_color_dict))
        taxa_layouts.append(taxon_layouts.LayoutEvolEvents(name='Taxa Evolutionary events', prop="evoltype",
            speciation_color="blue", 
            duplication_color="red", node_size = 2,
            legend=True))
        layouts = layouts + taxa_layouts
        level += 1
        total_color_dict.append(taxon_color_dict)
        
    #### prune at the last step in case of losing leaves information
    # prune tree by rank
    if args.rank_limit:
        tree = taxatree_prune(tree, rank_limit=args.rank_limit)
        
    # prune tree by condition 
    if args.pruned_by: # need to be wrap with quotes
        condition_strings = args.pruned_by
        tree = conditional_prune(tree, condition_strings, prop2type)

    #### Output #####
    if args.out_colordict:
        wrtie_color(total_color_dict)
    if args.plot:
        get_image(tree, layouts, args.port, os.path.abspath(args.plot))
    else:  
        tree.explore(daemon=False, compress=False, quiet=args.verbose, layouts=layouts, port=args.port, include_props=sorted(popup_prop_keys))
    
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

def get_label_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        color_dict = {} # key = value, value = color id
        if prop2type and prop2type.get(prop) == list:
            leaf_values = list(map(list,set(map(tuple,children_prop_array(tree, prop)))))
            
            prop_values = [val for sublist in leaf_values for val in sublist]
        else:
            prop_values = sorted(list(set(children_prop_array(tree, prop))))

        # normal text prop
        nvals = len(prop_values)            
        for i in range(0, nvals):
            if nvals <= 14:
                color_dict[prop_values[i]] = paried_color[i]
            else:
                color_dict[prop_values[i]] = random_color(h=None)

        layout = text_layouts.LayoutText(name='Label_'+prop, column=level, 
        color_dict=color_dict, text_prop=prop, width=column_width, padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_colorbranch_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        color_dict = {} # key = value, value = color id
        if prop2type and prop2type.get(prop) == list:
            leaf_values = list(map(list,set(map(tuple,children_prop_array(tree, prop)))))    
            prop_values = [val for sublist in leaf_values for val in sublist]
        else:
            prop_values = sorted(list(set(children_prop_array(tree, prop))))
        
        # normal text prop
        nvals = len(prop_values)            
        for i in range(0, nvals):
            if nvals <= 14:
                color_dict[prop_values[i]] = paried_color[i]
            else:
                color_dict[prop_values[i]] = random_color(h=None)
        layout = text_layouts.LayoutColorbranch(name='Colorbranch_'+prop, column=level, \
            color_dict=color_dict, text_prop=prop, width=column_width, \
                padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_rectangle_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        color_dict = {} # key = value, value = color id
        if prop2type and prop2type.get(prop) == list:
            leaf_values = list(map(list,set(map(tuple,children_prop_array(tree, prop)))))    
            prop_values = [val for sublist in leaf_values for val in sublist]
        else:
            prop_values = sorted(list(set(children_prop_array(tree, prop))))
        
        # normal text prop
        nvals = len(prop_values)            
        for i in range(0, nvals):
            if nvals <= 14:
                color_dict[prop_values[i]] = paried_color[i]
            else:
                color_dict[prop_values[i]] = random_color(h=None)
        layout = text_layouts.LayoutRect(name='Rectangular_'+prop, column=level,
                    color_dict=color_dict, text_prop=prop,
                    width=column_width, padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_binary_layouts(tree, props, level, prop2type, column_width=70, reverse=False, padding_x=1, padding_y=0):
    prop_color_dict = {}
    layouts = []

    for prop in props:
        color_dict = {} # key = value, value = color id
        prop_values = sorted(list(set(children_prop_array(tree, prop))))
        nvals = len(prop_values)

        for i in range(0, nvals): # only positive, negative, NaN, three options
            color_dict[prop_values[i]] = paried_color[i]
        
        color = random_color(h=None)
        if not reverse:
            layout = conditional_layouts.LayoutBinary('Binary_'+prop, level, color, color_dict, prop, width=column_width, padding_x=padding_x, padding_y=padding_y, reverse=reverse)
        else:
            layout = conditional_layouts.LayoutBinary('ReverseBinary_'+prop, level, color, color_dict, prop, width=column_width, padding_x=padding_x, padding_y=0, reverse=reverse)
        
        internal_prop = prop + '_' + 'counter'
        prop_color_dict[internal_prop] = color_dict
        prop_color_dict[prop] = color
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_barplot_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0, internal_rep='avg'):
    prop_color_dict = {}
    layouts = []
    barplot_padding_x = padding_x * 10 
    for prop in props:
        
        color_dict = {} # key = value, value = color id
        prop_values = children_prop_array(tree, prop)
        if prop_values:
            size_prop = prop
        else:
            size_prop = prop+'_'+internal_rep

        if level > len(paried_color):
            barplot_color =  random_color(h=None)
        else:
            barplot_color = paried_color[level]
        
        layout =  staple_layouts.LayoutBarplot(name='Barplot_'+prop, prop=prop, 
                                    width=column_width, color=barplot_color, 
                                    size_prop=size_prop, column=level, 
                                    internal_rep=internal_rep, padding_x=barplot_padding_x
                                    )

        prop_color_dict[prop] = barplot_color
        layouts.append(layout)  
        level += 1

    return layouts, level, prop_color_dict

def get_heatmap_layouts(tree, props, level, column_width=70, padding_x=1, padding_y=0, internal_rep='avg'):
    layouts = []
    for prop in props:
        prop_values = np.array(list(set(children_prop_array(tree, prop)))).astype('float64')
        prop_values = prop_values[~np.isnan(prop_values)]
        minval, maxval = prop_values.min(), prop_values.max()
        layout =  staple_layouts.LayoutHeatmap(name='Heatmap_'+prop, column=level, 
                    width=column_width, padding_x=padding_x, padding_y=padding_y,internal_rep=internal_rep, 
                    prop=prop, maxval=maxval, minval=minval)
        layouts.append(layout)  
        level += 1

    return layouts, level

def get_prop2type(node):
    output = {}
    prop2value = node.props

    if '_speciesFunction' in prop2value:
        del prop2value['_speciesFunction']
    
    for prop, value in prop2value.items():
        if isinstance(value, numbers.Number):
            output[prop] = float
        elif type(value) == list:
            output[prop] = list
        else:
            output[prop] = str

    # for prop, value in prop2value.items():
    #     output[prop] = type(value)
    
    return output

def props2matrix(tree, profiling_props, dtype=float):
    aa = [
        'A', 'R', 'N',
        'D', 'C', 'Q',
        'E', 'H',
        'I', 'S', 'K',
        'M', 'F', 'P',
        'L', 'T', 'W',
        'Z', 'V', 'B',
        'Y', 'X'
    ]
    absence_color = 'G'
    gradients = [
    'a', 'b', 'c',
    'd', 'e', 'f',
    'g', 'h', 'i',
    'j', 'k', 'l',
    'm', 'n', 'o',
    'p', 'q', 'r', 
    's', 't'
    ] #blue to red

    leaf2matrix = {}
    for node in tree.traverse():
        if node.is_leaf:
            leaf2matrix[node.name] = []
            for profiling_prop in profiling_props:
                if node.props.get(profiling_prop):
                    if dtype == float:
                        val = float(node.props.get(profiling_prop))
                    elif dtype == str:
                        val = node.props.get(profiling_prop)
                    leaf2matrix[node.name].append(val)
                else:
                    leaf2matrix[node.name].append(None)
    
    # gain all values from metadata
    if dtype == float:
        all_values = list(set(flatten([sublist for sublist in leaf2matrix.values()])))
        all_values = sorted(list(filter(lambda x: x is not None and not math.isnan(x), all_values)))
        maxval = max(all_values)
        minval = min(all_values)
        num = len(gradients)
        values = np.linspace(minval, maxval, num)

        matrix = ''
        for leaf, prop in leaf2matrix.items():
            matrix += '\n'+'>'+leaf+'\n'
            for i in range(len(prop)):
                search_value = prop[i]
                if search_value:
                    # Find the index of the closest element to the search value
                    index = np.abs(values - search_value).argmin()
                    matrix += gradients[index]
                else:
                    matrix += '-'
        return matrix, maxval, minval
    
    elif dtype == str:
        value2color = {}
        all_values = sorted(list(set(flatten([sublist for sublist in leaf2matrix.values()]))))
        for i in range(len(all_values)):
            val = all_values[i]
            if val != 'NaN':
                value2color[val] = aa[i]
            else:
                value2color[val] = absence_color
        
        matrix = ''
        for leaf, prop in leaf2matrix.items():
            matrix += '\n'+'>'+leaf+'\n'
            for item in prop:
                matrix += value2color[item]
        
        return matrix, value2color
  
# def categorical2profile(tree, profiling_prop):
#     aa = [
#         'A', 'R', 'N',
#         'D', 'C', 'Q',
#         'E', 'H',
#         'I', 'S', 'K',
#         'M', 'F', 'P',
#         'L', 'T', 'W',
#         'Z', 'V', 'B',
#         'Y', 'X'
#     ]
#     absence_color = 'G'

#     leaf2matrix = {}
#     for node in tree.traverse():
#         if node.is_leaf:
#             leaf2matrix[node.name] = []
#             #for profiling_prop in profiling_props:
#             if node.props.get(profiling_prop):
#                 val = node.props.get(profiling_prop)
#                 leaf2matrix[node.name].append(val)
#             else:
#                 leaf2matrix[node.name].append(None)

#     value2color = {}
#     all_values = list(set(flatten([sublist for sublist in leaf2matrix.values()])))
#     for i in range(len(all_values)):
#         val = all_values[i]
#         if val != 'NaN':
#             value2color[val] = aa[i]
#         else:
#             value2color[val] = absence_color
    
#     matrix = ''
#     for leaf, prop in leaf2matrix.items():
#         matrix += '\n'+'>'+leaf+'\n'
#         for item in prop:
#             matrix += value2color[item]
#     return matrix, value2color

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

def single2profile(tree, profiling_prop):
    all_values = sorted(list(set(flatten(children_prop_array(tree, profiling_prop)))), key=lambda x: (x != 'NaN', x))
    presence = 'a' # #E60A0A red
    absence = '-' # #EBEBEB lightgrey
    matrix = ''
    for leaf in tree.leaves():
        matrix += '\n'+'>'+leaf.name+'\n'
        if leaf.props.get(profiling_prop):
            for val in all_values:
                if val == leaf.props.get(profiling_prop):
                    matrix += presence
                else:
                    matrix += absence
        else:
            matrix += absence * len(all_values) +'\n'
    return matrix, all_values
    
def multiple2profile(tree, profiling_prop):
    all_values = sorted(list(set(flatten(children_prop_array(tree, profiling_prop)))), key=lambda x: (x != 'NaN', x))
    presence = 'a' # #E60A0A red
    absence = '-' # #EBEBEB lightgrey
    matrix = ''
    for leaf in tree.leaves():
        matrix += '\n'+'>'+leaf.name+'\n'
        if leaf.props.get(profiling_prop):
            for val in all_values:
                if val in leaf.props.get(profiling_prop):
                    matrix += presence
                else:
                    matrix += absence
        else:
            matrix += absence * len(all_values) +'\n'
    return matrix, all_values
      
