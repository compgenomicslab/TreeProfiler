#!/usr/bin/env python
from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from layouts import (
    text_layouts, taxon_layouts, staple_layouts, 
    conditional_layouts, seq_layouts, profile_layouts)
from utils import (
    ete4_parse, taxatree_prune, conditional_prune,
    children_prop_array, children_prop_array_missing, 
    flatten, get_consensus_seq)
import b64pickle
from tree_image import get_image

from collections import defaultdict
from itertools import islice
from io import StringIO
import numpy as np
import numbers
import math
import colorsys
import random
import sys
import os

paried_color = ["red", "darkblue", "lightgreen", "sienna", "lightCoral", "violet", "mediumturquoise",   "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal", "darkyellow"]

DESC = "plot tree"

def poplulate_plot_args(plot_args_p):
    """
    Parse the input parameters
    Return the parsed arguments.
    """
    group = plot_args_p.add_argument_group(title='Conditional display arguments',
        description="Conditional display  parameters")
    
    group.add_argument('--internal_plot_measure',
        default='avg',
        type=str,
        required=False,
        help="statistic measures to be shown in numerical layout for internal nodes, [default: avg]")  

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
    group.add_argument('--column_width',
        type=int,
        default=70,
        help="customize column width of each layout."
    )
    group.add_argument('--barplot_width',
        type=int,
        default=200,
        help="customize barplot width of barplot layout."
    )
    group.add_argument('--profiling_width',
        type=int,
        default=1000,
        help="customize profiling width of each profiling layout."
    )
    group.add_argument('--binary_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as binary_layout")
    group.add_argument('--revbinary_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as revbinary_layout")
    group.add_argument('--colorbranch_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as Textlayouts")
    group.add_argument('--label_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as label_layout")
    group.add_argument('--rectangular_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be plot as rectangular_layout")
    group.add_argument('--heatmap_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as heatmap_layout")
    group.add_argument('--barplot_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as barplot_layouts")
    group.add_argument('--taxonclade_layout',
        default=False,
        action='store_true',
        help="activate taxonclade_layout")
    group.add_argument('--taxonrectangular_layout',
        default=False,
        action='store_true',
        help="activate taxonrectangular_layout")
    group.add_argument('--emapper_layout',
        default=False,
        action='store_true',
        help="activate emapper_layout") #domain_layout
    group.add_argument('--domain_layout',
        default=False,
        action='store_true',
        help="activate domain_layout") #domain_layout
    group.add_argument('--alignment_layout',
        default=False,
        action='store_true',
        help="provide alignment file as fasta format")
    group.add_argument('--profiling_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index which need to be plot as profiling_layout for categorical columns")
    group.add_argument('--multi_profiling_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index which need to be plot as multi_profiling_layout for multiple values column")
    group.add_argument('--numerical_profiling_layout',
        type=lambda s: [item for item in s.split(',')],
        required=False,
        help="<col1,col2> names, column index which need to be plot as numerical_profiling_layout for numerical values column")

    group = plot_args_p.add_argument_group(title='Output arguments',
        description="Output parameters")
    # group.add_argument('--interactive',
    #     default=False,
    #     action='store_true',
    #     help="run interactive session")
    group.add_argument('--port',
        type=str,
        default=5000,
        help="run interactive session on custom port")
    group.add_argument('--plot',
        type=str,
        required=False,
        help="output as pdf")
    group.add_argument('--out_colordict',
        action="store_true", 
        required=False,
        help="print color dictionary of each property")


### visualize tree
def run(args):
    global prop2type, columns, tree
    node_props=[]
    columns = {}
    rank2values = {}

    total_color_dict = []
    layouts = []
    level = 1 # level 1 is the leaf name

    #parse input tree
    if args.tree:
        if args.tree_type == 'newick':
            tree = ete4_parse(args.tree)
        elif args.tree_type == 'ete':
            with open(args.tree, 'r') as f:
                file_content = f.read()
                tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
    
    #rest_prop = []
    if args.prop2type:
        prop2type = {}
        with open(args.prop2type, 'r') as f:
            for line in f:
                line = line.rstrip()
                prop, value = line.split('\t')
                prop2type[prop] = value
        
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
                'named_lineage': str
                }
        popup_prop_keys = list(prop2type.keys()) 

        if args.tree_type == 'ete':
            leafa, _, leafb, _ = tree._get_farthest_and_closest_leaves()
            leaf_prop2type = get_prop2type(leafa)
            leaf_prop2type.update(get_prop2type(leafb))
            internal_node_prop2type = get_prop2type(tree)
            prop2type.update(leaf_prop2type)
            prop2type.update(internal_node_prop2type)
            
        # elif args.tree_type == 'newick':
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
    # if args.num_stat != 'all':
    #     internal_num_rep = args.num_stat
    # else:
    #     internal_num_rep = args.internal_plot_measure

    internal_num_rep = args.internal_plot_measure

    # Get the input arguments in order
    input_order = []
    for arg in sys.argv[1:]:
        if arg.startswith('-') and arg.endswith('layout'):
            input_order.append(arg[2:])
        else:
            continue
    
    for layout in input_order:
        if layout == 'heatmap_layout':
            column_width = args.column_width
            props = []
            for i in args.heatmap_layout:
                props.append(i)

            heatmap_layouts = []
            for prop in props:
                prop_values = np.array(list(set(children_prop_array(tree, prop)))).astype('float64')
                prop_values = prop_values[~np.isnan(prop_values)]
                minval, maxval = prop_values.min(), prop_values.max()
                layout =  staple_layouts.LayoutHeatmap(name='Heatmap_'+prop, column=level, 
                                    width=column_width, internal_rep=internal_num_rep, 
                                    prop=prop, maxval=maxval, minval=minval)
                heatmap_layouts.append(layout)
                level += 1

            layouts.extend(heatmap_layouts)

        if layout == 'label_layout':
            label_layouts, level, color_dict = get_label_layouts(args.label_layout, level, prop2type=prop2type, column_width=args.column_width)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)

        if layout == 'colorbranch_layout':
            colorbranch_layouts, level, color_dict = get_colorbranch_layouts(args.colorbranch_layout, level, prop2type=prop2type, column_width=args.column_width)
            layouts.extend(colorbranch_layouts)
            total_color_dict.append(color_dict)

        if layout == 'rectangular_layout':
            rectangular_layouts, level, color_dict = get_rectangular_layouts(args.rectangular_layout, level, prop2type=prop2type, column_width=args.column_width)
            layouts.extend(rectangular_layouts)
            total_color_dict.append(color_dict)

        if layout == 'binary_layout':
            label_layouts, level, color_dict = get_binary_layouts(args.binary_layout, level, prop2type=prop2type, column_width=args.column_width, reverse=False)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)

        if layout == 'revbinary_layout':
            label_layouts, level, color_dict = get_binary_layouts(args.revbinary_layout, level, prop2type=prop2type, column_width=args.column_width, reverse=True)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)
        
        if layout == 'barplot_layout':
            barplot_layouts, level,color_dict = get_barplot_layouts(args.barplot_layout, level, prop2type, column_width=args.barplot_width, internal_rep=internal_num_rep)
            layouts.extend(barplot_layouts)
            total_color_dict.append(color_dict)

        if layout == 'alignment_layout':
            #fasta_file = args.alignment_layout
            lengh = len(max(children_prop_array(tree, 'alignment'),key=len))
            aln_layout = seq_layouts.LayoutAlignment(name='Alignment_layout', 
                        alignment_prop='alignment', column=level, scale_range=lengh,
                        summarize_inner_nodes=False)
            layouts.append(aln_layout)

        if layout == 'domain_layout':
            domain_layout = seq_layouts.LayoutDomain(name="Domain_layout", prop='dom_arq')
            layouts.append(domain_layout)

        if layout == 'profiling_layout':
            profiling_props = args.profiling_layout
            matrix = props2matrix(tree, profiling_props, dtype=str)
            profile_layout = profile_layouts.LayoutProfile(name='profiling_layout', mode='simple',
                alignment=matrix, profiles=profiling_props, column=level, width=args.profiling_width)
            level += 1
            layouts.append(profile_layout)

        if layout == 'multi_profiling_layout':
            profiling_props = args.multi_profiling_layout
            for profiling_prop in profiling_props:
                matrix, all_values = multiple2profile(tree, profiling_prop)
                profile_layout = profile_layouts.LayoutProfile(name=profiling_prop, mode='multi',
                alignment=matrix, profiles=all_values, column=level, summarize_inner_nodes=False, width=args.profiling_width)
                level += 1
                layouts.append(profile_layout)
        
        if layout == 'numerical_profiling_layout':
            profiling_props = args.numerical_profiling_layout
            matrix, maxval, minval = props2matrix(tree, profiling_props)
            #profile_layout = TreeLayout(name='numerical_profiling_layout', ns=get_alnface(alignment, level), aligned_faces = True)
            profile_layout = profile_layouts.LayoutProfile(name='numerical_profiling_layout', mode='numerical', 
                alignment=matrix, seq_format='gradients', profiles=profiling_props, value_range=[minval, maxval], column=level, width=args.profiling_width)
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
            'EC',
            ]
        #label_layouts, level, _ = get_layouts(text_props, 'rectangular', level, 'counter', prop2type=prop2type)
        label_layouts, level, _ = get_rectangular_layouts(text_props, level, prop2type=prop2type, column_width=args.column_width)
            
        layouts.extend(label_layouts)
        
        num_props = [
            #'evalue',
            'score'
        ]
        #barplot_layouts, level, _ = get_layouts(num_props, 'barplot', level, internal_num_rep, prop2type=prop2type)
        barplot_layouts, level, _ = get_barplot_layouts(num_props, level, prop2type, column_width=args.barplot_width, internal_rep=internal_num_rep)   
        layouts.extend(barplot_layouts)
        
        multiple_text_props = [
            'eggNOG_OGs', #28PAR@1|root,2QVY3@2759|Eukaryota
            'GOs', #GO:0000002,GO:0000003
            'KEGG_ko', #ko:K04451,ko:K10148
            'KEGG_Pathway', #ko01522,ko01524
            'KEGG_Module', #M00118 
            'KEGG_Reaction', #R00497
            'KEGG_rclass', #RC00141

            # cannot use kegg_get()
            'BRITE', #ko00000,ko00001,ko03000
            'KEGG_TC', #3.A.1.133.1 

            # Domains
            'CAZy',
            'BiGG_Reaction',
            'PFAMs'
        ]
        
        for multiple_text_prop in multiple_text_props:
            matrix, all_values = multiple2profile(tree, multiple_text_prop)
            multiple_text_prop_layout = profile_layouts.LayoutProfile(name=multiple_text_prop, 
            alignment=matrix, profiles=all_values, column=level)
            level += 1
            layouts.append(multiple_text_prop_layout)
            
    # Taxa layouts
    if args.taxonclade_layout or args.taxonrectangular_layout:
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
        for rank, value in sorted(rank2values.items(),):
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

            if args.taxonrectangular_layout:
                taxa_layout = taxon_layouts.TaxaRectangular(name = "TaxaRect_"+rank, rank=rank ,color_dict=color_dict, column=level)
                taxa_layouts.append(taxa_layout)
                #level += 1
            taxon_color_dict[rank] = color_dict
            
        #taxa_layouts.append(taxon_layouts.TaxaRectangular(name = "Last Common Ancester", color_dict=taxon_color_dict, column=level))
        taxa_layouts.append(taxon_layouts.LayoutSciName(name = 'Taxa Scientific name', color_dict=taxon_color_dict))
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
        tree.explore(tree_name='example',layouts=layouts, port=args.port, include_props=sorted(popup_prop_keys))
    
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

def get_label_layouts(props, level, prop2type, column_width=70):
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
        color_dict=color_dict, text_prop=prop, width=column_width)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_colorbranch_layouts(props, level, prop2type, column_width=70):
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
            color_dict=color_dict, text_prop=prop, width=column_width)
        layouts.append(layout)
    return layouts, level, prop_color_dict

def get_rectangular_layouts(props, level, prop2type, column_width=70):
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
                    width=column_width)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_binary_layouts(props, level, prop2type, column_width=70, reverse=False):
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
            layout = conditional_layouts.LayoutBinary('Binary_'+prop, level, color, color_dict, prop, reverse=reverse)
        else:
            layout = conditional_layouts.LayoutBinary('ReverseBinary_'+prop, level, color, color_dict, prop, reverse=reverse)
        
        internal_prop = prop + '_' + 'counter'
        prop_color_dict[internal_prop] = color_dict
        prop_color_dict[prop] = color
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_barplot_layouts(props, level, prop2type, column_width=70, internal_rep='avg'):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        
        color_dict = {} # key = value, value = color id
        if prop in prop2type and prop2type.get(prop) == float:
            size_prop = prop+'_'+internal_rep # using internal prop to set the range in case rank_limit cut all the leaves
        else:
            size_prop = prop

        if level > len(paried_color):
            barplot_color =  random_color(h=None)
        else:
            barplot_color = paried_color[level]
        
        layout =  staple_layouts.LayoutBarplot(name='Barplot_'+prop, prop=prop, 
                                    width=column_width, color=barplot_color, 
                                    size_prop=size_prop, column=level, 
                                    internal_rep=internal_rep,
                                    )

        
        prop_color_dict[prop] = barplot_color
        layouts.append(layout)  
        level += 1

    return layouts, level, prop_color_dict

def get_heatmap_layouts():
    return

def get_layouts(argv_inputs, layout_name, level, internal_rep, prop2type=None, column_width=70): 
    props = argv_inputs
    layouts = []
    prop_color_dict = {} # key = value, value = color id

    # load layout for each prop
    for idx, prop in enumerate(props):
        
        color_dict = {} # key = value, value = color id

        # binary layout should be set width
        if layout_name in ['binary', 'revbinary']:
            
            if columns:
                prop_values = sorted(list(set(columns[prop])))
            else:
                prop_values = sorted(list(set(children_prop_array(tree, prop))))
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
        # elif layout_name == 'heatmap':
        #     layout =  staple_layouts.LayoutHeatmap('Heatmap_'+prop, level, internal_rep, prop)
       
        elif layout_name == 'barplot':
            if prop in prop2type and prop2type.get(prop) == float:
                size_prop = prop+'_'+internal_rep # using internal prop to set the range in case rank_limit cut all the leaves
            else:
                size_prop = prop

            if level > len(paried_color):
                barplot_color =  random_color(h=None)
            else:
                barplot_color = paried_color[level]
            
            layout =  staple_layouts.LayoutBarplot(name='Barplot_'+prop, prop=prop, 
                                        width=column_width, color=barplot_color, 
                                        size_prop=size_prop, column=level, 
                                        internal_rep=internal_rep,
                                        )

            prop_color_dict[prop] = barplot_color

        # categorical layouts should be set width
        elif layout_name in ['label','rectangular', 'colorbranch']:
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
            if layout_name == 'label':
                #longest_val =  len(max(prop_values, key = len))
                layout = text_layouts.LayoutText(name='Label_'+prop, column=level, \
                    color_dict=color_dict, text_prop=prop, width=column_width)
                #layout = TreeLayout(name=prop+'_'+layout_name, ns=text_layouts.text_layout(prop, level, color_dict, internal_rep))
            
            elif layout_name == 'rectangular':
                layout = text_layouts.LayoutRect(name='Rectangular_'+prop, column=level,  \
                    color_dict=color_dict, text_prop=prop,\
                    width=column_width)
            
            elif layout_name == 'colorbranch':
                layout = text_layouts.LayoutColorbranch(name='Colorbranch_'+prop, column=level, \
                    color_dict=color_dict, text_prop=prop, width=column_width)
            
            prop_color_dict[prop] = color_dict
        
        layouts.append(layout)
        level += 1
        
    return layouts, level, prop_color_dict

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
        if node.is_leaf():
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
        all_values = list(filter(lambda x: x is not None and not math.isnan(x), all_values))
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
        all_values = list(set(flatten([sublist for sublist in leaf2matrix.values()])))
        for i in range(len(all_values)):
            val = all_values[i]
            if val != 'NaN':
                value2color[val] = aa[i]
            else:
                value2color[val] = 'G'
        
        matrix = ''
        for leaf, prop in leaf2matrix.items():
            matrix += '\n'+'>'+leaf+'\n'
            for item in prop:
                matrix += value2color[item]
        return matrix
  
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

def multiple2profile(tree, profiling_prop):
    all_values = sorted(list(set(flatten(children_prop_array(tree, profiling_prop)))))
    presence = 'D' # #E60A0A red
    absence = 'G' # #EBEBEB lightgrey
    matrix = ''
    for leaf in tree.iter_leaves():
        matrix += '\n'+'>'+leaf.name+'\n'
        if leaf.props.get(profiling_prop):
            for val in all_values:
                if val != 'NaN' and val in leaf.props.get(profiling_prop):
                    matrix += presence
                else:
                    matrix += absence
            # for index in range(len(all_values)):
            #     val = all_values[index]
            #     if val != 'NaN' and val in leaf.props.get(profiling_prop):
            #         matrix += aa[index % len(aa)]
            #     else:
            #         matrix += '-'
        else:
            matrix += absence * len(all_values) +'\n'
    return matrix, all_values

def categorical2profile(tree, profiling_prop):
    all_values = sorted(list(set(flatten(children_prop_array(tree, profiling_prop)))))
    aa = [
        'A', 'R', 'N',
        'D', 'C', 'Q',
        'E', 'G', 'H',
        'I', 'S', 'K',
        'M', 'F', 'P',
        'L', 'T', 'W',
        'Z', 'V', 'B',
        'Y', 'X'
    ]
    matrix = ''
    for leaf in tree.iter_leaves():
        matrix += '\n'+'>'+leaf.name+'\n'
        for i in range(len(all_values)):
            leaf_prop = leaf.props.get(profiling_prop)
            if leaf_prop and leaf_prop != 'NaN':
                matrix += aa[i % len(aa)]
            else:
                matrix += 'G'
    return matrix              
