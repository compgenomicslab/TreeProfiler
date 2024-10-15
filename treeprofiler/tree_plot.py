#!/usr/bin/env python
import numbers
import math
import sys
import os
import argparse
import csv
import logging

from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from itertools import islice
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np

from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from treeprofiler.tree_image import get_image
from treeprofiler.layouts import (
    text_layouts, taxon_layouts, staple_layouts, 
    conditional_layouts, seq_layouts, profile_layouts, phylosignal_layouts)

import treeprofiler.src.utils as utils
from treeprofiler.tree_annotate import can_convert_to_bool

paired_color = [
    '#9a312f', '#9b57d0', '#f8ce9a', '#f16017', '#28fef9', '#53707a',
    '#213b07', '#b5e5ac', '#9640b2', '#a9bd10', '#69e42b', '#b44d67',
    '#b110c1', '#0b08a3', '#d07671', '#29e23b', '#3f2bf4', '#9b2a08',
    '#b42b94', '#77566a', '#2dfee7', '#046904', '#e2835d', '#53db2b',
    '#0b97e9', '#e0f6e9', '#ba46d1', '#4aba53', '#d4d6db', '#7a5d7c',
    '#4b100e', '#9e6373', '#5f4945', '#7e057a', '#f8e372', '#209f87',
    '#383f59', '#9d59e9', '#40c9fb', '#4cfc8b', '#d94769', '#20feba',
    '#c53238', '#068b02', '#6b4c93', '#f1968e', '#86d720', '#076fa6',
    '#0dbcfe', '#4d74b2', '#7b3dd2', '#286d26', '#a0faca', '#97505d',
    '#159e7a', '#fc05df', '#5df454', '#9160e1', '#c2eb5e', '#304fce',
    '#033379', '#54770f', '#271211', '#ab8479', '#37d9a0', '#f12205',
    '#cdd7e2', '#578f56', '#5ad9be', '#8596e9', '#c999ee', '#5f6b8a',
    '#f5c3a1', '#8e0603', '#cc21cf', '#65e7d0', '#97b3b6', '#d6220c',
    '#29c1e1', '#a30139', '#c9a619', '#a19410', '#da874f', '#64246d',
    '#66f35d', '#b8366c', '#116c95', '#bd851a', '#27f7cb', '#512ca4',
    '#60e72e', '#d1941c', '#1045a8', '#c1b03a', '#0c62a5', '#7ac9b2',
    '#6bb9bd', '#cb30eb', '#26bad0', '#d9e557'
]

DESC = "plot tree"

# Set up the logger with INFO level by default
logger = logging.getLogger(__name__)

def setup_logger(level=logging.ERROR):
    """Sets up logging configuration."""
    logger.setLevel(level)
    
    # Create a StreamHandler to output to sys.stdout
    handler = logging.StreamHandler(sys.stdout)
    
    # Define the formatter with a custom format for errors
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s - [in %(filename)s:%(lineno)d]')
    
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def string_or_file(value):
    if os.path.isfile(value):
        if os.path.exists(value):
            return value
    else:
        return value  

def poplulate_plot_args(plot_args_p):
    """
    Parse the input parameters
    Return the parsed arguments.
    """
    group = plot_args_p.add_argument_group(title='Conditional display arguments',
        description="Conditional display  parameters")
    group.add_argument('--internal-plot-measure',
        default='avg',
        choices=['sum', 'avg', 'max', 'min', 'std', 'none'],
        type=str,
        required=False,
        help="statistic measures to be shown in numerical layout for internal nodes, [default: avg]")  

    group.add_argument('--collapsed-by', 
        type=string_or_file,
        required=False,
        action='append',
        help='target tree nodes collapsed by customized conditions')  
    group.add_argument('--highlighted-by', 
        type=string_or_file,
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
        type=float,
        default=200,
        help="customize barplot width of barplot layout.[default: 200]"
    )
    group.add_argument('--barplot-scale',
        type=str,
        default=None,
        help="find the barplot column as scale anchor.[default: None]"
    )
    group.add_argument('--barplot-colorby',
        type=str,
        default=None,
        help="Set the color of barplot by the a categorical property.[default: None]"
    )
    group.add_argument('--color-config',
        type=argparse.FileType('r'),
        default=None,
        help="Path to the file to find the color for each variables. [default: None]"
    )
    group.add_argument('-s', '--config-sep', default='\t',
        help="column separator of color table [default: \\t]")
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
    group.add_argument('--acr-discrete-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties which need to be plot as acr-discrete-layout")
    group.add_argument('--acr-continuous-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties which need to be plot as acr-continuous-layout")
    group.add_argument('--ls-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties which need to be plot as ls-layout")
    group.add_argument('--binary-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties which need to be plot as binary-layout which highlights the postives")
    group.add_argument('--binary-aggregate-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties which need to be plot as binary-aggregate-layout which highlights the postives")
    group.add_argument('--binary-unicolor-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties which need to be plot as binary-layout which highlights the postives")
    group.add_argument('--binary-unicolor-aggregate-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties which need to be plot as binary-aggregate-layout which highlights the postives")
    # group.add_argument('--revbinary-layout',
    #     nargs='+',
    #     required=False,
    #     help="<prop1> <prop2> names of properties which need to be plot as revbinary-layout which highlights the negatives")
    # group.add_argument('--revbinary-unicolor-layout',
    #     nargs='+',
    #     required=False,
    #     help="<prop1> <prop2> names of properties which need to be plot as revbinary-layout which highlights the negatives")
    group.add_argument('--colorbranch-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties where branches will be colored based on different values.")
    group.add_argument('--label-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties where values will be displayed on the aligned panel.")
    group.add_argument('--rectangle-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names  of properties where values will be label as rectangular color block on the aligned panel.")
    group.add_argument('--bubble-layout',
        nargs='+',
        required=False,
        help="(experimental) <prop1> <prop2> names of properties which need to be plot as bubble-layout")
    group.add_argument('--background-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties where values will be label as rectangular color block on the aligned panel.")
    group.add_argument('--piechart-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of properties whose internal nodes need to be plot as piechart-layout.")
    group.add_argument('--heatmap-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of numerical properties which need to be read as heatmap-layout")
    group.add_argument('--heatmap-mean-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of numerical properties which need to be read as heatmap-layout")
    group.add_argument('--heatmap-zscore-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of numerical properties which need to be read as heatmap-layout")
    group.add_argument('--barplot-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names of numerical properties which need to be read as barplot_layouts")
    # group.add_argument('--branchscore-layout',
    #     nargs='+',
    #     required=False,
    #     help="<prop1> <prop2> names of numerical properties which need to be read as branchscore_layouts")  
    group.add_argument('--taxonclade-layout',
        default=False,
        action='store_true',
        help="Activate taxonclade_layout which clades will be colored based on taxonomy of each node.")
    group.add_argument('--taxonrectangle-layout',
        default=False,
        action='store_true',
        help="Activate taxonrectangle-layout which taxonomy of each node will be display as rectangular blocks in aligned panel.")
    group.add_argument('--taxoncollapse-layout',
        default=False,
        action='store_true',
        help="Activate taxoncollapse-layout which taxonomy of each node will be display as rectangular blocks in aligned panel.")
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
        help="<prop1> <prop2> names of properties which need to be convert to presence-absence profiling matrix of each value")
    # group.add_argument('--multi-profiling-layout',
    #     nargs='+',
    #     required=False,
    #     help="<prop1> <prop2> names of properties containing values as list which need to be convert to presence-absence profiling matrix")
    group.add_argument('--categorical-matrix-layout',
        nargs='+',
        required=False,
        help="<prop1> <prop2> names which need to be plot as categorical_matrix_layout for categorical values")
    group.add_argument('--numerical-matrix-layout',
        nargs='+',
        required=False,
        help="numerical matrix that take into account ALL values into gradient from white to red. <prop1> <prop2> names which need to be plot as numerical_matrix_layout for numerical values ")
    group.add_argument('--binary-matrix-layout',
        nargs='+',
        required=False,
        help="names of properties which need to be plot as binary-matrix which highlights the postives")

    group = plot_args_p.add_argument_group(title='Visualizing output arguments',
        description="Visualizing output parameters")
    # group.add_argument('--interactive',
    #     default=False,
    #     action='store_true',
    #     help="run interactive session")

    group.add_argument('--hide-leaf-name', action='store_false',
        help='Hide the leaf names in the tree view.')
    group.add_argument('--hide-branch-support', action='store_false',
        help='Hide the branch support values in the tree view.')
    group.add_argument('--hide-branch-distance', action='store_false',
        help='Hide the branch distances in the tree view.')

    group.add_argument('--verbose',
        action="store_false", 
        required=False,
        help="show detail on prompt when visualizing taget tree.")
    group.add_argument('--port',
        type=str,
        default=5000,
        help="run interactive session on custom port.[default: 5000]")
    group.add_argument('--render',
        #type=str,
        action="store_true", 
        required=False,
        help="output as svg")
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

    # Setup logger with error level
    setup_logger()

    # parsing tree
    try:
        tree, eteformat_flag = utils.validate_tree(args.tree, args.input_type, args.internal)
    except utils.TreeFormatError as e:
        print(e)
        sys.exit(1)

    # resolve polytomy
    if args.resolve_polytomy:
        tree.resolve_polytomy()

    #rest_prop = []
    if args.prop2type:
        if eteformat_flag:
            logger.error("prop2type is not supported for ete format tree, only for newick tree.")
            sys.exit(1)
        else:
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
                'lca':str
                }
        popup_prop_keys = list(prop2type.keys()) 

        if eteformat_flag:
            for path, node in tree.iter_prepostorder():
                prop2type.update(get_prop2type(node))
                
    if args.collapsed_by: # need to be wrap with quotes
        condition_strings = args.collapsed_by
        for idx, condition in enumerate(condition_strings):
            if os.path.isfile(condition):
                color2conditions = build_color2conditions(condition, args.config_sep)
            else:
                syntax_sep = ','
                condition_list = condition.split(syntax_sep)
                color2conditions  = {}
                color2conditions[paired_color[idx]] = condition_list
            c_layout = conditional_layouts.LayoutCollapse(name='Collapsed_by_'+condition, color2conditions=color2conditions, column=level, prop2type = prop2type)
            layouts.append(c_layout)

    # label node by condition
    if args.highlighted_by: # need to be wrap with quotes
        condition_strings = args.highlighted_by
        for idx, condition in enumerate(condition_strings):
            if os.path.isfile(condition):
                color2conditions = build_color2conditions(condition, args.config_sep)
            else:
                syntax_sep = ','
                condition_list = condition.split(syntax_sep)
                color2conditions  = {}
                color2conditions[paired_color[idx]] = condition_list
            s_layout = conditional_layouts.LayoutHighlight(name='Highlighted_by_'+condition, color2conditions=color2conditions, column=level, prop2type = prop2type)
            layouts.append(s_layout)
    
    #### Layouts settings ####
    # numerical representative mearsure 
    internal_num_rep = args.internal_plot_measure

    # color configuration
    color_config = {}
    if args.color_config:
        color_config = read_config_to_dict(args.color_config, delimiter=args.config_sep)

    # Get the input arguments in order
    input_order = []
    for arg in sys.argv[1:]:
        if arg.startswith('-') and arg.endswith('layout'):
            input_order.append(arg[2:])
        else:
            continue
    
    visualized_props = []
    for layout in input_order:
        if layout == 'acr-discrete-layout':
            acr_discrete_layouts, level, color_dict = get_acr_discrete_layouts(tree, args.acr_discrete_layout, level, prop2type=prop2type, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config)
            layouts.extend(acr_discrete_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.acr_discrete_layout)

            #delta statistic 
            for suffix in ['delta', 'pval']:
                visualized_props.extend([utils.add_suffix(prop, suffix) for prop in args.acr_discrete_layout])
           
        if layout == 'acr-continuous-layout':
            acr_continuous_layouts = get_acr_continuous_layouts(tree, args.acr_continuous_layout, level, prop2type=prop2type, padding_x=args.padding_x, padding_y=args.padding_y)
            layouts.extend(acr_continuous_layouts)
            visualized_props.extend(args.acr_continuous_layout)

        if layout == 'ls-layout':
            ls_layouts, ls_props = get_ls_layouts(tree, args.ls_layout, level, prop2type=prop2type, padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config)
            layouts.extend(ls_layouts)
            visualized_props.extend(args.ls_layout)
            visualized_props.extend(ls_props)

        # if layout == 'heatmap-layout':
        #     numerical_props = args.heatmap_layout
        #     norm_method = 'min-max'
        #     heatmap_layouts, level = get_heatmap_matrix_layouts(layout, numerical_props, norm_method, internal_num_rep, color_config, args, level)
        #     layouts.extend(heatmap_layouts)

        # if layout == 'heatmap-mean-layout':
        #     numerical_props = args.heatmap_mean_layout
        #     norm_method = 'mean'
        #     heatmap_layouts, level = get_heatmap_matrix_layouts(layout, numerical_props, norm_method, internal_num_rep, color_config, args, level)
        #     layouts.extend(heatmap_layouts)

        # if layout == 'heatmap-zscore-layout':
        #     numerical_props = args.heatmap_zscore_layout
        #     norm_method = 'zscore'
        #     heatmap_layouts, level = get_heatmap_matrix_layouts(layout, numerical_props, norm_method, internal_num_rep, color_config, args, level)
        #     layouts.extend(heatmap_layouts)

        if layout == 'heatmap-layout':
            heatmap_layouts, level = get_heatmap_layouts(tree, args.heatmap_layout, level, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y, internal_rep=internal_num_rep, color_config=color_config, norm_method='min-max')
            layouts.extend(heatmap_layouts)
            visualized_props.extend(args.heatmap_layout)

        if layout == 'heatmap-mean-layout':
            heatmap_mean_layouts, level = get_heatmap_layouts(tree, args.heatmap_mean_layout, level, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y, internal_rep=internal_num_rep, color_config=color_config, norm_method='mean')
            layouts.extend(heatmap_mean_layouts)
            visualized_props.extend(args.heatmap_mean_layout)

        if layout == 'heatmap-zscore-layout':
            heatmap_zscore_layouts, level = get_heatmap_layouts(tree, args.heatmap_zscore_layout, level, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y, internal_rep=internal_num_rep, color_config=color_config, norm_method='zscore')
            layouts.extend(heatmap_zscore_layouts)
            visualized_props.extend(args.heatmap_zscore_layout)
            
        if layout == 'label-layout':
            label_layouts, level, color_dict = get_label_layouts(tree, args.label_layout, level, prop2type=prop2type, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.label_layout)

        if layout == 'colorbranch-layout':
            categorical_props = [prop for prop in args.colorbranch_layout if prop2type.get(prop) in [str, list, bool, None]]
            if categorical_props:
                colorbranch_layouts, level, color_dict = get_colorbranch_layouts(tree, categorical_props, level, prop2type=prop2type, column_width=args.column_width, padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config)
                layouts.extend(colorbranch_layouts)
                total_color_dict.append(color_dict)
                visualized_props.extend(categorical_props)
                #visualized_props.extend([utils.add_suffix(prop, 'counter') for prop in args.piechart_layout])

            numerical_props = [prop for prop in args.colorbranch_layout if prop2type.get(prop) in [float, int]]
            if numerical_props:
                branchscore_layouts = get_branchscore_layouts(tree, numerical_props, 
                prop2type, padding_x=args.padding_x, padding_y=args.padding_y, 
                internal_rep=internal_num_rep, color_config=color_config)
                layouts.extend(branchscore_layouts)
                visualized_props.extend(numerical_props)
        
        if layout == 'bubble-layout':
            categorical_props = [prop for prop in args.bubble_layout if prop2type.get(prop) in [str, list, bool, None]]
            if categorical_props:
                bubble_layouts, level, color_dict = get_categorical_bubble_layouts(tree, categorical_props, 
                level=level, prop2type=prop2type, 
                padding_x=args.padding_x, padding_y=args.padding_y, 
                color_config=color_config)
                layouts.extend(bubble_layouts)
                total_color_dict.append(color_dict)
                visualized_props.extend(categorical_props)
                #visualized_props.extend([utils.add_suffix(prop, 'counter') for prop in args.piechart_layout

            numerical_props = [prop for prop in args.bubble_layout if prop2type.get(prop) in [float, int]]
            if numerical_props:
                bubble_layouts, level, color_dict = get_numerical_bubble_layouts(tree, numerical_props, 
                level=level, prop2type=prop2type, 
                padding_x=args.padding_x, padding_y=args.padding_y, 
                internal_rep=internal_num_rep, color_config=color_config)
                layouts.extend(bubble_layouts)
                visualized_props.extend(numerical_props)

        if layout == "piechart-layout":
            piechart_layouts = get_piechart_layouts(tree, args.piechart_layout, 
            prop2type=prop2type, 
            padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config)
            layouts.extend(piechart_layouts)
            visualized_props.extend(args.piechart_layout)
            visualized_props.extend([utils.add_suffix(prop, 'counter') for prop in args.piechart_layout])

        if layout == 'rectangle-layout':
            rectangle_layouts, level, color_dict = get_rectangle_layouts(tree, args.rectangle_layout, 
            level, prop2type=prop2type, column_width=args.column_width, 
            padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config)
            layouts.extend(rectangle_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.rectangle_layout)
            visualized_props.extend([utils.add_suffix(prop, 'counter') for prop in args.rectangle_layout])

        if layout == 'background-layout':
            background_layouts, level, color_dict = get_background_layouts(tree, args.background_layout, 
            level, prop2type=prop2type, column_width=args.column_width, 
            padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config)
            layouts.extend(background_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.background_layout)
            visualized_props.extend([utils.add_suffix(prop, 'counter') for prop in args.background_layout])

        if layout == 'binary-layout':
            binary_layouts, level, color_dict = get_binary_layouts(tree, args.binary_layout, level, 
            prop2type=prop2type, column_width=args.column_width, reverse=False, 
            padding_x=args.padding_x, padding_y=args.padding_y, 
            color_config=color_config, same_color=False, aggregate=False)
            layouts.extend(binary_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.binary_layout)

        if layout == 'binary-aggregate-layout':
            binary_aggregate_layouts, level, color_dict = get_binary_layouts(tree, args.binary_aggregate_layout, level, prop2type=prop2type, column_width=args.column_width, reverse=False, padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config, same_color=False, aggregate=True)
            layouts.extend(binary_aggregate_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.binary_aggregate_layout)

        if layout == 'binary-unicolor-layout':
            binary2_layouts, level, color_dict = get_binary_layouts(tree, args.binary_unicolor_layout, level, prop2type=prop2type, column_width=args.column_width, reverse=False, padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config, same_color=True, aggregate=False)
            layouts.extend(binary2_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.binary_unicolor_layout)

        if layout == 'binary-unicolor-aggregate-layout':
            binary2_aggregate_layouts, level, color_dict = get_binary_layouts(tree, args.binary_unicolor_aggregate_layout, level, prop2type=prop2type, column_width=args.column_width, reverse=False, padding_x=args.padding_x, padding_y=args.padding_y, color_config=color_config, same_color=True, aggregate=True)
            layouts.extend(binary2_aggregate_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.binary_unicolor_aggregate_layout)

        if layout == 'barplot-layout':
            barplot_layouts, level, color_dict = get_barplot_layouts(tree, args.barplot_layout, level, 
            prop2type, column_width=args.barplot_width, padding_x=args.padding_x, padding_y=args.padding_y, 
            internal_rep=internal_num_rep, anchor_column=args.barplot_scale, color_config=color_config, barplot_colorby=args.barplot_colorby)
            layouts.extend(barplot_layouts)
            total_color_dict.append(color_dict)
            visualized_props.extend(args.barplot_layout)

        if layout == "branchscore-layout":
            branchscore_layouts = get_branchscore_layouts(tree, args.branchscore_layout, prop2type, padding_x=args.padding_x, padding_y=args.padding_y, internal_rep='avg')
            layouts.extend(branchscore_layouts)
            visualized_props.extend(args.branchscore_layout)

        if layout == 'alignment-layout':
            lengh = len(max(utils.tree_prop_array(tree, 'alignment'),key=len))
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
                matrix, value2color, all_profiling_values = multiple2matrix(tree, profiling_prop, prop2type=prop2type, color_config=color_config, eteformat_flag=eteformat_flag)
                matrix_layout = profile_layouts.LayoutPropsMatrixBinary(name=f"Profiling_{all_profiling_values}",
                matrix=matrix, matrix_props=all_profiling_values, value_range=[0,1],
                value_color=value2color, column=level, poswidth=args.column_width)
                level += 1
                layouts.append(matrix_layout)

        # categorical matrix
        if layout == 'categorical-matrix-layout':
            categorical_props = args.categorical_matrix_layout

            # drawing as array in matrix
            matrix, value2color = categorical2matrix(tree, categorical_props, color_config=color_config)
            matrix_layout = profile_layouts.LayoutPropsMatrixOld(name=f"Categorical_matrix_{categorical_props}",
                matrix=matrix, matrix_type='categorical', matrix_props=categorical_props,
                value_color=value2color, column=level, poswidth=args.column_width)
            
            level += 1
            layouts.append(matrix_layout)

        # numerical matrix
        if layout == 'numerical-matrix-layout':
            numerical_props = args.numerical_matrix_layout
            # matrix, value2color = float2matrix(tree, numerical_props, count_negative=True)
            # all_values = list(value2color.keys())
            # min_val, max_val = min(all_values), max(all_values)
            # matrix_layout = profile_layouts.LayoutPropsMatrix(name=f'Numerical_matrix_{numerical_props}', 
            #     matrix_type='numerical', alignment=matrix, matrix_props=numerical_props, 
            #     profiles=all_values, column=level, summarize_inner_nodes=False, 
            #     value_range = [min_val, max_val], value_color=value2color,
            #     poswidth=args.column_width)
            
            # if is list, it should provide more than one matrix
            matrix, minval, maxval, value2color, results_list, list_props, single_props = numerical2matrix(tree, 
                numerical_props, count_negative=True, internal_num_rep=internal_num_rep, 
                color_config=color_config, norm_method='min-max')
            
            if list_props:
                index_map = {value: idx for idx, value in enumerate(numerical_props)}
                sorted_list_props = sorted(list_props, key=lambda x: index_map[x])
                for list_prop in sorted_list_props:
                    matrix, minval, maxval, value2color = results_list[list_prop]
                    matrix_layout = profile_layouts.LayoutPropsMatrixOld(name=f"Numerical_matrix_{list_prop}", 
                        matrix=matrix, matrix_type='numerical', matrix_props=[list_prop], is_list=True, 
                        value_color=value2color, value_range=[minval, maxval], column=level,
                        poswidth=args.column_width)

                    level += 1
                    layouts.append(matrix_layout)

            if single_props:
                index_map = {value: idx for idx, value in enumerate(numerical_props)}
                sorted_single_props = sorted(single_props, key=lambda x: index_map[x])
                matrix_layout = profile_layouts.LayoutPropsMatrixOld(name=f"Numerical_matrix_{sorted_single_props}", 
                    matrix=matrix, matrix_type='numerical', matrix_props=sorted_single_props, is_list=False, 
                    value_color=value2color, value_range=[minval, maxval], column=level,
                    poswidth=args.column_width)

                level += 1
                layouts.append(matrix_layout)


        if layout == 'binary-matrix-layout':
            binary_props = args.binary_matrix_layout
            matrix, value2color, is_list = binary2matrix(tree, binary_props, color_config=color_config)
            all_values = list(value2color.keys())

            matrix_layout = profile_layouts.LayoutPropsMatrixBinary(name=f"Binary_matrix_{binary_props}",
                matrix=matrix, matrix_props=binary_props, value_range=[0,1],
                value_color=value2color, column=level, poswidth=args.column_width)

            level += 1
            layouts.append(matrix_layout)

        if layout == "taxoncollapse-layout":
            taxon_color_dict = {}
            taxa_layouts = []

            # generate a rank2values dict for pre taxonomic annotated tree
            if not rank2values:
                rank2values = defaultdict(list)
                for n in tree.traverse():
                    if n.props.get('lca'):
                        lca_dict = utils.string_to_dict(n.props.get('lca'))
                        for rank, sci_name in lca_dict.items():
                            rank2values[rank].append(sci_name)

                    current_rank = n.props.get('rank')
                    if current_rank and current_rank != 'Unknown':
                        rank2values[current_rank].append(n.props.get('sci_name',''))
            else:       
                pass
            
            # assign color for each value of each rank
            for rank, value in sorted(rank2values.items()):
                value = list(set(value))
                color_dict = utils.assign_color_to_values(value, paired_color)
                taxa_layout = taxon_layouts.TaxaCollapse(name = "TaxaCollapse_"+rank, rank=rank, rect_width=args.column_width, color_dict=color_dict, column=level)
                taxa_layouts.append(taxa_layout)
                
            layouts = layouts + taxa_layouts
            level += 1

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
            matrix, value2color, all_profiling_values = multiple2matrix(tree, multiple_text_prop)
            multiple_text_prop_layout = profile_layouts.LayoutPropsMatrixBinary(name=f"Profiling_{multiple_text_prop}",
            matrix=matrix, matrix_props=all_profiling_values, value_range=[0,1],
            active=False,
            value_color=value2color, column=level, poswidth=args.column_width)

            # matrix, all_values = multiple2profile(tree, multiple_text_prop)
            # multiple_text_prop_layout = profile_layouts.LayoutProfile(
            #     name="Profiling_"+multiple_text_prop, 
            #     mode='profiles', 
            #     alignment=matrix, 
            #     profiles=all_values, 
            #     active=False,
            #     column=level)

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
                    rank = n.props.get('rank')
                    rank2values[rank].append(n.props.get('sci_name',''))
        else:       
            pass

        
        # assign color for each value of each rank
        for rank, value in sorted(rank2values.items()):
            value = list(set(value))
            color_dict = utils.assign_color_to_values(value, paired_color)
            if args.taxonclade_layout:
                taxa_layout = taxon_layouts.TaxaClade(name='TaxaClade_'+rank, level=level, rank = rank, color_dict=color_dict)
                taxa_layouts.append(taxa_layout)

            if args.taxonrectangle_layout:
                taxa_layout = taxon_layouts.TaxaRectangular(name = "TaxaRect_"+rank, rank=rank, rect_width=args.column_width, color_dict=color_dict, column=level)
                taxa_layouts.append(taxa_layout)
                #level += 1

            # if args.taxoncollapse_layout:
            #     taxa_layout = taxon_layouts.TaxaCollapse(name = "TaxaCollapse_"+rank, rank=rank, rect_width=args.column_width, color_dict=color_dict, column=level)
            #     taxa_layouts.append(taxa_layout)

            taxon_color_dict[rank] = color_dict
            
        #taxa_layouts.append(taxon_layouts.TaxaRectangular(name = "Last Common Ancester", color_dict=taxon_color_dict, column=level))
        taxa_layouts.append(taxon_layouts.LayoutSciName(name = 'Taxa Scientific name', color_dict=taxon_color_dict))
        taxa_layouts.append(taxon_layouts.LayoutEvolEvents(name='Taxa Evolutionary events', prop="evoltype",
            speciation_color="blue", 
            duplication_color="red", node_size = 3,
            legend=True))
        layouts = layouts + taxa_layouts
        level += 1
        total_color_dict.append(taxon_color_dict)
        
    #### prune at the last step in case of losing leaves information
    # prune tree by rank
    if args.rank_limit:
        tree = utils.taxatree_prune(tree, rank_limit=args.rank_limit)
        
    # prune tree by condition 
    if args.pruned_by: # need to be wrap with quotes
        condition_strings = args.pruned_by
        tree = utils.conditional_prune(tree, condition_strings, prop2type)

    #### Output #####
    popup_prop_keys.extend(list(set(visualized_props)))
    popup_prop_keys = tuple(popup_prop_keys)
    
    if args.out_colordict:
        wrtie_color(total_color_dict)
    if args.render:
        file_path = "tree-1.svg"
        get_image(tree, layouts, args.port, os.path.abspath(file_path))
    else:
        tree.explore(keep_server=True, compress=False, quiet=args.verbose, 
        layouts=layouts, port=args.port, include_props=sorted(popup_prop_keys),
        show_leaf_name=args.hide_leaf_name, show_branch_support=args.hide_branch_support,
        show_branch_length=args.hide_branch_distance)

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

def read_config_to_dict(file_obj, delimiter):
    """
    Reads a configuration file to a dictionary.

    The configuration file should have the format:
    prop,value,color
    random_type,low,green
    ...
    :param filename: Path to the file.
    :param delimiter: Delimiter used in the configuration file.
    :return: A dictionary with (prop, value) tuple as keys and color as values.
    """
    config_dict = {}
    # Reset file pointer to start, in case it's been accessed before
    file_obj.seek(0)
    lines = file_obj.readlines()
    filtered_lines = [line for line in lines if not line.strip().startswith('#') and line.strip()]

    #reader = csv.DictReader(file_obj, delimiter=delimiter)
    reader = csv.reader(filtered_lines, delimiter=delimiter)
    headers = next(reader)  # Get the headers
    for row in reader:
        # Map the headers to the row values
        row_dict = dict(zip(headers, row))

        prop = row_dict['PROP']
        value = row_dict['VALUE']
        color = row_dict['COLOR']
        detail = row_dict.get('CONDITION')

        # Initialize property if not present
        if prop not in config_dict:
            config_dict[prop] = {"value2color": {}, "detail2color": {}}

        # Assign colors based on presence of detail or value
        if detail:
            if detail.lower() not in config_dict[prop]["detail2color"]:
                config_dict[prop]["detail2color"][detail.lower()] = (color, value)
            else:
                config_dict[prop]["detail2color"][detail.lower()] = []
                config_dict[prop]["detail2color"][detail.lower()].append((color, value))
        if value:
            config_dict[prop]["value2color"][value] = color

    return config_dict

def process_common_ancestors(color_dict, tree, common_ancestor_separator='||'):
    for key in list(color_dict.keys()):
        if common_ancestor_separator in key:
            children = key.split(common_ancestor_separator)
            ancestor = tree.common_ancestor(children)
            if ancestor:
                if ancestor.name:
                    # If the ancestor has a name, update its color in the dictionary
                    color_dict[ancestor.name] = color_dict[key]
                    del color_dict[key]
                else:
                    # If the ancestor has no name, assign the key as its name
                    ancestor.name = key
    return color_dict

def build_color2conditions(condition_file, config_sep):
    color2conditions = {}
    
    with open(condition_file, 'r') as f:
        for line in f:
            line = line.rstrip()  # Remove trailing whitespace
            if line and not line.startswith('#'):  # Ignore empty lines and comments
                if not line.startswith('PROP'):  # Ignore lines starting with 'PROP'
                    sep = config_sep
                    parts = line.split(sep)
                    
                    # Ensure there are at least 4 parts after splitting
                    if len(parts) >= 4:
                        left = parts[0]
                        right = parts[1]
                        color = parts[2]
                        operator = parts[3]
                        
                        # Construct the condition string
                        condition_string = ''.join([left, operator, right])
                        
                        # Add the condition string to the dictionary under the corresponding color
                        if color not in color2conditions:
                            color2conditions[color] = []
                        color2conditions[color].append(condition_string)
                    else:
                        logger.error(f"Invalid line: {line}")
                        sys.exit(1)
    return color2conditions

def get_acr_discrete_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0, color_config=None):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        if prop2type and prop2type.get(prop) == list:
            leaf_values = list(map(list,set(map(tuple, utils.tree_prop_array(tree, prop)))))    
            prop_values = [val for sublist in leaf_values for val in sublist]
        else:
            prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))
        
        color_dict = {} # key = value, value = color id
        if color_config and color_config.get(prop):
            if color_config.get(prop).get('value2color'):
                color_dict = color_config.get(prop).get('value2color')

            # Check if all property values have an assigned color
            existing_values = set(color_dict.keys())
            additional_values = set(prop_values) - existing_values
            if additional_values:
                # Fetch new colors for the additional values
                additional_colors = utils.assign_color_to_values(sorted(additional_values), paired_color)
                color_dict.update(additional_colors)
        else:
            # normal text prop
            color_dict = utils.assign_color_to_values(prop_values, paired_color)

        layout = phylosignal_layouts.LayoutACRDiscrete(name='acr_'+prop, column=level, \
            color_dict=color_dict, acr_prop=prop, width=column_width, \
                padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_acr_continuous_layouts(tree, props, level, prop2type, padding_x=1, padding_y=0):
    gradientscolor = utils.build_color_gradient(20, colormap_name='jet')
    layouts = []
    for prop in props:
        try:
            all_values = np.array(sorted(list(set(utils.tree_prop_array(tree, prop, numeric=True))))).astype('float64')
        except ValueError:
            logger.error(f"Property {prop} is not numeric. Please check the property type.")
            sys.exit(1)
        all_values = all_values[~np.isnan(all_values)]
        minval, maxval = all_values.min(), all_values.max()
        num = len(gradientscolor)
        index_values = np.linspace(minval, maxval, num)
        value2color = {}
        for search_value in all_values:
            index = np.abs(index_values - search_value).argmin()+1
            value2color[search_value] = gradientscolor[index]
        layout = phylosignal_layouts.LayoutACRContinuous(name='acr_'+prop, column=level, \
            color_dict=value2color, score_prop=prop, value_range=[minval, maxval], \
            color_range=[gradientscolor[20], gradientscolor[10], gradientscolor[1]])
        layouts.append(layout)
    return layouts

def get_ls_layouts(tree, props, level, prop2type, padding_x=1, padding_y=0, color_config=None):
    precision_suffix = "prec"
    sensitivity_suffix = "sens"
    f1_suffix  = "f1"
    ls_clade_suffix = "ls_clade"
    ls_clade_props = [utils.add_suffix(prop, ls_clade_suffix) for prop in props]
    lsprop2color = utils.assign_color_to_values(ls_clade_props, paired_color)

    layouts = []
    ls_props = []
    for prop in props:
        value2color = {}
        if color_config and color_config.get(prop) is not None:
            prop_config = color_config[prop]
            
            color_dict = {}
            # First, try to use value2color mappings if they exist and are applicable
            if 'value2color' in prop_config and prop_config['value2color']:
                color_dict = prop_config['value2color']
                sorted_color_dict = {float(key): value for key, value in color_dict.items()}
                gradientscolor = sorted_color_dict.values()
            elif 'detail2color' in prop_config and prop_config['detail2color']:
                min_color = prop_config['detail2color'].get('color_min', 'white')
                max_color = prop_config['detail2color'].get('color_max', 'red')
                mid_color = prop_config['detail2color'].get('color_mid', None)
                gradientscolor = utils.build_custom_gradient(20, min_color, max_color, mid_color)
        else:
            gradientscolor = utils.build_color_gradient(20, colormap_name='bwr')

        for suffix in [precision_suffix, sensitivity_suffix, f1_suffix]:
            
            ls_prop = utils.add_suffix(prop, suffix)
            minval, maxval = 0, 1
            
            # get value
            internalnode_all_values = np.array(sorted(list(set(utils.tree_prop_array(tree, ls_prop, numeric=True))))).astype('float64')
            all_values = internalnode_all_values[~np.isnan(internalnode_all_values)]
            num = len(gradientscolor)
            index_values = np.linspace(minval, maxval, num)
            if all_values:
                for search_value in all_values:
                    if search_value not in value2color:
                        index = np.abs(index_values - search_value).argmin()+1
                        value2color[search_value] = gradientscolor[index]
            else:
                logger.error(f"Property {ls_prop} is empty. Please check annotation.")
                sys.exit(1)
            # layout = staple_layouts.LayoutBranchScore(name='BranchScore_'+prop, \
            # color_dict=gradientscolor, score_prop=prop, internal_rep=internal_rep, \
            # value_range=[minval, maxval], \
            # color_range=[gradientscolor[20], gradientscolor[10], gradientscolor[1]])
            if suffix != "f1":
                layout = staple_layouts.LayoutBranchScore(name='ls_'+ls_prop, \
                    color_dict=value2color, score_prop=ls_prop, value_range=[minval, maxval], \
                    color_range=[gradientscolor[20], gradientscolor[10], gradientscolor[1]], 
                    show_score=True, active=False)
            else:
                layout = staple_layouts.LayoutBranchScore(name='ls_'+ls_prop, \
                    color_dict=value2color, score_prop=ls_prop, value_range=[minval, maxval], \
                    color_range=[gradientscolor[20], gradientscolor[10], gradientscolor[1]], 
                    show_score=True)
            
            layouts.append(layout)
            ls_props.append(ls_prop)
    
        ls_clade_prop = utils.add_suffix(prop, ls_clade_suffix)
        ls_clade_layout = phylosignal_layouts.LayoutLineageSpecific(name=f'Linear Specific Clade {prop}', \
            ls_prop=ls_clade_prop, color=lsprop2color[ls_clade_prop])

        layouts.append(ls_clade_layout)
        ls_props.append(ls_clade_prop)

    return layouts, ls_props

def get_piechart_layouts(tree, props, prop2type, padding_x=1, padding_y=0, radius=20, color_config=None):
    layouts = []
    for prop in props:
        color_dict = {}
        if color_config and color_config.get(prop):
            if color_config.get(prop).get('value2color'):
                color_dict = color_config.get(prop).get('value2color')
        else: 
            if prop2type and prop2type.get(prop) == list:
                leaf_values = list(map(list, set(map(tuple, utils.tree_prop_array(tree, prop)))))
                prop_values = [val for sublist in leaf_values for val in sublist]
            else:
                prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))

            if not prop_values:
                logger.error(f"Property {prop} is empty. Please check annotation.")
                sys.exit(1)
            color_dict = utils.assign_color_to_values(prop_values, paired_color)
        layout = text_layouts.LayoutPiechart(name='Piechart_'+prop, color_dict=color_dict, text_prop=prop, radius=radius)
        layouts.append(layout)
    return layouts

def get_label_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0, color_config=None):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        color_dict = {}
        if color_config and color_config.get(prop):
            if color_config.get(prop).get('value2color'):
                color_dict = color_config.get(prop).get('value2color')
        else: 
            if prop2type and prop2type.get(prop) == list:
                leaf_values = list(map(list,set(map(tuple,utils.tree_prop_array(tree, prop)))))
                prop_values = [val for sublist in leaf_values for val in sublist]
            else:
                prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))

            if not prop_values:
                logger.error(f"Property {prop} is empty. Please check annotation.")
                sys.exit(1)

            color_dict = utils.assign_color_to_values(prop_values, paired_color)

        layout = text_layouts.LayoutText(name='Label_'+prop, column=level, 
        color_dict=color_dict, text_prop=prop, width=column_width, padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_colorbranch_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0, color_config=None):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        color_dict = {} # key = value, value = color id
        if color_config and color_config.get(prop):
            if color_config.get(prop).get('value2color'):
                color_dict = color_config.get(prop).get('value2color')
            
            if prop == 'name':
                color_dict = process_common_ancestors(color_dict, tree, common_ancestor_separator='||')

            # Check if all property values have an assigned color
            # prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))
            # existing_values = set(color_dict.keys())
            # additional_values = set(prop_values) - existing_values
            # if additional_values:
            #     # Fetch new colors for the additional values
            #     additional_colors = utils.assign_color_to_values(sorted(additional_values), paired_color)
            #     color_dict.update(additional_colors)

        else: 
            if prop2type and prop2type.get(prop) == list:
                leaf_values = list(map(list,set(map(tuple,utils.tree_prop_array(tree, prop)))))    
                prop_values = [val for sublist in leaf_values for val in sublist]
            else:
                prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))
            
            # normal text prop
            color_dict = utils.assign_color_to_values(prop_values, paired_color)

        layout = text_layouts.LayoutColorbranch(name='Colorbranch_'+prop, column=level, \
            color_dict=color_dict, text_prop=prop, width=column_width, \
                padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_rectangle_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0, color_config=None):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        color_dict = {} # key = value, value = color id
        if color_config and color_config.get(prop):
            if color_config.get(prop).get('value2color'):
                color_dict = color_config.get(prop).get('value2color')
        else:
            if prop2type and prop2type.get(prop) == list:
                leaf_values = list(map(list,set(map(tuple,utils.tree_prop_array(tree, prop)))))    
                prop_values = [val for sublist in leaf_values for val in sublist]
            else:
                prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))
            
            if not prop_values:
                logger.error(f"Property {prop} is empty. Please check annotation.")
                sys.exit(1)

            # normal text prop
            color_dict = utils.assign_color_to_values(prop_values, paired_color)

        layout = text_layouts.LayoutRect(name='Rectangular_'+prop, column=level,
                    color_dict=color_dict, text_prop=prop,
                    width=column_width, padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_background_layouts(tree, props, level, prop2type, column_width, padding_x=1, padding_y=0, color_config=None):
    prop_color_dict = {}
    layouts = []
    for prop in props:
        color_dict = {} # key = value, value = color id
        if color_config and color_config.get(prop):
            if color_config.get(prop).get('value2color'):
                color_dict = color_config.get(prop).get('value2color')

            if prop == 'name':
                color_dict = process_common_ancestors(color_dict, tree, common_ancestor_separator='||')
        else:
            if prop2type and prop2type.get(prop) == list:
                leaf_values = list(map(list,set(map(tuple,utils.tree_prop_array(tree, prop)))))    
                prop_values = [val for sublist in leaf_values for val in sublist]
            else:
                prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))
            
            if not prop_values:
                logger.error(f"Property {prop} is empty. Please check annotation.")
                sys.exit(1)
                
            # normal text prop
            color_dict = utils.assign_color_to_values(prop_values, paired_color)
            
        layout = text_layouts.LayoutBackground(name='Background_'+prop, 
                    column=level, width=column_width,
                    color_dict=color_dict, text_prop=prop,
                    padding_x=padding_x, padding_y=padding_y)

        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_binary_layouts(tree, props, level, prop2type, column_width=70, reverse=False, padding_x=1, padding_y=0, color_config=None, same_color=False, aggregate=False):
    prop_color_dict = {}
    layouts = []
    
    for prop in props:
        #prop_values = sorted(list(set(utils.tree_prop_array(tree, prop, leaf_only=True))))
        prop_values = utils.tree_prop_array(tree, prop, leaf_only=True)

        if not reverse:
            max_count = utils.find_bool_representations(prop_values)
        else:
            max_count = utils.find_bool_representations(prop_values, rep=False)

        # If you still need a sorted list of unique property values
        prop_values = sorted(set(prop_values))

        if can_convert_to_bool(prop_values):
            if color_config and color_config.get(prop):
                if color_config.get(prop).get('value2color'):
                    color_dict = color_config.get(prop).get('value2color')

                    if can_convert_to_bool(color_dict.keys()):
                        color_dict = {eval(k): v for k, v in color_dict.items()}
                        color = color_dict.get(True, "#ff0000") #get true color
                        if color_dict.get(False):
                            print("Warning: False color is not supported in binary layout. Ignored. ")
            else:
                if same_color:
                    color = "#ff0000"
                else:
                    if level >= len(paired_color):
                        color =  utils.random_color(h=None)
                    else:
                        color = paired_color[level]

            if not reverse:
                layout = conditional_layouts.LayoutBinary('Binary_'+prop, level, bool_prop=prop, color=color, width=column_width, padding_x=padding_x, padding_y=padding_y, reverse=reverse, aggregate=aggregate, max_count=max_count)
            else:
                layout = conditional_layouts.LayoutBinary('ReverseBinary_'+prop, level, bool_prop=prop, width=column_width, padding_x=padding_x, padding_y=0, reverse=reverse, aggregate=aggregate, max_count=max_count)
            
            internal_prop = utils.add_suffix(prop, 'counter')

            layouts.append(layout)
            level += 1
        else:
            logger.error(f"Property {prop} is not binary trait.")
            sys.exit(1)
    return layouts, level, prop_color_dict

def get_branchscore_layouts(tree, props, prop2type, padding_x=1, padding_y=0, internal_rep='avg', color_config=None):
    """
    Output dictionary of each score prop and corresponding color.
    """

    def parse_color_config(prop, color_config, minval, maxval):
        max_color = 'red'
        min_color = 'white'
        mid_color = None
        value2color = {}

        prop_config = color_config.get(prop, {})
        color_dict = prop_config.get('value2color', {})

        if color_dict:
            value2color = {float(key): value for key, value in color_dict.items()}

        detail2color = prop_config.get('detail2color', {})

        temp_min_color, temp_min_val = detail2color.get('color_min', (None, None))
        temp_max_color, temp_max_val = detail2color.get('color_max', (None, None))
        temp_mid_color, temp_mid_val = detail2color.get('color_mid', (None, None))

        if temp_max_color:
            max_color = temp_max_color
        if temp_min_color:
            min_color = temp_min_color
        if temp_mid_color:
            mid_color = temp_mid_color

        if temp_min_val:
            minval = float(temp_min_val)
        if temp_max_val:
            maxval = float(temp_max_val)

        gradientscolor = utils.build_custom_gradient(20, min_color, max_color, mid_color)

        return gradientscolor, value2color, minval, maxval

    layouts = []

    for prop in props:
        # Get leaf values of each prop
        leaf_all_values = np.array(sorted(list(set(utils.tree_prop_array(tree, prop, numeric=True))))).astype('float64')

        # Get internal values of each prop
        internal_prop = utils.add_suffix(prop, internal_rep)
        internalnode_all_values = np.array(sorted(list(set(utils.tree_prop_array(tree, internal_prop, numeric=True))))).astype('float64')
        all_values = np.concatenate((leaf_all_values, internalnode_all_values))
        all_values = all_values[~np.isnan(all_values)]
        value2color = {}
        minval, maxval = all_values.min(), all_values.max()

        if color_config and color_config.get(prop) is not None:
            gradientscolor, value2color, minval, maxval = parse_color_config(prop, color_config, minval, maxval)
        else:
            gradientscolor = utils.build_color_gradient(20, colormap_name='jet')

        # Preload corresponding gradient color of each value
        num = len(gradientscolor)
        index_values = np.linspace(minval, maxval, num)

        for search_value in all_values:
            if search_value not in value2color:
                index = np.abs(index_values - search_value).argmin() + 1
                value2color[search_value] = gradientscolor[index]

        # Get corresponding gradient color on the fly of visualization
        layout = staple_layouts.LayoutBranchScore(
            name='BranchScore_' + prop,
            color_dict=value2color,
            score_prop=prop,
            internal_rep=internal_rep,
            value_range=[minval, maxval],
            color_range=[gradientscolor[20], gradientscolor[10], gradientscolor[1]]
        )
        layouts.append(layout)

    return layouts

def get_barplot_layouts(tree, props, level, prop2type, column_width=70, padding_x=1, padding_y=0, internal_rep='avg', anchor_column=None, color_config=None, barplot_colorby=None):
    def get_barplot_color(level):
        global paired_color
        """Determines the color for the barplot based on the level and available paired colors."""
        if level > len(paired_color):
            return utils.random_color(h=None)
        else:
            return paired_color[level]

    def process_prop_values(tree, prop):
        """Extracts and processes property values, excluding NaNs."""
        prop_values = np.array(list(set(utils.tree_prop_array(tree, prop)))).astype('float64')
        prop_values = prop_values[~np.isnan(prop_values)]
        if prop_values.size != 0:
            return prop_values
        else:
            logger.error(f"Tree doesn't have '{prop}' property")
            sys.exit(1)

    def calculate_column_width(prop_values, anchormax=None):
        """Calculates new column width based on property values and optional anchormax."""
        if anchormax is not None:
            minval, maxval = prop_values.min(), prop_values.max()
            return maxval / (anchormax / column_width)
        return column_width

    def configure_layout(prop, new_column_width, color_dict, color_prop, size_prop, barplot_color=None, size_range=[]):
        """Configures and returns the layout for the current property."""
        layout_params = {
            'name': f'Barplot_{prop}',
            'prop': prop,
            'width': new_column_width,
            'color': None if color_dict else barplot_color,
            'colors': color_dict,
            'color_prop': color_prop,
            'size_prop': size_prop,
            'column': level,
            'internal_rep': internal_rep,
            'padding_x': padding_x * 10,
            'size_range': size_range,
        }
        if color_dict is None:
            del layout_params['colors']
            del layout_params['color_prop']
        else:
            del layout_params['color']
        return staple_layouts.LayoutBarplot(**layout_params)

    prop_color_dict = {}
    layouts = []
    barplot_minval = 0

    # Initialize anchor column values if provided
    anchormax = None
    if anchor_column:
        anchor_column_values = process_prop_values(tree, anchor_column)
        anchormax = anchor_column_values.max()

    for prop in props:
        prop_values = process_prop_values(tree, prop)    
        maxval = prop_values.max()
        size_prop = prop if prop_values.any() else f"{prop}_{internal_rep}"
        new_column_width = calculate_column_width(prop_values, anchormax)
        

        # Determine color configuration if available
        if color_config:
            for key, value in color_config.items():
                if 'barplot_colorby' in list(color_config.get(key).get('detail2color').keys()):
                    color_dict = color_config.get(key).get('value2color')
                    color_prop = key
                    barplot_color = None
                else:        
                    color_dict = color_config.get(prop, {}).get('value2color', None)
                    color_prop = None
                    barplot_color = None
        else:
            # Apply default color logic
            color_dict = None
            color_prop = None
            barplot_color = None
            if barplot_colorby:
                prop_values = sorted(list(set(utils.tree_prop_array(tree, barplot_colorby))))
                color_dict = utils.assign_color_to_values(prop_values, paired_color)
                color_prop = barplot_colorby
            else:
                barplot_color = get_barplot_color(level)
                prop_color_dict[prop] = barplot_color

        # Configure and add layout
        if maxval and maxval > barplot_minval:
            size_range = [barplot_minval, maxval]
        else:
            size_range = []
        layout = configure_layout(prop, new_column_width, color_dict, color_prop, size_prop, barplot_color, size_range)
        layouts.append(layout)
        level += 1
    
    return layouts, level, prop_color_dict

def get_categorical_bubble_layouts(tree, props, level, prop2type, column_width=70, padding_x=0, padding_y=0, color_config=None):
    prop_color_dict = {}
    layouts = []
    max_radius = 15
    for prop in props:
        color_dict = {} # key = value, value = color id
        if color_config and color_config.get(prop):
            if color_config.get(prop).get('value2color'):
                color_dict = color_config.get(prop).get('value2color')
        else:
            if prop2type and prop2type.get(prop) == list:
                leaf_values = list(map(list,set(map(tuple,utils.tree_prop_array(tree, prop)))))    
                prop_values = [val for sublist in leaf_values for val in sublist]
            else:
                prop_values = sorted(list(set(utils.tree_prop_array(tree, prop))))
            
            if not prop_values:
                logger.error(f"Property {prop} is empty. Please check annotation.")
                sys.exit(1)

            # normal text prop
            color_dict = utils.assign_color_to_values(prop_values, paired_color)

        # layout = text_layouts.LayoutRect(name='Rectangular_'+prop, column=level,
        #             color_dict=color_dict, text_prop=prop,
        #             width=column_width, padding_x=padding_x, padding_y=padding_y)
        # Configure and add layout
        layout = text_layouts.LayoutBubbleCategorical(name=f'Bubble_{prop}', column=level, 
        prop=prop, color_dict=color_dict, 
        max_radius=max_radius, padding_x=padding_x, padding_y=padding_y)

        layouts.append(layout)
        level += 1
    return layouts, level, prop_color_dict

def get_numerical_bubble_layouts(tree, props, level, prop2type, padding_x=0, padding_y=0, internal_rep='avg', color_config=None):
    def process_prop_values(tree, prop):
        """Extracts and processes property values, excluding NaNs."""
        prop_values = np.array(list(set(utils.tree_prop_array(tree, prop)))).astype('float64')
        return prop_values[~np.isnan(prop_values)]

    prop_color_dict = {}
    layouts = []
    max_radius = 15

    for prop in props:
        prop_values = process_prop_values(tree, prop)

        #minval, maxval = all_prop_values.min(), all_prop_values.max()
        abs_maxval = np.abs(prop_values).max()
        size_prop = prop if prop_values.any() else f"{prop}_{internal_rep}"

        # Configure and add layout
        layout = staple_layouts.LayoutBubbleNumerical(name=f'Bubble_{prop}', column=level, 
        prop=prop, max_radius=max_radius, abs_maxval=abs_maxval, 
        padding_x=padding_x, padding_y=padding_y)
        layouts.append(layout)
        level += 1

    return layouts, level, prop_color_dict

def get_heatmap_layouts(tree, props, level, column_width=70, padding_x=1, padding_y=0, internal_rep='avg', color_config=None, norm_method='min-max'):
    def min_max_normalize(value, minval, maxval):
        if maxval - minval == 0:
            return 0
        else:
            return (value - minval) / (maxval - minval)

    def mean_normalize(value, mean_val, minval, maxval):
        if maxval - minval == 0:
            return 0
        else:
            return (value - mean_val) / (maxval - minval)

    def z_score_normalize(value, mean_val, std_val):
        if std_val == 0:
            return 0
        else:
            return (value - mean_val) / std_val
    
    def check_list_type(lst):
        """
        Check if the input is:
        - A list of strings or floats
        - A list of lists of strings or floats
        """
        # Check if the first element is a list
        if isinstance(lst[0], list):
            logger.error("Array are not supported. please use make sure you have use 'treeprofiler annotate ... --data-matrix <data.array>' and now use --numerical-matrix-layout <data.array> instead.")
            sys.exit(1)
        else:
            # Check if all elements are strings or floats
            if all(isinstance(el, float) for el in lst):
                return True
            else:
                False


    def parse_color_config(prop, color_config, minval, maxval):
        max_color = 'red'
        min_color = 'white'
        mid_color = None
        nan_color = '#EBEBEB'
        value2color = {}

        prop_config = color_config.get(prop, {})
        color_dict = prop_config.get('value2color', {})

        if color_dict:
            value2color = {float(key): value for key, value in color_dict.items()}

        detail2color = prop_config.get('detail2color', {})

        temp_min_color, temp_min_val = detail2color.get('color_min', (None, None))
        temp_max_color, temp_max_val = detail2color.get('color_max', (None, None))
        temp_mid_color, temp_mid_val = detail2color.get('color_mid', (None, None))
        temp_none_color, _ = detail2color.get('color_nan', (None, None))

        if temp_max_color:
            max_color = temp_max_color
        if temp_min_color:
            min_color = temp_min_color
        if temp_mid_color:
            mid_color = temp_mid_color
        if temp_none_color:
            nan_color = temp_none_color
        if temp_min_val:
            minval = float(temp_min_val)
        if temp_max_val:
            maxval = float(temp_max_val)

        gradientscolor = utils.build_custom_gradient(20, min_color, max_color, mid_color)

        return gradientscolor, value2color, minval, maxval, nan_color

    layouts = []
    all_values = []

    for prop in props:
        value2color = {}
        leaf_all_values_raw = utils.tree_prop_array(tree, prop, numeric=True)
        leaf_all_values = np.array(sorted(list(set(leaf_all_values_raw)))).astype('float64')
        internal_prop = utils.add_suffix(prop, internal_rep)
        internalnode_all_values = np.array(sorted(list(set(utils.tree_prop_array(tree, internal_prop, numeric=True))))).astype('float64')
        prop_all_values = np.concatenate((leaf_all_values, internalnode_all_values))
        prop_all_values = prop_all_values[~np.isnan(prop_all_values)]

        minval, maxval = np.min(prop_all_values), np.max(prop_all_values)
        mean_val = np.mean(prop_all_values)
        std_val = np.std(prop_all_values)
        
        if color_config and color_config.get(prop) is not None:
            gradientscolor, value2color, minval, maxval, nan_color = parse_color_config(prop, color_config, minval, maxval)
        else:
            gradientscolor = None
            nan_color = '#EBEBEB'

        if not gradientscolor:
            if norm_method == 'min-max':
                gradientscolor = utils.build_color_gradient(20, colormap_name="Reds")
            else: # "mean" "zscore"
                gradientscolor = utils.build_color_gradient(20, colormap_name="coolwarm")
        
        num = len(gradientscolor)
        for search_value in prop_all_values:
            if search_value is None or math.isnan(search_value):
                value2color[search_value] = nan_color
            #value2color[search_value] = _get_color(search_value, gradientscolor, norm_method)
            else:
                search_value = float(search_value)
                if search_value not in value2color:
                    if norm_method == "min-max":
                        normalized_value = min_max_normalize(search_value, minval, maxval)
                        index_values = np.linspace(0, 1, num)
                    elif norm_method == "mean":
                        normalized_value = mean_normalize(search_value, mean_val, minval, maxval)
                        index_values = np.linspace(-1, 1, num)
                    elif norm_method == "zscore":
                        normalized_value = z_score_normalize(search_value, mean_val, std_val)
                        index_values = np.linspace(-3, 3, num)
                    else:
                        logger.error("Unsupported normalization method.")
                        sys.exit(1)
                    index = np.abs(index_values - normalized_value).argmin() + 1
                    value2color[search_value] = gradientscolor.get(index, "")

        layout = staple_layouts.LayoutHeatmap(name=f'Heatmap_{prop}_{norm_method}', column=level,
                    width=column_width, padding_x=padding_x, padding_y=padding_y, \
                    internal_rep=internal_rep, heatmap_prop=prop, maxval=maxval, minval=minval,\
                    value_color=value2color, value_range=[minval, maxval], color_range=gradientscolor,
                    absence_color=nan_color)

        layouts.append(layout)
        level += 1

    return layouts, level
    
def get_heatmap_matrix_layouts(layout_name, numerical_props, norm_method, internal_num_rep, color_config, args, level):
    layouts = []
    matrix, minval, maxval, value2color, results_list, list_props = numerical2matrix(tree,
                                                            numerical_props,
                                                            count_negative=True,
                                                            internal_num_rep=internal_num_rep,
                                                            color_config=color_config,
                                                            norm_method=norm_method)

    if not list_props:
        matrix_layout = profile_layouts.LayoutPropsMatrixOld(name=f"Heatmap_{numerical_props}",
                                                            matrix=matrix,
                                                            matrix_type='numerical',
                                                            matrix_props=numerical_props,
                                                            is_list=False,
                                                            value_color=value2color,
                                                            value_range=[minval, maxval],
                                                            column=level,
                                                            poswidth=args.column_width)

        level += 1
        layouts.append(matrix_layout)
    else:
        list_props = list(list_props)
        for list_prop in list_props:
            matrix, minval, maxval, value2color = results_list[list_prop]
            matrix_layout = profile_layouts.LayoutPropsMatrixOld(name=f"Heatmap_{list_prop}",
                                                            matrix=matrix,
                                                            matrix_type='numerical',
                                                            matrix_props=list_prop,
                                                            is_list=True,
                                                            value_color=value2color,
                                                            value_range=[minval, maxval],
                                                            column=level,
                                                            poswidth=args.column_width)

            level += 1
            layouts.append(matrix_layout)
    return layouts, level

def get_prop2type(node):
    output = {}
    prop2value = node.props
    if '_speciesFunction' in prop2value:
        del prop2value['_speciesFunction']
    
    for prop, value in prop2value.items():
        if value != 'NaN':
            if isinstance(value, numbers.Number):
                output[prop] = float
            elif type(value) == list:
                output[prop] = list
            else:
                output[prop] = str    
    return output

def categorical2matrix(tree, profiling_props, dtype=str, color_config=None):
    """
    Input:
    tree: A tree structure with nodes, each having properties.
    profiling_props: A list of property names to be processed for each leaf in the tree.
    
    Output:
    A dictionary of matrix representation of the tree leaves and their properties.
    A sorted dictionary mapping property values to their corresponding colors.
    """
    absence_value = "NaN"
    absence_color = "#EBEBEB"
    leaf2matrix = {}
    for node in tree.traverse():
        if node.is_leaf:
            leaf2matrix[node.name] = []
            for profiling_prop in profiling_props:
                if node.props.get(profiling_prop) is not None:
                    if dtype == str:
                        val = node.props.get(profiling_prop)
                    leaf2matrix[node.name].append(val)
                else:
                    leaf2matrix[node.name].append(absence_value)

    # get color
    value2color = {} # key = value, value = color id
    matrix_prop = '*' # special case for matrix
    
    if color_config and color_config.get(matrix_prop):
        if color_config.get(matrix_prop).get('value2color'):
            value2color = color_config.get(matrix_prop).get('value2color')
            value2color[absence_value] = absence_color
    else:
        all_values = sorted(list(set(utils.flatten([sublist for sublist in leaf2matrix.values()]))))
        value2color = utils.assign_color_to_values(all_values, paired_color)
        if absence_value in value2color:
            value2color[absence_value] = absence_color
    
    return leaf2matrix, value2color

def numerical2matrix(tree, profiling_props, count_negative=True, internal_num_rep=None, color_config=None, norm_method='min-max'):
    """
    Input:
    tree: A tree structure with nodes, each having properties.
    profiling_props: A list of property names to be processed for each leaf in the tree.

    Output:
    A dictionary of matrix representation of the tree leaves and their properties.
    A sorted dictionary mapping property values to their corresponding colors.
    """
    def flatten(l):
        return [item for sublist in l for item in sublist]

    def min_max_normalize(value, minval, maxval):
        if maxval - minval == 0:
            return 0
        else:
            return (value - minval) / (maxval - minval)

    def mean_normalize(value, mean_val, minval, maxval):
        if maxval - minval == 0:
            return 0
        else:
            return (value - mean_val) / (maxval - minval)

    def z_score_normalize(value, mean_val, std_val):
        if std_val == 0:
            return 0
        else:
            return (value - mean_val) / std_val

    def _get_color(search_value, color_dict, norm_method='min-max'):
        num = len(color_dict)
        search_value = float(search_value)
        if norm_method == "min-max":
            normalized_value = min_max_normalize(search_value)
            index_values = np.linspace(0, 1, num)
        elif norm_method == "mean":
            normalized_value = mean_normalize(search_value)
            index_values = np.linspace(-1, 1, num)
        elif norm_method == "zscore":
            normalized_value = z_score_normalize(search_value)
            index_values = np.linspace(-3, 3, num)
        else:
            logger.error("Unsupported normalization method.")
            sys.exit(1)
        index = np.abs(index_values - normalized_value).argmin() + 1
        #index = np.abs(index_values - search_value).argmin() + 1
        return color_dict.get(index, "")

    def parse_color_config(color_config, profiling_props, all_props_wildcard, minval, maxval):
        gradientscolor = None
        nan_color = '#EBEBEB'
        max_color = 'red'
        min_color = 'white'
        mid_color = None
        value2color = {}
        
        if color_config.get(all_props_wildcard) is not None:
            prop_config = color_config[all_props_wildcard]
            if 'value2color' in prop_config and prop_config['value2color']:
                value2color = prop_config['value2color']
                value2color = {float(key): value for key, value in value2color.items()}
            if 'detail2color' in prop_config and prop_config['detail2color']:
                detail2color = prop_config.get('detail2color', {})
                temp_min_color, temp_min_val = detail2color.get('color_min', (None, None))
                temp_max_color, temp_max_val = detail2color.get('color_max', (None, None))
                temp_mid_color, temp_mid_val = detail2color.get('color_mid', (None, None))
                temp_none_color, _ = detail2color.get('color_nan', (None, None))

                if temp_max_color:
                    max_color = temp_max_color
                if temp_min_color:
                    min_color = temp_min_color
                if temp_mid_color:
                    mid_color = temp_mid_color
                if temp_none_color:
                    nan_color = temp_none_color
                if temp_min_val:
                    minval = float(temp_min_val)
                if temp_max_val:
                    maxval = float(temp_max_val)

                gradientscolor = utils.build_custom_gradient(20, min_color, max_color, mid_color)
        
        if profiling_props:
            for profiling_prop in profiling_props:
                if color_config.get(profiling_prop) is not None:
                    prop_config = color_config[profiling_prop]
                    if 'value2color' in prop_config and prop_config['value2color']:
                        value2color = prop_config['value2color']
                        value2color = {float(key): value for key, value in value2color.items()}
                    if 'detail2color' in prop_config and prop_config['detail2color']:
                        detail2color = prop_config.get('detail2color', {})
                        temp_min_color, temp_min_val = detail2color.get('color_min', (None, None))
                        temp_max_color, temp_max_val = detail2color.get('color_max', (None, None))
                        temp_mid_color, temp_mid_val = detail2color.get('color_mid', (None, None))

                        if temp_max_color:
                            max_color = temp_max_color
                        if temp_min_color:
                            min_color = temp_min_color
                        if temp_mid_color:
                            mid_color = temp_mid_color
                        if temp_min_val:
                            minval = float(temp_min_val)
                        if temp_max_val:
                            maxval = float(temp_max_val)
                        gradientscolor = utils.build_custom_gradient(20, min_color, max_color, mid_color)
        
        return value2color, gradientscolor, minval, maxval, nan_color

    def process_color_configuration(node2matrix, profiling_props=None):
        
        value2color = {}
        all_props_wildcard = '*'

        # Get color configuration
        all_values_raw = list(set(flatten([sublist for sublist in node2matrix.values()])))
        all_values = sorted(list(filter(lambda x: x is not None and not math.isnan(x), all_values_raw)))

        if not count_negative:
            positive_values = sorted(list(filter(lambda x: x is not None and not math.isnan(x) and x >= 0, all_values)))
            minval, maxval = min(positive_values), max(positive_values)
            mean_val = np.mean(positive_values)
            std_val = np.std(positive_values)
        else:
            minval, maxval = min(all_values), max(all_values)
            mean_val = np.mean(all_values)
            std_val = np.std(all_values)

        if color_config:
            value2color, gradientscolor, minval, maxval, nan_color = parse_color_config(color_config, profiling_props, all_props_wildcard, minval, maxval)
        else:
            gradientscolor = None
            nan_color = '#EBEBEB'

        if not gradientscolor:
            if norm_method == 'min-max':
                gradientscolor = utils.build_color_gradient(20, colormap_name="Reds")
            else: # "mean" "zscore"
                gradientscolor = utils.build_color_gradient(20, colormap_name="coolwarm")
        
        num = len(gradientscolor)
        for search_value in all_values_raw:
            if search_value is None:
                value2color[search_value] = nan_color
            elif math.isnan(search_value):
                value2color[search_value] = nan_color
            #value2color[search_value] = _get_color(search_value, gradientscolor, norm_method)
            else:
                search_value = float(search_value)
                if search_value not in value2color:
                    if not count_negative and search_value < 0:
                        value2color[search_value] = nan_color
                    else:
                        if norm_method == "min-max":
                            normalized_value = min_max_normalize(search_value, minval, maxval)
                            index_values = np.linspace(0, 1, num)
                        elif norm_method == "mean":
                            normalized_value = mean_normalize(search_value, mean_val, minval, maxval)
                            index_values = np.linspace(-1, 1, num)
                        elif norm_method == "zscore":
                            normalized_value = z_score_normalize(search_value, mean_val, std_val)
                            index_values = np.linspace(-3, 3, num)
                        else:
                            logger.error("Unsupported normalization method.")
                            sys.exit(1)
                        index = np.abs(index_values - normalized_value).argmin() + 1
                        value2color[search_value] = gradientscolor.get(index, "")
        return minval, maxval, value2color

    node2matrix_single = {}
    node2matrix_list = {prop: {} for prop in profiling_props}

    single_props = set()
    list_props = set()
    
    for node in tree.traverse():
        node2matrix_single[node.name] = []
        for profiling_prop in profiling_props:
            if node.is_leaf:
                prop_value = node.props.get(profiling_prop)
                if prop_value is not None:
                    if isinstance(prop_value, list):
                        list_props.add(profiling_prop)
                        prop_value = list(map(float, prop_value))
                        if node.name not in node2matrix_list[profiling_prop]:
                            node2matrix_list[profiling_prop][node.name] = []
                        node2matrix_list[profiling_prop][node.name] = prop_value
                    else:
                        single_props.add(profiling_prop)
                        node2matrix_single[node.name].append(float(prop_value))
                else:
                    node2matrix_single[node.name].append(None)
                    if node.name not in node2matrix_list[profiling_prop]:
                        node2matrix_list[profiling_prop][node.name] = []
                    node2matrix_list[profiling_prop][node.name].append(None)
            else:
                if internal_num_rep != 'none':
                    representative_prop = utils.add_suffix(profiling_prop, internal_num_rep)
                    prop_value = node.props.get(representative_prop)
                    if prop_value is not None:
                        if isinstance(prop_value, list):
                            list_props.add(profiling_prop)
                            prop_value = list(map(float, prop_value))
                            if node.name not in node2matrix_list[profiling_prop]:
                                node2matrix_list[profiling_prop][node.name] = []
                            node2matrix_list[profiling_prop][node.name] = prop_value
                        else:
                            node2matrix_single[node.name].append(float(prop_value))
                    else:
                        node2matrix_single[node.name].append(None)
                        if node.name not in node2matrix_list[profiling_prop]:
                            node2matrix_list[profiling_prop][node.name] = []
                        node2matrix_list[profiling_prop][node.name].append(None)

    # Process single values
    if single_props:
        minval_single, maxval_single, value2color_single = process_color_configuration(node2matrix_single, profiling_props)
    else:
        minval_single, maxval_single, value2color_single = None, None, None
    
    if list_props:
        # Process list values for each profiling_prop
        results_list = {}
        for prop in profiling_props:
            if prop in list_props:
                minval_list, maxval_list, value2color_list = process_color_configuration(node2matrix_list[prop])
                results_list[prop] = (node2matrix_list[prop], minval_list, maxval_list, value2color_list)
            else:
                results_list[prop] = (None, None, None, None)
    else:
        results_list = None

    return node2matrix_single, minval_single, maxval_single, value2color_single, results_list, list_props, single_props

def binary2matrix(tree, profiling_props, color_config=None):
    """
    Input:
    tree: A tree structure with nodes, each having properties.
    profiling_props: A list of property names to be processed for each leaf in the tree.

    Output:
    A dictionary of matrix representation of the tree leaves and their properties.
    A sorted dictionary mapping property values to their corresponding colors.
    """
    is_list = False
    binary2color = {True: 1, False: 0}
    node2matrix = {}
    counter_separator = '||'
    all_props_wildcard = '*'
    value2color = {}

    for node in tree.traverse():
        node2matrix[node.name] = []
        for profiling_prop in profiling_props:
            if node.is_leaf:
                prop_value = node.props.get(profiling_prop)  
                if prop_value is not None:  
                    if isinstance(prop_value, list):  # Check if the property value is a list
                        is_list = True  # Set is_array to True upon finding the first list
                        for array_element in prop_value:
                            node2matrix[node.name].append(binary2color.get(utils.str2bool(array_element)))
                    else:  # If not a list, directly handle the single value case
                        node2matrix[node.name].append(binary2color.get(utils.str2bool(prop_value)))
                else:  # If prop_value is None, append None
                    node2matrix[node.name].append(None)
            
            else: # for internal nodes parse counter of True/Total percentage 
                representative_prop = utils.add_suffix(profiling_prop, "counter")
                if node.props.get(representative_prop):
                    ratio = utils.counter2ratio(node, representative_prop)
                    node2matrix[node.name].append(ratio)

    if color_config:
        if color_config.get(all_props_wildcard) is not None:
            prop_config = color_config[all_props_wildcard]
            if 'value2color' in prop_config and prop_config['value2color']:
                value2color = prop_config['value2color']
                value2color = {float(key): value for key, value in value2color.items()}
        if profiling_props:
            for profiling_prop in profiling_props:
                if color_config.get(profiling_prop) is not None:
                    prop_config = color_config[profiling_prop]
                    if 'value2color' in prop_config and prop_config['value2color']:
                        value2color = prop_config['value2color']
                        value2color = {float(key): value for key, value in value2color.items()}
    
    gradientscolor = utils.build_color_gradient(20, colormap_name='Reds')
    
    # get color for binary value 0 to 1
    all_values_raw = list(set(utils.flatten([sublist for sublist in node2matrix.values()])))
    all_values = sorted(list(filter(lambda x: x is not None and not math.isnan(x), all_values_raw)))
    num = len(gradientscolor)
    index_values = np.linspace(0, 1, num) # binary value 0 to 1
    for search_value in all_values:
        if search_value not in value2color:
            index = np.abs(index_values - search_value).argmin() + 1
            value2color[search_value] = gradientscolor[index]
    
    return node2matrix, value2color, is_list

def float2matrix(tree, profiling_props, count_negative=True):
    """
    Input:
    tree: A tree structure with nodes, each having properties.
    profiling_props: A list of property names to be processed for each leaf in the tree.
    
    Output:
    A fasta text output of matrix representation of the tree leaves and their properties in color codes.
    A sorted dictionary mapping property values to their corresponding colors in the gradient.
    """
    def process_value(value):
        """Process a single value or list of values to float."""
        if isinstance(value, list):
            return [float(v) if v is not None else None for v in value]
        else:
            return float(value) if value is not None else None

    gradients = [
    'a', 'b', 'c',
    'd', 'e', 'f',
    'g', 'h', 'i',
    'j', 'k', 'l',
    'm', 'n', 'o',
    'p', 'q', 'r', 
    's', 't'
    ] #white to red
    
    absence_color = '-'
    negative_color = 'x'
    leaf2matrix = {}
    for node in tree.traverse():
        if node.is_leaf:
            #leaf2matrix[node.name] = [float(node.props.get(prop)) if node.props.get(prop) is not None else None for prop in profiling_props]
            leaf2matrix[node.name] = []
            for prop in profiling_props:
                value = node.props.get(prop)
                processed_value = process_value(value)
                leaf2matrix[node.name].extend(processed_value if isinstance(processed_value, list) else [processed_value])

    value2color = {}
    all_values = list(set(utils.flatten(leaf2matrix.values())))
    if count_negative:
        all_values = sorted([x for x in all_values if x is not None and not math.isnan(x)])
    else:
        all_values = sorted([x for x in all_values if x is not None and not math.isnan(x) and x >= 0])
    
    maxval, minval = max(all_values), min(all_values)
    num = len(gradients)
    values = np.linspace(minval, maxval, num)

    matrix = ''
    for leaf, prop in leaf2matrix.items():
        matrix += '\n' + '>' + leaf + '\n'
        for search_value in prop:
            if search_value is not None:
                if not count_negative and search_value < 0:
                    matrix += negative_color
                    value2color[search_value] = negative_color
                else:
                    index = np.abs(values - search_value).argmin()
                    matrix += gradients[index]
                    value2color[search_value] = gradients[index]
            else:
                matrix += '-'
                value2color[search_value] = absence_color

    sorted_value2color = OrderedDict(sorted((k, v) for k, v in value2color.items() if k is not None))
    return matrix, sorted_value2color

def str2matrix(tree, profiling_props):
    """
    Input:
    tree: A tree structure with nodes, each having properties.
    profiling_props: A list of property names to be processed for each leaf in the tree.
    
    Output:
    A fasta text output of matrix representation of the tree leaves and their properties in color codes.
    A sorted dictionary mapping property values to their corresponding colors in the 20 amino acid color codes.
    """
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
    leaf2matrix = {}
    for node in tree.traverse():
        if node.is_leaf:
            leaf2matrix[node.name] = [node.props.get(prop) if node.props.get(prop) is not None else None for prop in profiling_props]

    value2color = {}
    all_values = sorted(set(utils.flatten(leaf2matrix.values())))
    for i, val in enumerate(all_values):
        if val != 'NaN':
            value2color[val] = aa[i % len(aa)]  # Use modulo to avoid out-of-range errors
        else:
            value2color[val] = absence_color

    matrix = ''
    for leaf, prop in leaf2matrix.items():
        matrix += '\n' + '>' + leaf + '\n' + ''.join([value2color.get(item, '-') for item in prop])

    return matrix, value2color

def single2matrix(tree, profiling_prop):
    precence_color = '#E60A0A' # #E60A0A red
    absence_color = '#EBEBEB' # #EBEBEB lightgrey
    all_categorical_values = sorted(list(set(utils.flatten(utils.tree_prop_array(tree, profiling_prop)))), key=lambda x: (x != 'NaN', x))
    
    node2matrix = {}
    for node in tree.traverse():
        if node.is_leaf:
            # Leaf node processing for presence/absence
            if node.props.get(profiling_prop):
                node2matrix[node.name] = [1 if val == node.props.get(profiling_prop) else 0 for val in all_categorical_values]
        else:
            # Internal node processing to add a counter of True/Total percentage
            representative_prop = utils.add_suffix(profiling_prop, "counter")
            if node.props.get(representative_prop):
                ratios = utils.categorical2ratio(node, representative_prop, all_categorical_values)
                node2matrix[node.name] = ratios  # Extend the ratio across all possible values
    
    #get color
    gradientscolor = utils.build_color_gradient(20, colormap_name='Reds')
    #value2color = {}
    value2color = {1: precence_color, 0: absence_color}
    # get color for binary value 0 to 1
    all_values_raw = list(set(utils.flatten([sublist for sublist in node2matrix.values()])))
    all_values = sorted(list(filter(lambda x: x is not None and not math.isnan(x), all_values_raw)))
    num = len(gradientscolor)
    index_values = np.linspace(0, 1, num) # binary value 0 to 1
    for search_value in all_values:
        if search_value not in value2color:
            index = np.abs(index_values - search_value).argmin() + 1
            value2color[search_value] = gradientscolor[index]
    
    #value2color = {1: precence_color, 0: absence_color}
    return node2matrix, value2color, all_categorical_values

# def single2matrix(tree, profiling_prop):
#     precence_color = '#E60A0A' # #E60A0A red
#     absence_color = '#EBEBEB' # #EBEBEB lightgrey
#     all_values = sorted(list(set(utils.flatten(utils.tree_prop_array(tree, profiling_prop)))), key=lambda x: (x != 'NaN', x))
#     leaf2matrix = {}
#     for leaf in tree.leaves():
#         leaf2matrix[leaf.name] = []
#         for val in all_values:
#             if val == leaf.props.get(profiling_prop):
#                 leaf2matrix[leaf.name].append(1)
#             else:
#                 leaf2matrix[leaf.name].append(0)
#     value2color = {1: precence_color, 0: absence_color}
#     return leaf2matrix, value2color, all_values

def single2profile(tree, profiling_prop):
    all_values = sorted(list(set(utils.flatten(utils.tree_prop_array(tree, profiling_prop)))), key=lambda x: (x != 'NaN', x))
    presence = 'p' # #E60A0A red
    absence = 'z' # #EBEBEB lightgrey
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

def multiple2matrix(tree, profiling_prop, prop2type=None, color_config=None, eteformat_flag=False):
    precence_color = '#E60A0A'  # red
    absence_color = '#EBEBEB'   # grey

    # Determine the data type of the profiling property
    data_type = prop2type.get(profiling_prop)

    # Get all categorical values based on whether data_type is a list and eteformat_flag
    if data_type and data_type == list:
        tree_prop_array = utils.tree_prop_array(tree, profiling_prop, leaf_only=True, list_type=not eteformat_flag)
    else:
        tree_prop_array = utils.tree_prop_array(tree, profiling_prop, leaf_only=True)

    all_categorical_values = sorted(list(set(utils.flatten(tree_prop_array))), key=lambda x: (x != 'NaN', x))

    # Create node to matrix mappings
    node2matrix = {}
    for node in tree.traverse():
        node_prop = node.props.get(profiling_prop)
        if node.is_leaf and node_prop:
            node2matrix[node.name] = [1 if val in node_prop else 0 for val in all_categorical_values]
        else:
            representative_prop = utils.add_suffix(profiling_prop, "counter")
            if node.props.get(representative_prop):
                ratios = utils.categorical2ratio(node, representative_prop, all_categorical_values)
                node2matrix[node.name] = ratios

    # Build a color gradient for binary values
    gradientscolor = utils.build_color_gradient(20, colormap_name='Reds')
    value2color = {1: precence_color, 0: absence_color}

    # Get unique values from node2matrix and sort non-NaN values
    all_values_raw = list(set(utils.flatten(node2matrix.values())))
    all_values = sorted([x for x in all_values_raw if x is not None and not math.isnan(x)])

    # Map non-binary values to gradient colors
    num = len(gradientscolor)
    index_values = np.linspace(0, 1, num)
    for search_value in all_values:
        if search_value not in value2color:
            index = np.abs(index_values - search_value).argmin() + 1
            value2color[search_value] = gradientscolor[index]

    return node2matrix, value2color, all_categorical_values


def multiple2profile(tree, profiling_prop):
    all_values = sorted(list(set(utils.flatten(utils.tree_prop_array(tree, profiling_prop)))), key=lambda x: (x != 'NaN', x))
    presence = 'p' # #E60A0A red
    absence = 'z' # #EBEBEB lightgrey
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