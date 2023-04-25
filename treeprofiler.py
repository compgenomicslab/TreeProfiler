#!/usr/bin/env python

from argparse import ArgumentParser
import argparse
import sys
import tree_annotate
import tree_plot


__author__ = 'Ziqi DENG'
__license__ = "GPL v2"
__email__ = 'dengziqi1234@gmail.com'
__version__ = '0.0.1'
__date__ = '01-11-2022'
__description__ = ('A program for profiling metadata on target '
                    'tree and conduct summary analysis')

def populate_main_args(main_args_p):
    """
    Parse the input parameters
    Return the parsed arguments.
    """
    # input parameters group
    group = main_args_p.add_argument_group(title='SOURCE TREE INPUT',
        description="Source tree input parameters")
    group.add_argument('-t', '--tree',
        type=str,
        required=False,
        help="Input tree, .nw file, customized tree input")
    group.add_argument('--annotated_tree',
        default=False,
        action='store_true',
        required=False,
        help="input tree already annotated by treeprofiler")
    group.add_argument('--tree_type',
        type=str,
        default='newick',
        required=False,
        help="statistic calculation to perform for numerical data in internal nodes, [newick, ete]")
    
    group.add_argument('--prop2type',
        type=str,
        required=False,
        help="config tsv file where determine the datatype of target properties, if your input tree type is .ete, it's note necessary")
    
    group = main_args_p.add_argument_group(title='Pruning parameters',
        description="Auto pruning parameters")

    group.add_argument('--rank_limit',
        type=str,
        required=False,
        help="TAXONOMIC_LEVEL prune annotate tree by rank limit")
    group.add_argument('--pruned_by', 
        type=str,
        required=False,
        action='append',
        help='target tree pruned by customized conditions')
    
    # group = main_args_p.add_argument_group(title='OUTPUT options',
    #     description="")
    # group.add_argument('--ete4out',
    #     default=False,
    #     action='store_true',
    #     help="export intermediate tree in ete4")
    # group.add_argument('-o', '--outdir',
    #     type=str,
    #     required=False,
    #     help="output annotated tree")
    # group.add_argument('--outtsv',
    #     type=str,
    #     required=False,
    #     help="output annotated tsv file")
    # group.add_argument('--out_colordict',
    #     action="store_true", 
    #     required=False,
    #     help="print color dictionary of each property")

    #args = parser.parse_args()
    #return args

def populate_annotate_args(annotate_args_p):
    group = annotate_args_p.add_argument_group(title='METADATA TABLE parameters',
        description="Input parameters of METADATA")
    # group.add_argument('-t', '--tree',
    #     type=str,
    #     required=False,
    #     help="Input tree, .nw file, customized tree input")
    group.add_argument('-d', '--metadata',
        required=False,
        help="<metadata.csv> .csv, .tsv. mandatory input",
        type=lambda s: [item for item in s.split(',')])
    group.add_argument('--no_colnames',
        default=False,
        action='store_true',
        required=False,
        help="metadata table doesn't contain columns name")
    group.add_argument('--text_prop',
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as categorical data",
        type=lambda s: [item for item in s.split(',')])
    group.add_argument('--multiple_text_prop',
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as categorical data which contains more than one value and seperate by ',' such as GO:0000003,GO:0000902,GO:0000904,GO:0003006",
        type=lambda s: [item for item in s.split(',')])
    group.add_argument('--num_prop',
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as numerical data",
        type=lambda s: [item for item in s.split(',')])
    group.add_argument('--bool_prop',
        required=False,
        help="<col1,col2> names, column index or index range of columns which need to be read as boolean data",
        type=lambda s: [item for item in s.split(',')])
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
    
    group.add_argument('--emapper_annotations',
        type=str,
        required=False,
        help="out.emapper.annotations")
    group.add_argument('--emapper_pfam',
        type=str,
        required=False,
        help="out.emapper.pfams")
    group.add_argument('--emapper_smart',
        type=str,
        required=False,
        help="out.emapper.smart")
    group.add_argument('--alignment',
        type=str,
        required=False,
        help="Sequence alignment, .fasta format")


    group = annotate_args_p.add_argument_group(title='Annotation arguments',
        description="Annotation parameters")
    group.add_argument('--taxonomic_profile',
        default=False,
        action='store_true',
        required=False,
        help="Determine if you need taxonomic annotation on tree")
    group.add_argument('--num_stat',
        default='all',
        type=str,
        required=False,
        help="statistic calculation to perform for numerical data in internal nodes, [all, sum, avg, max, min, std] ")  
    group.add_argument('--counter_stat',
        default='raw',
        type=str,
        required=False,
        help="statistic calculation to perform for categorical data in internal nodes, raw count or in percentage [raw, relative] ")  
    
    group = annotate_args_p.add_argument_group(title='OUTPUT options',
        description="")
    group.add_argument('--ete4out',
        default=False,
        action='store_true',
        help="export intermediate tree in ete4")
    group.add_argument('-o', '--outdir',
        type=str,
        required=False,
        help="output annotated tree")
    group.add_argument('--outtsv',
        type=str,
        required=False,
        help="output annotated tsv file")

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

def main():
    main_args_p = argparse.ArgumentParser(description=
        "treeprofiler.py (ver. "+__version__+
        " of "+__date__+")." + __description__+ " Authors: "+
        __author__+" ("+__email__+")",
        formatter_class=argparse.RawTextHelpFormatter, add_help=False)

    populate_main_args(main_args_p)

    parser = argparse.ArgumentParser(description="this is tree profiler ",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser = parser.add_subparsers(title="AVAILABLE PROGRAMS")
    
    ## - ANNOTATE - 
    annotate_args_p = subparser.add_parser('annotate', parents=[main_args_p],
                                            description='annotate tree')
    populate_annotate_args(annotate_args_p)
    annotate_args_p.set_defaults(func=tree_annotate.run)

    ## - PLOT - 
    plot_args_p = subparser.add_parser('plot', parents=[main_args_p],
                                            description='annotate plot')
    poplulate_plot_args(plot_args_p)
    plot_args_p.set_defaults(func=tree_plot.run )

    ## - RUN -
    args = parser.parse_args(sys.argv[1:])
    args.func(args)
    
if __name__ == '__main__':
    main()