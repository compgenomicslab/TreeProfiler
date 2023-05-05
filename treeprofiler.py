#!/usr/bin/env python

import sys
import argparse as ap

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
        help="Input tree, .nw file, customized tree input")
    group.add_argument('--annotated_tree',
        default=False,
        action='store_true',
        help="input tree already annotated by treeprofiler")
    group.add_argument('--tree_type',
        type=str,
        default='newick',
        help="statistic calculation to perform for numerical data in internal nodes, [newick, ete]")

    group.add_argument('--prop2type',
        type=str,
        help="config tsv file where determine the datatype of target properties, if your input tree type is .ete, it's note necessary")

    group = main_args_p.add_argument_group(title='Pruning parameters',
        description="Auto pruning parameters")

    group.add_argument('--rank_limit',
        type=str,
        help="TAXONOMIC_LEVEL prune annotate tree by rank limit")
    group.add_argument('--pruned_by',
        type=str,
        action='append',
        help='target tree pruned by customized conditions')

    # NOTE(JBC): I would not use the "required=False" since it is the default
    # and makes things more verbose.

    # NOTE(JBC): I would remove all the commented code.

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

def main():
    main_args_p = ap.ArgumentParser(description=
        "treeprofiler.py (ver. "+__version__+
        " of "+__date__+")." + __description__+ " Authors: "+
        __author__+" ("+__email__+")",
        formatter_class=ap.RawTextHelpFormatter, add_help=False)

    populate_main_args(main_args_p)

    parser = ap.ArgumentParser(description="this is tree profiler ",
                               formatter_class=ap.RawDescriptionHelpFormatter)
    subparser = parser.add_subparsers(title="AVAILABLE PROGRAMS")

    ## - ANNOTATE -
    annotate_args_p = subparser.add_parser('annotate', parents=[main_args_p],
                                            description='annotate tree')
    tree_annotate.populate_annotate_args(annotate_args_p)
    annotate_args_p.set_defaults(func=tree_annotate.run)

    ## - PLOT -
    plot_args_p = subparser.add_parser('plot', parents=[main_args_p],
                                       description='annotate plot')
    tree_plot.poplulate_plot_args(plot_args_p)
    plot_args_p.set_defaults(func=tree_plot.run)

    ## - RUN -
    args = parser.parse_args(sys.argv[1:])
    args.func(args)

if __name__ == '__main__':
    main()
