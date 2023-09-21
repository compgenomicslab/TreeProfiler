#!/usr/bin/env python

import sys
import argparse as ap

from treeprofiler import tree_annotate
from treeprofiler import tree_plot


__author__ = 'anonymous'
__license__ = "GPL v3"
__email__ = 'anonymous@gmail.com'
__version__ = '1.1.0'
__date__ = '14-07-2023'
__description__ = ('A program for profiling metadata on target '
                    'tree and conduct summary analysis')

def populate_main_args(main_args_p):
    """
    Parse the input parameters
    """

    # input parameters group
    group = main_args_p.add_argument_group(title='SOURCE TREE INPUT',
        description="Source tree input parameters")
    group.add_argument('-t', '--tree',
        type=str,
        required=True,
        help="Input tree, .nw file, customized tree input")
    group.add_argument('--annotated-tree',
        default=False,
        action='store_true',
        help="input tree already annotated by treeprofiler if you want to skip the annotate part.")
    group.add_argument('--internal-parser',
        default="name",
        choices=["name", "support"],
        type=str,
        required=False,
        help="To specify how to interpret internal nodes in newick format. [default: name]")
    group.add_argument('--input-type',
        type=str,
        default="ete",
        choices=["newick", "ete"],
        help="Specify input tree format. [newick, ete]. [default: ete]")
    group.add_argument('--prop2type',
        type=str,
        help="config tsv file where determine the datatype of target properties, if your input tree type is .ete, it's note necessary")

    group = main_args_p.add_argument_group(title='Pruning parameters',
        description="Auto pruning parameters")
    group.add_argument('--rank-limit',
        type=str,
        help="TAXONOMIC_LEVEL prune annotate tree by rank limit")
    group.add_argument('--pruned-by',
        type=str,
        action='append',
        help='target tree pruned by customized conditions, such as --pruned_by "name contains FALPE"')

def main():
    desc = (f'treeprofiler.py (ver. {__version__} of {__date__}).'
            f'{__description__} Authors: {__author__} ({__email__})')
    main_args_p = ap.ArgumentParser(description=desc,
                                    formatter_class=ap.RawTextHelpFormatter,
                                    add_help=False)

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
    if len(sys.argv[1:]) < 1:
        print(parser.print_usage())
        return
    
    args = parser.parse_args(sys.argv[1:])
    args.func(args)

if __name__ == '__main__':
    main()
