#!/usr/bin/env python3

import os, math, re
import logging
import sys
import time
import random
import csv
import tarfile

from io import StringIO
from itertools import islice
from collections import defaultdict, Counter
import numpy as np

from scipy import stats

from ete4.parser.newick import NewickError
from ete4.core.seqgroup import SeqGroup
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from treeprofiler.src.utils import (
    validate_tree, TreeFormatError, get_internal_parser,
    taxatree_prune, conditional_prune,
    children_prop_array, children_prop_array_missing, 
    flatten, get_consensus_seq, add_suffix, clear_extra_features)
from treeprofiler.src.phylosignal import run_acr_discrete, run_delta
from treeprofiler.src.ls import run_ls
from treeprofiler.src import b64pickle

from multiprocessing.pool import ThreadPool
from multiprocessing import Pool

DESC = "annotate tree"

TAXONOMICDICT = {# start with leaf name
                'rank': str,
                'sci_name': str,
                'taxid': str,
                'lineage':str,
                'named_lineage': str,
                'evoltype': str,
                'dup_sp': str,
                'dup_percent': float,
                }

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def populate_annotate_args(parser):
    gmeta = parser.add_argument_group(
        title='METADATA TABLE parameters',
        description="Input parameters of METADATA")
    add = gmeta.add_argument
    add('-d', '--metadata', nargs='+',
        help="<metadata.csv> .csv, .tsv. mandatory input")
    add('-sep', '--metadata_sep', default='\t',
        help="column separator of metadata table [default: \\t]")
    add('--no_colnames', action='store_true',
        help="metadata table doesn't contain columns name")
    add('--aggregate-duplicate', action='store_true',
        help="treeprofiler will aggregate duplicated metadata to a list as a property if metadata contains duplicated row")
    add('--text-prop', nargs='+',
        help=("<col1> <col2> names, column index or index range of columns which "
              "need to be read as categorical data"))
    add('--multiple-text-prop', nargs='+',
        help=("<col1> <col2> names, column index or index range of columns which "
              "need to be read as categorical data which contains more than one"
              " value and seperate by ',' such "
              "as GO:0000003,GO:0000902,GO:0000904,GO:0003006"))
    add('--num-prop', nargs='+',
        help=("<col1> <col2> names, column index or index range of columns which "
              "need to be read as numerical data"))
    add('--bool-prop', nargs='+',
        help=("<col1> <col2> names, column index or index range of columns which "
              "need to be read as boolean data"))
    add('--text-prop-idx', nargs='+',
        help="1 2 3 or [1-5] index of columns which need to be read as categorical data")
    add('--num-prop-idx', nargs='+',
        help="1 2 3 or [1-5] index columns which need to be read as numerical data")
    add('--bool-prop-idx', nargs='+',
        help="1 2 3 or [1-5] index columns which need to be read as boolean data")
    add('--acr-discrete-columns', nargs='+',
        help=("<col1> <col2> names to perform acr analysis for discrete traits"))
    # add('--acr-continuous-columns', nargs='+',
    #     help=("<col1> <col2> names to perform acr analysis for continuous traits"))
    add('--ls-columns', nargs='+',
        help=("<col1> <col2> names to perform lineage specificity analysis"))
    # add('--taxatree',
    #     help=("<kingdom|phylum|class|order|family|genus|species|subspecies> "
    #           "reference tree from taxonomic database"))
    add('--taxadb', type=str.upper,
        choices=['NCBI', 'GTDB'],
        help="<NCBI|GTDB> for taxonomic profiling or fetch taxatree")
    add('--taxa-dump', type=str,
        help='Path to taxonomic database dump file for specific version, such as https://github.com/etetoolkit/ete-data/raw/main/gtdb_taxonomy/gtdblatest/gtdb_latest_dump.tar.gz')
    add('--taxon-column',
        help="Activate taxonomic annotation using <col1> name of columns which need to be read as taxon data. \
            Unless taxon data in leaf name, please use 'name' as input such as --taxon-column name")
    add('--taxon-delimiter', default=None,
        help="delimiter of taxa columns. [default: None]")
    add('--taxa-field', type=int, default=0,
        help="field of taxa name after delimiter. [default: 0]")
    add('--emapper-annotations',
        help="attach eggNOG-mapper output out.emapper.annotations")
    add('--emapper-pfam',
        help="attach eggNOG-mapper pfam output out.emapper.pfams")
    add('--emapper-smart',
        help="attach eggNOG-mapper smart output out.emapper.smart")
    add('--alignment',
        help="Sequence alignment, .fasta format")

    annotation_group = parser.add_argument_group(title='Internal nodes annotation arguments',
        description="Annotation parameters")
    annotation_group.add_argument('--column-summary-method', 
        nargs='+',
        required=False,
        help="Specify summary method for individual columns in the format ColumnName=Method")
    annotation_group.add_argument('--num-stat',
        default='all',
        choices=['all', 'sum', 'avg', 'max', 'min', 'std', 'none'],
        type=str,
        required=False,
        help="statistic calculation to perform for numerical data in internal nodes, [all, sum, avg, max, min, std, none]. If 'none' was chosen, numerical properties won't be summarized nor annotated in internal nodes. [default: all]")  
    annotation_group.add_argument('--counter-stat',
        default='raw',
        choices=['raw', 'relative', 'none'],
        type=str,
        required=False,
        help="statistic calculation to perform for categorical data in internal nodes, raw count or in percentage [raw, relative, none]. If 'none' was chosen, categorical and boolean properties won't be summarized nor annotated in internal nodes [default: raw]")  
    
    acr_group = parser.add_argument_group(title='Ancestral Character Reconstruction arguments',
        description="ACR parameters")
    acr_group.add_argument('--prediction-method',
        default='MPPA',
        choices=['MPPA','MAP','JOINT','DOWNPASS','ACCTRAN','DELTRAN','COPY','ALL','ML','MP'],
        type=str,
        required=False,
        help="prediction method for ACR discrete analysis [default: MPPA]"
        )
    acr_group.add_argument('--model',
        default='F81',
        choices=['JC','F81','EFT','HKY','JTT','CUSTOM_RATES'],
        type=str,
        required=False,
        help="Evolutionary model for ML methods in ACR discrete analysis [default: F81]"
        )
    acr_group.add_argument('--threads',
        default=4,
        type=int,
        required=False,
        help="Number of threads to use for annotation [default: 4]")
    delta_group = parser.add_argument_group(title='Ancestral Character Reconstruction arguments',
        description="Delta statistic parameters")
    delta_group.add_argument('--delta-stats',
        action='store_true',
        required=False,
        help="Calculate delta statistic for discrete traits in ACR analysis, ONLY for MPPA or MAP prediction method.[default: False]"
    )
    delta_group.add_argument('--ent_type',
        default='SE',
        choices=['LSE', 'SE', 'GINI'],
        type=str,
        required=False,
        help="Entropy method to measure the degree of phylogenetic signal between discrete trati and phylogeny. \
            [default: SE] for Shannon Entropy, other options are GINI for Gini impurity and LSE for Linear Shannon Entropy."
    )
    delta_group.add_argument('--iteration',
        default=10000,
        type=int,
        required=False,
        help="Number of iterations for delta statistic calculation. [default: 100]"
    )
    delta_group.add_argument('--lambda0', 
        type=float, 
        default=0.1, 
        help='Rate parameter of the delta statistic calculation.'
    )
    delta_group.add_argument('--se', 
        type=float, 
        default=0.5, 
        help='Standard deviation of the delta statistic calculation.')
    delta_group.add_argument('--thin', 
        type=int, 
        default=10, 
        help='Keep only each xth iterate.')
    delta_group.add_argument('--burn', 
        type=int, 
        default=100, 
        help='Burned-in iterates.')
    ls_group = parser.add_argument_group(title='Lineage Specificity Analysis arguments',
        description="ls parameters")
    ls_group.add_argument('--prec-cutoff',
        default=0.95,
        type=float,
        required=False,
        help="Precision cutoff for lineage specificity analysis [default: 0.95]")
    ls_group.add_argument('--sens-cutoff',
        default=0.95,
        type=float,
        required=False,
        help="Sensitivity threshold for lineage specificity analysis [default: 0.95]")
    
    group = parser.add_argument_group(title='OUTPUT options',
        description="")
    group.add_argument('-o', '--outdir',
        type=str,
        required=True,
        help="Directory for annotated outputs.")

def run_tree_annotate(tree, input_annotated_tree=False,
        metadata_dict={}, node_props=[], columns={}, prop2type={},
        emapper_annotations=None,
        text_prop=[], text_prop_idx=[], multiple_text_prop=[], num_prop=[], num_prop_idx=[],
        bool_prop=[], bool_prop_idx=[], prop2type_file=None, alignment=None, emapper_pfam=None,
        emapper_smart=None, counter_stat='raw', num_stat='all', column2method={},
        taxadb='GTDB', taxa_dump=None, taxon_column='name',
        taxon_delimiter='', taxa_field=0, rank_limit=None, pruned_by=None, 
        acr_discrete_columns=None, prediction_method="MPPA", model="F81", 
        delta_stats=False, ent_type="SE", 
        iteration=100, lambda0=0.1, se=0.5, thin=10, burn=100, 
        ls_columns=None, prec_cutoff=0.95, sens_cutoff=0.95, 
        threads=1, outdir='./'):

    total_color_dict = []
    layouts = []
    level = 1 # level 1 is the leaf name

    if emapper_annotations:
        emapper_metadata_dict, emapper_node_props, emapper_columns = parse_emapper_annotations(emapper_annotations)
        metadata_dict.update(emapper_metadata_dict)
        node_props.extend(emapper_node_props)
        columns.update(emapper_columns)
    
        prop2type.update({
            'name': str,
            'dist': float,
            'support': float,
            'seed_ortholog': str,
            'evalue': float,
            'score': float,
            'eggNOG_OGs': list,
            'max_annot_lvl': str,
            'COG_category': str,
            'Description': str,
            'Preferred_name': str,
            'GOs': list,
            'EC':str,
            'KEGG_ko': list,
            'KEGG_Pathway': list,
            'KEGG_Module': list,
            'KEGG_Reaction':list,
            'KEGG_rclass':list,
            'BRITE':list,
            'KEGG_TC':list,
            'CAZy':list,
            'BiGG_Reaction':list,
            'PFAMs':list
        })

    if text_prop:
        text_prop = text_prop
    else:
        text_prop = []

    if multiple_text_prop:
        multiple_text_prop = multiple_text_prop
    else:
        multiple_text_prop = []

    if num_prop:
        num_prop = num_prop
    else:
        num_prop = []

    if bool_prop:
        bool_prop = bool_prop
    else:
        bool_prop = []

    if text_prop_idx:

        index_list = []
        for i in text_prop_idx:

            if i[0] == '[' and i[-1] == ']':
                text_prop_start, text_prop_end = get_range(i)
                for j in range(text_prop_start, text_prop_end+1):
                    index_list.append(j)
            else:
                index_list.append(int(i))

        text_prop = [node_props[index-1] for index in index_list]

    if num_prop_idx:
        index_list = []
        for i in num_prop_idx:
            if i[0] == '[' and i[-1] == ']':
                num_prop_start, num_prop_end = get_range(i)
                for j in range(num_prop_start, num_prop_end+1):
                    index_list.append(j)
            else:
                index_list.append(int(i))

        num_prop = [node_props[index-1] for index in index_list]

    if bool_prop_idx:
        index_list = []
        for i in bool_prop_idx:
            if i[0] == '[' and i[-1] == ']':
                bool_prop_start, bool_prop_end = get_range(i)
                for j in range(bool_prop_start, bool_prop_end+1):
                    index_list.append(j)
            else:
                index_list.append(int(i))

        bool_prop = [node_props[index-1] for index in index_list]

    #rest_prop = []
    if prop2type_file:
        prop2type = {}
        with open(prop2type_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                prop, value = line.split('\t')
                prop2type[prop] = eval(value)
    else:
        # output datatype of each property of each tree node including internal nodes
        if prop2type:
            for key, dtype in prop2type.items():
                if key in text_prop+multiple_text_prop+num_prop+bool_prop:
                    pass
                
                # taxon prop wouldn be process as numerical/text/bool/list value
                elif (taxon_column and key in taxon_column):
                    pass

                else:
                    if dtype == list:
                        multiple_text_prop.append(key)
                    if dtype == str:
                        if key not in multiple_text_prop:
                            text_prop.append(key)
                        else:
                            pass
                    if dtype == float:
                        num_prop.append(key)
                    if dtype == bool:
                        bool_prop.append(key)

        # paramemters can over write the default
        if emapper_annotations:
            text_prop.extend([
                'seed_ortholog',
                'max_annot_lvl',
                'COG_category',
                'EC'
            ])
            num_prop.extend([
                'evalue',
                'score'
            ])
            multiple_text_prop.extend([
                'eggNOG_OGs', 'GOs', 'KEGG_ko', 'KEGG_Pathway',
                'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass',
                'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'])

        
        for prop in text_prop:
            prop2type[prop] = str
  

        for prop in bool_prop:
            prop2type[prop] = bool


        for prop in multiple_text_prop:
            prop2type[prop] = list


        for prop in num_prop:
            prop2type[prop] = float


        prop2type.update({# start with leaf name
                'name':str,
                'dist':float,
                'support':float,
                })

    # load annotations to leaves
    start = time.time()

    # alignment annotation
    if alignment:
        alignment_prop = 'alignment'
        name2seq = parse_fasta(alignment)
        for leaf in tree.leaves():
            leaf.add_prop(alignment_prop, name2seq.get(leaf.name,''))
        prop2type.update({
            alignment_prop:str
            })

    # domain annotation before other annotation
    if emapper_pfam:
        domain_prop = 'dom_arq'
        if not alignment:
            raise ValueError("Please provide alignment file using '--alignment' for pfam annotation.")
        annot_tree_pfam_table(tree, emapper_pfam, alignment, domain_prop=domain_prop)
        prop2type.update({
            domain_prop:str
            })
    if emapper_smart:
        domain_prop = 'dom_arq'
        if not alignment:
            raise ValueError("Please provide alignment file using '--alignment' for smart annotation.")
        annot_tree_smart_table(tree, emapper_smart, alignment, domain_prop=domain_prop)
        prop2type.update({
            domain_prop:str
            })

    # load all metadata to leaf nodes

    # input_annotated_tree determines if input tree is already annotated, if annotated, no longer need metadata
    
    if not input_annotated_tree:
        if taxon_column: # to identify taxon column as taxa property from metadata
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type, taxon_column=taxon_column, taxon_delimiter=taxon_delimiter, taxa_field=taxa_field)
        else:
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type)
    else:
        annotated_tree = tree

    end = time.time()
    print('Time for load_metadata_to_tree to run: ', end - start)

    
    # Ancestor Character Reconstruction analysis
    # data preparation
    if acr_discrete_columns:
        logging.info(f"Performing ACR analysis with Character {acr_discrete_columns} via {prediction_method} method with {model} model.......\n")
        # need to be discrete traits
        discrete_traits = text_prop + bool_prop
        for k in acr_discrete_columns:
            if k not in discrete_traits:
                raise ValueError(f"Character {k} is not discrete trait, please check your input.")

        #############################
        start = time.time()
        acr_discrete_columns_dict = {k: v for k, v in columns.items() if k in acr_discrete_columns}
        acr_results, annotated_tree = run_acr_discrete(annotated_tree, acr_discrete_columns_dict, \
        prediction_method=prediction_method, model=model, threads=threads, outdir=outdir)
        
        # Clear extra features
        clear_extra_features([annotated_tree], prop2type.keys())
        
        # get observed delta
        # only MPPA,MAP method has marginal probabilities to calculate delta
        if delta_stats:
            if prediction_method in ['MPPA', 'MAP']:
                logging.info(f"Performing Delta Statistic analysis with Character {acr_discrete_columns}...\n")
                prop2delta = run_delta(acr_results, annotated_tree, ent_type=ent_type, 
                lambda0=lambda0, se=se, sim=iteration, burn=burn, thin=thin, 
                threads=threads)

                for prop, delta_result in prop2delta.items():
                    logging.info(f"Delta statistic of {prop} is: {delta_result}")
                    tree.add_prop(add_suffix(prop, "delta"), delta_result)

                # start calculating p_value
                logging.info(f"Calculating p_value for delta statistic...")
                # get a copy of the tree
                dump_tree = annotated_tree.copy()
                clear_extra_features([dump_tree], ["name", "dist", "support"])
                
                prop2array = {}
                for prop in columns.keys():
                    prop2array.update(convert_to_prop_array(metadata_dict, prop))
                
                prop2delta_array = get_pval(prop2array, dump_tree, acr_discrete_columns_dict, \
                    iteration=100, prediction_method=prediction_method, model=model,
                    ent_type=ent_type, lambda0=lambda0, se=se, sim=iteration, burn=burn, thin=thin, 
                    threads=threads)

                for prop, delta_array in prop2delta_array.items():
                    p_value = np.sum(np.array(delta_array) > prop2delta[prop]) / len(delta_array)
                    logging.info(f"p_value of {prop} is {p_value}")
                    tree.add_prop(add_suffix(prop, "pval"), p_value)
                    prop2type.update({
                        add_suffix(prop, "pval"): float
                    })

                for prop in acr_discrete_columns:
                    prop2type.update({
                        add_suffix(prop, "delta"): float
                    })
            else:
                logging.warning(f"Delta statistic analysis only support MPPA and MAP prediction method, {prediction_method} is not supported.")

        end = time.time()
        print('Time for acr to run: ', end - start)

    # lineage specificity analysis
    if ls_columns:
        logging.info(f"Performing Lineage Specificity analysis with Character {ls_columns}...\n")
        if all(column in bool_prop for column in ls_columns):
            best_node, qualified_nodes = run_ls(annotated_tree, props=ls_columns, 
            precision_cutoff=prec_cutoff, sensitivity_cutoff=sens_cutoff)
            for prop in ls_columns:
                prop2type.update({
                    add_suffix(prop, "prec"): float,
                    add_suffix(prop, "sens"): float,
                    add_suffix(prop, "f1"): float
                })
        else:
            logging.warning(f"Lineage specificity analysis only support boolean properties, {ls_columns} is not boolean property.")

    # statistic method
    counter_stat = counter_stat #'raw' or 'relative'
    num_stat = num_stat

    # merge annotations depends on the column datatype
    start = time.time()
    
    # choose summary method based on datatype
    for prop in text_prop+multiple_text_prop+bool_prop:
        if not prop in column2method:
            column2method[prop] = counter_stat
        if column2method[prop] != 'none':
            prop2type[add_suffix(prop, "counter")] = str

    for prop in num_prop:
        if not prop in column2method:
            column2method[prop] = num_stat
        if column2method[prop] == 'all':
            prop2type[add_suffix(prop, "avg")] = float
            prop2type[add_suffix(prop, "sum")] = float
            prop2type[add_suffix(prop, "max")] = float
            prop2type[add_suffix(prop, "min")] = float
            prop2type[add_suffix(prop, "std")] = float
        elif column2method[prop] == 'none':
            pass
        else:
            prop2type[add_suffix(prop, column2method[prop])] = float

    if not input_annotated_tree:
        node2leaves = annotated_tree.get_cached_content()

        # Prepare data for all nodes
        nodes_data = []
        for node in annotated_tree.traverse("postorder"):
            if not node.is_leaf:
                node_data = (node, node2leaves[node], text_prop, multiple_text_prop, bool_prop, num_prop, column2method, alignment, name2seq)
                nodes_data.append(node_data)

        # Process nodes in parallel if more than one thread is specified
        if threads > 1:
            with Pool(threads) as pool:
                results = pool.map(process_node, nodes_data)
        else:
            # For single-threaded execution, process nodes sequentially
            results = map(process_node, nodes_data)

        # Integrate the results back into your tree
        for result in results:
            node, internal_props, consensus_seq = result
            for key, value in internal_props.items():
                node.add_prop(key, value)
            if consensus_seq:
                node.add_prop(alignment_prop, consensus_seq)

        # #pre load node2leaves to save time
        # node2leaves = annotated_tree.get_cached_content()
        
        # for i, node in enumerate(annotated_tree.traverse("postorder")):
        #     internal_props = {}
        #     if not node.is_leaf:
        #         if text_prop:
        #             internal_props_text = merge_text_annotations(node2leaves[node], text_prop, column2method)
        #             internal_props.update(internal_props_text)

        #         if multiple_text_prop:
        #             internal_props_multi = merge_multitext_annotations(node2leaves[node], multiple_text_prop, column2method)
        #             internal_props.update(internal_props_multi)

        #         if bool_prop:
        #             internal_props_bool = merge_text_annotations(node2leaves[node], bool_prop, column2method)
        #             internal_props.update(internal_props_bool)

        #         if num_prop:
        #             internal_props_num = merge_num_annotations(node2leaves[node], num_prop, column2method)                        
        #             if internal_props_num:
        #                 internal_props.update(internal_props_num)
                
        #         for key,value in internal_props.items():
        #             node.add_prop(key, value)

        #         if alignment:
        #             # matrix = ''
        #             # for leaf in node.leaves():
        #             #     if name2seq.get(leaf.name):
        #             #         matrix += ">"+leaf.name+"\n"
        #             #         matrix += name2seq.get(leaf.name)+"\n"
        #             # consensus_seq = get_consensus_seq(StringIO(matrix), 0.7)
        #             # node.add_prop(alignment_prop, consensus_seq)
                    
        #             #matrix_string = build_matrix_string(node, name2seq)
        #             #consensus_seq = get_consensus_seq(matrix_string)

        #             def _consensus_node(node, name2seq, threshold=0.7):
        #                 matrix_string = build_matrix_string(node, name2seq)
        #                 consensus_seq = get_consensus_seq(matrix_string, threshold=threshold)
        #                 return consensus_seq
                    
        #             if threads > 1:
        #                 with Pool(threads - 1) as pool:
        #                     consensus_seq = pool.map(get_consensus_seq, matrix_string)
        #             else:
        #                 consensus_seq = _consensus_node(node, name2seq)
        #             node.add_prop(alignment_prop, consensus_seq)

                
    else:
        pass
        
    end = time.time()
    print('Time for merge annotations to run: ', end - start)


    # taxa annotations
    start = time.time()
    if taxon_column:
        if not taxadb:
            raise Exception('Please specify which taxa db using --taxadb <GTDB|NCBI>')
        else:
            
            if taxa_dump and taxadb == 'GTDB':
                logging.info(f"Loading GTDB database dump file {taxa_dump}...")
                GTDBTaxa().update_taxonomy_database(taxa_dump)
            elif taxa_dump and taxadb == 'NCBI':
                logging.info(f"Loading NCBI database dump file {taxa_dump}...")
                NCBITaxa().update_taxonomy_database(taxa_dump)
                
            annotated_tree, rank2values = annotate_taxa(annotated_tree, db=taxadb, \
                taxid_attr=taxon_column, sp_delimiter=taxon_delimiter, sp_field=taxa_field)
                
        # evolutionary events annotation
        annotated_tree = annotate_evol_events(annotated_tree, sp_delimiter=taxon_delimiter, sp_field=taxa_field)
        prop2type.update(TAXONOMICDICT)
    else:
        rank2values = {}

    end = time.time()
    print('Time for annotate_taxa to run: ', end - start)
    
    # prune tree by rank
    if rank_limit:
        annotated_tree = taxatree_prune(annotated_tree, rank_limit=rank_limit)

    # prune tree by condition
    if pruned_by: # need to be wrap with quotes
        condition_strings = pruned_by
        annotated_tree = conditional_prune(annotated_tree, condition_strings, prop2type)
    
    # name internal nodes
    annotated_tree = name_nodes(annotated_tree)
    return annotated_tree, prop2type

def run(args):
    total_color_dict = []
    layouts = []
    level = 1 # level 1 is the leaf name
    prop2type = {}
    metadata_dict = {}
    column2method = {}

    # checking file and output exists
    if not os.path.exists(args.tree):
        raise FileNotFoundError(f"Input tree {args.tree} does not exist.") 
    
    if args.metadata:
        for metadata_file in args.metadata:
            if not os.path.exists(metadata_file):
                raise FileNotFoundError(f"Metadata {metadata_file} does not exist.") 

    if not os.path.exists(args.outdir):
        raise FileNotFoundError(f"Output directory {args.outdir} does not exist.") 
        

    # parsing tree
    try:
        tree, eteformat_flag = validate_tree(args.tree, args.input_type, args.internal_parser)
    except TreeFormatError as e:
        print(e)
        sys.exit(1)

    # resolve polytomy
    if args.resolve_polytomy:
        tree.resolve_polytomy()
        
    # parse csv to metadata table
    start = time.time()
    print("start parsing...")
    # parsing metadata
    if args.metadata: # make a series aof metadatas
        metadata_dict, node_props, columns, prop2type = parse_csv(args.metadata, delimiter=args.metadata_sep, \
        no_colnames=args.no_colnames, aggregate_duplicate=args.aggregate_duplicate)
    else: # annotated_tree
        node_props=[]
        columns = {}
    end = time.time()
    print('Time for parse_csv to run: ', end - start)

    if args.emapper_annotations:
        emapper_metadata_dict, emapper_node_props, emapper_columns = parse_emapper_annotations(args.emapper_annotations)
        metadata_dict.update(emapper_metadata_dict)
        node_props.extend(emapper_node_props)
        columns.update(emapper_columns)
        prop2type.update({
            'name': str,
            'dist': float,
            'support': float,
            'seed_ortholog': str,
            'evalue': float,
            'score': float,
            'eggNOG_OGs': list,
            'max_annot_lvl': str,
            'COG_category': str,
            'Description': str,
            'Preferred_name': str,
            'GOs': list,
            'EC':str,
            'KEGG_ko': list,
            'KEGG_Pathway': list,
            'KEGG_Module': list,
            'KEGG_Reaction':list,
            'KEGG_rclass':list,
            'BRITE':list,
            'KEGG_TC':list,
            'CAZy':list,
            'BiGG_Reaction':list,
            'PFAMs':list
        })

    # start annotation
    if args.column_summary_method:
        column2method = process_column_summary_methods(args.column_summary_method)
    
    annotated_tree, prop2type = run_tree_annotate(tree, input_annotated_tree=args.annotated_tree,
        metadata_dict=metadata_dict, node_props=node_props, columns=columns,
        prop2type=prop2type,
        text_prop=args.text_prop, text_prop_idx=args.text_prop_idx,
        multiple_text_prop=args.multiple_text_prop, num_prop=args.num_prop, num_prop_idx=args.num_prop_idx,
        bool_prop=args.bool_prop, bool_prop_idx=args.bool_prop_idx,
        prop2type_file=args.prop2type, alignment=args.alignment,
        emapper_pfam=args.emapper_pfam, emapper_smart=args.emapper_smart, 
        counter_stat=args.counter_stat, num_stat=args.num_stat, column2method=column2method, 
        taxadb=args.taxadb, taxa_dump=args.taxa_dump, taxon_column=args.taxon_column,
        taxon_delimiter=args.taxon_delimiter, taxa_field=args.taxa_field,
        rank_limit=args.rank_limit, pruned_by=args.pruned_by, 
        acr_discrete_columns=args.acr_discrete_columns, 
        prediction_method=args.prediction_method, model=args.model, 
        delta_stats=args.delta_stats, ent_type=args.ent_type, 
        iteration=args.iteration, lambda0=args.lambda0, se=args.se,
        thin=args.thin, burn=args.burn,
        ls_columns=args.ls_columns, prec_cutoff=args.prec_cutoff, sens_cutoff=args.sens_cutoff, 
        threads=args.threads, outdir=args.outdir)

    if args.outdir:
        base=os.path.splitext(os.path.basename(args.tree))[0]
        out_newick = base + '_annotated.nw'
        out_prop2tpye = base + '_prop2type.txt'
        out_ete = base+'_annotated.ete'
        out_tsv = base+'_annotated.tsv'

        ### out newick
        annotated_tree.write(outfile=os.path.join(args.outdir, out_newick), props=None, 
                    parser=get_internal_parser(args.internal_parser), format_root_node=True)
        
        ### output prop2type
        with open(os.path.join(args.outdir, base+'_prop2type.txt'), "w") as f:
            #f.write(first_line + "\n")
            for key, value in prop2type.items():
                f.write("{}\t{}\n".format(key, value.__name__))
        ### out ete
        with open(os.path.join(args.outdir, base+'_annotated.ete'), 'w') as f:
            f.write(b64pickle.dumps(annotated_tree, encoder='pickle', pack=False))

        ### out tsv
        prop_keys = list(prop2type.keys())
        if args.taxon_column:
            prop_keys.extend(list(TAXONOMICDICT.keys()))
        if args.annotated_tree:
            tree2table(annotated_tree, internal_node=True, props=None, outfile=os.path.join(args.outdir, out_tsv))
        else:
            tree2table(annotated_tree, internal_node=True, props=prop_keys, outfile=os.path.join(args.outdir, out_tsv))

    # if args.outtsv:
    #     tree2table(annotated_tree, internal_node=True, outfile=args.outtsv)
    return

def check_missing(input_string):
    """
    define missing:
    1) One or more non-word characters at the beginning of the string.
    2) The exact strings "none", "None", "null", or "NaN".
    3) An empty string (zero characters).
    """
    pattern = r'^(?:\W+|none|None|null|Null|NaN|)$'
    
    if input_string is None:
        return True
    elif re.match(pattern, input_string):
        #print("Input contains only non-alphanumeric characters, 'none', a missing value, or an empty value.")
        return True
    else:
        return False


def check_tar_gz(file_path):
    try:
        with tarfile.open(file_path, 'r:gz') as tar:
            return True
    except tarfile.ReadError:
        return False

def parse_csv(input_files, delimiter='\t', no_colnames=False, aggregate_duplicate=False):
    """
    Takes tsv table as input
    Return
    metadata, as dictionary of dictionaries for each node's metadata
    node_props, a list of property names(column names of metadata table)
    columns, dictionary of property name and it's values
    """
    metadata = {}
    columns = defaultdict(list)
    prop2type = {}
    def update_metadata(reader, node_header):
        for row in reader:
            nodename = row[node_header]
            del row[node_header]
            #row = {k: 'NaN' if (not v or v.lower() == 'none') else v for k, v in row.items() } ## replace empty to NaN
            for k, v in row.items(): # replace missing value
                if check_missing(v):
                    row[k] = 'NaN'
                else:
                    row[k] = v

            if nodename in metadata.keys():
                for prop, value in row.items():
                    if aggregate_duplicate:
                        if prop in metadata[nodename]:
                            exisiting_value = metadata[nodename][prop]
                            new_value = ','.join([exisiting_value,value])
                            metadata[nodename][prop] = new_value
                            columns[prop].append(new_value)
                        else:
                            metadata[nodename][prop] = value
                            columns[prop].append(value)
                    else:
                        metadata[nodename][prop] = value
                        columns[prop].append(value)
            else:
                metadata[nodename] = dict(row)
                for (prop, value) in row.items(): # go over each column name and value
                    columns[prop].append(value) # append the value into the appropriate list
                                    # based on column name k

    def update_prop2type(node_props):
        for prop in node_props:
            if set(columns[prop])=={'NaN'}:
                #prop2type[prop] = np.str_
                prop2type[prop] = str
            else:
                dtype = infer_dtype(columns[prop])
                prop2type[prop] = dtype # get_type_convert(dtype)
    
    for input_file in input_files:
        # check file
        if check_tar_gz(input_file):
            with tarfile.open(input_file, 'r:gz') as tar:
                for member in tar.getmembers():
                    if member.isfile() and member.name.endswith('.tsv'):
                        with tar.extractfile(member) as tsv_file:
                            tsv_text = tsv_file.read().decode('utf-8').splitlines()
                            if no_colnames:
                                fields_len = len(tsv_text[0].split(delimiter))
                                headers = ['col'+str(i) for i in range(fields_len)]
                                reader = csv.DictReader(tsv_text, delimiter=delimiter,fieldnames=headers)
                            else:
                                reader = csv.DictReader(tsv_text, delimiter=delimiter)
                                headers = reader.fieldnames
                            node_header, node_props = headers[0], headers[1:]
                            update_metadata(reader, node_header)
                        
                        update_prop2type(node_props)
                            
        else:          
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

                    #row = {k: 'NaN' if (not v or v.lower() == 'none') else v for k, v in row.items() } ## replace empty to NaN

                    for k, v in row.items(): # replace missing value
                        if check_missing(v):
                            row[k] = 'NaN'
                        else:
                            row[k] = v

                    if nodename in metadata.keys():
                        for prop, value in row.items():
                            if aggregate_duplicate:
                                if prop in metadata[nodename]:
                                    exisiting_value = metadata[nodename][prop]
                                    new_value = ','.join([exisiting_value,value])
                                    metadata[nodename][prop] = new_value
                                    columns[prop].append(new_value)
                                else:
                                    metadata[nodename][prop] = value
                                    columns[prop].append(value)
                            else:
                                metadata[nodename][prop] = value
                                columns[prop].append(value)
                    else:
                        metadata[nodename] = dict(row)
                        for (prop, value) in row.items(): # go over each column name and value
                            columns[prop].append(value) # append the value into the appropriate list
                                            # based on column name k
            update_prop2type(node_props)

    return metadata, node_props, columns, prop2type

def process_column_summary_methods(column_summary_methods):
    column_methods = {}
    if column_summary_methods:
        for entry in column_summary_methods:
            try:
                column, method = entry.split('=')
                column_methods[column] = method
            except ValueError:
                raise ValueError(f"Invalid format for --column-summary-method: '{entry}'. Expected format: ColumnName=Method")
    return column_methods

def get_comma_separated_values(lst):
    for item in lst:
        if isinstance(item, str) and any(',' in x for x in item.split()):
            return True
    return False

def can_convert_to_bool(column):
    true_values = {'true', 't', 'yes', 'y', '1'}
    false_values = {'false', 'f', 'no', 'n', '0'}
    ignore_values = {'nan', 'none', ''}  # Add other representations of NaN as needed

    # Initialize sets to hold the representations of true and false values
    true_representations = set()
    false_representations = set()

    for value in column:
        str_val = str(value).strip()  # Preserving the original capitalization
        if str_val.lower() in ignore_values:
            continue  # Skip this value
        if str_val.lower() in true_values:
            true_representations.add(str_val)
        elif str_val.lower() in false_values:
            false_representations.add(str_val)
        else:
            return False

    # Check that all true values and all false values have exactly one representation
    return len(true_representations) <= 1 and len(false_representations) <= 1


def convert_column_data(column, np_dtype):
    #np_dtype = np.dtype(dtype).type
    try:
        data = np.array(column).astype(np_dtype)
        return np_dtype
    except ValueError as e:
        return None

def convert_to_prop_array(metadata_dict, prop):
    """
    Convert a dictionary of metadata to a structured array format.

    Parameters:
    metadata_dict (dict): The original dictionary containing metadata.
    prop (str): The property to extract from the metadata.

    Returns:
    dict: A dictionary with keys as properties and values as structured arrays.
    
    # {"leaf1":{"prop":"value1"},{"leaf2":{"prop":"value2"}}} 
    # to 
    # {"prop":[["leaf1", "leaf2"],["value1", "value2"]]}
    """
    prop_array = {prop: [[], []]}
    for leaf, value in metadata_dict.items():
        prop_array[prop][0].append(leaf)  # Append key
        prop_array[prop][1].append(value.get(prop, None))  # Append property value, handle missing values

    return prop_array

def convert_back_to_original(prop2array):
    """
    Convert the structured array format back to the original dictionary format.

    Parameters:
    prop2array (dict): The structured array format dictionary.

    Returns:
    dict: The original format of the data.
    """
    metadata_dict = {}
    for key in prop2array:
        identifiers, values = prop2array[key]
        for identifier, value in zip(identifiers, values):
            if identifier not in metadata_dict:
                metadata_dict[identifier] = {}
            metadata_dict[identifier][key] = value
    return metadata_dict

def infer_dtype(column):
    if get_comma_separated_values(column):
        return list
    elif can_convert_to_bool(column):
        return bool
    else:
        dtype_dict = {
            float:np.float64,
            str:np.str_
            }
        #dtype_order = ['float64', 'str']
        for dtype, np_dtype in dtype_dict.items():
            result = convert_column_data(column, np_dtype)
            if result is not None:
                # Successful inference, exit from the loop
                return dtype
        return None

def load_metadata_to_tree(tree, metadata_dict, prop2type={}, taxon_column=None, taxon_delimiter='', taxa_field=0):
    #name2leaf = {}
    multi_text_seperator = ','

    name2leaf = defaultdict(list)
    # preload all leaves to save time instead of search in tree
    for leaf in tree.leaves():
        name2leaf[leaf.name].append(leaf)
    
    # load all metadata to leaf nodes
    for node, props in metadata_dict.items():
        if node in name2leaf.keys():
            target_nodes = name2leaf[node]
            for target_node in target_nodes:
                for key,value in props.items():
                    # taxa
                    if key == taxon_column:
                        if taxon_delimiter:
                            taxon_prop = value.split(taxon_delimiter)[taxa_field]
                        else:
                            taxon_prop = value
                        target_node.add_prop(key, taxon_prop)
                    
                    # numerical
                    elif key in prop2type and prop2type[key]==float:
                        try:
                            flot_value = float(value)
                            if math.isnan(flot_value):
                                target_node.add_prop(key, 'NaN')
                            else:
                                target_node.add_prop(key, flot_value)
                        except (ValueError,TypeError):
                            target_node.add_prop(key, 'NaN')

                    # categorical
                    # list
                    elif key in prop2type and prop2type[key]==list:
                        value_list = value.split(multi_text_seperator)
                        target_node.add_prop(key, value_list)
                    # str
                    else:
                        target_node.add_prop(key, value)
        else:
            pass

        # hits = tree.get_leaves_by_name(node)
        # if hits:
        #     for target_node in hits:
        #         for key,value in props.items():
        #             if key == taxon_column:
        #                 taxon_prop = value.split(taxon_delimiter)[-1]
        #                 target_node.add_prop(key, taxon_prop)
        #             elif key in prop2type and prop2type[key]=='num':
        #                 if math.isnan(float(value)):
        #                     target_node.add_prop(key, value)
        #                 else:
        #                     target_node.add_prop(key, float(value))
        #             else:
        #                 target_node.add_prop(key, value)
        # else:
        #     pass
        #hits = tree.search_nodes(name=node) # including internal nodes

    return tree

def process_node(node_data):
    node, node_leaves, text_prop, multiple_text_prop, bool_prop, num_prop, column2method, alignment, name2seq = node_data
    internal_props = {}

    # Process text, multitext, bool, and num properties
    if text_prop:
        internal_props_text = merge_text_annotations(node_leaves, text_prop, column2method)
        internal_props.update(internal_props_text)

    if multiple_text_prop:
        internal_props_multi = merge_multitext_annotations(node_leaves, multiple_text_prop, column2method)
        internal_props.update(internal_props_multi)

    if bool_prop:
        internal_props_bool = merge_text_annotations(node_leaves, bool_prop, column2method)
        internal_props.update(internal_props_bool)

    if num_prop:
        internal_props_num = merge_num_annotations(node_leaves, num_prop, column2method)
        if internal_props_num:
            internal_props.update(internal_props_num)

    # Generate consensus sequence
    consensus_seq = None
    if alignment:  # Assuming 'alignment' is a condition to check
        matrix_string = build_matrix_string(node, name2seq)  # Assuming 'name2seq' is accessible here
        consensus_seq = get_consensus_seq(matrix_string, threshold=0.7)

    return node, internal_props, consensus_seq

def merge_text_annotations(nodes, target_props, column2method):
    pair_seperator = "--"
    item_seperator = "||"
    internal_props = {}
    for target_prop in target_props:
        counter_stat = column2method.get(target_prop, "raw")
        if counter_stat == 'raw':
            prop_list = children_prop_array_missing(nodes, target_prop)
            internal_props[add_suffix(target_prop, 'counter')] = item_seperator.join([add_suffix(str(key), value, pair_seperator) for key, value in sorted(dict(Counter(prop_list)).items())])

        elif counter_stat == 'relative':
            prop_list = children_prop_array_missing(nodes, target_prop)
            counter_line = []

            total = sum(dict(Counter(prop_list)).values())

            for key, value in sorted(dict(Counter(prop_list)).items()):

                rel_val = '{0:.2f}'.format(float(value)/total)
                counter_line.append(add_suffix(key, rel_val, pair_seperator))
            internal_props[add_suffix(target_prop, 'counter')] = item_seperator.join(counter_line)
            #internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])

        else:
            #print('Invalid stat method')
            pass

    return internal_props

def merge_multitext_annotations(nodes, target_props, column2method):
    #seperator of multiple text 'GO:0000003,GO:0000902,GO:0000904'
    multi_text_seperator = ','
    pair_seperator = "--"
    item_seperator = "||"

    internal_props = {}
    for target_prop in target_props:
        counter_stat = column2method.get(target_prop, "raw")
        if counter_stat == 'raw':
            prop_list = children_prop_array(nodes, target_prop)
            multi_prop_list = []

            for elements in prop_list:
                for j in elements:
                    multi_prop_list.append(j)
            internal_props[add_suffix(target_prop, 'counter')] = item_seperator.join([add_suffix(str(key), value, pair_seperator) for key, value in sorted(dict(Counter(multi_prop_list)).items())])

        elif counter_stat == 'relative':
            prop_list = children_prop_array(nodes, target_prop)
            multi_prop_list = []

            for elements in prop_list:
                for j in elements:
                    multi_prop_list.append(j)

            counter_line = []

            total = sum(dict(Counter(multi_prop_list)).values())

            for key, value in sorted(dict(Counter(multi_prop_list)).items()):
                rel_val = '{0:.2f}'.format(float(value)/total)
                counter_line.append(add_suffix(key, rel_val, pair_seperator))
            internal_props[add_suffix(target_prop, 'counter')] = item_seperator.join(counter_line)
            #internal_props[add_suffix(target_prop, 'counter')] = '||'.join([add_suffix(key, value, '--') for key, value in dict(Counter(prop_list)).items()])
        else:
            #print('Invalid stat method')
            pass

    return internal_props

def merge_num_annotations(nodes, target_props, column2method):
    internal_props = {}
    for target_prop in target_props:
        num_stat = column2method.get(target_prop, None)
        if num_stat != 'none':
            if target_prop != 'dist' and target_prop != 'support':
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
                        #print(target_prop)
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
                        #print('Invalid stat method')
                        pass
                else:
                    pass

    if internal_props:
        return internal_props
    else:
        return None

def name_nodes(tree):
    for i, node in enumerate(tree.traverse("postorder")):
        if not node.name:
            if not node.is_root:
                node.name = 'N'+str(i)
            else:
                node.name = 'Root'
    return tree

def gtdb_accession_to_taxid(accession):
        """Given a GTDB accession number, returns its complete accession"""
        if accession.startswith('GCA'):
            prefix = 'GB_'
            return prefix+accession
        elif accession.startswith('GCF'):
            prefix = 'RS_'
            return prefix+accession
        else:
            return accession

def annotate_taxa(tree, db="GTDB", taxid_attr="name", sp_delimiter='.', sp_field=0):
    global rank2values
    logging.info(f"\n==============Annotating tree with {db} taxonomic database============")
    
    def return_spcode_ncbi(leaf):
        try:
            return leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field]
        except (IndexError, ValueError):
            return leaf.props.get(taxid_attr)

    def return_spcode_gtdb(leaf):
        try:
            if sp_delimiter:
                species_attribute = leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field]
                return gtdb_accession_to_taxid(species_attribute)
            else:
                return gtdb_accession_to_taxid(leaf.props.get(taxid_attr))
        except (IndexError, ValueError):
            return gtdb_accession_to_taxid(leaf.props.get(taxid_attr))

    if db == "GTDB":
        gtdb = GTDBTaxa()
        tree.set_species_naming_function(return_spcode_gtdb)
        gtdb.annotate_tree(tree,  taxid_attr="species")

    elif db == "NCBI":
        ncbi = NCBITaxa()
        # extract sp codes from leaf names
        tree.set_species_naming_function(return_spcode_ncbi)
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

def annotate_evol_events(tree, sp_delimiter='.', sp_field=0):
    def return_spcode(leaf):
        try:
            return leaf.name.split(sp_delimiter)[sp_field]
        except (IndexError, ValueError):
            return leaf.name

    tree.set_species_naming_function(return_spcode)

    node2species = tree.get_cached_content('species')
    for n in tree.traverse():
        n.props['species'] = node2species[n]
        if len(n.children) == 2:
            dup_sp = node2species[n.children[0]] & node2species[n.children[1]]
            if dup_sp:
                n.props['evoltype'] = 'D'
                n.props['dup_sp'] = ','.join(dup_sp)
                n.props['dup_percent'] = round(len(dup_sp)/len(node2species[n]), 3) * 100
            else:
                n.props['evoltype'] = 'S'
        n.del_prop('_speciesFunction')
    return tree

def get_range(input_range):
    column_range = input_range[input_range.find("[")+1:input_range.find("]")]
    column_start, column_end = [int(i) for i in column_range.split('-')]
    #column_list_idx = [i for i in range(column_start, column_end+1)]
    return column_start, column_end

def parse_emapper_annotations(input_file, delimiter='\t', no_colnames=False):
    metadata = {}
    columns = defaultdict(list)
    prop2type = {}
    headers = ["#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
               "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs",
               "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
               "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]

    with open(input_file, 'r') as f:
        # Skip lines starting with '##'
        filtered_lines = (line for line in f if not line.startswith('##'))

        if no_colnames:
            reader = csv.DictReader(filtered_lines, delimiter=delimiter, fieldnames=headers)
        else:
            reader = csv.DictReader(filtered_lines, delimiter=delimiter)

        node_header, node_props = headers[0], headers[1:]
        for row in reader:
            nodename = row[node_header]
            del row[node_header]

            for k, v in row.items():  # Replace missing value
                row[k] = 'NaN' if check_missing(v) else v
            metadata[nodename] = dict(row)
            for k, v in row.items():  # Go over each column name and value
                columns[k].append(v)  # Append the value into the appropriate list based on column name k

    return metadata, node_props, columns

def annot_tree_pfam_table(post_tree, pfam_table, alg_fasta, domain_prop='dom_arq'):
    pair_delimiter = "@"
    item_seperator = "||"
    fasta = SeqGroup(alg_fasta) # aligned_fasta
    raw2alg = defaultdict(dict)

    for num, (name, seq, _) in enumerate(fasta):
        p_raw = 1
        for p_alg, (a) in enumerate(seq, 1):
            if a != '-':
                raw2alg[name][p_raw] = p_alg
                p_raw +=1
    
    seq2doms = defaultdict(list)
    with open(pfam_table) as f_in:
        for line in f_in:
            if not line.startswith('#'):
                info = line.strip().split('\t')
                seq_name = info[0]
                dom_name = info[1]
                dom_start = int(info[7])
                dom_end = int(info[8])
                if raw2alg.get(seq_name):
                    trans_dom_start = raw2alg[seq_name][dom_start]
                    trans_dom_end = raw2alg[seq_name][dom_end]
                    dom_info_string = pair_delimiter.join([dom_name, str(trans_dom_start), str(trans_dom_end)])
                    seq2doms[seq_name].append(dom_info_string)

    for l in post_tree:
        if l.name in seq2doms.keys():
            domains = seq2doms[l.name]
            domains_string = item_seperator.join(domains)
            l.add_prop(domain_prop, domains_string)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_node_domains = n.get_closest_leaf()[0].props.get(domain_prop, 'none@none@none')
            n.add_prop(domain_prop, random_node_domains)

    # for n in post_tree.traverse():
    #     print(n.name, n.props.get('dom_arq'))

def annot_tree_smart_table(post_tree, smart_table, alg_fasta, domain_prop='dom_arq'):
    pair_delimiter = "@"
    item_seperator = "||"
    fasta = SeqGroup(alg_fasta) # aligned_fasta
    raw2alg = defaultdict(dict)
    for num, (name, seq, _) in enumerate(fasta):
        p_raw = 1
        for p_alg, (a) in enumerate(seq, 1):
            if a != '-':
                raw2alg[name][p_raw] = p_alg
                p_raw +=1

    seq2doms = defaultdict(list)
    with open(smart_table) as f_in:
        for line in f_in:
            if not line.startswith('#'):
                info = line.strip().split('\t')
                seq_name = info[0]
                dom_name = info[1]
                dom_start = int(info[2])
                dom_end = int(info[3])
                if raw2alg.get(seq_name):
                    trans_dom_start = raw2alg[seq_name][dom_start]
                    trans_dom_end = raw2alg[seq_name][dom_end]

                dom_info_string = pair_delimiter.join([dom_name, str(trans_dom_start), str(trans_dom_end)])
                seq2doms[seq_name].append(dom_info_string)

    for l in post_tree:
        if l.name in seq2doms.keys():
            domains = seq2doms[l.name]
            domains_string = item_seperator.join(domains)
            l.add_prop(domain_prop, domains_string)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_node_domains = n.get_closest_leaf()[0].props.get(domain_prop, 'none@none@none')
            n.add_prop(domain_prop, random_node_domains)

    # for n in post_tree.traverse():
    #     print(n.name, n.props.get('dom_arq'))

def parse_fasta(fastafile):
    fasta_dict = {}
    with open(fastafile,'r') as f:
        head = ''
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq != '':
                    fasta_dict[head] = seq
                    seq = ''
                    head = line[1:]
                else:
                    head = line[1:]
            else:
                seq += line
    fasta_dict[head] = seq
    return fasta_dict

# def goslim_annotation(gos_input, relative=True):
#     """
#     deprecated
#     """
#     output_dict = {}
#     all_golsims_dict = {}
#     goslim_script = os.path.join(os.path.dirname(__file__), 'goslim_list.R')
#     output = subprocess.check_output([goslim_script,gos_input])
#     #output_list = [line.split(' \t ') for line in output.decode('utf-8').split('\n') if line ]
#     for line in output.decode('utf-8').split('\n'):
#         if line:
#             name, entries, desc, count = line.split(' \t ')
#             if entries != '-':
#                 entries = entries.split(',')
#                 desc = desc.split('|')
#                 count = np.array(count.split('|')).astype(int)
#                 if relative:
#                     count = [float(i)/sum(count) for i in count]
#             output_dict[name] = [entries, desc, count]
#             for i in range(len(entries)):
#                 entry = entries[i]
#                 single_desc = desc[i]
#                 if entry not in all_golsims_dict:
#                     all_golsims_dict[entry] = single_desc
#     return output_dict, all_golsims_dict

def get_pval(prop2array, dump_tree, acr_discrete_columns_dict, iteration=100, 
            prediction_method="MPPA", model="F81", ent_type='SE', 
            lambda0=0.1, se=0.5, sim=10000, burn=100, thin=10, threads=1):
    prop2delta_array = {}
    for _ in range(iteration):
        shuffled_dict = {}
        for column, trait in acr_discrete_columns_dict.items():
            trait = acr_discrete_columns_dict[column]
            #shuffle traits
            shuffled_trait = np.random.choice(trait, len(trait), replace=False)
            prop2array[column][1] = list(shuffled_trait)
            shuffled_dict[column] = list(shuffled_trait)

        # Converting back to the original dictionary format
        # # annotate new metadata to leaf
        new_metadata_dict = convert_back_to_original(prop2array)
        dump_tree = load_metadata_to_tree(dump_tree, new_metadata_dict)
        
        # # run acr
        random_acr_results, dump_tree = run_acr_discrete(dump_tree, shuffled_dict, \
        prediction_method="MPPA", model="F81", threads=threads, outdir=None)
        random_delta = run_delta(random_acr_results, dump_tree, ent_type=ent_type, 
                lambda0=lambda0, se=se, sim=sim, burn=burn, thin=thin, 
                threads=threads)

        for prop, delta_result in random_delta.items():
            if prop in prop2delta_array:
                prop2delta_array[prop].append(delta_result)
            else:
                prop2delta_array[prop] = [delta_result]
        clear_extra_features([dump_tree], ["name", "dist", "support"])

    return prop2delta_array

# Function to build the matrix string for a node
def build_matrix_string(node, name2seq):
    matrix = ''
    for leaf in node.leaves():
        if name2seq.get(leaf.name):
            matrix += f">{leaf.name}\n{name2seq.get(leaf.name)}\n"
    return matrix

def tree2table(tree, internal_node=True, props=None, outfile='tree2table.csv'):
    node2leaves = {}
    leaf2annotations = {}
    if not props:
        props = set()
        for node in tree.traverse():
            props |= node.props.keys()
        props = [ p for p in props if not p.startswith("_") ]

    with open(outfile, 'w', newline='') as csvfile:
        if '_speciesFunction' in props:
            props.remove('_speciesFunction')
        fieldnames = ['name', 'dist', 'support']
        fieldnames.extend(x for x in sorted(props) if x not in fieldnames)

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for node in tree.traverse():
            if node.name:
                # if '_speciesFunction' in node.props:
                #     node.del_prop('_speciesFunction')

                if node.is_leaf:
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
