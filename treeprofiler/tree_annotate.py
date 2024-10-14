#!/usr/bin/env python3

import os, math, re
import logging
import sys
import time
import random
import csv
import tarfile

from collections import defaultdict, Counter
import numpy as np
from scipy import stats
import requests

from ete4.parser.newick import NewickError
from ete4 import SeqGroup
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa

from treeprofiler.src import utils
from treeprofiler.src.phylosignal import run_acr_discrete, run_acr_continuous, run_delta
from treeprofiler.src.ls import run_ls
from treeprofiler.src import b64pickle

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

# Global Variable for emapper headers
EMAPPER_HEADERS = ["#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
                   "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs",
                   "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
                   "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]

# Available methods and models for ACR
# Discrete traits
DISCRETE_METHODS = ['MPPA', 'MAP', 'JOINT', 'DOWNPASS', 'ACCTRAN', 'DELTRAN', 'COPY', 'ALL', 'ML', 'MP']
DISCRETE_MODELS = ['JC', 'F81', 'EFT', 'HKY', 'JTT']

# Continuous traits
CONTINUOUS_METHODS = ['ML', 'BAYESIAN']
CONTINUOUS_MODELS = ['BM', 'OU']

# Set up the logger with INFO level by default
logger = logging.getLogger(__name__)
def setup_logger():
    """Sets up logging configuration."""
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

def populate_annotate_args(parser):
    gmeta = parser.add_argument_group(
        title='METADATA TABLE parameters',
        description="Input parameters of METADATA")
    add = gmeta.add_argument
    add('-m', '--metadata', nargs='+',
        help="<metadata.csv> .csv, .tsv. mandatory input")
    # add('--data-matrix', nargs='+',
    #     help="<metadata.csv> .csv, .tsv. optional input")
    add('--data-matrix',  nargs='+',
        help="<datamatrix.csv> .csv, .tsv. matrix data metadata table as array to tree, please do not provide column headers in this file")
    add('-s', '--metadata-sep', default='\t',
        help="column separator of metadata table [default: \\t]")
    add('--no-headers', action='store_true',
        help="metadata table doesn't contain columns name, namespace col+index will be assigned as the key of property such as col1.")
    add('--duplicate', action='store_true',
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
    add('--ls-columns', nargs='+',
        help=("<col1> <col2> names to perform lineage specificity analysis"))
    # add('--taxatree',
    #     help=("<kingdom|phylum|class|order|family|genus|species|subspecies> "
    #           "reference tree from taxonomic database"))
    add('--taxadb', type=str.upper,
        choices=['NCBI', 'GTDB', 'customdb'],
        help="<NCBI|GTDB> for taxonomic annotation or fetch taxatree")
    add('--gtdb-version', type=int,
        choices=[95, 202, 207, 214, 220],
        help='GTDB version for taxonomic annotation, such as 220. If it is not provided, the latest version will be used.')
    add('--taxa-dump', type=str,
        help='Path to taxonomic database dump file for specific version, such as gtdb taxadump https://github.com/etetoolkit/ete-data/raw/main/gtdb_taxonomy/gtdblatest/gtdb_latest_dump.tar.gz or NCBI taxadump https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz')
    add('--taxon-column',
        help="Activate taxonomic annotation using <col1> name of columns which need to be read as taxon data. \
            Unless taxon data in leaf name, please use 'name' as input such as --taxon-column name")
    add('--taxon-delimiter', default=None,
        help="delimiter of taxa columns. [default: None]")
    add('--taxa-field', type=int, default=0,
        help="field of taxa name after delimiter. [default: 0]")
    add('--ignore-unclassified', action='store_true',
        help="Ignore unclassified taxa in taxonomic annotation")
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
        help="Specify summary method for individual columns in the format COL=METHOD. Method option can be seen in --counter-stat and --num-stat.")
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
    # ACR for discrete traits columns
    acr_group.add_argument('--acr-discrete-columns', nargs='+',
        help=("List of column names (e.g., <col1> <col2>) to perform ACR analysis for discrete traits."))
    # ACR for continuous traits columns
    acr_group.add_argument('--acr-continuous-columns', nargs='+',
        help=("List of column names (e.g., <col1> <col2>) to perform ACR analysis for continuous traits."))
    acr_group.add_argument('--prediction-method',
        default='MPPA',
        choices=DISCRETE_METHODS + CONTINUOUS_METHODS,
        type=str,
        required=False,
        help=("Prediction method for ACR analysis. "
              f"For discrete traits: {', '.join(DISCRETE_METHODS)}. "
              f"For continuous traits: {', '.join(CONTINUOUS_METHODS)}. [default: MPPA]"))
    acr_group.add_argument('--model',
        default='F81',
        choices=DISCRETE_MODELS + CONTINUOUS_MODELS,
        type=str,
        required=False,
        help=("Evolutionary model for ML methods in ACR analysis. "
              f"For discrete traits: {', '.join(DISCRETE_MODELS)}. "
              f"For continuous traits: {', '.join(CONTINUOUS_MODELS)}. [default: F81]"))
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
    delta_group.add_argument('--ent-type',
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
    group.add_argument('--quiet',
        default=False,
        action='store_true',
        help="Suppress logging messages")
    group.add_argument('--stdout',
        default=False,
        action='store_true',
        help="Print the output to stdout")
    group.add_argument('-o', '--outdir',
        type=str,
        required=False,
        help="Directory for annotated outputs.")

def run_tree_annotate(tree, input_annotated_tree=False,
        metadata_dict={}, node_props=[], columns={}, prop2type={},
        text_prop=[], text_prop_idx=[], multiple_text_prop=[], num_prop=[], num_prop_idx=[],
        bool_prop=[], bool_prop_idx=[], prop2type_file=None, alignment=None, emapper_mode=False, emapper_pfam=None,
        emapper_smart=None, counter_stat='raw', num_stat='all', column2method={},
        taxadb='GTDB', gtdb_version=None, taxa_dump=None, taxon_column=None,
        taxon_delimiter='', taxa_field=0, ignore_unclassified=False,
        rank_limit=None, pruned_by=None, 
        acr_discrete_columns=None, acr_continuous_columns=None, prediction_method="MPPA", model="F81", 
        delta_stats=False, ent_type="SE", 
        iteration=100, lambda0=0.1, se=0.5, thin=10, burn=100, 
        ls_columns=None, prec_cutoff=0.95, sens_cutoff=0.95, 
        threads=1, outdir='./'):

    total_color_dict = []
    layouts = []
    level = 1 # level 1 is the leaf name
    
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
            logger.error("Please provide alignment file using '--alignment' for pfam annotation.")
            sys.exit(1)
        annot_tree_pfam_table(tree, emapper_pfam, alignment, domain_prop=domain_prop)
        prop2type.update({
            domain_prop:str
            })
    if emapper_smart:
        domain_prop = 'dom_arq'
        if not alignment:
            logger.error("Please provide alignment file using '--alignment' for smart annotation.")
            sys.exit(1)
        annot_tree_smart_table(tree, emapper_smart, alignment, domain_prop=domain_prop)
        prop2type.update({
            domain_prop:str
            })

    # load all metadata to leaf nodes

    # input_annotated_tree determines if input tree is already annotated, if annotated, no longer need metadata
    
    if not input_annotated_tree:
        if taxon_column: # to identify taxon column as taxa property from metadata
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type, taxon_column=taxon_column, taxon_delimiter=taxon_delimiter, taxa_field=taxa_field, ignore_unclassified=ignore_unclassified)
        else:
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type)
    else:
        annotated_tree = tree

    end = time.time()
    logger.info(f'Time for load_metadata_to_tree to run: {end - start}')

    
    # Ancestor Character Reconstruction analysis
    # discrete data preparation
    if acr_discrete_columns:
        logger.info(f"Performing ACR analysis with discrete traits {acr_discrete_columns} via {prediction_method} method with {model} model.......\n")
        # need to be discrete traits
        discrete_traits = text_prop + bool_prop
        for k in acr_discrete_columns:
            if k:
                if k not in discrete_traits:
                    logger.error(f"Character {k} is not discrete trait, please check your input.")
                    sys.exit(1)

        if prediction_method in DISCRETE_METHODS:
            if model in DISCRETE_MODELS:
                pass
            else:
                logger.error(f"Model {model} is not supported for discrete traits, please check your input.")
                sys.exit(1)
        else:
            logger.error(f"Prediction method {prediction_method} is not supported for discrete traits, please check your input.")
            sys.exit(1)
        #############################
        start = time.time()
        acr_discrete_columns_dict = {k: v for k, v in columns.items() if k in acr_discrete_columns}
        acr_results, annotated_tree = run_acr_discrete(annotated_tree, acr_discrete_columns_dict, \
        prediction_method=prediction_method, model=model, threads=threads, outdir=outdir)
        
        # Clear extra features
        utils.clear_extra_features([annotated_tree], prop2type.keys())
        
        # get observed delta
        # only MPPA,MAP method has marginal probabilities to calculate delta
        if delta_stats:
            if prediction_method in ['MPPA', 'MAP']:
                logger.info(f"Performing Delta Statistic analysis with Character {acr_discrete_columns}...\n")
                prop2delta = run_delta(acr_results, annotated_tree, ent_type=ent_type, 
                lambda0=lambda0, se=se, sim=iteration, burn=burn, thin=thin, 
                threads=threads)

                for prop, delta_result in prop2delta.items():
                    logger.info(f"Delta statistic of {prop} is: {delta_result}")
                    tree.add_prop(utils.add_suffix(prop, "delta"), delta_result)

                # start calculating p_value
                logger.info(f"Calculating p_value for delta statistic...")
                # get a copy of the tree
                dump_tree = annotated_tree.copy()
                utils.clear_extra_features([dump_tree], ["name", "dist", "support"])
                
                prop2array = {}
                for prop in columns.keys():
                    prop2array.update(convert_to_prop_array(metadata_dict, prop))
                
                prop2delta_array = get_pval(prop2array, dump_tree, acr_discrete_columns_dict, \
                    iteration=100, prediction_method=prediction_method, model=model,
                    ent_type=ent_type, lambda0=lambda0, se=se, sim=iteration, burn=burn, thin=thin, 
                    threads=threads)

                for prop, delta_array in prop2delta_array.items():
                    p_value = np.sum(np.array(delta_array) > prop2delta[prop]) / len(delta_array)
                    logger.info(f"p_value of {prop} is {p_value}")
                    tree.add_prop(utils.add_suffix(prop, "pval"), p_value)
                    prop2type.update({
                        utils.add_suffix(prop, "pval"): float
                    })

                for prop in acr_discrete_columns:
                    prop2type.update({
                        utils.add_suffix(prop, "delta"): float
                    })
            else:
                logger.warning(f"Delta statistic analysis only support MPPA and MAP prediction method, {prediction_method} is not supported.")

        end = time.time()
        logger.info(f'Time for acr to run: {end - start}')

    # continuous data preparation
    if acr_continuous_columns:
        logger.info(f"Performing ACR analysis with continuous traits {acr_continuous_columns} via {prediction_method} method with {model} model.......\n")
        # need to be discrete traits
        continuous_traits = num_prop
        for k in acr_continuous_columns:
            if k:
                if k not in continuous_traits:
                    logger.error(f"Character {k} is not continuous trait, please check your input.")
                    sys.exit(1)
        if prediction_method in CONTINUOUS_METHODS:
            if model in CONTINUOUS_MODELS:
                pass
            else:
                logger.error(f"Model {model} is not supported for continuous traits, please check your input.")
                sys.exit(1)
        else:
            logger.error(f"Prediction method {prediction_method} is not supported for continuous traits, please check your input.")
            sys.exit(1)
        # convert metadata to observed traits
        transformed_dict = {key: {} for key in acr_continuous_columns}
        for leaf, props in metadata_dict.items():
            for prop in acr_continuous_columns:
                transformed_dict[prop][leaf] = float(props[prop])

        start = time.time()
        acr_results, tree = run_acr_continuous(annotated_tree, transformed_dict, model=model, prediction_method=prediction_method, threads=threads, outdir=outdir)
        end = time.time()
        logger.info(f'Time for acr to run: {end - start}')

    # lineage specificity analysis
    if ls_columns:
        logger.info(f"Performing Lineage Specificity analysis with Character {ls_columns}...\n")
        if all(column in bool_prop for column in ls_columns):
            best_node, qualified_nodes = run_ls(annotated_tree, props=ls_columns, 
            precision_cutoff=prec_cutoff, sensitivity_cutoff=sens_cutoff)
            for prop in ls_columns:
                prop2type.update({
                    utils.add_suffix(prop, "prec"): float,
                    utils.add_suffix(prop, "sens"): float,
                    utils.add_suffix(prop, "f1"): float
                })
        else:
            logger.warning(f"Lineage specificity analysis only support boolean properties, {ls_columns} is not boolean property.")

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
            prop2type[utils.add_suffix(prop, "counter")] = str

    for prop in num_prop:
        if not prop in column2method:
            column2method[prop] = num_stat
        if column2method[prop] == 'all':
            prop2type[utils.add_suffix(prop, "avg")] = float
            prop2type[utils.add_suffix(prop, "sum")] = float
            prop2type[utils.add_suffix(prop, "max")] = float
            prop2type[utils.add_suffix(prop, "min")] = float
            prop2type[utils.add_suffix(prop, "std")] = float
        elif column2method[prop] == 'none':
            pass
        else:
            prop2type[utils.add_suffix(prop, column2method[prop])] = float

    if not input_annotated_tree:
        node2leaves = annotated_tree.get_cached_content()

        # Prepare data for all nodes
        nodes_data = []
        nodes = []

        for node in annotated_tree.traverse("postorder"):
            if not node.is_leaf:
                nodes.append(node)
                node_data = (node, node2leaves[node], text_prop, multiple_text_prop, bool_prop, num_prop, column2method, alignment if 'alignment' in locals() else None, name2seq if 'name2seq' in locals() else None, emapper_mode)
                nodes_data.append(node_data)
        
        # Process nodes in parallel if more than one thread is specified
        if threads > 1:
            with Pool(threads) as pool:
                results = pool.map(process_node, nodes_data)
        else:
            # For single-threaded execution, process nodes sequentially
            results = map(process_node, nodes_data)

        # Integrate the results back into tree
        for node, result in zip(nodes, results):
            internal_props, consensus_seq = result
            for key, value in internal_props.items():
                node.add_prop(key, value)
            if consensus_seq:
                node.add_prop(alignment_prop, consensus_seq)

    else:
        pass
        
    end = time.time()
    logger.info(f'Time for merge annotations to run: {end - start}')

    # taxa annotations
    start = time.time()
    if taxon_column:
        if not taxadb:
            logger.error('Please specify which taxa db using --taxadb <GTDB|NCBI>')
            sys.exit(1)
        else:
            if taxadb == 'GTDB':
                if gtdb_version and taxa_dump:
                    logger.error('Please specify either GTDB version or taxa dump file, not both.')
                    sys.exit(1)
                if gtdb_version:
                    # get taxadump from ete-data
                    gtdbtaxadump = get_gtdbtaxadump(gtdb_version)
                    logger.info(f"Loading GTDB database dump file {gtdbtaxadump}...")
                    GTDBTaxa().update_taxonomy_database(gtdbtaxadump)
                elif taxa_dump:
                    logger.info(f"Loading GTDB database dump file {taxa_dump}...")
                    GTDBTaxa().update_taxonomy_database(taxa_dump)
                else:
                    logger.info("No specific version or dump file provided; using latest GTDB data...")
                    GTDBTaxa().update_taxonomy_database()
            elif taxadb == 'NCBI':
                if taxa_dump:
                    logger.info(f"Loading NCBI database dump file {taxa_dump}...")
                    NCBITaxa().update_taxonomy_database(taxa_dump)
                # else:
                #     NCBITaxa().update_taxonomy_database()
                
            annotated_tree, rank2values = annotate_taxa(annotated_tree, db=taxadb, \
                    taxid_attr=taxon_column, sp_delimiter=taxon_delimiter, sp_field=taxa_field, \
                    ignore_unclassified=ignore_unclassified)
                
        # evolutionary events annotation
        annotated_tree = annotate_evol_events(annotated_tree, sp_delimiter=taxon_delimiter, sp_field=taxa_field)
        prop2type.update(TAXONOMICDICT)
    else:
        rank2values = {}
    utils.clear_specific_features(annotated_tree, ['species'], leaf_only=False, internal_only=True)
    end = time.time()
    logger.info(f'Time for annotate_taxa to run: {end - start}')
    
    # prune tree by rank
    if rank_limit:
        annotated_tree = utils.taxatree_prune(annotated_tree, rank_limit=rank_limit)

    # prune tree by condition
    if pruned_by: # need to be wrap with quotes
        condition_strings = pruned_by
        annotated_tree = utils.conditional_prune(annotated_tree, condition_strings, prop2type)
    
    # name internal nodes
    annotated_tree = name_nodes(annotated_tree)
    return annotated_tree, prop2type


def run_array_annotate(tree, array_dict, num_stat='none', column2method={}):
    matrix_props = list(array_dict.keys())
    # annotate to the leaves
    for node in tree.traverse():
        if node.is_leaf:
            for filename, array in array_dict.items():
                if array.get(node.name):
                    node.add_prop(filename, array.get(node.name))


    # merge annotations to internal nodes
    for node in tree.traverse():
        if not node.is_leaf:
            for prop in matrix_props:
                # get the array from the children leaf nodes
                arrays = [child.get_prop(prop) for child in node.leaves() if child.get_prop(prop) is not None]
                
                if column2method.get(prop) is not None:
                    num_stat = column2method.get(prop)

                stats = compute_matrix_statistics(arrays, num_stat=num_stat)
                if stats:
                    for stat, value in stats.items():
                        node.add_prop(utils.add_suffix(prop, stat), value.tolist())
                        #prop2type[utils.add_suffix(prop, stat)] = float
    return tree


def run(args):
    total_color_dict = []
    layouts = []
    level = 1 # level 1 is the leaf name
    prop2type = {}
    metadata_dict = {}
    column2method = {}

    setup_logger()

    if args.metadata:
        for metadata_file in args.metadata:
            if not os.path.exists(metadata_file):
                logger.error(f"Metadata {metadata_file} does not exist.") 
                sys.exit(1)

    # Validation: Ensure at least one of --outdir or --stdout is selected
    if not args.outdir and not args.stdout:
        parser.error("You must specify either --outdir or --stdout to output results.")
    
    if args.outdir:
        if not os.path.exists(args.outdir):
            logger.error(f"Output directory {args.outdir} does not exist.") 
            sys.exit(1)
        

    # parsing tree
    try:
        tree, eteformat_flag = utils.validate_tree(args.tree, args.input_type, args.internal)
    except utils.TreeFormatError as e:
        logger.error(e)
        sys.exit(1)

    # resolve polytomy
    if args.resolve_polytomy:
        tree.resolve_polytomy()
    
    # set logger level
    if args.quiet:
        logger.setLevel(logging.CRITICAL)  # Mute all log levels below CRITICA

    logger.info(f'Loaded tree: {args.tree} \n{tree.describe()}')

    # parse csv to metadata table
    start = time.time()
    logger.info(f'start parsing...')
    # parsing metadata
    if args.metadata: # make a series of metadatas
        metadata_dict, node_props, columns, prop2type = parse_csv(args.metadata, delimiter=args.metadata_sep, \
        no_headers=args.no_headers, duplicate=args.duplicate)
    else: # annotated_tree
        node_props=[]
        columns = {}
    
    if args.data_matrix:
        array_dict = parse_tsv_to_array(args.data_matrix, delimiter=args.metadata_sep)
    end = time.time()
    logger.info(f'Time for parse_csv to run: {end - start}')
    
    if args.emapper_annotations:
        emapper_mode = True
        emapper_metadata_dict, emapper_node_props, emapper_columns = parse_emapper_annotations(args.emapper_annotations)
        metadata_dict = utils.merge_dictionaries(metadata_dict, emapper_metadata_dict)
        node_props.extend(emapper_node_props)
        columns.update(emapper_columns)
        prop2type.update({
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
    else:
        emapper_mode = False

    # start annotation
    if args.column_summary_method:
        column2method = process_column_summary_methods(args.column_summary_method)
    
    # Group metadata-related arguments
    metadata_options = {
        "metadata_dict": metadata_dict,
        "node_props": node_props,
        "columns": columns,
        "prop2type": prop2type,
        "text_prop": args.text_prop,
        "text_prop_idx": args.text_prop_idx,
        "multiple_text_prop": args.multiple_text_prop,
        "num_prop": args.num_prop,
        "num_prop_idx": args.num_prop_idx,
        "bool_prop": args.bool_prop,
        "bool_prop_idx": args.bool_prop_idx,
        "prop2type_file": args.prop2type,
    }

    # Group analysis-related arguments (ACR and Lineage Specificity options)
    analytic_options = {
        "acr_discrete_columns": args.acr_discrete_columns,
        "acr_continuous_columns": args.acr_continuous_columns,
        "prediction_method": args.prediction_method,
        "model": args.model,
        "delta_stats": args.delta_stats,
        "ent_type": args.ent_type,
        "iteration": args.iteration,
        "lambda0": args.lambda0,
        "se": args.se,
        "thin": args.thin,
        "burn": args.burn,
        "ls_columns": args.ls_columns,
        "prec_cutoff": args.prec_cutoff,
        "sens_cutoff": args.sens_cutoff,
    }

    # Group taxonomic-related arguments
    taxonomic_options = {
        "taxadb": args.taxadb,
        "gtdb_version": args.gtdb_version,
        "taxa_dump": args.taxa_dump,
        "taxon_column": args.taxon_column,
        "taxon_delimiter": args.taxon_delimiter,
        "taxa_field": args.taxa_field,
        "ignore_unclassified": args.ignore_unclassified,
    }

    # Group emapper-related arguments
    emapper_options = {
        "emapper_mode": emapper_mode,
        "emapper_pfam": args.emapper_pfam,
        "emapper_smart": args.emapper_smart,
    }

    # Group output and miscellaneous options
    output_options = {
        "rank_limit": args.rank_limit,
        "pruned_by": args.pruned_by,
        "threads": args.threads,
        "outdir": args.outdir,
    }

    # Simplified function call with grouped arguments
    annotated_tree, prop2type = run_tree_annotate(
        tree,
        input_annotated_tree=args.annotated_tree,
        **metadata_options,
        alignment=args.alignment,
        counter_stat=args.counter_stat,
        num_stat=args.num_stat,
        column2method=column2method,
        **taxonomic_options,
        **analytic_options,
        **emapper_options,
        **output_options
    )

    if args.data_matrix:
        annotated_tree = run_array_annotate(annotated_tree, array_dict, num_stat=args.num_stat, column2method=column2method)

    
    if args.outdir:
        base=os.path.splitext(os.path.basename(args.tree))[0]
        out_newick = base + '_annotated.nw'
        out_prop2tpye = base + '_prop2type.txt'
        out_ete = base+'_annotated.ete'
        out_tsv = base+'_annotated.tsv'

        
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

        ### out newick
        ## need to correct wrong symbols in the newick tree, such as ',' -> '||'
        # Find all keys where the value is of type list
        list_keys = [key for key, value in prop2type.items() if value == list]
        # Replace all commas in the tree with '||'
        list_sep = '||'
        for node in annotated_tree.leaves():
            for key in list_keys:
                if node.props.get(key):
                    list2str = list_sep.join(node.props.get(key))
                    node.add_prop(key, list2str)
        annotated_tree.write(outfile=os.path.join(args.outdir, out_newick), props=None, 
                    parser=utils.get_internal_parser(args.internal), format_root_node=True)
    
    if args.stdout:
        print(annotated_tree.write(props=None, parser=utils.get_internal_parser(args.internal), format_root_node=True))

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

def parse_csv(input_files, delimiter='\t', no_headers=False, duplicate=False):
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
                    if duplicate:
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
                            if no_headers:
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
                # Read the first line to determine the number of fields
                first_line = next(f)
                fields_len = len(first_line.split(delimiter))

                

                # Reset the file pointer to the beginning
                f.seek(0)

                if no_headers:
                    # Generate header names
                    headers = ['col'+str(i) for i in range(fields_len)]
                    # Create a CSV reader with the generated headers
                    reader = csv.DictReader(f, delimiter=delimiter, fieldnames=headers)
                else:
                    # Use the existing headers in the file
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
                            if duplicate:
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

def parse_tsv_to_array(input_files, delimiter='\t', no_headers=True):
    """
    Parses a TSV file into a dictionary with the first item of each row as the key
    and the rest of the items in the row as a list in the value.

    :param filename: Path to the TSV file to be parsed.
    :return: A dictionary with keys as the first item of each row and values as lists of the remaining items.
    """
    is_float = True
    matrix2array = {}
    
    for input_file in input_files:
        leaf2array = {}
        prefix = os.path.basename(input_file)
        with open(input_file, 'r') as file:
            for line in file:
                # Split each line by tab, strip removes trailing newline
                row = line.strip().split(delimiter)
                node = row[0]  
                value = row[1:]  # The rest of the items as value
                # Replace empty string with np.nan
                value_list = [np.nan if x == '' else x for x in value]
                try:
                    np_array = np.array(value_list).astype(np.float64)
                    leaf2array[node] = np_array.tolist()
                except ValueError:
                    # Handle the case where conversion fails
                    logger.warning(f"Warning: Non-numeric data found in {prefix} for node {node}. Skipping.")
                    leaf2array[node] = None
                    is_float = False

        matrix2array[prefix] = leaf2array        
    return matrix2array

def process_column_summary_methods(column_summary_methods):
    column_methods = {}
    if column_summary_methods:
        for entry in column_summary_methods:
            try:
                column, method = entry.split('=')
                column_methods[column] = method
            except ValueError:
                logger.error(f"Invalid format for --column-summary-method: '{entry}'. Expected format: ColumnName=Method")
                sys.exit(1)
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
        str_val = str(value).strip()  
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

def load_metadata_to_tree(tree, metadata_dict, prop2type={}, taxon_column=None, taxon_delimiter='', taxa_field=0, ignore_unclassified=False):
    #name2leaf = {}
    multi_text_seperator = ','
    common_ancestor_seperator = '||'

    name2node = defaultdict(list)
    # preload all leaves to save time instead of search in tree
    for node in tree.traverse():
        if node.name:
            name2node[node.name].append(node)

    # load all metadata to leaf nodes
    for node, props in metadata_dict.items():
        if node in name2node.keys():
            target_nodes = name2node[node]
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
            if common_ancestor_seperator in node:
                # get the common ancestor
                children = node.split(common_ancestor_seperator)
                target_node = tree.common_ancestor(children)
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
    node, node_leaves, text_prop, multiple_text_prop, bool_prop, num_prop, column2method, alignment, name2seq, emapper_mode = node_data
    internal_props = {}

    # Process text, multitext, bool, and num properties
    if text_prop:
        internal_props_text = merge_text_annotations(node_leaves, text_prop, column2method, emapper_mode)
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
    if alignment and name2seq is not None:  # Check alignment and name2seq together
        aln_sum = column2method.get('alignment')
        if aln_sum is None or aln_sum != 'none':
            matrix_string = build_matrix_string(node, name2seq)  # Assuming 'name2seq' is accessible here
            consensus_seq = utils.get_consensus_seq(matrix_string, threshold=0.7)
        
    return internal_props, consensus_seq

def merge_text_annotations(nodes, target_props, column2method, emapper_mode=False):
    pair_seperator = "--"
    item_seperator = "||"
    internal_props = {}
    counters = {}
    
    for target_prop in target_props:
        counter_stat = column2method.get(target_prop, "raw")
        prop_list = utils.children_prop_array_missing(nodes, target_prop)
        counter = dict(Counter(prop_list))  # Store the counter
        if 'NaN' in counter:
            del counter['NaN']
        counters[target_prop] = counter  # Add the counter to the counters dictionary

        if counter_stat == 'raw':
            # Find the key with the highest count
            if emapper_mode and counter:
                most_common_key = max(counter, key=counter.get)
                internal_props[target_prop] = most_common_key

            # Add the raw counts to internal_props
            internal_props[utils.add_suffix(target_prop, 'counter')] = item_seperator.join(
                [utils.add_suffix(str(key), value, pair_seperator) for key, value in sorted(counter.items())]
            )

        elif counter_stat == 'relative':
            # Find the key with the highest count
            if emapper_mode and counter:
                most_common_key = max(counter, key=counter.get)
                internal_props[target_prop] = most_common_key

            total = sum(counter.values())

            # Add the relative counts to internal_props
            internal_props[utils.add_suffix(target_prop, 'counter')] = item_seperator.join(
                [utils.add_suffix(str(key), '{0:.2f}'.format(float(value)/total), pair_seperator) for key, value in sorted(counter.items())]
            )
        elif counter_stat == 'none':
            pass
        else:
            logger.error("invalid counter_stat")
            sys.exit(1)
    return internal_props

def merge_multitext_annotations(nodes, target_props, column2method):
    # Seperator of multiple text 'GO:0000003,GO:0000902,GO:0000904'
    multi_text_seperator = ','
    pair_seperator = "--"
    item_seperator = "||"

    internal_props = {}
    counters = {}

    for target_prop in target_props:
        counter_stat = column2method.get(target_prop, "raw")
        prop_list = utils.children_prop_array(nodes, target_prop)
        
        # Flatten the list of lists into a single list
        multi_prop_list = [item for sublist in prop_list for item in sublist]
        counter = dict(Counter(multi_prop_list))  # Store the counter
        counters[target_prop] = counter  # Add the counter to the counters dictionary

        if counter_stat == 'raw':
            # Add the raw counts to internal_props
            internal_props[utils.add_suffix(target_prop, 'counter')] = item_seperator.join(
                [utils.add_suffix(str(key), value, pair_seperator) for key, value in sorted(counter.items())]
            )

        elif counter_stat == 'relative':
            total = sum(counter.values())

            # Add the relative counts to internal_props
            internal_props[utils.add_suffix(target_prop, 'counter')] = item_seperator.join(
                [utils.add_suffix(str(key), '{0:.2f}'.format(float(value) / total), pair_seperator) for key, value in sorted(counter.items())]
            )

        else:
            # Handle invalid counter_stat, if necessary
            pass

    return internal_props


def merge_num_annotations(nodes, target_props, column2method):
    internal_props = {}
    for target_prop in target_props:
        num_stat = column2method.get(target_prop, None)
        if num_stat != 'none':
            if target_prop != 'dist' and target_prop != 'support':
                prop_array = np.array(utils.children_prop_array(nodes, target_prop),dtype=np.float64)
                prop_array = prop_array[~np.isnan(prop_array)] # remove nan data
                
                
                if prop_array.any():
                    n, (smin, smax), sm, sv, ss, sk = stats.describe(prop_array)

                    if num_stat == 'all':
                        internal_props[utils.add_suffix(target_prop, 'avg')] = sm
                        internal_props[utils.add_suffix(target_prop, 'sum')] = np.sum(prop_array)
                        internal_props[utils.add_suffix(target_prop, 'max')] = smax
                        internal_props[utils.add_suffix(target_prop, 'min')] = smin
                        if math.isnan(sv) == False:
                            internal_props[utils.add_suffix(target_prop, 'std')] = sv
                        else:
                            internal_props[utils.add_suffix(target_prop, 'std')] = 0

                    elif num_stat == 'avg':
                        internal_props[utils.add_suffix(target_prop, 'avg')] = sm
                    elif num_stat == 'sum':
                        internal_props[utils.add_suffix(target_prop, 'sum')] = np.sum(prop_array)
                    elif num_stat == 'max':
                        internal_props[utils.add_suffix(target_prop, 'max')] = smax
                    elif num_stat == 'min':
                        internal_props[utils.add_suffix(target_prop, 'min')] = smin
                    elif num_stat == 'std':
                        if math.isnan(sv) == False:
                            internal_props[utils.add_suffix(target_prop, 'std')] = sv
                        else:
                            internal_props[utils.add_suffix(target_prop, 'std')] = 0
                    else:
                        #print('Invalid stat method')
                        pass
                else:
                    pass

    if internal_props:
        return internal_props
    else:
        return None

def compute_matrix_statistics(matrix, num_stat=None):
    """
    Computes specified statistics for the given matrix based on the num_stat parameter.
    
    :param matrix: A list of lists representing the matrix.
    :param num_stat: Specifies which statistics to compute. Can be "avg", "max", "min", "sum", "std", "all", or None.
    :return: A dictionary with the requested statistics or an empty dict/message.
    """
    
    stats = {}

    if num_stat == 'none':
        return stats  # Return an empty dictionary if no statistics are requested
    
    # Replace None with np.nan or another appropriate value before creating the array
    if matrix is not None:
        cleaned_matrix = [[0 if x is None else x for x in row] for row in matrix]
        np_matrix = np.array(cleaned_matrix, dtype=np.float64)
    else:
        return {}  # Return an empty dictionary if the matrix is empty

    if np_matrix.size == 0:
        return {}  # Return an empty dictionary if the matrix is empty

 
    if np_matrix.ndim == 2 and np_matrix.shape[1] > 0:
        available_stats = {
            'avg': np_matrix.mean(axis=0),
            'max': np_matrix.max(axis=0),
            'min': np_matrix.min(axis=0),
            'sum': np_matrix.sum(axis=0),
            'std': np_matrix.std(axis=0)
        }

        if num_stat == "all":
            return available_stats
        elif num_stat in available_stats:
            stats[num_stat] = available_stats[num_stat]
        else:
            logger.error(f"Unsupported stat '{num_stat}'. Supported stats are 'avg', 'max', 'min', 'sum', 'std', or 'all'.")
            sys.exit(1)
    return stats

def name_nodes(tree):
    for i, node in enumerate(tree.traverse("postorder")):
        if not node.name or node.name == 'None':
            if not node.is_root:
                node.name = 'N'+str(i)
            else:
                node.name = 'Root'
    return tree

def gtdb_accession_to_taxid(accession):
        """Given a GTDB accession number, returns its complete accession"""
        if accession.startswith('GCA'):
            prefix = 'GB_'
            return prefix+accessionac
        elif accession.startswith('GCF'):
            prefix = 'RS_'
            return prefix+accession
        else:
            return accession

def get_gtdbtaxadump(version):
    """
    Download GTDB taxonomy dump
    """
    url = f"https://github.com/etetoolkit/ete-data/raw/main/gtdb_taxonomy/gtdb{version}/gtdb{version}dump.tar.gz"
    fname = f"gtdb{version}dump.tar.gz"
    logger.info(f'Downloading GTDB taxa dump fname from {url} ...')
    with open(fname, 'wb') as f:
        f.write(requests.get(url).content)
    return fname

def annotate_taxa(tree, db="GTDB", taxid_attr="name", sp_delimiter='.', sp_field=0, ignore_unclassified=False):
    global rank2values
    logger.info(f"\n==============Annotating tree with {db} taxonomic database============")
    
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

    def merge_dictionaries(dict_ranks, dict_names):
        """
        Merges two dictionaries into one where the key is the rank from dict_ranks 
        and the value is the corresponding name from dict_names.

        :param dict_ranks: Dictionary where the key is a numeric id and the value is a rank.
        :param dict_names: Dictionary where the key is the same numeric id and the value is a name.
        :return: A new dictionary where the rank is the key and the name is the value.
        """
        merged_dict = {}
        for key, rank in dict_ranks.items():
            if key in dict_names:  # Ensure the key exists in both dictionaries
                if rank not in merged_dict or rank == 'no rank':  # Handle 'no rank' by not overwriting existing entries unless it's the first encounter
                    merged_dict[rank] = dict_names[key]

        return merged_dict


    if db == "GTDB":
        gtdb = GTDBTaxa()
        tree.set_species_naming_function(return_spcode_gtdb)
        gtdb.annotate_tree(tree,  taxid_attr="species", ignore_unclassified=ignore_unclassified)
        suffix_to_rank_dict = {
            'd__': 'superkingdom',  # Domain or Superkingdom
            'p__': 'phylum',
            'c__': 'class',
            'o__': 'order',
            'f__': 'family',
            'g__': 'genus',
            's__': 'species'
        }
        gtdb_re = r'^(GB_GCA_[0-9]+\.[0-9]+|RS_GCF_[0-9]+\.[0-9]+)'
        for n in tree.traverse():
            # in case miss something
            if n.props.get('named_lineage'):
                lca_dict = {}
                for taxa in n.props.get("named_lineage"):
                    if re.match(gtdb_re, taxa):
                        potential_rank = 'subspecies'
                        lca_dict[potential_rank] = n.props.get("sci_name")
                    else:
                        potential_rank = suffix_to_rank_dict.get(taxa[:3], None)
                        if potential_rank:
                            lca_dict[potential_rank] = taxa
                n.add_prop("lca", utils.dict_to_string(lca_dict))

    elif db == "NCBI":
        ncbi = NCBITaxa()
        # extract sp codes from leaf names
        tree.set_species_naming_function(return_spcode_ncbi)
        ncbi.annotate_tree(tree, taxid_attr="species", ignore_unclassified=ignore_unclassified)
        for n in tree.traverse():
            if n.props.get('lineage'):
                lca_dict = {}
                #for taxa in n.props.get("lineage"):
                lineage2rank = ncbi.get_rank(n.props.get("lineage"))
                taxid2name = ncbi.get_taxid_translator(n.props.get("lineage"))
                lca_dict = merge_dictionaries(lineage2rank, taxid2name)
                n.add_prop("named_lineage", list(taxid2name.values()))
                n.add_prop("lca", utils.dict_to_string(lca_dict))

    # tree.annotate_gtdb_taxa(taxid_attr='name')
    # assign internal node as sci_name
    rank2values = defaultdict(list)
    for n in tree.traverse():
        if db == 'NCBI':
            n.del_prop('_speciesFunction')
        if n.props.get('rank') and n.props.get('rank') != 'Unknown':
            rank2values[n.props.get('rank')].append(n.props.get('sci_name',''))

        # if n.name:
        #     pass
        # else:
        #     n.name = n.props.get("sci_name", "")
        
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

def parse_emapper_annotations(input_file, delimiter='\t', no_headers=False):
    metadata = {}
    columns = defaultdict(list)
    prop2type = {}
    # EMAPPER_HEADERS = ["#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
    #            "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs",
    #            "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
    #            "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]

    with open(input_file, 'r') as f:
        # Skip lines starting with '##'
        filtered_lines = (line for line in f if not line.startswith('##'))

        if no_headers:
            reader = csv.DictReader(filtered_lines, delimiter=delimiter, fieldnames=headers)
        else:
            reader = csv.DictReader(filtered_lines, delimiter=delimiter)

        node_header, node_props = EMAPPER_HEADERS[0], EMAPPER_HEADERS[1:]
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
                    try:
                        trans_dom_start = raw2alg[seq_name][dom_start]
                        trans_dom_end = raw2alg[seq_name][dom_end]
                        dom_info_string = pair_delimiter.join([dom_name, str(trans_dom_start), str(trans_dom_end)])
                        seq2doms[seq_name].append(dom_info_string)
                    except KeyError:
                        logger.error(f"Cannot find {dom_start} or {dom_end} in {seq_name}")
                        sys.exit(1)

    for l in post_tree:
        if l.name in seq2doms.keys():
            domains = seq2doms[l.name]
            domains_string = item_seperator.join(domains)
            l.add_prop(domain_prop, domains_string)

    for n in post_tree.traverse():
        # get the most common domain
        if not n.is_leaf:
            prop_list = utils.children_prop_array(n, domain_prop)
            counter = dict(Counter(prop_list))
            most_common_key = max(counter, key=counter.get)
            n.add_prop(domain_prop, most_common_key)

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
        # get the most common domain
        if not n.is_leaf:
            prop_list = utils.children_prop_array(n, domain_prop)
            counter = dict(Counter(prop_list))
            most_common_key = max(counter, key=counter.get)
            n.add_prop(domain_prop, most_common_key)

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

def _worker_function(iteration_data):
    # Unpack the necessary data for one iteration
    prop2array, dump_tree, acr_discrete_columns_dict, prediction_method, model, ent_type, lambda0, se, sim, burn, thin, threads = iteration_data

    shuffled_dict = {}
    for column, trait in acr_discrete_columns_dict.items():
        # Shuffle traits
        shuffled_trait = np.random.choice(trait, len(trait), replace=False)
        prop2array[column][1] = list(shuffled_trait)
        shuffled_dict[column] = list(shuffled_trait)

    # Converting back to the original dictionary format
    new_metadata_dict = convert_back_to_original(prop2array)
    updated_tree = load_metadata_to_tree(dump_tree, new_metadata_dict)

    # Run ACR
    random_acr_results, updated_tree = run_acr_discrete(updated_tree, shuffled_dict, 
                                                        prediction_method=prediction_method, 
                                                        model=model, threads=threads, outdir=None)
    random_delta = run_delta(random_acr_results, updated_tree, ent_type=ent_type, 
                             lambda0=lambda0, se=se, sim=sim, burn=burn, thin=thin, 
                             threads=threads)

    # Clear extra features from the tree
    utils.clear_extra_features([updated_tree], ["name", "dist", "support"])
    return random_delta
    
def get_pval(prop2array, dump_tree, acr_discrete_columns_dict, iteration=100, 
             prediction_method="MPPA", model="F81", ent_type='SE', 
             lambda0=0.1, se=0.5, sim=10000, burn=100, thin=10, threads=1):
    prop2delta_array = {}

    # Prepare data for each iteration
    iteration_data = [(prop2array, dump_tree, acr_discrete_columns_dict, prediction_method, model, ent_type, lambda0, se, sim, burn, thin, threads) for _ in range(iteration)]

    # Use multiprocessing pool
    if threads > 1:
        with Pool(threads) as pool:
            results = pool.map(_worker_function, iteration_data)
    else:
        results = map(_worker_function, iteration_data)

    # Aggregate results
    for delta_result in results:
        for prop, result in delta_result.items():
            if prop in prop2delta_array:
                prop2delta_array[prop].append(result)
            else:
                prop2delta_array[prop] = [result]

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
