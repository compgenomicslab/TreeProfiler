#!/usr/bin/env python3

import sys
import os
import time
import random
import math
import re
import csv
import tarfile
from io import StringIO
from itertools import islice
from collections import defaultdict, Counter
import numpy as np
from scipy import stats


from ete4.parser.newick import NewickError
from ete4.coretype.seqgroup import SeqGroup
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from treeprofiler.src import b64pickle
from treeprofiler.src.utils import (
    ete4_parse, taxatree_prune, conditional_prune,
    children_prop_array, children_prop_array_missing,
    flatten, get_consensus_seq)

DESC = "annotate tree"

def populate_annotate_args(parser):
    gmeta = parser.add_argument_group(
        title='METADATA TABLE parameters',
        description="Input parameters of METADATA")

    add = gmeta.add_argument
    add('-d', '--metadata', nargs='+',
        help="<metadata.csv> .csv, .tsv. mandatory input")
    add('--no_colnames', action='store_true',
        help="metadata table doesn't contain columns name")
    add('--aggregate-duplicate', action='store_true',
        help="treeprofiler will aggregate duplicated metadata to a list as a property if metadata contains duplicated row")
    add('--text-prop', nargs='+',
        help=("<col1,col2> names, column index or index range of columns which "
              "need to be read as categorical data"))
    add('--multiple-text-prop', nargs='+',
        help=("<col1,col2> names, column index or index range of columns which "
              "need to be read as categorical data which contains more than one"
              " value and seperate by ',' such "
              "as GO:0000003,GO:0000902,GO:0000904,GO:0003006"))
    add('--num-prop', nargs='+',
        help=("<col1,col2> names, column index or index range of columns which "
              "need to be read as numerical data"))
    add('--bool-prop', nargs='+',
        help=("<col1,col2> names, column index or index range of columns which "
              "need to be read as boolean data"))
    add('--text-prop-idx',
        help="1,2,3 or [1-5] index of columns which need to be read as categorical data")
    add('--num-prop-idx',
        help="1,2,3 or [1-5] index columns which need to be read as numerical data")
    add('--bool-prop-idx',
        help="1,2,3 or [1-5] index columns which need to be read as boolean data")
    # add('--taxatree',
    #     help=("<kingdom|phylum|class|order|family|genus|species|subspecies> "
    #           "reference tree from taxonomic database"))
    add('--taxadb', default='GTDB',
        help="<NCBI|GTDB> for taxonomic profiling or fetch taxatree [default: GTDB]")
    add('--taxon-column',
        help="<col1> name of columns which need to be read as taxon data")
    add('--taxon-delimiter', default='',
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

    group = parser.add_argument_group(title='Annotation arguments',
        description="Annotation parameters")
    group.add_argument('--taxonomic-profile',
        default=False,
        action='store_true',
        required=False,
        help="Activate taxonomic annotation on tree")
    group.add_argument('--num-stat',
        default='all',
        choices=['all', 'sum', 'avg', 'max', 'min', 'std', 'none'],
        type=str,
        required=False,
        help="statistic calculation to perform for numerical data in internal nodes, [all, sum, avg, max, min, std, none]. If 'none' was chosen, numerical properties won't be summarized nor annotated in internal nodes. [default: all]")  

    group.add_argument('--counter-stat',
        default='raw',
        choices=['raw', 'relative', 'none'],
        type=str,
        required=False,
        help="statistic calculation to perform for categorical data in internal nodes, raw count or in percentage [raw, relative, none]. If 'none' was chosen, categorical and boolean properties won't be summarized nor annotated in internal nodes [default: raw]")  
    
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
        emapper_smart=None, counter_stat='raw', num_stat='all',
        taxonomic_profile=False, taxadb='GTDB', taxon_column='name',
        taxon_delimiter='', taxa_field=0, rank_limit=None, pruned_by=None,
        outdir='./'):

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
                    # if dtype == bool:
                    #     bool_prop.append(key)

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

        for prop in text_prop+bool_prop:
            prop2type[prop] = str
            prop2type[prop+'_counter'] = str

        for prop in multiple_text_prop:
            prop2type[prop] = list
            prop2type[prop+'_counter'] = str

        for prop in num_prop:
            prop2type[prop] = float
            prop2type[prop+'_avg'] = float
            prop2type[prop+'_sum'] = float
            prop2type[prop+'_max'] = float
            prop2type[prop+'_min'] = float
            prop2type[prop+'_std'] = float

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

    # domain annotation before other annotation
    if emapper_pfam:
        if not alignment:
            raise ValueError("Please provide alignment file using '--alignment' for pfam annotation.")
        annot_tree_pfam_table(tree, emapper_pfam, alignment)

    if emapper_smart:
        if not alignment:
            raise ValueError("Please provide alignment file using '--alignment' for smart annotation.")
        annot_tree_smart_table(tree, emapper_smart, alignment)


    # load all metadata to leaf nodes
    #taxon_column = []

    # input_annotated_tree determines if input tree is already annotated, if annotated, no longer need metadata
    if not input_annotated_tree:
        if taxon_column: # to identify taxon column as taxa property from metadata
            #taxon_column.append(taxon_column)
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type, taxon_column=taxon_column, taxon_delimiter=taxon_delimiter, taxa_field=taxa_field)
        else:
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type)
    else:
        annotated_tree = tree

    end = time.time()
    print('Time for load_metadata_to_tree to run: ', end - start)

    # stat method
    counter_stat = counter_stat #'raw' or 'relative'
    num_stat = num_stat

    # merge annotations depends on the column datatype
    start = time.time()
    if not input_annotated_tree:

        #pre load node2leaves to save time
        node2leaves = annotated_tree.get_cached_content()
        for i, node in enumerate(annotated_tree.traverse("postorder")):
            internal_props = {}
            if not node.is_leaf:
                if counter_stat != 'none':
                    if text_prop:
                        internal_props_text = merge_text_annotations(node2leaves[node], text_prop, counter_stat=counter_stat)
                        internal_props.update(internal_props_text)

                    if multiple_text_prop:
                        internal_props_multi = merge_multitext_annotations(node2leaves[node], multiple_text_prop, counter_stat=counter_stat)
                        internal_props.update(internal_props_multi)

                    if bool_prop:
                        internal_props_bool = merge_text_annotations(node2leaves[node], bool_prop, counter_stat=counter_stat)
                        internal_props.update(internal_props_bool)

                if num_stat != 'none':
                    if num_prop:
                        internal_props_num = merge_num_annotations(node2leaves[node], num_prop, num_stat=num_stat)                        
                        if internal_props_num:
                            internal_props.update(internal_props_num)

                # deprecated
                # if rest_column:
                #     internal_props_rest = merge_text_annotations(node2leaves[node], rest_column, counter_stat=counter_stat)
                #     internal_props.update(internal_props_rest)

                #internal_props = {**internal_props_text, **internal_props_num, **internal_props_rest}
                #print(internal_props.keys())

                for key,value in internal_props.items():
                    node.add_prop(key, value)

                if alignment:
                    matrix = ''
                    for leaf in node.leaves():
                        if name2seq.get(leaf.name):
                            matrix += ">"+leaf.name+"\n"
                            matrix += name2seq.get(leaf.name)+"\n"
                    consensus_seq = get_consensus_seq(StringIO(matrix), 0.7)
                    node.add_prop(alignment_prop, consensus_seq)

    else:
        pass
    end = time.time()
    print('Time for merge annotations to run: ', end - start)


    # taxa annotations
    start = time.time()
    if taxonomic_profile:
        # taxonomic annotation
        if not taxadb:
            print('Please specify which taxa db using --taxadb <GTDB|NCBI>')
        else:
            if taxadb == 'GTDB':
                if taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=taxadb, taxid_attr=taxon_column)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=taxadb, taxid_attr="name")
            elif taxadb == 'NCBI':
                if taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=taxadb, taxid_attr=taxon_column, sp_delimiter=taxon_delimiter, sp_field=taxa_field)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=taxadb, taxid_attr="name", sp_delimiter=taxon_delimiter, sp_field=taxa_field)
        
        # evolutionary events annotation
        annotated_tree = annotate_evol_events(annotated_tree, sp_delimiter=taxon_delimiter, sp_field=taxa_field)
        
        prop2type.update({# start with leaf name
                'rank': str,
                'sci_name': str,
                'taxid': str,
                'lineage':str,
                'named_lineage': str,
                'evoltype': str,
                'dup_sp': str,
                'dup_percent': float,
                })
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
    if args.tree:
        if args.input_type == 'newick':
            try:
                tree = ete4_parse(open(args.tree), internal_parser=args.internal_parser)
            except Exception as e:
                print(e)
                sys.exit(1)
        elif args.input_type == 'ete':
            try:
                with open(args.tree, 'r') as f:
                    file_content = f.read()
                    tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
            except ValueError as e:
                print(e)
                print("In valid ete format.")
                sys.exit(1)
    # if refer tree from taxadb, input tree will be ignored
    elif taxatree and taxadb:
        tree = ''
    else:
        sys.exit('empty input')

    # parse csv to metadata table
    start = time.time()
    print("start parsing...")
    # parsing metadata
    if args.metadata: # make a series aof metadatas
        metadata_dict, node_props, columns, prop2type = parse_csv(args.metadata, no_colnames=args.no_colnames, aggregate_duplicate=args.aggregate_duplicate)
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
    annotated_tree, prop2type = run_tree_annotate(tree, input_annotated_tree=args.annotated_tree,
        metadata_dict=metadata_dict, node_props=node_props, columns=columns,
        prop2type=prop2type,
        text_prop=args.text_prop, text_prop_idx=args.text_prop_idx,
        multiple_text_prop=args.multiple_text_prop, num_prop=args.num_prop, num_prop_idx=args.num_prop_idx,
        bool_prop=args.bool_prop, bool_prop_idx=args.bool_prop_idx,
        prop2type_file=args.prop2type, alignment=args.alignment,
        emapper_pfam=args.emapper_pfam, emapper_smart=args.emapper_smart,
        counter_stat=args.counter_stat, num_stat=args.num_stat,
        taxonomic_profile=args.taxonomic_profile, taxadb=args.taxadb, taxon_column=args.taxon_column,
        taxon_delimiter=args.taxon_delimiter, taxa_field=args.taxa_field,
        rank_limit=args.rank_limit, pruned_by=args.pruned_by, outdir=args.outdir)

    if args.outdir:
        base=os.path.splitext(os.path.basename(args.tree))[0]
        out_newick = base + '_annotated.nw'
        out_prop2tpye = base + '_prop2type.txt'
        out_ete = base+'_annotated.ete'
        out_tsv = base+'_annotated.tsv'

        ### out newick
        annotated_tree.write(outfile=os.path.join(args.outdir, out_newick), props=None, format_root_node=True)
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
        if args.taxonomic_profile:
            prop_keys.extend([
                'rank',
                'sci_name',
                'taxid',
                'lineage',
                'named_lineage'
            ])
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

    if re.match(pattern, input_string):
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

# def get_type_convert(np_type):
#     """
#     convert np_type to python type
#     """
#     convert_type = type(np.zeros(1,np_type).tolist()[0])
#     return (np_type, convert_type)

def get_comma_separated_values(lst):
    for item in lst:
        if isinstance(item, str) and any(',' in x for x in item.split()):
            return True
    return False

def convert_column_data(column, np_dtype):
    #np_dtype = np.dtype(dtype).type
    try:
        data = np.array(column).astype(np_dtype)
        return np_dtype
    except ValueError:
        return None
    #data.astype(np.float)

def infer_dtype(column):
    if get_comma_separated_values(column):
        return list
    else:
        dtype_dict = {
            float:np.float64,
            bool:np.bool_,
            str:np.str_
            }
        #dtype_order = ['float64', 'bool', 'str']
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
                    elif key in prop2type and prop2type[key]==list:
                        value_list = value.split(multi_text_seperator)
                        target_node.add_prop(key, value_list)
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

def merge_text_annotations(nodes, target_props, counter_stat='raw'):
    pair_seperator = "--"
    item_seperator = "||"
    internal_props = {}
    for target_prop in target_props:
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
            print('Invalid stat method')
            break

    return internal_props

def merge_multitext_annotations(nodes, target_props, counter_stat='raw'):
    #seperator of multiple text 'GO:0000003,GO:0000902,GO:0000904'
    multi_text_seperator = ','
    pair_seperator = "--"
    item_seperator = "||"

    internal_props = {}
    for target_prop in target_props:
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
                    print('Invalid stat method')
                    pass
            else:
                pass

    if internal_props:
        return internal_props
    else:
        return None


def add_suffix(name, suffix, delimiter='_'):
    return str(name) + delimiter + str(suffix)

def name_nodes(tree):
    for i, node in enumerate(tree.traverse("postorder")):
        if not node.name:
            if not node.is_root:
                node.name = 'N'+str(i)
            else:
                node.name = 'Root'
    return tree

def annotate_taxa(tree, db="GTDB", taxid_attr="name", sp_delimiter='.', sp_field=0):
    global rank2values
    def return_spcode(leaf):
        #print(leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field])
        try:
            return leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field]
        except (IndexError, ValueError):
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

def annotate_evol_events(tree, sp_delimiter='.', sp_field=0):
    def return_spcode(leaf):
        try:
            return leaf.name.split(sp_delimiter)[sp_field]
        except (IndexError, ValueError):
            return leaf.name

    tree.set_species_naming_function(return_spcode)

    node2species = tree.get_cached_content(store_attr='species')
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
    headers = ["#query", "seed_ortholog", "evalue",	"score","eggNOG_OGs",
            "max_annot_lvl","COG_category","Description","Preferred_name","GOs",
            "EC","KEGG_ko","KEGG_Pathway",	"KEGG_Module", "KEGG_Reaction",	"KEGG_rclass",
            "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]

    lines_count = len(open(input_file).readlines())

    with open(input_file, 'r') as f:

        if no_colnames:
            reader = csv.DictReader(f, delimiter=delimiter, fieldnames=headers)
        else:
            lines_count = len(f.readlines())
            skip_header = 4
            skip_footer = 3
            f.seek(0) # using f twice
            reader = csv.DictReader(islice(f,skip_header,lines_count-skip_footer), delimiter=delimiter)

        node_header, node_props = headers[0], headers[1:]
        for row in reader:
            nodename = row[node_header]
            del row[node_header]

            for k, v in row.items(): # replace missing value

                if check_missing(v):
                    row[k] = 'NaN'
                else:
                    row[k] = v
            metadata[nodename] = dict(row)
            for (k,v) in row.items(): # go over each column name and value
                columns[k].append(v) # append the value into the appropriate list
                                    # based on column name k

    return metadata, node_props, columns

def annot_tree_pfam_table(post_tree, pfam_table, alg_fasta):
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
            l.add_prop('dom_arq', domains_string)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_node_domains = n.get_closest_leaf()[0].props.get('dom_arq', 'none@none@none')
            n.add_prop('dom_arq', random_node_domains)

    # for n in post_tree.traverse():
    #     print(n.name, n.props.get('dom_arq'))

def annot_tree_smart_table(post_tree, smart_table, alg_fasta):
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
            l.add_prop('dom_arq', domains_string)

    for n in post_tree.traverse():
        if not n.is_leaf:
            random_node_domains = n.get_closest_leaf()[0].props.get('dom_arq', 'none@none@none')
            n.add_prop('dom_arq', random_node_domains)

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
    return
