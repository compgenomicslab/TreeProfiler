#!/usr/bin/env python
from ete4.parser.newick import NewickError
from ete4.coretype.seqgroup import SeqGroup
from ete4 import Tree, PhyloTree
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from utils import (
    ete4_parse, taxatree_prune, conditional_prune,
    children_prop_array, children_prop_array_missing, 
    flatten, get_consensus_seq)
import b64pickle

from collections import defaultdict
from collections import Counter
import numpy as np
from io import StringIO
from scipy import stats
from itertools import islice
import random
import math
import csv
import time
import re
import os

DESC = "annotate tree"

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

def run(args):
    total_color_dict = []
    layouts = []
    level = 2 # level 1 is the leaf name
    prop2type = {}
    metadata_dict = {}

    # TODO: I want to change this into something like this:
    #     tree, metada = args
    #     run_tree_annotate(tree, metadata)
    
    # parse csv to metadata table
    start = time.time()
    print("start parsing...")
    if args.metadata: # make a series aof metadatas
        if args.no_colnames:
            # property key will be named col1, col2, col3, ... if without headers
            metadata_dict, node_props, columns, prop2type = parse_csv(args.metadata, no_colnames=args.no_colnames)
        else:
            metadata_dict, node_props, columns, prop2type = parse_csv(args.metadata)
    else: # annotated_tree
        node_props=[]
        columns = {}
        
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

    end = time.time()
    print('Time for parse_csv to run: ', end - start)
    
    # parse input tree
    if args.tree:
        if args.tree_type == 'newick':
            tree = ete4_parse(args.tree)
        elif args.tree_type == 'ete':
            with open(args.tree, 'r') as f:
                file_content = f.read()
                tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
    # if refer tree from taxadb, input tree will be ignored
    elif args.taxatree and args.taxadb:
        tree = ''
    else:
        sys.exit('empty input')

    if args.text_prop:
        text_prop = args.text_prop
    else:
        text_prop = []

    if args.multiple_text_prop:
        multiple_text_prop = args.multiple_text_prop
    else:
        multiple_text_prop = []

    if args.num_prop:
        num_prop = args.num_prop
    else:
        num_prop = []
    
    if args.bool_prop:
        bool_prop = args.bool_prop
    else:
        bool_prop = []

    if args.text_prop_idx:
        text_prop_idx = []
        for i in args.text_prop_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                text_prop_start, text_prop_end = get_range(i)
                for j in range(text_prop_start, text_prop_end+1):
                    text_prop_idx.append(j)
            else:
                text_prop_idx.append(int(i))

        text_prop = [node_props[index-1] for index in text_prop_idx]
    
    if args.num_prop_idx:
        num_prop_idx = []
        for i in args.num_prop_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                num_prop_start, num_prop_end = get_range(i)
                for j in range(num_prop_start, num_prop_end+1):
                    num_prop_idx.append(j)
            else:
                num_prop_idx.append(int(i))

        num_prop = [node_props[index-1] for index in num_prop_idx]

    if args.bool_prop_idx:
        bool_prop_idx = []
        for i in args.bool_prop_idx.split(','):
            if i[0] == '[' and i[-1] == ']':
                bool_prop_start, bool_prop_end = get_range(i)
                for j in range(bool_prop_start, bool_prop_end+1):
                    bool_prop_idx.append(j)
            else:
                bool_prop_idx.append(int(i))

        bool_prop = [node_props[index-1] for index in bool_prop_idx]

    #rest_prop = []
    if args.prop2type:
        prop2type = {}
        with open(args.prop2type, 'r') as f:
            for line in f:
                line = line.rstrip()
                prop, value = line.split('\t')
                prop2type[prop] = value 
    else:
        # output datatype of each property of each tree node including internal nodes
        if prop2type:
            for key, dtype in prop2type.items():
                if key in text_prop+multiple_text_prop+num_prop+bool_prop:
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
        if args.emapper_annotations:
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
    if args.alignment:
        alignment_prop = 'alignment'
        name2seq = parse_fasta(args.alignment)
        for leaf in tree.iter_leaves():
            leaf.add_prop(alignment_prop, name2seq.get(leaf.name,''))
    
    # domain annotation before other annotation
    if args.emapper_pfam:
        annot_tree_pfam_table(tree, args.emapper_pfam, args.alignment)
    
    if args.emapper_smart:
        annot_tree_smart_table(tree, args.emapper_smart, args.alignment)
    

    # load all metadata to leaf nodes
    taxon_column = []
    if not args.annotated_tree:
        if args.taxon_column: # to identify taxon column as taxa property from metadata
            taxon_column.append(args.taxon_column)
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type, taxon_column=args.taxon_column, taxon_delimiter=args.taxon_delimiter)
        else:
            annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type)
    else:
        annotated_tree = tree

    end = time.time()
    print('Time for load_metadata_to_tree to run: ', end - start)
    
    # stat method
    counter_stat = args.counter_stat #'raw' or 'relative'
    num_stat = args.num_stat

    # merge annotations depends on the column datatype
    start = time.time()
    if not args.annotated_tree:
        #annotate_tree(t)
        
        #pre load node2leaves to save time
        node2leaves = annotated_tree.get_cached_content()
        count = 0
        for node in annotated_tree.traverse("postorder"):
            internal_props = {}
            if not node.is_leaf():
                if text_prop:
                    internal_props_text = merge_text_annotations(node2leaves[node], text_prop, counter_stat=counter_stat)
                    internal_props.update(internal_props_text)

                if multiple_text_prop:
                    internal_props_multi = merge_multitext_annotations(node2leaves[node], multiple_text_prop, counter_stat=counter_stat)
                    internal_props.update(internal_props_multi)

                if num_prop:
                    internal_props_num = merge_num_annotations(node2leaves[node], num_prop, num_stat=num_stat)                        
                    if internal_props_num:
                        internal_props.update(internal_props_num)

                if bool_prop:
                    internal_props_bool = merge_text_annotations(node2leaves[node], bool_prop, counter_stat=counter_stat)
                    internal_props.update(internal_props_bool)

                
                # deprecated
                # if rest_column:
                #     internal_props_rest = merge_text_annotations(node2leaves[node], rest_column, counter_stat=counter_stat)
                #     internal_props.update(internal_props_rest)
                
                #internal_props = {**internal_props_text, **internal_props_num, **internal_props_rest}
                #print(internal_props.keys())

                for key,value in internal_props.items():
                    node.add_prop(key, value)

                if args.alignment:
                    matrix = ''
                    for leaf in node.iter_leaves():
                        matrix += ">"+leaf.name+"\n"
                        matrix += name2seq.get(leaf.name, '')+"\n"
                    consensus_seq = get_consensus_seq(StringIO(matrix), 0.7)
                    node.add_prop(alignment_prop, consensus_seq)

    else:
        pass
    end = time.time()
    print('Time for merge annotations to run: ', end - start)
    
    
    # taxa annotations
    start = time.time()
    if args.taxonomic_profile:
        if not args.taxadb:
            print('Please specify which taxa db using --taxadb <GTDB|NCBI>')
        else:
            if args.taxadb == 'GTDB':
                if args.taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr=args.taxon_column)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name")
            elif args.taxadb == 'NCBI':
                if args.taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr=args.taxon_column, sp_delimiter=args.taxon_delimiter, sp_field=args.taxa_field)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name", sp_delimiter=args.taxon_delimiter, sp_field=args.taxa_field)
    else:
        rank2values = {}
    
    end = time.time()
    print('Time for annotate_taxa to run: ', end - start)

    # prune tree by rank
    if args.rank_limit:
        annotated_tree = taxatree_prune(annotated_tree, rank_limit=args.rank_limit)

    # prune tree by condition 
    if args.pruned_by: # need to be wrap with quotes
        condition_strings = args.pruned_by
        annotated_tree = conditional_prune(annotated_tree, condition_strings, prop2type)

    # output tree
    # if args.out_colordict:
    #     wrtie_color(total_color_dict)
    
    if args.ete4out:
        base=os.path.splitext(os.path.basename(args.tree))[0]
        with open(base+'_annotated.ete', 'w') as f:
            f.write(b64pickle.dumps(annotated_tree, encoder='pickle', pack=False))
    
    if args.outdir:
        base=os.path.splitext(os.path.basename(args.tree))[0]
        out_newick = base + '_annotated.nw'
        out_prop2tpye = base + '_prop2type.txt'
        out_ete = base+'_annotated.ete'
        out_tsv = base+'_annotated.tsv'

        ### out newick
        annotated_tree.write(outfile=os.path.join(args.outdir, out_newick), properties = [], format=1)
        ### output prop2type
        with open(os.path.join(args.outdir, base+'_prop2type.txt'), "w") as f:
            #f.write(first_line + "\n")
            for key, value in prop2type.items():
                f.write("{}\t{}\n".format(key, value))
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
        tree2table(annotated_tree, internal_node=True, props=prop_keys, outfile=os.path.join(args.outdir, out_tsv))
                
    # if args.outtsv:
    #     tree2table(annotated_tree, internal_node=True, outfile=args.outtsv)
    return 

def check_missing(input_string):

    pattern = r'^(?:\W+|none|None|null|NaN|)$'

    if re.match(pattern, input_string):
        #print("Input contains only non-alphanumeric characters, 'none', a missing value, or an empty value.")
        return True
    else:
        return False
    
def parse_csv(input_files, delimiter='\t', no_colnames=False):
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
    for input_file in input_files:
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
                        if prop in metadata[nodename]:
                            exisiting_value = metadata[nodename][prop]
                            new_value = ','.join([exisiting_value,value])
                            metadata[nodename][prop] = new_value
                            columns[prop].append(new_value)
                        else:
                            metadata[nodename][prop] = value
                            columns[prop].append(value)
                else:
                    metadata[nodename] = dict(row)
                    for (prop, value) in row.items(): # go over each column name and value             
                        columns[prop].append(value) # append the value into the appropriate list
                                        # based on column name k

        for prop in node_props:
            if set(columns[prop])=={'NaN'}:
                #prop2type[prop] = np.str_
                prop2type[prop] = str
            else:
                dtype = infer_dtype(columns[prop])
                prop2type[prop] = dtype # get_type_convert(dtype)
    return metadata, node_props, columns, prop2type

def get_type_convert(np_type):
    """
    convert np_type to python type
    """
    convert_type = type(np.zeros(1,np_type).tolist()[0])
    return (np_type, convert_type)

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

# def multiple_text_profile(tree, profiling_prop):
#     all_gos = children_prop_array(tree, profiling_prop)
#     all_gos = flatten(all_gos)
    
#     for go in all_gos:
#         for n in tree.iter_leaves():
#             print(n.props.get(profiling_prop))
#     return

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

def load_metadata_to_tree(tree, metadata_dict, prop2type={}, taxon_column=None, taxon_delimiter=';'):
    #name2leaf = {}
    multi_text_seperator = ','

    name2leaf = defaultdict(list)
    # preload all leaves to save time instead of search in tree
    for leaf in tree.iter_leaves():
        name2leaf[leaf.name].append(leaf)

    # load all metadata to leaf nodes
    for node, props in metadata_dict.items():
        if node in name2leaf.keys():
            target_nodes = name2leaf[node]
            for target_node in target_nodes:
                for key,value in props.items():
                    # taxa
                    if key == taxon_column:
                        taxon_prop = value.split(taxon_delimiter)[-1]
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

# def children_prop_array(nodes, prop):
#     #array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
#     array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
#     return array

# def children_prop_array_missing(nodes, prop):
#     array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
#     #array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
#     return array

def annotate_taxa(tree, db="GTDB", taxid_attr="name", sp_delimiter='.', sp_field=0):
    global rank2values
    def return_spcode(leaf):
        #print(leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field])
        try:
            return leaf.props.get(taxid_attr).split(sp_delimiter)[sp_field]
        except IndexError:
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

def get_range(input_range):
    column_range = input_range[input_range.find("[")+1:input_range.find("]")]
    column_start, column_end = [int(i) for i in column_range.split('-')]
    #column_list_idx = [i for i in range(column_start, column_end+1)]
    return column_start, column_end

### emapper annotate tree
def tree_emapper_annotate(args):
    print("start mapping emapper annotation")
    prop2type = {}
    metadata_dict = {}
    #parse input tree
    if args.tree:
        if args.tree_type == 'newick':
            tree = ete4_parse(args.tree)
        elif args.tree_type == 'ete':
            with open(args.tree, 'r') as f:
                file_content = f.read()
                tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
    
    # parse emapper annotation
    if args.emapper_annotations:
        metadata_dict, node_props, columns = parse_emapper_annotations(args.emapper_annotations)
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
    else: 
        # annotated_tree
        node_props=[]
        columns = {}
    
    if args.emapper_pfam:
        annot_tree_pfam_table(tree, args.emapper_pfam, args.alignment)
        pass
    
    if args.emapper_smart:
        annot_tree_smart_table(tree, args.emapper_smart, args.alignment)
        pass
    
    if args.alignment:
        #name2seq = parse_fasta(args.alignment)
        # for leaf in tree.iter_leaves():
        #     leaf.add_prop('seq', name2seq.get(leaf.name,''))
        pass

    popup_prop_keys = list(prop2type.keys())
    if args.taxonomic_profile:
        popup_prop_keys.extend([
            'rank',
            'sci_name',
            'taxid',
            'lineage',
            'named_lineage'
        ])
        
    # load metadata to leaf nodes
    if args.taxon_column: # to identify taxon column as taxa property from metadata
        taxon_column.append(args.taxon_column)
        annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type, taxon_column=args.taxon_column, taxon_delimiter=args.taxon_delimiter)
    else:
        annotated_tree = load_metadata_to_tree(tree, metadata_dict, prop2type=prop2type)
    
    # merge annotations depends on the column datatype
    start = time.time()
    if args.emapper_annotations:
        text_prop = ['seed_ortholog', 'max_annot_lvl', 'COG_category', 'EC', ]
        num_prop = ['evalue', 'score']
        multiple_text_prop = ['eggNOG_OGs', 'GOs', 'KEGG_ko', 'KEGG_Pathway', 
                            'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass',
                            'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'] # Pfams

        counter_stat = 'raw'
        num_stat = 'all'

        #pre load node2leaves to save time
        node2leaves = annotated_tree.get_cached_content()
        for node in annotated_tree.traverse("postorder"):
            internal_props = {}
            if node.is_leaf():
                pass
            else:
                if text_prop:
                    internal_props_text = merge_text_annotations(node2leaves[node], text_prop, counter_stat=counter_stat)
                    internal_props.update(internal_props_text)
                if multiple_text_prop:
                    internal_props_multi = merge_multitext_annotations(node2leaves[node], multiple_text_prop, counter_stat=counter_stat)
                    internal_props.update(internal_props_multi)

                if num_prop:
                    internal_props_num = merge_num_annotations(node2leaves[node], num_prop, num_stat=num_stat)                        
                    if internal_props_num:
                        internal_props.update(internal_props_num)

                for key,value in internal_props.items():
                    node.add_prop(key, value)

    end = time.time()
    print('Time for merge annotations to run: ', end - start)

    # taxa annotations
    if args.taxonomic_profile:
        if not args.taxadb:
            print('Please specify which taxa db using --taxadb <GTDB|NCBI>')
        else:
            if args.taxadb == 'GTDB':
                if args.taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr=args.taxon_column)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name")
            elif args.taxadb == 'NCBI':
                if args.taxon_column:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr=args.taxon_column, sp_delimiter=args.taxon_delimiter, sp_field=args.taxa_field)
                else:
                    annotated_tree, rank2values = annotate_taxa(annotated_tree, db=args.taxadb, taxid_attr="name", sp_delimiter=args.taxon_delimiter, sp_field=args.taxa_field)
    else:
        rank2values = {}

    # prune tree by rank
    if args.rank_limit:
        annotated_tree = taxatree_prune(annotated_tree, rank_limit=args.rank_limit)

    # prune tree by condition 
    if args.pruned_by: # need to be wrap with quotes
        condition_strings = args.pruned_by
        annotated_tree = conditional_prune(annotated_tree, condition_strings, prop2type)

    if args.outdir:
        base=os.path.splitext(os.path.basename(args.tree))[0]
        out_newick = base + '_annotated.nw'
        out_prop2tpye = base + '_prop2type.txt'
        out_ete = base+'_annotated.ete'
        out_tsv = base+'_annotated.tsv'

        ### out newick
        annotated_tree.write(outfile=os.path.join(args.outdir, out_newick), properties = [], format=1)
        ### output prop2type
        with open(os.path.join(args.outdir, base+'_prop2type.txt'), "w") as f:
            #f.write(first_line + "\n")
            for key, value in prop2type.items():
                f.write("{}\t{}\n".format(key, value))
        ### out ete
        with open(os.path.join(args.outdir, base+'_annotated.ete'), 'w') as f:
            f.write(b64pickle.dumps(annotated_tree, encoder='pickle', pack=False))
        
        ### out tsv
        tree2table(annotated_tree, internal_node=True, props=popup_prop_keys, outfile=os.path.join(args.outdir, out_tsv))
    return 


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
        if not n.is_leaf():
            random_seq_name = random.choice(n.get_leaf_names())
            random_node = post_tree.search_nodes(name=random_seq_name)[0]
            random_node_domains = random_node.props.get('dom_arq', 'none@none@none')
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
                print(info)
                seq_name = info[0]
                dom_name = info[1]
                dom_start = int(info[2])
                dom_end = int(info[3])

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
        if not n.is_leaf():
            random_seq_name = random.choice(n.get_leaf_names())
            random_node = post_tree.search_nodes(name=random_seq_name)[0]
            random_node_domains = random_node.props.get('dom_arq', 'none@none@none')
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

def goslim_annotation(gos_input, relative=True):
    """
    deprecated
    """
    output_dict = {}
    all_golsims_dict = {}
    goslim_script = os.path.join(os.path.dirname(__file__), 'goslim_list.R')
    output = subprocess.check_output([goslim_script,gos_input])
    #output_list = [line.split(' \t ') for line in output.decode('utf-8').split('\n') if line ]
    for line in output.decode('utf-8').split('\n'):
        if line:
            name, entries, desc, count = line.split(' \t ')
            if entries != '-':
                entries = entries.split(',')
                desc = desc.split('|') 
                count = np.array(count.split('|')).astype(int)
                if relative:
                    count = [float(i)/sum(count) for i in count]
            output_dict[name] = [entries, desc, count]
            for i in range(len(entries)):
                entry = entries[i]
                single_desc = desc[i]
                if entry not in all_golsims_dict:
                    all_golsims_dict[entry] = single_desc
    return output_dict, all_golsims_dict 

def tree2table(tree, internal_node=True, props=[], outfile='tree2table.csv'):
    node2leaves = {}
    leaf2annotations = {}
    if not props:
        leaf = tree.get_farthest_leaf()[0]
        props = list(leaf.props)
        
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

                if node.is_leaf():
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