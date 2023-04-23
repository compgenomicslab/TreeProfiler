#!/usr/bin/env python


from ete4.parser.newick import NewickError
from ete4 import Tree, PhyloTree
from ete4.coretype.seqgroup import SeqGroup
from ete4 import GTDBTaxa
from ete4 import NCBITaxa
from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview.renderer.gardening import remove
#from ete4.smartview  import SeqFace, SeqMotifFace, AlignmentFace
from layouts import (text_layouts, taxon_layouts, staple_layouts, 
                    conditional_layouts, seq_layouts, profile_layouts)
from tree_plot import get_image
from utils import get_consensus_seq

from argparse import ArgumentParser
import argparse
from collections import defaultdict
from collections import Counter
from itertools import islice
from scipy import stats
from io import StringIO
import colorsys
import random
import b64pickle
import itertools
import math
import numpy as np
import csv
import sys
import time
import os

__author__ = 'Ziqi DENG'
__license__ = "GPL v2"
__email__ = 'dengziqi1234@gmail.com'
__version__ = '0.0.1'
__date__ = '01-11-2022'
__description__ = ('A program for profiling metadata on target '
                    'tree and conduct summary analysis')


#colors_50 = ["#E41A1C","#C72A35","#AB3A4E","#8F4A68","#735B81","#566B9B","#3A7BB4","#3A85A8","#3D8D96","#419584","#449D72","#48A460","#4CAD4E","#56A354","#629363","#6E8371","#7A7380","#87638F","#93539D","#A25392","#B35A77","#C4625D","#D46A42","#E57227","#F67A0D","#FF8904","#FF9E0C","#FFB314","#FFC81D","#FFDD25","#FFF12D","#F9F432","#EBD930","#DCBD2E","#CDA12C","#BF862B","#B06A29","#A9572E","#B65E46","#C3655F","#D06C78","#DE7390","#EB7AA9","#F581BE","#E585B8","#D689B1","#C78DAB","#B791A5","#A8959F","#999999"]
paried_color = ["red", "darkblue", "lightgreen", "sienna", "lightCoral", "violet", "mediumturquoise",   "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal", "darkyellow"]

### annotate tree ####
def tree_annotate(args):
    total_color_dict = []
    layouts = []
    level = 2 # level 1 is the leaf name
    prop2type = {}
    metadata_dict = {}

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

    

    #code goes here
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
    
    ### decide popup keys
    # if args.annotated_tree:
    if args.tree_type == 'ete':
        leafa, _, leafb, _ = tree._get_farthest_and_closest_leaves()
        leaf_prop2type = get_prop2type(leafa)
        leaf_prop2type.update(get_prop2type(leafb))
        
        internal_node_prop2type = get_prop2type(tree)
        prop2type.update(leaf_prop2type)
        prop2type.update(internal_node_prop2type)
        
        # exisiting props in internal node
        existing_internal_props = list(tree.props.keys())
        # exisiting props in leaf node
        existing_leaf_props = list(leaf_prop2type.keys()) 
        popup_prop_keys = list(set(existing_internal_props+existing_leaf_props))
    elif args.tree_type == 'newick':
        leafa, _, leafb, _ = tree._get_farthest_and_closest_leaves()
        # props which add in the arguments
        required_internal_props = list(prop2type.keys()) 
        # exisiting prop in leaf node
        existing_leaf_props = list(leafa.props.keys()) + list(leafb.props.keys())
        popup_prop_keys = list(set(required_internal_props + existing_leaf_props))
    # else:
    #     # all the metadata to the leaves, no internal
    #     popup_prop_keys = list(prop2type.keys())
    popup_prop_keys = list(prop2type.keys())
    if args.taxonomic_profile:
        popup_prop_keys.extend([
            'rank',
            'sci_name',
            'taxid',
            'lineage',
            'named_lineage'
        ])

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
    if args.counter_stat:
        counter_stat = args.counter_stat

    if args.num_stat:
        num_stat = args.num_stat

    # merge annotations depends on the column datatype
    start = time.time()

    if not args.annotated_tree:
        #pre load node2leaves to save time
        node2leaves = annotated_tree.get_cached_content()
        count = 0
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
        tree2table(annotated_tree, internal_node=True, props=popup_prop_keys, outfile=os.path.join(args.outdir, out_tsv))
                
    # if args.outtsv:
    #     tree2table(annotated_tree, internal_node=True, outfile=args.outtsv)
    
    return 

def ete4_parse(newick):
    #parse tree via ete4
    try:
        tree = PhyloTree(newick)
    except NewickError:
        try:
            tree = PhyloTree(newick, format=1)            
        except NewickError:
            tree = PhyloTree(newick, format=1, quoted_node_names=True)

    # Correct 0-dist trees
    has_dist = False
    for n in tree.traverse(): 
        if float(n.dist) > 0: 
            has_dist = True
            break
    if not has_dist: 
        for n in tree.iter_descendants(): 
            n.dist = 1

    return tree

def check_missing(input_string):
    import re

    pattern = r'^(?:\W+|none|None|null|NaN|)$'

    if re.match(pattern, input_string):
        #print("Input contains only non-alphanumeric characters, 'none', a missing value, or an empty value.")
        return True
    else:
        return False
    
# 
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

def flatten(nasted_list):
    """
    input: nasted_list - this contain any number of nested lists.
    ------------------------
    output: list_of_lists - one list contain all the items.
    """

    list_of_lists = []
    for item in nasted_list:
        if type(item) == list:
            list_of_lists.extend(item)
        else:
            list_of_lists.extend(nasted_list)
    return list_of_lists

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

def children_prop_array(nodes, prop):
    #array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
    array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
    return array

def children_prop_array_missing(nodes, prop):
    array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
    #array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
    return array

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

import numbers

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

# def parse_pfam_annotations(input_file, delimiter='\t', no_colnames=False):
#     return

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

import subprocess
def goslim_annotation(gos_input, relative=True):
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
                

### visualize tree
def tree_plot(args):
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
                popup_prop_keys.append(prop)
                popup_prop_keys.append(prop+"_"+internal_num_rep)

            #heatmap_layouts, level, _ = staple_layouts.LayoutHeatmap('Heatmap_'+args.heatmap_layout, level, internal_num_rep, args.heatmap_layout)
            layouts.extend(heatmap_layouts)

        if layout == 'colorbranch_layout':
            colorbranch_layouts, level, color_dict = get_layouts(args.colorbranch_layout, 'colorbranch',  
                                                                level, 'counter', prop2type=prop2type, 
                                                                column_width=args.column_width)
            layouts.extend(colorbranch_layouts)
            total_color_dict.append(color_dict)

        if layout == 'label_layout':
            label_layouts, level, color_dict = get_layouts(args.label_layout, 'label', level, 'counter', prop2type=prop2type, column_width=args.column_width)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)

        if layout == 'rectangular_layout':
            rectangular_layouts, level, color_dict = get_layouts(args.rectangular_layout, 'rectangular', level, 'counter', prop2type=prop2type, column_width=args.column_width)
            layouts.extend(rectangular_layouts)
            total_color_dict.append(color_dict)

        if layout == 'binary_layout':
            label_layouts, level, color_dict = get_layouts(args.binary_layout, 'binary', level, 'counter', column_width=args.column_width)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)

        if layout == 'revbinary_layout':
            label_layouts, level, color_dict = get_layouts(args.revbinary_layout, 'revbinary', level, 'counter', column_width=args.column_width)
            layouts.extend(label_layouts)
            total_color_dict.append(color_dict)
        
        if layout == 'barplot_layout':
            barplot_layouts, level,color_dict = get_layouts(args.barplot_layout, 'barplot', level, internal_num_rep, prop2type=prop2type, column_width=args.barplot_width)
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
        label_layouts, level, _ = get_layouts(text_props, 'rectangular', level, 'counter', prop2type=prop2type)
        layouts.extend(label_layouts)
        
        num_props = [
            #'evalue',
            'score'
        ]
        barplot_layouts, level, _ = get_layouts(num_props, 'barplot', level, internal_num_rep, prop2type=prop2type)
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
            # if multiple_text_prop == 'GOs':
            #     pair_seperator = "--"
            #     item_seperator = "||"
            #     target_prop = 'GOslims'
            #     gos_input = os.path.join(os.path.dirname(__file__) + 'gos_input.tsv')
            #     node2leaves = tree.get_cached_content()

            #     # run goslim_list.r to retrieve goslims
            #     with open(gos_input, 'w') as f:
            #         for leaf in tree.iter_leaves():
            #             if leaf.props.get(multiple_text_prop):
            #                 go_prop = ','.join(leaf.props.get(multiple_text_prop))
            #                 line = leaf.name + "\t" + go_prop + "\n"
            #                 f.write(line)
                
            #     output, all_golsims = goslim_annotation(gos_input)

            #     # load to leaves
            #     for leaf in tree.iter_leaves():
            #         leaf_goslim = output.get(leaf.name,'')
            #         if leaf_goslim:
            #             leaf.add_prop(target_prop, leaf_goslim[0])

            #     # sum to parent nodes
            #     for node in tree.traverse("postorder"):
            #         if node.is_leaf():
            #             pass
            #         else:
            #             prop_list = children_prop_array(node2leaves[node], target_prop)
            #             multi_prop_list = []
            #             for elements in prop_list:
            #                 for j in elements:
            #                     multi_prop_list.append(j)
            #             node.add_prop(add_suffix(target_prop, 'counter'), item_seperator.join([add_suffix(str(key), value, pair_seperator) for key, value in sorted(dict(Counter(multi_prop_list)).items())]))

            #     # ouput to layouts 
            #     for entry, desc in all_golsims.items():
            #         if entry != '-':
            #             golayout = profile_layouts.LayoutGOslim(name=f'GOslims:{desc}({entry})', column=level,
            #                                 go_propfile=[entry, desc], goslim_prop=target_prop, padding_x=2, 
            #                                 padding_y=2, legend=True)
            #             level+=1
            #             layouts.append(golayout)
            #     popup_prop_keys.append('GOslims')
            #     popup_prop_keys.append(add_suffix(target_prop, 'counter'))

            # else:
            #     matrix, all_values = multiple2profile(tree, multiple_text_prop)
            #     multiple_text_prop_layout = profile_layouts.LayoutProfile(name=multiple_text_prop, 
            #     alignment=matrix, profiles=all_values, column=level)
            #     level += 1
            #     layouts.append(multiple_text_prop_layout)
            
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
    
    # if args.interactive:
    #     tree.explore(tree_name='example',layouts=layouts, port=args.port, popup_prop_keys=sorted(popup_prop_keys))
    # elif args.plot:
    #     plot(tree, layouts, args.port, args.plot)
    if args.plot:
        get_image(tree, layouts, args.port, os.path.abspath(args.plot))
    else:
        tree.explore(tree_name='example',layouts=layouts, port=args.port, include_props=sorted(popup_prop_keys))
    
    return

import re
def taxatree_prune(tree, rank_limit='subspecies'):
    ranks = ['domain','superkingdom','kingdom','subkingdom','infrakingdom','superphylum','phylum','division','subphylum','subdivision','infradivision','superclass','class','subclass','infraclass','subterclass','parvclass','megacohort','supercohort','cohort','subcohort','infracohort','superorder','order','suborder','infraorder','parvorder','superfamily','family','subfamily','supertribe','tribe','subtribe','genus','subgenus','section','subsection','species group','series','species subgroup','species','infraspecies','subspecies','forma specialis','variety','varietas','subvariety','race','stirp','form','forma','morph','subform','biotype','isolate','pathogroup','serogroup','serotype','strain','aberration']
    no_ranks = ['clade','unspecified','no rank','unranked','Unknown']
    
    # rank_limit = rank_limit.lower()
    
    # ex = False
    # while not ex:
    #     ex = True
    #     for n in tree.traverse('preorder'):
    #         if not n.is_root():
    #             rank_prop = n.props.get('rank')
    #             if rank_prop in ranks: 
    #                 rank_idx = ranks.index(rank_prop)
    #                 limit_rank_idx = ranks.index(rank_limit)
    #                 if rank_idx >= limit_rank_idx:
    #                     for child in n.get_children():
    #                         child.detach()
    #                         ex = False
                    
    # ex = False
    # while not ex:
    #     ex = True
    #     for n in tree.iter_leaves():
    #         if n.props.get('rank') != rank_limit:
    #             n.detach()
    #             ex = False
    for node in tree.traverse("preorder"):
        if node.props.get('rank') == rank_limit:
            children = node.children.copy()
            for ch in children:
                print("prune", ch.name)
                remove(ch)
    return tree

from utils import to_code, call, counter_call
def conditional_prune(tree, conditions_input, prop2type):
    conditional_output = []
    for line in conditions_input:
        single_one = to_code(line)
        conditional_output.append(single_one)

    ex = False
    while not ex:
        ex = True
        for n in tree.traverse():
            if not n.is_root():
                final_call = False
                for or_condition in conditional_output:
                    for condition in or_condition:
                        op = condition[1]
                        if op == 'in':
                            value = condition[0]
                            prop = condition[2]
                            datatype = prop2type.get(prop)
                            final_call = call(n, prop, datatype, op, value)
                        elif ":" in condition[0]:
                            internal_prop, leaf_prop = condition[0].split(':')
                            value = condition[2]
                            datatype = prop2type[internal_prop]
                            final_call = counter_call(n, internal_prop, leaf_prop, datatype, op, value)
                        else:
                            prop = condition[0]
                            value = condition[2]
                            prop = condition[0]
                            value = condition[2]
                            datatype = prop2type.get(prop)
                            final_call = call(n, prop, datatype, op, value)
                        if final_call == False:
                            break
                        else:
                            continue
                    if final_call:
                        #print('cut', n.name)
                        n.detach()
                        ex = False
                    else:
                        pass
            else:
                if n.dist == 0: 
                    n.dist = 1
    return tree

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

def get_layouts(argv_inputs, layout_name, level, internal_rep, prop2type=None, column_width=70): 
    props = []
    layouts = []
    prop_color_dict = {} # key = value, value = color id
    # identify range [1-5], index 1,2,3 and column names
    for i in argv_inputs:
        if i[0] == '[' and i[-1] == ']':
            column_start, column_end = get_range(i)
            for j in range(column_start, column_end+1):
                props.append(node_props[j-1])
        else:
            try:
                i = int(i)
                props.append(node_props[i-1])
            except ValueError:
                props.append(i)

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

def populate_emapper_annotate_args(emapper_annotate_args_p):
    group = emapper_annotate_args_p.add_argument_group(title='emapper annotate parameters',
        description="Input parameters of emapper annotate ")
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
    
    group.add_argument('--taxonomic_profile',
        default=False,
        action='store_true',
        required=False,
        help="Determine if you need taxonomic annotation on tree")
    
    group = emapper_annotate_args_p.add_argument_group(title='OUTPUT options',
        description="")
    group.add_argument('-o', '--outdir',
        type=str,
        required=False,
        help="output annotated tree")

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
    _main(sys.argv)

def _main(arguments):
    # global prop2type
    # global text_prop, num_prop, bool_prop
    # global annotated_tree, node_props, columns

    # CREATE REUSABLE PARSER OPTIONS
    # main args
    # args = read_args() ## old

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
    annotate_args_p.set_defaults(func=tree_annotate)

    ## - EMAPPER ANNOTATE - 
    # emapper_annotate_args_p = subparser.add_parser('emapper-annotate', parents=[main_args_p],
    #                                         description='annotate tree with eggnog mapper annotation data')
    # populate_emapper_annotate_args(emapper_annotate_args_p)
    # emapper_annotate_args_p.set_defaults(func=tree_emapper_annotate)
    
    ## - PLOT - 
    plot_args_p = subparser.add_parser('plot', parents=[main_args_p],
                                            description='annotate plot')
    poplulate_plot_args(plot_args_p)
    plot_args_p.set_defaults(func=tree_plot)

    ## - RUN -
    args = parser.parse_args(arguments[1:])
    args.func(args)
    
if __name__ == '__main__':
    main()