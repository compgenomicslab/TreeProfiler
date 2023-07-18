from __future__ import annotations
from ete4.parser.newick import NewickError
from ete4.smartview.renderer.gardening import remove
from ete4 import Tree, PhyloTree
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from itertools import chain
import operator
import math
import Bio
import re

# conditional syntax calling
operator_dict = {
                '<':operator.lt,
                '<=':operator.le,
                '=':operator.eq,
                '!=':operator.ne,
                '>':operator.gt,
                '>=':operator.ge,
                }

def check_nan(value):
    try:
        return math.isnan (float(value))
    except ValueError:
        return False

def counter_call(node, internal_prop, leaf_prop, datatype, operator_string, right_value):
    pair_delimiter = "--"
    item_seperator = "||"
    if datatype == str:
        counter_props = node.props.get(internal_prop)
            
        if counter_props:
            counter_datas = counter_props.split(item_seperator)
            for counter_data in counter_datas:
                k, v = counter_data.split(pair_delimiter)
                if k == leaf_prop:
                    left_value = float(v)
                    return operator_dict[operator_string](left_value, float(right_value))
                else:
                    pass
        else:    
            return False
        
    else:
        return False

def call(node, prop, datatype, operator_string, right_value):
    num_operators = [ '<', '<=', '>', '>=' ] 
    if datatype == str:
        if operator_string in num_operators:
            return False
        elif operator_string == 'contains':
            left_value = node.props.get(prop)
            if left_value:
                return right_value in left_value 
        elif operator_string == 'in':
            left_value = right_value
            right_value = node.props.get(prop)
            if right_value:
                return left_value in right_value 
        else:
            left_value = node.props.get(prop)
            
            if left_value:
                return operator_dict[operator_string](left_value, right_value)
    
    elif datatype == float:
        left_value = node.props.get(prop)
        if left_value:
            return operator_dict[operator_string](float(left_value), float(right_value))
        else:
            return False
    
    elif datatype == list:
        if operator_string in num_operators:
            return False
        elif operator_string == 'contains':
            left_value = node.props.get(prop)
            if left_value:
                return right_value in left_value 
        

def to_code(string):
    conditional_output = []
    operators = [ '<', '<=', '>', '>=', '=', '!=', 'in', 'contains'] 
    
    r = re.compile( '|'.join( '(?:{})'.format(re.escape(o)) for o in sorted(operators, reverse=True, key=len)) )

    #codes = string.split(',')
    # code = code.replace(",", " and ") # ',' means and
    # code = code.replace(";", " or ") # ';' means or
    # code = code.replace("=", " = ")
    # code = code.replace('>', " > ")
    # code = code.replace('>=', ' >= ')

    condition_strings = string.split(',')
    for condition_string in condition_strings:
        ops = r.findall(condition_string)
        for op in ops:
            condition_string = re.sub(op, ' '+op+' ', condition_string)
            left_value, op, right_value = condition_string.split(None,2)
            conditional_output.append([left_value, op, right_value])

    return conditional_output

SeqRecord = Bio.SeqRecord.SeqRecord
def get_consensus_seq(filename: Path | str, threshold=0.7) -> SeqRecord:
    #https://stackoverflow.com/questions/73702044/how-to-get-a-consensus-of-multiple-sequence-alignments-using-biopython
    common_alignment = MultipleSeqAlignment(
        chain(*AlignIO.parse(filename, "fasta"))
    )
    summary = SummaryInfo(common_alignment)
    consensus = summary.dumb_consensus(threshold, "-")
    return consensus


# parse ete4 Tree
def ete4_parse(newick, internal_parser="name"):
    if internal_parser == "name":
        tree = PhyloTree(newick, parser=1)
    elif internal_parser == "support":
        tree = PhyloTree(newick, parser=0)
        
    # Correct 0-dist trees
    has_dist = False
    for n in tree.traverse(): 
        if n.dist and float(n.dist) > 0: 
            has_dist = True
            break
    if not has_dist: 
        for n in tree.descendants(): 
            n.dist = 1
    return tree

# pruning
def taxatree_prune(tree, rank_limit='subspecies'):
    for node in tree.traverse("preorder"):
        if node.props.get('rank') == rank_limit:
            children = node.children.copy()
            for ch in children:
                print("prune", ch.name)
                remove(ch)
    return tree

def conditional_prune(tree, conditions_input, prop2type):
    conditional_output = []
    for line in conditions_input:
        single_one = to_code(line)
        
        conditional_output.append(single_one)

    ex = False
    while not ex:
        ex = True
        for n in tree.traverse():
            if not n.is_root:
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
                        n.detach()
                        ex = False
                    else:
                        pass
            # else:
            #     if n.dist == 0: 
            #         n.dist = 1
    return tree


    #array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
    array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
    return array

def children_prop_array(nodes, prop):
    #array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
    array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
    return array

def children_prop_array_missing(nodes, prop):
    """replace empty to missing value 'NaN' """ 
    array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
    #array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
    return array

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