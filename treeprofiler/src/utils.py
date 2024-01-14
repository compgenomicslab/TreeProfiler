from __future__ import annotations
from ete4.parser.newick import NewickError
from ete4.core.operations import remove
from ete4 import Tree, PhyloTree
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from itertools import chain
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import random
import colorsys
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
def get_internal_parser(internal_parser="name"):
    if internal_parser == "name":
        return 1
    elif internal_parser == "support":
        return 0
         
def ete4_parse(newick, internal_parser="name"):
    tree = PhyloTree(newick, parser=get_internal_parser(internal_parser))
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
    array = []
    for n in nodes:
        prop_value = n.props.get(prop)
        if prop_value is not None:
            # Check if the property value is a set
            if isinstance(prop_value, set):
                # Extract elements from the set
                array.extend(prop_value)
            else:
                # Directly append the property value
                array.append(prop_value)
    return array

def children_prop_array_missing(nodes, prop):
    """replace empty to missing value 'NaN' """ 
    array = [n.props.get(prop) if n.props.get(prop) else 'NaN' for n in nodes] 
    #array = [n.props.get(prop) for n in nodes if n.props.get(prop) ] 
    return array

def convert_to_int_or_float(column):
    """
    Convert a column to integer if possible, otherwise convert to float64.
    
    Args:
    column (list): The input column data.
    
    Returns:
    np.ndarray: Array converted to integer or float64.
    """
    np_column = np.array(column)

    # Try converting to integer
    try:
        return np_column.astype(np.int64)
    except ValueError:
        # If conversion to integer fails, convert to float64
        return np_column.astype(np.float64)
        
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

def random_color(h=None, l=None, s=None, num=None, sep=None, seed=None):
    """Return the RGB code of a random color.

    Hue (h), Lightness (l) and Saturation (s) of the generated color
    can be specified as arguments.
    """
    def rgb2hex(rgb):
        return '#%02x%02x%02x' % rgb

    def hls2hex(h, l, s):
        return rgb2hex( tuple([int(x*255) for x in colorsys.hls_to_rgb(h, l, s)]))

    if not h:
        if seed:
            random.seed(seed)
        color = 1.0 / random.randint(1, 360)
    else:
        color = h

    if not num:
        n = 1
        sep = 1
    else:
        n = num

    if not sep:
        n = num
        sep = (1.0/n)

    evenly_separated_colors =  [color + (sep*n) for n in range(n)]

    rcolors = []
    for h in evenly_separated_colors:
        if not s:
            s = 0.5
        if not l:
            l = 0.5
        rcolors.append(hls2hex(h, l, s))

    if num:
        return rcolors
    else:
        return rcolors[0]

def assign_color_to_values(values, paired_colors):
    """Assigns colors to values, either from a predefined list or generates new ones."""
    color_dict = {}

    if len(values) <= len(paired_colors):
        # Use predefined colors if enough are available
        color_dict = {val: paired_colors[i] for i, val in enumerate(values)}
    else:
        # Use the assign_colors function to generate colors if not enough predefined colors
        color_dict = assign_colors(values)

    return dict(sorted(color_dict.items()))

def rgba_to_hex(rgba):
    """Convert RGBA to Hexadecimal."""
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

def assign_colors(variables, cmap_name='tab20'):
    """Assigns colors to variables using a matplotlib colormap."""
    cmap = plt.cm.get_cmap(cmap_name, len(variables))  # Get the colormap
    colors = [rgba_to_hex(cmap(i)) for i in range(cmap.N)]  # Generate colors in hex format
    return dict(zip(variables, colors))


def build_color_gradient(n_colors, colormap_name="viridis"):
    """
    Build a color gradient based on the specified matplotlib colormap.

    Parameters:
    n_colors (int): Number of distinct colors to include in the gradient.
    colormap_name (str): Name of the matplotlib colormap to use. "viridis"  # Replace with "plasma", "inferno", "magma", etc., as needed

    Returns:
    dict: A dictionary mapping indices to colors in the specified colormap.
    """
    cmap = plt.get_cmap(colormap_name)
    indices = np.linspace(0, 1, n_colors)
    color_gradient = {i: mcolors.rgb2hex(cmap(idx)) for i, idx in enumerate(indices, 1)}
    return color_gradient

def clear_extra_features(forest, features):
    features = set(features) | {'name', 'dist', 'support'}
    for tree in forest:
        for n in tree.traverse():
            for f in set(n.props) - features:
                if f not in features:
                    n.del_prop(f)
            
            for key, value in n.props.items():
                # Check if the value is a set
                if isinstance(value, set):
                    # Convert the set to a string representation
                    # You can customize the string conversion as needed
                    n.props[key] = ','.join(map(str, value))

def add_suffix(name, suffix, delimiter='_'):
    return str(name) + delimiter + str(suffix)