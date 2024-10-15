from __future__ import annotations
from treeprofiler.src import b64pickle
from ete4.parser.newick import NewickError
from ete4.core.operations import remove
from ete4 import Tree, PhyloTree
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from itertools import chain
from distutils.util import strtobool
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
from scipy import stats
import random
import colorsys
import operator
import math
import Bio
import re
import sys, os
from io import StringIO

# conditional syntax calling
operator_dict = {
                '<':operator.lt,
                '<=':operator.le,
                '=':operator.eq,
                '!=':operator.ne,
                '>':operator.gt,
                '>=':operator.ge,
                }

_true_set = {'yes', 'true', 't', 'y', '1'}
_false_set = {'no', 'false', 'f', 'n', '0'}

def str2bool(value, raise_exc=False):
    if isinstance(value, str) or sys.version_info[0] < 3 and isinstance(value, basestring):
        value = value.lower()
        if value in _true_set:
            return True
        if value in _false_set:
            return False

    if raise_exc:
        raise ValueError('Expected "%s"' % '", "'.join(_true_set | _false_set))
    return None

def str2bool_exc(value):
    return str2bool(value, raise_exc=True)

def check_float_array(array):
    """
    Checks if all elements in the array can be treated as floats
    
    :param array: The array to check.
    :return: True if all elements can be converted to floats, False otherwise.
    """
    try:
        # Attempt to convert the list to a NumPy array of floats
        np_array = np.array(array, dtype=np.float64)
        return not np.isnan(np_array).any()  # Check if the conversion resulted in any NaNs
    except ValueError:
        return False  # Conversion to float failed, indicating non-numerical data

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
    if datatype == str or datatype is None:
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
        

def to_code(condition_strings):
    conditional_output = []
    operators = [ '<', '<=', '>', '>=', '=', '!=', 'in', 'contains'] 
    
    r = re.compile( '|'.join( '(?:{})'.format(re.escape(o)) for o in sorted(operators, reverse=True, key=len)) )

    #codes = string.split(',')
    # code = code.replace(",", " and ") # ',' means and
    # code = code.replace(";", " or ") # ';' means or
    # code = code.replace("=", " = ")
    # code = code.replace('>', " > ")
    # code = code.replace('>=', ' >= ')

    #condition_strings = strings.split(',')
    for condition_string in condition_strings:
        ops = r.findall(condition_string)
        for op in ops:
            condition_string = re.sub(op, ' '+op+' ', condition_string)
            left_value, op, right_value = condition_string.split(None,2)
            conditional_output.append([left_value, op, right_value])

    return conditional_output

SeqRecord = Bio.SeqRecord.SeqRecord
def get_consensus_seq(matrix_string: Path | str, threshold=0.7) -> SeqRecord:
    #https://stackoverflow.com/questions/73702044/how-to-get-a-consensus-of-multiple-sequence-alignments-using-biopython
    common_alignment = MultipleSeqAlignment(
        chain(*AlignIO.parse(StringIO(matrix_string), "fasta"))
    )
    summary = SummaryInfo(common_alignment)
    consensus = summary.dumb_consensus(threshold, "-")
    return consensus

def counter2ratio(node, prop, minimum=0.05):
    counter_separator = '||'
    items_separator = '--'
    count_missing = True
    total = 0
    positive = 0

    counter_props = node.props.get(prop).split(counter_separator)
    for counter_prop in counter_props:
        k, v = counter_prop.split('--')
        if count_missing:
            if not check_nan(k):
                if strtobool(k):
                    positive = float(v)
            total += float(v) # here consider missing data in total
        else:
            if not check_nan(k):
                total += float(v) # here doesn't consider missing data in total
                if strtobool(k):
                    positive = float(v)
            
    total = int(total)
    if total != 0:
        ratio = positive / total
    else:
        ratio = 0
    
    if ratio < minimum and ratio != 0: # show minimum color for too low
        ratio = 0.05
    
    return ratio

def categorical2ratio(node, prop, all_values, minimum=0.05):
    counter_separator = '||'
    items_separator = '--'
    count_missing = True
    total = 0
    positive = 0
    ratios = []

    counter_props = node.props.get(prop).split(counter_separator)
    counter_dict = {k: v for k, v in [counter_prop.split('--') for counter_prop in counter_props]}
    total = sum([int(v) for v in counter_dict.values()])
    for value in all_values:
        if value in counter_dict:
            positive = int(counter_dict[value])
        else:
            positive = 0
        ratio = positive / total
        if ratio < minimum and ratio != 0: # show minimum color for too low
            ratio = 0.05
        ratios.append(ratio)
    
    return ratios

def dict_to_string(d, pair_seperator="--", item_seperator="||"):
    return item_seperator.join([f"{key}{pair_seperator}{value}" for key, value in d.items()])

def string_to_dict(s, pair_seperator="--", item_seperator="||"):
    return {item.split(pair_seperator)[0]: item.split(pair_seperator)[1] for item in s.split(item_seperator)}

def merge_dictionaries(dict1, dict2):
    for key, value in dict2.items():
        if key in dict1:
            dict1[key].update(value)
        else:
            dict1[key] = value
    return dict1
    
# validate tree format
class TreeFormatError(Exception):
    pass

def validate_tree(tree_path, input_type, internal_parser=None):
    tree = None  # Initialize tree to None
    eteformat_flag = False
    if input_type in ['ete', 'auto']:
        try:
            with open(tree_path, 'r') as f:
                file_content = f.read()
            tree = b64pickle.loads(file_content, encoder='pickle', unpack=False)
            eteformat_flag = True
        except Exception as e:
            if input_type == 'ete':
                raise TreeFormatError(f"Error loading tree in 'ete' format: {e}")

    if input_type in ['newick', 'auto'] and tree is None:
        #try:
        if tree_path == '-':
            tree = ete4_parse(sys.stdin, internal_parser=internal_parser)
        else:
            # checking file and output exists
            if not os.path.exists(tree_path):
                raise FileNotFoundError(f"Input tree {tree_path} does not exist.")
            
            tree = ete4_parse(open(tree_path), internal_parser=internal_parser)
        #except Exception as e:
        #    raise TreeFormatError(f"Error loading tree in 'newick' format: {e}\n"
        #                          "Please try using the correct parser with --internal-parser option, or check the newick format.")

    # if tree is None:
    #     raise TreeFormatError("Failed to load the tree in either 'ete' or 'newick' format.")

    return tree, eteformat_flag

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
        rank = node.props.get('rank')
        if rank == rank_limit:
            children = node.children.copy()
            for ch in children:
                print("prune", ch.name)
                remove(ch)
        lca_string = node.props.get('lca')
        if lca_string:
            lca_dict = string_to_dict(lca_string)
            if lca_dict:
                lca = lca_dict.get(rank_limit, None)
                if lca:
                    node.name = lca
                    children = node.children.copy()
                    for ch in children:
                        print("prune", ch.name)
                        remove(ch)
    return tree

def conditional_prune(tree, conditions_input, prop2type):
    conditional_output = []
    single_one = to_code(conditions_input)
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

# def _tree_prop_array(node, prop, leaf_only=False, numeric=False, list_type=False):
#     array = []
#     sep = '||'
#     if leaf_only:
#         for n in node.leaves():
#             prop_value = n.props.get(prop)
#             if prop_value is not None:
#                 if list_type:
#                     prop_value = prop_value.split(sep)
#                     if numeric:
#                         try:
#                             prop_value = [float(p) if p else np.nan for p in prop_value]
#                         except TypeError:
#                             raise TypeError(f"Cannot treat value '{prop_value}' as a number. Please check data type.")
#                     array.append(prop_value)
#                 else:
#                     # Check if the property value is a set
#                     if isinstance(prop_value, set):
#                         # Extract elements from the set
#                         array.extend(prop_value)
#                     else:
#                         if numeric:
#                             if prop_value == 'NaN':
#                                 array.append(np.nan)
#                             else:
#                                 try:
#                                     array.append(float(prop_value))
#                                 except TypeError:
#                                     raise TypeError(f"Cannot treat value '{prop_value}' as a number. Please check data type or use --numerical-matrix-layout")
#                         else:
#                             array.append(prop_value)
#     else:
#         for n in node.traverse():
#             prop_value = n.props.get(prop)
#             if prop_value is not None:
#                 if list_type:
#                     prop_value = prop_value.split(sep)
#                     if numeric:
#                         try:
#                             prop_value = [float(p) if p else np.nan for p in prop_value]
#                         except TypeError:
#                             raise TypeError(f"Cannot treat value '{prop_value}' as a number. Please check data type.")
#                     array.append(prop_value)
#                 else:
#                     # Check if the property value is a set
#                     if isinstance(prop_value, set):
#                         # Extract elements from the set
#                         array.extend(prop_value)
#                     else:
#                         if numeric:
#                             if prop_value == 'NaN':
#                                 array.append(np.nan)
#                             else:
#                                 try:
#                                     array.append(float(prop_value))
#                                 except TypeError:
#                                     raise TypeError(f"Cannot treat value '{prop_value}' as a number. Please check data type or use --numerical-matrix-layout")
#                         else:
#                             # Directly append the property value
#                             array.append(prop_value)
#     return array

def tree_prop_array(node, prop, leaf_only=False, numeric=False, list_type=False):
    array = []
    sep = '||'
    
    # Decide whether to traverse all nodes or only leaves
    nodes = node.leaves() if leaf_only else node.traverse()
    
    # Iterate over the selected nodes
    for n in nodes:
        prop_value = n.props.get(prop)
        if prop_value is not None:
            
            # Handle list-type property values
            if list_type:
                prop_value = prop_value.split(sep)
                if numeric:
                    try:
                        # Convert elements to floats, replace empty with NaN
                        prop_value = [float(p) if p else np.nan for p in prop_value]
                    except ValueError:
                        raise TypeError(f"Cannot treat value '{prop_value}' as a number. Please check data type.")
                array.append(prop_value)
            
            # Handle non-list property values
            else:
                # Handle sets (extend array with set elements)
                if isinstance(prop_value, set):
                    array.extend(prop_value)
                
                # Handle numeric values
                elif numeric:
                    if prop_value == 'NaN':
                        array.append(np.nan)
                    else:
                        try:
                            array.append(float(prop_value))
                        except ValueError:
                            raise TypeError(f"Cannot treat value '{prop_value}' as a number. Please check data type or use --numerical-matrix-layout")
                
                # Handle non-numeric, non-list values
                else:
                    array.append(prop_value)

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
        color_dict = assign_colors(values, cmap_name='terrain')
        
    return dict(sorted(color_dict.items()))

def rgba_to_hex(rgba):
    """Convert RGBA to Hexadecimal."""
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

def assign_colors(variables, cmap_name='tab20'):
    """Assigns colors to variables using a matplotlib colormap."""
    cmap = plt.cm.get_cmap(cmap_name, len(variables))  # Get the colormap
    colors = [rgba_to_hex(cmap(i)) for i in range(cmap.N)]  # Generate colors in hex format
    #random.shuffle(colors)
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

def build_custom_gradient(n_colors, min_color, max_color, mid_color=None):
    """
    Build a color gradient between two specified colors.

    Parameters:
    n_colors (int): Number of distinct colors to include in the gradient.
    min_color (str): Hex code or named color for the start of the gradient.
    max_color (str): Hex code or named color for the end of the gradient.

    Returns:
    dict: A dictionary mapping indices to colors in the generated gradient.
    """
    # Convert min and max colors to RGB
    min_rgb = mcolors.to_rgb(min_color)
    max_rgb = mcolors.to_rgb(max_color)
    mid_rgb = mcolors.to_rgb(mid_color) if mid_color else None
    
    color_gradient = {}
    # Determine if we're using a mid_color and split the range accordingly
    if mid_color:
        # Halfway point for the gradient transition
        mid_point = n_colors // 2
        
        for i in range(1, n_colors + 1):
            if i <= mid_point:
                # Transition from min_color to mid_color
                interpolated_rgb = [(mid_c - min_c) * (i - 1) / (mid_point - 1) + min_c for min_c, mid_c in zip(min_rgb, mid_rgb)]
            else:
                # Transition from mid_color to max_color
                interpolated_rgb = [(max_c - mid_c) * (i - mid_point - 1) / (n_colors - mid_point - 1) + mid_c for mid_c, max_c in zip(mid_rgb, max_rgb)]
            color_gradient[i] = mcolors.to_hex(interpolated_rgb)
    else:
        # If no mid_color, interpolate between min_color and max_color directly
        for i in range(1, n_colors + 1):
            interpolated_rgb = [(max_c - min_c) * (i - 1) / (n_colors - 1) + min_c for min_c, max_c in zip(min_rgb, max_rgb)]
            color_gradient[i] = mcolors.to_hex(interpolated_rgb)
    
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

def clear_specific_features(tree, features, leaf_only=False, internal_only=False):
    if leaf_only and not internal_only:
        for n in tree.leaves():
            for f in features:
                if f in n.props:
                    n.del_prop(f)
    elif internal_only and not leaf_only:
        for n in tree.traverse():
            if not n.is_leaf:
                for f in features:
                    if f in n.props:
                        n.del_prop(f)
    if leaf_only and internal_only:
        for n in tree.traverse():
            for f in features:
                if f in n.props:
                    n.del_prop(f)

def add_suffix(name, suffix, delimiter='_'):
    return str(name) + delimiter + str(suffix)

def normalize_values(values, normalization_method="min-max"):
    """
    Normalizes a list of numeric values using the specified method.
    
    Parameters:
    - values: List of elements to be normalized.
    - normalization_method: String indicating the normalization method. 
      Options are "min-max", "mean-norm", and "z-score".
      
    Returns:
    - A numpy array of normalized values.
    """
    
    def try_convert_to_float(values):
        """Attempts to convert values to float, raises error for non-convertible values."""
        converted = []
        for v in values:
            try:
                if v.lower() != 'nan':  # Assuming 'nan' is used to denote missing values
                    converted.append(float(v))
                else:
                    converted.append(np.nan)  # Convert 'nan' string to numpy NaN for consistency
            except ValueError:
                raise ValueError(f"Cannot treat value '{v}' as a number.")
        return np.array(converted)
    
    numeric_values = try_convert_to_float(values)
    valid_values = numeric_values[~np.isnan(numeric_values)]
    
    if normalization_method == "min-max":
        normalized = (valid_values - valid_values.min()) / (valid_values.max() - valid_values.min())
    elif normalization_method == "mean-norm":
        normalized = (valid_values - valid_values.mean()) / (valid_values.max() - valid_values.min())
    elif normalization_method == "z-score":
        normalized = stats.zscore(valid_values)
    else:
        raise ValueError("Unsupported normalization method.")
    
    return normalized

def find_bool_representations(column, rep=True):
    true_values = {'true', 't', 'yes', 'y', '1'}
    false_values = {'false', 'f', 'no', 'n', '0'}
    ignore_values = {'nan', 'none', ''}  # Add other representations of NaN as needed

    # Initialize sets to hold the representations of true and false values
    count = 0

    for value in column:
        str_val = str(value).strip().lower()  # Normalize the string value
        if str_val in ignore_values:
            continue  # Skip this value
        if rep:
            if str_val in true_values:
                count += 1
        else:
            if str_val in false_values:
               count += 1

    return count 


def color_gradient(c1, c2, mix=0):
    """ Fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1) """
    # https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def make_color_darker_log(hex_color, total, base=10):
    """Darkens the hex color based on a logarithmic scale of the total."""
    # Calculate darkening factor using a logarithmic scale
    darkening_factor = math.log(1 + total, base) / 50  # Adjust base and divisor as needed
    return make_color_darker(hex_color, darkening_factor)
    
def make_color_darker(hex_color, darkening_factor):
    """Darkens the hex color by a factor. Simplified version for illustration."""
    # Simple darkening logic for demonstration
    c = mcolors.hex2color(hex_color)  # Convert hex to RGB
    darker_c = [max(0, x - darkening_factor) for x in c]  # Darken color
    return mcolors.to_hex(darker_c)

def make_color_darker_scaled(hex_color, positive, maximum, base=10, scale_factor=10, min_darkness=0.6):
    """
    Darkens the hex color based on the positive count, maximum count, and a scaling factor.
    
    :param hex_color: The original color in hex format.
    :param positive: The current count.
    :param maximum: The maximum count achievable, corresponding to the darkest color.
    :param base: The base for the logarithmic calculation, affecting darkening speed.
    :param scale_factor: Factor indicating how much darker the color can get at the maximum count.
    :param min_darkness: The minimum darkness level allowed.
    :return: The darkened hex color.
    """
    if positive > maximum:
        raise ValueError("Positive count cannot exceed the maximum specified.")
    
    # Calculate the normalized position of 'positive' between 0 and 'maximum'
    normalized_position = positive / maximum if maximum != 0 else 0
    
    # Calculate the logarithmic scale position
    log_position = math.log(1 + normalized_position * (scale_factor - 1), base) / math.log(scale_factor, base)
    
    # Ensure the log_position respects the min_darkness threshold
    if log_position >= min_darkness:
        log_position = min_darkness

    # Convert hex to RGB
    rgb = mcolors.hex2color(hex_color)
    
    # Apply the darkening based on log_position
    darkened_rgb = [(1 - log_position) * channel for channel in rgb]
    
    return mcolors.to_hex(darkened_rgb)

# def transform_columns(columns, treat_as_whole=True, normalization_method="min-max"):
#     transformed = defaultdict(dict)
    
#     def normalize(values, method):
#         if method == "min-max":
#             return (values - values.min()) / (values.max() - values.min())
#         elif method == "mean-norm":
#             return (values - values.mean()) / (values.max() - values.min())
#         elif method == "z-score":
#             return stats.zscore(values)
#         else:
#             raise ValueError("Unsupported normalization method.")

#     def try_convert_to_float(values):
#         converted = []
#         for v in values:
#             try:
#                 if v != 'NaN':  # Assuming 'NaN' is used to denote missing values
#                     converted.append(float(v))
#                 else:
#                     converted.append(np.nan)  # Convert 'NaN' string to numpy NaN for consistency
#             except ValueError:
#                 raise ValueError(f"Cannot treat value '{v}' as a number.")
#         return np.array(converted)
    
#     if treat_as_whole:
#         # Attempt to concatenate all numeric columns into one array
#         all_numeric_values = np.concatenate([
#             try_convert_to_float(values) for values in columns.values()
#         ])
        
#         normalized_values = normalize(all_numeric_values[~np.isnan(all_numeric_values)], normalization_method)
        
#         start_idx = 0
#         for prop, values in columns.items():
#             num_values = len([v for v in values if v != 'NaN'])
#             transformed[prop][normalization_method] = normalized_values[start_idx:start_idx+num_values]
#             start_idx += num_values
#     else:
#         for prop, values in columns.items():
#             numeric_values = try_convert_to_float(values)
#             # Normalize only the non-NaN values
#             transformed[prop][normalization_method] = normalize(numeric_values[~np.isnan(numeric_values)], normalization_method)

#     return transformed
