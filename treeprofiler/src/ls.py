#!/usr/bin/env python3
try:
    from distutils.util import strtobool
except ImportError:
    from treeprofiler.src.utils import strtobool
    
from treeprofiler.src.utils import add_suffix

# Lineage specificity analysis
# Function to calculate precision, sensitivity, and F1 score
def calculate_metrics(node, total_with_trait, prop):
    if not node.is_leaf:
        clade_with_trait = sum(1 for child in node.leaves() if bool_checker(child, prop))
        clade_total = len([leave for leave in node.leaves()])
        precision = clade_with_trait / clade_total if clade_total else 0
        sensitivity = clade_with_trait / total_with_trait if total_with_trait else 0
        f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) else 0
        return precision, sensitivity, f1
    return 0, 0, 0

# Total number of nodes with the trait
def get_total_trait(tree, prop):
    return sum(1 for node in tree.leaves() if bool_checker(node, prop))

def bool_checker(node, prop):
    """
    Check if the property of a node can be interpreted as a boolean 'True'.

    :param node: The node whose property is to be checked.
    :param prop: The property name to check.
    :return: True if the property exists and is a boolean 'True', False otherwise.
    """
    prop_value = node.props.get(prop)
    if prop_value is not None:
        try:
            return bool(strtobool(str(prop_value)))
        except ValueError:
            return False
    return False

###### start lineage specificity analysis ######
def run_ls(tree, props, precision_cutoff=0.95, sensitivity_cutoff=0.95):
    best_node = None
    qualified_nodes = []
    best_f1 = -1
    for prop in props:
        total_with_trait = get_total_trait(tree, prop)
        # Calculating metrics for each clade
        for node in tree.traverse("postorder"):
            if not node.is_leaf:
                #node.add_prop(trait=int(node.name[-1]) if node.is_leaf else 0)
                precision, sensitivity, f1 = calculate_metrics(node, total_with_trait, prop)
                node.add_prop(add_suffix(prop, "prec"), precision)
                node.add_prop(add_suffix(prop, "sens"), sensitivity)
                node.add_prop(add_suffix(prop, "f1"), f1)
                #node.add_prop(precision=precision, sensitivity=sensitivity, f1_score=f1)
                #print(f"Node: {node.name} , Precision: {precision}, Sensitivity: {sensitivity}, F1 Score: {f1}")

                # Check if the node meets the lineage-specific criteria
                if not node.is_root:
                    if precision >= precision_cutoff and sensitivity >= sensitivity_cutoff:
                        node.add_prop(add_suffix(prop, "ls_clade"), True)
                        qualified_nodes.append(node)
                        if f1 > best_f1:
                            best_f1 = f1
                            best_node = node
        # if best_node:
        #     print(f"Root of Lineage-Specific Clade of {prop} Trait with f1 score {best_node.props.get(add_suffix(prop, 'f1'))}")
        #     best_node.add_prop(add_suffix(prop, "ls_clade"), True)

    return best_node, qualified_nodes

# #### find lineage-specific clades ####
# def find_lineage_specific_root(tree):
#     best_node = None
#     best_f1 = -1
#     for node in tree.traverse("postorder"):
#         if not node.is_leaf:
#             precision, sensitivity, f1 = calculate_metrics(node, total_with_trait)
#             node.add_props(precision=precision, sensitivity=sensitivity, f1_score=f1)
            
#             # Check if the node meets the lineage-specific criteria
#             if precision >= 0.5 and sensitivity >= 0.5 and f1 > best_f1:
#                 best_f1 = f1
#                 best_node = node
#     return best_node
