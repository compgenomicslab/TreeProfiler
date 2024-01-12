#!/usr/bin/env python3
from ete4 import Tree

# Lineage specificity analysis

# Function to calculate precision, sensitivity, and F1 score
def calculate_metrics(node, total_with_trait):
    if not node.is_leaf:
        clade_with_trait = sum(1 for child in node.leaves() if int(child.props.get("trait")) == 1)
        clade_total = len([leave for leave in node.leaves()])
        precision = clade_with_trait / clade_total if clade_total else 0
        sensitivity = clade_with_trait / total_with_trait if total_with_trait else 0
        f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) else 0
        return precision, sensitivity, f1
    return 0, 0, 0

