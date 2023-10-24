#!/bin/bash

## annotate tree 
treeprofiler annotate \
--tree ./basic_example1.nw \
--input-type newick \
--metadata ./basic_example1_metadata1.tsv \
--bool-prop bool_type bool_type2 \
--counter-stat relative \
-o ./

# Conditional pruning, prune leaf node whose name contain "FALPE"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--pruned-by "name contains FALPE"