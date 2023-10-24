#!/bin/bash

## annotate tree 
treeprofiler annotate \
--tree ./basic_example1.nw \
--input-type newick \
--metadata ./basic_example1_metadata1.tsv \
--bool-prop bool_type bool_type2 \
--counter-stat relative \
-o ./

# select tree internal node where sample1_avg feature < 0.50
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--heatmap-layout sample1 \
--collapsed-by "sample1_avg < 0.50"

# collapse tree internal nodes, where `high` relative counter > 0.35 in random_type_counter property
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--rectangle-layout random_type \
--collapsed-by "random_type_counter:high > 0.35" \
--column-width 70