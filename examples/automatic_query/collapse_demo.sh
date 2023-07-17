#!/bin/bash

## annotate tree 
treeprofiler annotate \
--tree ./basic_example1.nw \
--input_type newick \
--metadata ./basic_example1_metadata1.tsv \
--bool_prop bool_type,bool_type2 \
--counter_stat relative \
-o ./

# select tree internal node where sample1_avg feature < 0.50
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--heatmap_layout sample1 \
--collapsed_by "sample1_avg < 0.50"

# collapse tree internal nodes, where `high` relative counter > 0.35 in random_type_counter property
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--rectangle_layout random_type \
--collapsed_by "random_type_counter:high > 0.35" \
--column_width 70