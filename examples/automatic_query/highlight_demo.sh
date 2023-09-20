#!/bin/bash

## annotate tree 
treeprofiler annotate \
--tree ./basic_example1.nw \
--input-type newick \
--metadata ./basic_example1_metadata1.tsv \
--bool-prop bool_type bool_type2 \
--counter-stat relative \
-o ./


# select tree node whose name contains `FALPE` character
treeprofiler plot \
--tree ./basic_example1_annotated.nw \
--input-type newick \
--highlighted-by "name contains FALPE"


# select tree node whose sample1 feature > 0.50, here we using ete format which can resume the datatype 
treeprofiler plot \
--tree ./basic_example1_annotated.ete \
--input-type ete \
--highlighted-by "sample1 > 0.50" \
--heatmap-layout sample1

# if use tree in newick format, we need to attach the prop2type file which can resume the datatype
treeprofiler plot \
--tree basic_example1_annotated.nw \
--input-type newick \
--prop2type basic_example1_prop2type.txt \
--highlighted-by "sample1 > 0.50" \
--heatmap-layout sample1


# select tree  node where sample1 feature > 0.50 AND sample2 < 0.2
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--heatmap-layout sample1 sample2 sample3 sample4 sample5 \
--highlighted-by "sample1>0.50,sample2<0.2" 

# select tree node where sample1 feature > 0.50 OR sample2 < 0.2
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--heatmap-layout sample1 sample2 sample3 sample4 sample5 \
--highlighted-by "sample1>0.50" \
--highlighted-by "sample2<0.2" 