#!/bin/bash

# *Annotate multiple metadatas into tree, seperated by ,
echo "Annotate example tree with two metadata tables"
treeprofiler.py annotate \
--tree basic_example1.nw \
--metadata basic_example1.tsv,basic_example1_metadata2.tsv \
--bool_prop bool_type \
-o ./

# **Visualize categorical properties 'random_type' with differen layouts
echo "visualize annotated example tree by showing categorical property random_type with label_layout, rectangular_layout and colorbranch_layout."
treeprofiler.py plot \
--tree basic_example1_annotated.ete \
--tree_type ete \
--rectangular_layout random_type \
--colorbranch_layout random_type \
--label_layout random_type \

# Convert ategorical properties random_type into presence-absence profiling matrix using --profiling_layout
echo "visualize random_type into presence-absence profiling matrix"
treeprofiler.py plot \
--tree basic_example1_annotated.ete \
--tree_type ete \
--profiling_layout random_type


# Visualize numerical data 
echo "visualize annotated example tree by showing numerical property sample[1-5] with heatmap_layout, and abs_data with barplot_layout"
treeprofiler.py plot \
--tree basic_example1_annotated.ete \
--tree_type ete \
--heatmap_layout sample1,sample2,sample3,sample4,sample5 \
--barplot_layout abs_data

# Visualize boolean data
echo "visualize annotated example tree by showing boolean property bool_type with binary_layout and bool_type2 with revbinary_layout."
treeprofiler.py plot \
--tree basic_example1_annotated.ete \
--tree_type ete \
--binary_layout bool_type \
--revbinary_layout bool_type2 \

# Visualize list data
echo "visualize annotated example tree by showing presence-absence matrix of composition of each element of property list_data with multi_profiling_layout"
treeprofiler.py plot \
--tree basic_example1_annotated.ete \
--tree_type ete \
--multi_profiling_layout list_data

#
#- *if boolean value is 1 or 0, treeprofiler will infer it as numerical data, hence we determine it as boolean value by using --bool_prop arguments
#
#- *here we use ete format because ete format can resume the original datatype of each property, especially when we need to visualize such as list, boolean or numerical data.
#