#!/bin/bash

# *Annotate multiple metadatas into tree, seperated by ,
echo "Annotate example tree with two metadata tables"
treeprofiler annotate \
--tree basic_example1.nw \
--input_type newick \
--metadata basic_example1_metadata1.tsv,basic_example1_metadata2.tsv \
--bool_prop bool_type \
-o ./

echo "Visualize properties categorical data random_type in rectangle_layout, numerical data sample1, sample2 in heatmap_layout and barplot_layout."
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--rectangle_layout random_type \
--binary_layout bool_type \
--heatmap_layout sample1,sample2,sample3 \
--barplot_layout sample4,sample5 \
--profiling_layout random_type \
--column_width 40 \
--padding_x 3

# **Visualize categorical properties 'random_type' with differen layouts
echo "Visualize annotated example tree by showing categorical property random_type with label_layout, rectangle_layout and colorbranch_layout."
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--rectangle_layout random_type \
--colorbranch_layout random_type \
--label_layout random_type \

# Convert ategorical properties random_type into presence-absence profiling matrix using --profiling_layout
echo "Visualize random_type into presence-absence profiling matrix"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--profiling_layout random_type


# Visualize numerical data 
echo "Visualize annotated example tree by showing numerical property sample[1-5] with heatmap_layout, and abs_data with barplot_layout"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--heatmap_layout sample1,sample2,sample3,sample4,sample5 \
--barplot_layout abs_data

# Visualize boolean data
echo "Visualize annotated example tree by showing boolean property bool_type with binary_layout and bool_type2 with revbinary_layout."
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--binary_layout bool_type \
--revbinary_layout bool_type2 \

# Visualize list data
echo "Visualize annotated example tree by showing presence-absence matrix of composition of each element of property list_data with multi_profiling_layout"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input_type ete \
--multi_profiling_layout list_data

#
#- *if boolean value is 1 or 0, treeprofiler will infer it as numerical data, hence we determine it as boolean value by using --bool_prop arguments
#
#- *here we use ete format because ete format can resume the original datatype of each property, especially when we need to Visualize such as list, boolean or numerical data.
#