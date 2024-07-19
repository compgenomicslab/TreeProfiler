#!/bin/bash

# *Annotate multiple metadatas into tree, seperated by ,
echo "Annotate example tree with two metadata tables"
treeprofiler annotate \
--tree basic_example1.nw \
--input-type newick \
--metadata basic_example1_metadata1.tsv basic_example1_metadata2.tsv \
-o ./

echo "Visualize properties categorical data random_type in rectangle_layout, numerical data sample1, sample2 in heatmap_layout and barplot_layout."
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--rectangle-layout random_type \
--binary-layout bool_type \
--heatmap-layout sample1 sample2 sample3 \
--barplot-layout sample4 sample5 \
--profiling-layout random_type \
--column-width 40 \
--padding-x 3

# **Visualize categorical properties 'random_type' with differen layouts
echo "Visualize annotated example tree by showing categorical property random_type with label_layout, rectangle_layout and colorbranch_layout."
treeprofiler plot \
--tree basic_example1_annotated.ete \
--rectangle-layout random_type 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--colorbranch-layout random_type 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--label-layout random_type 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--background-layout random_type 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--piechart-layout random_type

# Convert ategorical properties random_type into presence-absence profiling matrix using --profiling_layout
echo "Visualize random_type into presence-absence profiling matrix"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--profiling-layout random_type

# Visualize list data
echo "Visualize annotated example tree by showing presence-absence matrix of composition of each element of property list_data with multi_profiling_layout"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--profiling-layout list_data

# Visualize binary data
echo "Visualize annotated example tree by showing boolean property bool_type with binary_layout and bool_type2 with revbinary_layout."
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--binary-layout bool_type bool_type2

treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--binary-aggregate-layout bool_type bool_type2

treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--binary-unicolor-layout bool_type bool_type2

treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--binary-unicolor-aggregate-layout bool_type bool_type2

treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--binary-matrix-layout bool_type bool_type2

# Visualize numerical data 
echo "Visualize annotated example tree by showing numerical property sample[1-5] with heatmap_layout, and abs_data with barplot_layout"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--colorbranch-layout dist sample1

treeprofiler plot \
--tree basic_example1_annotated.ete \
--heatmap-layout sample1 sample2 sample3 sample4 sample5 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--heatmap-mean-layout sample1 sample2 sample3 sample4 sample5 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--heatmap-zscore-layout sample1 sample2 sample3 sample4 sample5 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--bubble-layout sample1 sample2 sample3 sample4 sample5 

# Visualize numerical data with barplot
echo "Visualize annotated example tree by showing numerical property abs_data[1-2] with barplot_layout"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--barplot-layout abs_data abs_data2 \
--barplot-width 250 \
--barplot-scale abs_data

# Visualize numerical data with different heatmap normalization methods
echo "Visualize annotated example tree by showing numerical property sample[1-5] with heatmap_layout"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--heatmap-layout sample1 sample2 sample3 sample4 sample5 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--heatmap-mean-layout sample1 sample2 sample3 sample4 sample5 

treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--heatmap-zscore-layout sample1 sample2 sample3 sample4 sample5 

# Visualize numerical data with numerial matrix layout
echo "Visualize annotated example tree by showing numerical property sample[1-5] with numerical matrix"
treeprofiler plot \
--tree basic_example1_annotated.ete \
--input-type ete \
--numerical-matrix-layout sample1 sample2 sample3 sample4 sample5 

