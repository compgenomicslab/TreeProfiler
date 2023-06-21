#!/bin/bash

# annotate diauxic tree with diauxic numerical matrix metadata
echo "annotate diauxic tree with diauxic numerical matrix metadata"
treeprofiler.py annotate \
--tree diauxic.nw \
--metadata examples/basic_example2/diauxic.array \
--outdir ./

# visualize annotated diauxic tree by showing numerical data from col1-col7 with numerical matrix
echo "visualize annotated diauxic tree by showing numerical data from col1-col7 with numerical matrix"
treeprofiler.py plot \
--tree diauxic_annotated.ete \
--tree_type ete\
--numerical_matrix_layout col1,col2,col3,col4,col5,col6,col7


# annotate MCC_FluA_H3 tree with MCC_FluA_H3_Genotype metadata and MCC_FluA_H3 multiple sequence alignment
echo "annotate MCC_FluA_H3 tree with MCC_FluA_H3_Genotype metadata and MCC_FluA_H3 multiple sequence alignment"
treeprofiler.py annotate \
--tree MCC_FluA_H3.nw \
--metadata MCC_FluA_H3_Genotype.txt \
--alignment FluA_H3_AA.fas \
--outdir ./

# Visualize annotated MCC_FluA_H3 tree all feature with retangular block in aligned panel using --categorical_matrix_layout
echo "visualize annotated MCC_FluA_H3 tree all feature with retangular block"
treeprofiler.py plot \
--tree MCC_FluA_H3_annotated.nw \
--categorical_matrix_layout PB2,PB1,PA,HA,NP,NA,M,NS

# visualize MCC_FluA_H3 MSA 
echo "visualize MCC_FluA_H3 tree with MSA "
treeprofiler.py plot \
--tree MCC_FluA_H3_annotated.nw \
--alignment_layout