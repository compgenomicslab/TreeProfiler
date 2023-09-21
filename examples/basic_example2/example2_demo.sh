#!/bin/bash

# annotate diauxic tree with diauxic numerical matrix metadata
echo "annotate diauxic tree with diauxic numerical matrix metadata"
treeprofiler annotate \
--tree diauxic.nw \
--input-type newick \
--metadata diauxic.array \
--outdir ./

# visualize annotated diauxic tree by showing numerical data from col1-col7 with numerical matrix
echo "visualize annotated diauxic tree by showing numerical data from col1-col7 with numerical matrix"
treeprofiler plot \
--tree diauxic_annotated.ete \
--input-type ete \
--numerical-matrix-layout col1 col2 col3 col4 col5 col6 col7


# annotate MCC_FluA_H3 tree with MCC_FluA_H3_Genotype metadata and MCC_FluA_H3 multiple sequence alignment
echo "annotate MCC_FluA_H3 tree with MCC_FluA_H3_Genotype metadata and MCC_FluA_H3 multiple sequence alignment"
treeprofiler annotate \
--tree MCC_FluA_H3.nw \
--input-type newick \
--internal-parser support \
--metadata MCC_FluA_H3_Genotype.txt \
--alignment FluA_H3_AA.fas \
--outdir ./

# Visualize annotated MCC_FluA_H3 tree all feature with retangular block in aligned panel using --categorical_matrix_layout
echo "visualize annotated MCC_FluA_H3 tree all feature with retangular block"
treeprofiler plot \
--tree MCC_FluA_H3_annotated.nw \
--input-type newick \
--categorical-matrix-layout PB2 PB1 PA HA NP NA M NS

# visualize MCC_FluA_H3 MSA 
echo "visualize MCC_FluA_H3 tree with MSA "
treeprofiler plot \
--tree MCC_FluA_H3_annotated.nw \
--input-type newick \
--alignment-layout
