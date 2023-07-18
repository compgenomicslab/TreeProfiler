#!/bin/bash

# Annotate tree with GTDB taxonomic annotation
treeprofiler annotate \
--tree gtdb_example1.nw \
--input_type newick \
--metadata gtdb_example1.tsv \
--taxon_column name \
--taxonomic_profile \
--taxadb GTDB \
--outdir ./

# Visualize tree with colored taxonclade
treeprofiler plot \
--tree gtdb_example1_annotated.nw \
--input_type newick \
--taxonrectangle_layout \
--taxonclade_layout