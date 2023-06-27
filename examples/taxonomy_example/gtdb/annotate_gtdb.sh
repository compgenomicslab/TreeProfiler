#!/bin/bash

# Annotate tree with GTDB taxonomic annotation
treeprofiler annotate \
--tree gtdb_example1.nw \
--metadata gtdb_example1.tsv \
--taxon_column name \
--taxonomic_profile \
--taxadb GTDB \
--outdir ./

# Visualize tree with colored taxonclade
treeprofiler plot \
--tree gtdb_example1_annotated.nw \
--taxonrectangular_layout \
--taxonclade_layout