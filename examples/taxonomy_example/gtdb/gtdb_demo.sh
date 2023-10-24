#!/bin/bash

# Annotate tree with GTDB taxonomic annotation
treeprofiler annotate \
--tree gtdb_example1.nw \
--input-type newick \
--metadata gtdb_example1.tsv \
--taxon-column name \
--taxonomic-profile \
--taxadb GTDB \
--outdir ./

# Visualize tree with colored taxonclade
treeprofiler plot \
--tree gtdb_example1_annotated.nw \
--input-type newick \
--taxonrectangle-layout \
--taxonclade-layout