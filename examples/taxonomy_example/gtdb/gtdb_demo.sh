#!/bin/bash

# Annotate tree with GTDB taxonomic annotation
treeprofiler annotate \
--tree gtdb_example1.nw \
--input-type newick \
--metadata gtdb_example1.tsv \
--taxon-column name \
--taxadb GTDB \
--taxa-dump gtdb202dump.tar.gz \
--outdir ./

# Visualize tree with colored taxonclade
treeprofiler plot \
--tree gtdb_example1_annotated.nw \
--input-type newick \
--taxonclade-layout

# Visualize tree with colored rectangle
treeprofiler plot \
--tree gtdb_example1_annotated.nw \
--input-type newick \
--taxonrectangle-layout

# Visualize tree with colored rectangle with automatic collapse
treeprofiler plot \
--tree gtdb_example1_annotated.ete \
--input-type ete \
--taxoncollapse-layout
