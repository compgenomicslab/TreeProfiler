#!/bin/bash

# Annotate tree with NCBI taxonomic annotation
treeprofiler annotate \
--tree ncbi_example.nw \
--input-type newick \
--metadata ncbi_example.tsv \
--taxonomic-profile \
--taxadb NCBI \
--taxon-column name \
--taxon-delimiter . \
--taxa-field 0  \
--outdir ./

# Visualize tree with colored taxonclade
treeprofiler plot \
--tree ncbi_example_annotated.nw \
--input-type newick \
--taxonclade-layout
