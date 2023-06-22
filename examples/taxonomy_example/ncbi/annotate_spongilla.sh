#!/bin/bash

# Annotate tree with NCBI taxonomic annotation
treeprofiler.py annotate \
--tree spongilla_example.nw \
--metadata spongilla_example.tsv \
--taxonomic_profile \
--taxadb NCBI \
--taxon_column name \
--taxon_delimiter . \
--taxa_field 0  \
--outdir ./

# Visualize tree with colored taxonclade
treeprofiler.py plot \
--tree spongilla_example_annotated.nw \
--taxonclade_layout