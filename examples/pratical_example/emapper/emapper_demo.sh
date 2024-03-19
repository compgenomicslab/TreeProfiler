#!/bin/bash

# annotate tree of nifH with eggNOG-mapper annotations and taxonomic annotations
echo "Start mapping tree with annotation metadata"
treeprofiler annotate \
--tree  nifH.nw \
--input-type newick \
--emapper-annotation nifH.out.emapper.annotations  \
--emapper-pfam nifH.out.emapper.pfam \
--alignment nifH.faa.aln \
--taxon-column name \
--taxadb NCBI \
--taxon-delimiter '.' \
--taxa-field 0 \
-o ./


# visualize annotated trees with all eggnog mapper features
echo "Visualizing annotated tree with all eggnog mapper features......"
treeprofiler plot \
--tree nifH_annotated.ete \
--input-type ete \
--emapper-layout \
--domain-layout \
--taxonclade-layout \
--column-width 70
