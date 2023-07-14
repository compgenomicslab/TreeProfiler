#!/bin/bash

# annotate tree of nifH with eggNOG-mapper annotations and taxonomic annotations
echo "Start mapping tree with annotation metadata"
treeprofiler annotate \
--tree  nifH.nw \
--input_type newick \
--emapper_annotation nifH.out.emapper.annotations  \
--emapper_pfam nifH.out.emapper.pfam \
--alignment nifH.faa.aln \
--taxonomic_profile \
--taxadb NCBI \
--taxon_delimiter . \
--taxa_field 0 \
-o ./


# visualize annotated trees with all eggnog mapper features
echo "Visualizing annotated tree with all eggnog mapper features......"
treeprofiler plot \
--tree nifH_annotated.ete \
--input_type ete \
--emapper_layout \
--domain_layout \
--taxonclade_layout \
--column_width 70