#!/bin/bash

# annotate tree of COG1348 with eggNOG-mapper annotations and taxonomic annotations
echo "Start mapping tree with annotation metadata"
treeprofiler.py annotate \
--tree  COG1348.tree \
--emapper_annotation COG1348.out.emapper.annotations  \
--emapper_pfam COG1348.out.emapper.pfam \
--alignment COG1348.faa.aln \
--taxonomic_profile \
--taxadb NCBI \
--taxon_delimiter . \
--taxa_field 0 \
-o ./


# visualize annotated trees with all eggnog mapper features
echo "Visualizing annotated tree with all eggnog mapper features......"
treeprofiler.py plot \
--tree COG1348_annotated.ete \
--tree_type ete \
--emapper_layout \
--domain_layout \
--taxonclade_layout \
