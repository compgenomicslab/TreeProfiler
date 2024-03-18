# annotate tree of progenome with habitat information and taxonomic annotation with NCBI database
echo "Annotate GTDB tree with progenome3 habitat information"
treeprofiler annotate \
--tree progenome3.nw \
--input-type newick \
--metadata progenome3.tsv \
--taxon-column name \
--taxadb NCBI \
--taxon-delimiter . \
--taxa-field 0 \
--outdir ./

# visualize annotated trees with all features
echo "Visualizing annotated GTDB tree with progenome3 habitat information aquatic_habitat, host_associated, soil_habitat"
treeprofiler plot \
--tree progenome3_annotated.ete \
--input-type ete \
--barplot-layout GC size \
--binary-layout aquatic_habitat host_associated soil_habitat \
--taxonclade-layout 

