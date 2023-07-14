# annotate tree of progenome with habitat information and taxonomic annotation with NCBI database
echo "Annotate GTDB tree with progenome3 habitat information"
treeprofiler annotate \
--tree progenome3.nw \
--input_type newick \
--metadata progenome3.tsv \
--taxonomic_profile \
--taxadb NCBI \
--taxon_delimiter . \
--taxa_field 0 \
--bool_prop aquatic_habitat,host_associated,soil_habitat \
--outdir ./

# visualize annotated trees with all features
echo "Visualizing annotated GTDB tree with progenome3 habitat information aquatic_habitat, host_associated, soil_habitat"
treeprofiler plot \
--tree progenome3_annotated.ete \
--input_type ete \
--barplot_layout GC,size \
--binary_layout aquatic_habitat,host_associated,soil_habitat \
--taxonclade_layout 

#if bool value is 1 or 0, treeprofiler will infer it as numerical data, hence we determine it as boolean value by using `--bool_prop` arguments