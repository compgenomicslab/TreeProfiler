
# annotate gtdb tree with ar122 metadata lite, bac120 metadata lite and progenome3 habitat information
echo "Annotate GTDB tree with ar122 metadata lite, bac120 metadata lite and progenome3 habitat information..."
treeprofiler annotate \
--tree gtdbv202.nw \
--metadata ar122_metadata_r202_lite.tar.gz,bac120_metadata_r202_lite.tar.gz,progenome3.tar.gz \
--taxonomic_profile \
--taxadb GTDB \
-o .

# visualize annotated tree with selected properties
echo "Visualizing annotated GTDB tree with GTDB metadata, which are genome_size, protein_count, gc_percentage, ncbi_assembly_level, ncbi_genome_category"
echo "And progenome3 habitat information aquatic_habitat, host_associated, soil_habitat..."
treeprofiler plot \
--tree gtdbv202_annotated.ete \
--tree_type ete \
--barplot_layout genome_size,protein_count \
--heatmap_layout gc_percentage \
--binary_layout aquatic_habitat,host_associated,soil_habitat \
--rectangular_layout ncbi_assembly_level,ncbi_genome_category \
--taxonclade_layout \
--column_width 70