# download archaea and bacteria trees from gtdb
# wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_r202.tree
# wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_r202.tree

#merge two trees into one as gtdb tree
# python merge_gtdbtree.py ar122_r202.tree bac120_r202.tree > gtdbv202.nw

# download metadata of archaea and bacteria
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_metadata_r202.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz

# annotate gtdb tree with ar122 metadata, bac120 metadata and progenome3 habitat information
treeprofiler.py annotate \
--tree gtdbv202.nw \
--metadata ar122_metadata_r202.tar.gz,bac120_metadata_r202.tar.gz,progenome3.tsv \
--bool_prop aquatic_habitat,host_associated,soil_habitat \
--taxonomic_profile \
--taxadb GTDB \
-o .

# visualize annotated tree with selected properties
treeprofiler.py plot \
--tree gtdbv202_annotated.ete \
--tree_type ete \
--barplot_layout genome_size,protein_count \
--heatmap_layout gc_percentage \
--binary_layout aquatic_habitat,host_associated,soil_habitat \
--rectangular_layout ncbi_assembly_level,ncbi_genome_category \
--taxonclade_layout


