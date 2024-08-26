# full gtdb annotation requires >6G disk space and >15G RAM memory
# download archaea and bacteria trees from gtdb
echo "Downloading  archaea and bacteria trees from GTDB database, version 202"
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_r202.tree
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_r202.tree

#merge two trees into one as gtdb tree
echo "Merging archaea and bacteria trees into one"
python merge_gtdbtree.py ar122_r202.tree bac120_r202.tree > gtdbv202.nw

# download metadata of archaea and bacteria
echo "Downloading metadata of archaea and bacteria from GTDB database, version 202"
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_metadata_r202.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz

# annotate gtdb tree with ar122 metadata, bac120 metadata and progenome3 habitat information
echo "Annotate GTDB tree with ar122 metadata, bac120 metadata and progenome3 habitat information..."
treeprofiler annotate \
--tree gtdbv202.nw \
--input-type newick \
--metadata ar122_metadata_r202.tar.gz bac120_metadata_r202.tar.gz progenome3.tsv \
--taxon-column name \
--taxadb GTDB \
--gtdb-version 202 \
-o .

# visualize annotated tree with selected properties
echo "Visualizing annotated GTDB tree with GTDB metadata, which are genome_size, protein_count, gc_percentage, ncbi_assembly_level, ncbi_genome_category"
echo "And progenome3 habitat information aquatic_habitat, host_associated, soil_habitat..."
treeprofiler plot \
--tree gtdbv202_annotated.ete \
--input-type ete \
--barplot-layout genome_size protein_count \
--heatmap-layout gc_percentage \
--binary-layout aquatic_habitat host_associated soil_habitat \
--rectangle-layout ncbi_assembly_level ncbi_genome_category \
--taxonclade-layout \
--column-width 70


