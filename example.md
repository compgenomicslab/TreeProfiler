# claudia's tree with binary layout
python treeprofiler.py --tree ../claudia_tree/concatenate_pg3_faa_ft.nw -d /home/deng/Projects/metatree_drawer/claudia_tree/annotated_concatenate_pg3_faa_ft.tsv --taxon_column GCA --taxonomic_profile --num_prop GC,size --bool_prop aquatic_habitat,host_associated,soil_habitat --BinaryLayout aquatic_habitat,host_associated,soil_habitat

# highlighted by 
python treeprofiler.py -t demo/p__Thermoproteota.nw -d /home/deng/Projects/metatree_drawer/metatreedrawer/demo/metadata_p__Thermoproteota_relative_random.txt --taxonomic_profile --num_prop sample1,sample2,sample3,sample4,sample5 --highlighted_by 'name!=GB_GCA_000494185.1,sample1>=0.5' --port 5002

python treeprofiler.py --tree ../claudia_tree/concatenate_pg3_faa_ft.nw --metadata ../claudia_tree/annotated_concatenate_pg3_faa_ft.tsv --num_prop GC,size --bool_prop aquatic_habitat,host_associated,soil_habitat --counter_stat relative --BinaryLayout aquatic_habitat,host_associated,soil_habitat --highlighted_by 'soil_habitat_counter: 1>0.5' --highlighted_by 'GC_avg>55'   --port 5001 

#Taxon 
python treeprofiler.py -t demo/p__Thermoproteota.nw -d /home/deng/Projects/metatree_drawer/metatreedrawer/demo/metadata_p__Thermoproteota_relative_random.txt --taxonomic_profile --num_prop sample1,sample2,sample3,sample4,sample5 --bool_prop bool_type,bool_type2 --counter_stat relative --TaxonLayout name --port 5003 --HeatmapLayout [1-5] --collapsed_by 'bool_type2_counter:True>0.5'


python treeprofiler.py -t demo/p__Thermoproteota.nw -d /home/deng/Projects/metatree_drawer/metatreedrawer/demo/metadata_p__Thermoproteota_relative_random.txt --taxonomic_profile --num_prop sample1,sample2,sample3,sample4,sample5 --text_prop random_type --bool_prop bool_type,bool_type2 --TaxonLayout name --BinaryLayout bool_type,bool_type2 --BarplotLayout sample1,sample2

python treeprofiler.py --tree examples/gtdb_example1.nw --metadata examples/gtdb_example1.tsv --num_prop sample1 --taxonomic_profile --rank_limit family --interactive --BarplotLayout sample1

#NCBI
python treeprofiler.py --tree examples/spongilla_example.nw --taxonomic_profile --annotated_tree --taxadb NCBI --taxon_delimiter . --taxa_field 0 --TaxonLayout --interactive


# progenome

python treeprofiler.py --tree examples/progenome3.nw --metadata examples/progenome3.tsv --taxon_column GCF --taxonomic_profile --num_prop GC,size --bool_prop aquatic_habitat,host_associated,soil_habitat --BarplotLayout GC,size --TaxonLayout --BinaryLayout aquatic_habitat,host_associated,soil_habitat --interactive --outtree examples/progenome3_annotated.nw

python treeprofiler.py --tree examples/progenome3_annotated.nw --metadata examples/progenome3.tsv --BarplotLayout GC,size --TaxonLayout --BinaryLayout aquatic_habitat,host_associated,soil_habitat --interactive

##ncbi
python treeprofiler.py --tree examples/progenome3.nw --metadata examples/progenome3.tsv --taxonomic_profile --taxadb NCBI --taxon_delimiter . --taxa_field 0 --num_prop GC,size --bool_prop aquatic_habitat,host_associated,soil_habitat --BarplotLayout GC,size --TaxonLayout --BinaryLayout aquatic_habitat,host_associated,soil_habitat --interactive

##bug 
python treeprofiler.py --tree examples/progenome3.nw --metadata examples/progenome3.tsv --taxonomic_profile --taxadb NCBI --taxon_delimiter . --taxa_field 0 --num_prop GC,size --bool_prop aquatic_habitat,host_associated,soil_habitat --BarplotLayout GC,size --TaxonLayout --BinaryLayout aquatic_habitat,host_associated,soil_habitat --interactive --rank_limit phylum


#
python treeprofiler.py --tree examples/progenome3_annotated.nw --annotated_tree --num_prop GC,size --bool_prop aquatic_habitat,host_associated,soil_habitat --interactive --BarplotLayout size --BinaryLayout aquatic_habitat,host_associated,soil_habitat


########## multi ############
# emapper_annotation 
treeprofiler.py annotate --tree examples/emapper/7955/7955.ENSDARP00000116736.fasta.final_tree.nw --emapper_annotations examples/emapper/7955/out.emapper.annotations --taxonomic_profile --taxadb NCBI --taxa_field 0 --taxon_delimiter . -o examples/emapper/7955/emapper_annotations/

treeprofiler.py plot --tree examples/emapper/7955/emapper_annotations/7955.ENSDARP00000116736.fasta.final_tree_annotated.ete --tree_type ete --emapper_layout

# emapper_pfam
treeprofiler.py annotate 
--tree examples/emapper/7955/7955.ENSDARP00000116736.fasta.final_tree.nw 
--emapper_pfam examples/emapper/7955/out.emapper.pfam 
--seq examples/emapper/7955/7955.ENSDARP00000116736.fasta.final_tree.fa 
--outdir examples/emapper/7955/emapper_pfams/

treeprofiler.py plot --tree examples/emapper/7955/emapper_pfams/7955.ENSDARP00000116736.fasta.final_tree_annotated.ete --tree_type ete --domain_layout

# alignment
treeprofiler.py plot --tree examples/emapper/7955/7955.ENSDARP00000116736.fasta.final_tree.nw  --alignment_layout examples/emapper/7955/7955.ENSDARP00000116736.fasta.final_tree.fa

# multi metadata input, auto detect list/str/num/bool
treeprofiler.py annotate --tree examples/emapper/7955/7955.ENSDARP00000116736.fasta.final_tree.nw --metadata examples/emapper/7955/out.emapper.annotations.clean,examples/emapper/7955/emapper_multi/GO.out --taxonomic_profile --taxadb NCBI --taxa_field 0 --taxon_delimiter . -o examples/emapper/7955/emapper_multi/

treeprofiler.py plot --tree examples/emapper/7955/emapper_multi/7955.ENSDARP00000116736.fasta.final_tree_annotated.ete --tree_type ete --label_layout ontology,eggNOG_OGs --rectangular_layout COG_category --barplot_layout score

treeprofiler.py plot --tree examples/emapper/7955/emapper_multi/7955.ENSDARP00000116736.fasta.final_tree_annotated.ete --tree_type ete --profiling_layout eggNOG_OGs

