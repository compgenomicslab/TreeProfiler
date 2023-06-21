# annotate tree of progenome with habitat information and taxonomic annotation with NCBI database
treeprofiler.py annotate \
--tree examples/progenome3/progenome3.nw \
--metadata examples/progenome3/progenome3.tsv \
--taxonomic_profile \
--taxadb NCBI \
--taxon_delimiter . \
--taxa_field 0 \
--bool_prop aquatic_habitat,host_associated,soil_habitat \
--outdir examples/progenome3/

# visualize annotated trees with all features
treeprofiler.py plot \
--tree examples/progenome3/progenome3_annotated.ete \
--tree_type ete \
--barplot_layout GC,size \
--binary_layout aquatic_habitat,host_associated,soil_habitat \
--taxonclade_layout 

#if bool value is 1 or 0, treeprofiler will infer it as numerical data, hence we determine it as boolean value by using `--bool_prop` arguments