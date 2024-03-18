# start ACR computing for Country trait and Lineage-Specificity analysis for is_Albania and is_Greece traits
treeprofiler annotate \
--tree Albanian.tree.152tax.nwk \
--internal-parser name \
-d metadata_tab_ls.csv \
--acr-discrete-columns Country \
--delta-stats \
--ls-columns is_Albania is_Greece \
--prec-cutoff 0.7 \
--sens-cutoff 0.7 \
--threads 6 \
-o ./

# Visualize ACR annotation
treeprofiler plot \
-t Albanian.tree.152tax_annotated.ete \
--acr-discrete-layout Country

# Visualize Lineage-Specificity annotation
treeprofiler plot -t Albanian.tree.152tax_annotated.ete \
--ls-layout is_Albania is_Greece