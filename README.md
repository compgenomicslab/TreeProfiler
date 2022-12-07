# MetaTreeProfiler 

```
metatreeprofiler.py <annotate|draw> -t tree.nw -d metadata.csv
```

params:

    Input arguments
    -t TREE, --tree TREE  Input tree, .nw file, customized tree input
    --taxatree TAXATREE   <kingdom|phylum|class|order|family|genus|species|subspecies> reference tree from taxonomic database
    --taxadb TAXADB       <NCBI|GTDB> for taxonomic profiling or fetch taxatree default [GTDB]
    --taxonomic_profile   taxonomic profiling tree, default False

    -d METADATA, --metadata METADATA
                            <metadata.csv> .csv, .tsv. mandatory input
    --text_column TEXT_COLUMN
                            <col1,col2> names of columns which need to be read as categorical data
    --num_column NUM_COLUMN
                            <col1,col2> names of columns which need to be read as numerical data
    --bool_column BOOL_COLUMN
                            <col1,col2> names of columns which need to be read as boolean data
    --taxon_column TAXON_COLUMN
                            <col1> name of columns which need to be read as taxon data
    --taxon_delimiter TAXON_DELIMITER
                            delimiter of taxa columns. default [;]
    --text_column_idx       1,2,3 or 1-5 index of columns which need to be read as categorical data
    --num_column_idx        1,2,3 or 1-5 index columns which need to be read as numerical data
    --bool_column_idx       1,2,3 or 1-5 index columns which need to be read as boolean data

    Analysis arguments:

    --num_stat NUM_STAT   statistic calculation to perform for numerical data in internal nodes, [all, sum, avg, max, min, std] 
    --internal_plot_measure INTERNAL_PLOT_MEASURE
                            statistic measures to be shown in numerical layout for internal nodes, [default: avg]
    --counter_stat COUNTER_STAT
                            statistic calculation to perform for categorical data in internal nodes, raw count or in percentage [raw, relative] 
    --rank_limit RANK_LIMIT
                            TAXONOMIC_LEVEL prune annotate tree by rank limit
    --pruned_by PRUNED_BY
                            target tree pruned by customized conditions
    --collapsed_by COLLAPSED_BY
                            target tree collapsed by customized conditions
    --highlighted_by HIGHLIGHTED_BY
                            target tree highlighted by customized conditions
    
    syntax
    `--pruned_by|collapsed_by|highlighted_by 'prop>value,prop contains value, value in prop'
    --pruned_by|collapsed_by|highlighted_by 'random_type:low<0.5'` 

    Treelayout parameters

    --drawer DRAWER       Circular or Rectangular
    --collapse_level COLLAPSE_LEVEL
                            default collapse level, default is 10
    --ultrametric         ultrametric tree
    
    Plot arguments
    --BinaryLayout BINARYLAYOUT
                        <col1,col2> names of columns which need to be plot as BinaryLayout
    --RevBinaryLayout REVBINARYLAYOUT
                            <col1,col2> names of columns which need to be plot as RevBinaryLayout
    --ColorbranchLayout COLORBRANCHLAYOUT
                            <col1,col2> names of columns which need to be plot as Textlayouts
    --LabelLayout LABELLAYOUT
                            <col1,col2> names of columns which need to be plot as LabelLayout
    --RectangularLayout RECTANGULARLAYOUT
                            <col1,col2> names of columns which need to be plot as RectangularLayout
    --HeatmapLayout HEATMAPLAYOUT
                            <col1,col2> names of columns which need to be read as HeatmapLayout
    --BarplotLayout BARPLOTLAYOUT
                            <col1,col2> names of columns which need to be read as BarplotLayouts
    --TaxonLayout TAXONLAYOUT
                            <col1,col2> names of columns which need to be read as TaxonLayouts
    
    Output arguments:

    --interactive         run interactive session
    --plot PLOT           output as pdf
    -o OUTTREE, --outtree OUTTREE
                            output annotated tree
    --outtsv OUTTSV       output annotated tsv file