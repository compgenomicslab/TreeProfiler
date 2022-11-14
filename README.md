# MetaTreeProfiler 

```
metatreeprofiler.py <annotate|draw> -t tree.nw -d metadata.csv
```

params:

    Input arguments
    -t TREE, --tree TREE  Input tree, .nw file, customized tree input
    --taxatree TAXATREE   <kingdom|phylum|class|order|family|genus|species|subspecies> reference tree from taxonomic database
    --taxadb TAXADB       <NCBI|GTDB> for taxonomic profiling or fetch taxatree default [GTDB]
    --taxonomic_profile   taxonomic profiling tree

    -d METADATA, --metadata METADATA
                            <metadata.csv> .csv, .tsv. mandatory input
    --text_column TEXT_COLUMN
                            <col1,col2> names of columns which need to be read as categorical data
    --num_column NUM_COLUMN
                            <col1,col2> names of columns which need to be read as numerical data
    --taxon_column TAXON_COLUMN
                            <col1> name of columns which need to be read as taxon data
    --taxon_delimiter TAXON_DELIMITER
                            delimiter of taxa columns. default [;]


    Collapse arguments:
    --rank_limit RANK_LIMIT
                            TAXONOMIC_LEVEL prune annotate tree by rank limit
    --collapse condition
    for numerical
    --stats 
    --min_value ''
    --max_value

    for categorical
    --tag 


    Plot arguments

    --TextLayout TEXTLAYOUT
                            <col1,col2> names of columns which need to be plot as Textlayouts
    --LabelLayout LABELLAYOUT
                            <col1,col2> names of columns which need to be plot as LabelLayout
    --RectangularLayout RECTANGULARLAYOUT
                            <col1,col2> names of columns which need to be plot as RectangularLayout
    --HeatmapLayout HEATMAPLAYOUT
                            <col1,col2> names of columns which need to be read as HeatmapLayout
    --BarplotLayout BARPLOTLAYOUT
                            <col1,col2> names of columns which need to be read as BarplotLayouts
    --TaxonLayout 
                            names of columns which need to be read as TaxonLayout
    
    Output arguments:

    --interactive         run interactive session
    --plot PLOT           output as pdf
    -o OUTTREE, --outtree OUTTREE
                            output annotated tree
    --outtsv OUTTSV       output annotated tsv file