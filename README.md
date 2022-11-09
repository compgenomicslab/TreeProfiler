# MetaTreeProfiler 

```
metatreeprofiler.py <annotate|draw|plot> -t tree.nw -d metadata.csv
```

params:

    Input arguments
    -t <tree.nw> .nw file, customized tree input
    --taxa <kingdom|phylum|class|order|family|genus|species|subspecies> reference tree from taxonomic database
    -r <NCBI|GTDB> default [GTDB]  
    -d --data <metadata.csv> .csv, .tsv. mandatory input 
    --sep seperator of metadata columns. default [tab]
    --text_column col1,col2
    --num_column col1,col2
    --taxon_column col1

    Analysis arguments
    --rank_limit TAXONOMIC_LEVEL prune annotate tree by rank limit

    Plot arguments
    --Textlayouts col1,col2
    --RectangularLayouts col1,col2
    --HeatmapLayouts col1,col2
    --BarplotLayouts col1,col2

    
    Output arguments
    --interactive 
    -o --outfile output file <annotate.nw>
    --plot output picture