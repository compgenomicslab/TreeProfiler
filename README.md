# MetaTreeProfiler Tutorial

## Overview
MetaTreeProfiler is command-line tool for profiling metadata table into phylogenetic tree with descriptive analysis and output visualization

## Installation & Basic Usage
MetaTreeProfiler requires to install ete4 toolkit
```
# install ete4 dependencies Cython
conda install -c anaconda cython
pip install Flask
pip install Flask-RESTful Flask-HTTPAuth Flask-Compress Flask-Cors
conda install six numpy scipy

# install ete4
git clone https://github.com/etetoolkit/ete.git
cd ete/
git branch checkout ete4
pip install -e .
```

### Input files
MetaTreeProfiler takes following file types as input 

| Input    |      Filetype  | 
|----------|-------------   |
| Tree     |      Newick    | 
| Metadata |      TSV       |

### basic usage

Choose your input tree and corresponding metadata, MetaTreeProfiler will map all the metadata into related tree node, and visualize it on the local server browser in interactive interface.

```
metatreeprofiler.py --tree tree.nw --metadata metadata.tsv --interactive
```

### Output arguments

#### Tree visuliaze interface 
add flag `--interactive`
```
metatreeprofiler.py --tree tree.nw --metadata metadata.tsv --interactive
```

#### Output as annotated tree 
```
metatreeprofiler.py --tree tree.nw --metadata metadata.tsv --outtree annotated_tree.nw
```

#### Output as tsv file
```
metatreeprofiler.py --tree tree.nw --metadata metadata.tsv --outtsv annotated_tree.tsv
```



# Using MetaTreeProfiler
In this Tutorial we will use MetaTreeProfiler and demostrate basic usage with data in examples/


```
cd examples/
ls 
basic_example1.nw  basic_example1.tsv

```

## Maping metadata into tree
### Simple map metadata into corresponding leaf nodes 
treeprofiler will start local server which metadata will be mapped to corresponding leaf nodes
```
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --interactive
```

It's recommended to save the annotated tree as outtree for re-running MetaTreeProfiler much faster, using `--outtree` or `-o`

```
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --outtree examples/annotated_basic_example1.nw
```
re-run annotated tree by adding flag `--annotated_tree`
```
python treeprofiler.py --tree examples/annotated_basic_example1.nw --metadata examples/basic_example1.tsv --annotated_tree --interactive
```

### Map metadata into tree and profile tree internal nodes annotations and analysis
At the above example, we only mapped metadata to leaf nodes, in this example, we will also profile **internal nodes** annotation and analysis of their children nodes.


### Mapping Categorical data
For categorical dataset, each internal node will count the selected feature of its children nodes as counter, as shown as ```<feature_name>_counter``` in internal node. To label categorical feature metadata, using following arguments

```
# label categorical data by column name(s) in metadata (for multiple columns, seperate by ","), using --text_column <header>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column random_type --interactive

# label categorical data by column index, using --text_column_idx <idx>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column_idx 12 --interactive  

# label column index by range, "[star_idx-end_idx]"
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column_idx [1-5] --interactive 
```

Categorical data will be process as counter in each internal node. Users can choose either counter is raw or relative count by using `--counter_stat`
```
# raw count, example internal_node shown as: ```random_type_counter: medium--3||high--2```
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column seed_ortholog,random_type --counter_stat raw --interactive

# relative count, example internal_node shown as: ```random_type_counter: medium--0.60||high--0.40```
 internal_node example shown as, random_type_counter: medium--3||high--2
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column seed_ortholog,random_type --counter_stat relative --interactive
```
### Mapping Boolean data
For Boolean dataset, each internal node will count the selected feature(s) of its children nodes as counter as categorical, as shown as `_counter` of suffix feature name(s) of internal node. To label Boolean feature metadata, using following arguments


```
# label boolean data by column name(s) in metadata (for multiple columns, seperate by ","), using --bool_column <header>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column bool_type --interactive

# label boolean data by column index, using --bool_column_idx <idx>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column_idx 13 --interactive  

# label column index by range, "[star_idx-end_idx]"
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column_idx [1-5] --interactive 
```

Boolean counter stats follows rule as categorical data

### Mapping Numerical data
For numerical dataset, each internal node will perform folwing descriptive statistic analysis of all of its children node of selected feature(s)

| internal_node properties  |      statistic method  | 
|----------|-------------   |
| `<feature name>`_avg      |      average    | 
| `<feature name>`_sum      |      sum    | 
| `<feature name>`_max      |      maximum    | 
| `<feature name>`_min      |      minimum    | 
| `<feature name>`_std      |      standard deviation    | 

To label numerical data
```
# label numerical data by column name(s) in metadata (for multiple columns, seperate by ","), using --num_column <header>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column evalue,score --interactive

# label numerical data by column index, using --num_column_idx <idx>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx 2,3 --interactive  

# label column index by range, "[star_idx-end_idx]"
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx [2-3] --interactive 
```

By default, numerical feature will be calculated all the descriptive statistic, but users can choose specific one to be calculated by using `--num_stat [all, sum, avg, max, min, std] `

--num_stat NUM_STAT   statistic calculation to perform for numerical data in internal nodes, [all, sum, avg, max, min, std] 
```
# by default
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column evalue,score --num_stat all --interactive

# only average calculation
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column evalue,score --num_stat avg --interactive
```

### Mapping Taxon data

## Visualizing annotated tree with layouts

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
                            
## Conditional query in annotated tree

### Query Syntax
### pruned_by
### collapsed_by
### highlighted_by

params:

    input parameters:
        Input parameters

    -t TREE, --tree TREE  Input tree, .nw file, customized tree input
    -d METADATA, --metadata METADATA
                            <metadata.csv> .csv, .tsv. mandatory input
    --annotated_tree      inputtree already annot                 ated by treeprofileer
    --no_colnames         metadata table doesn't contain columns name
    --text_column TEXT_COLUMN
                            <col1,col2> names of columns which need to be read as categorical data
    --num_column NUM_COLUMN
                            <col1,col2> names of columns         which need to be read as numerical data
    --bool_column BOOL_COLUMN
                            <col1,col2> names of columns which need to be read as boolean data
    --text_column_idx TEXT_COLUMN_IDX
                            1,2,3 or 1-5 index of columns which need to be read as categorical data
    --num_column_idx NUM_COLUMN_IDX
                            1,2,3 or 1-5 index columns which need to be read as numerical data
    --bool_column_idx BOOL_COLUMN_IDX
                            1,2,3 or 1-5 index columns which need to be read as boolean data
    --taxatree TAXATREE   <kingdom|phylum|class|order|family|genus|species|subspecies> reference tree from taxonomic database
    --taxadb TAXADB       <NCBI|GTDB> for taxonomic profiling or fetch taxatree default [GTDB]
    --taxon_column TAXON_COLUMN
                            <col1> name of columns which need to be read as taxon data
    --taxon_delimiter TAXON_DELIMITER
                            delimiter of taxa columns. default [;]
    --taxonomic_profile   Determine if you need taxonomic profile on tree

    Analysis arguments:
    Analysis parameters

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

    basic treelayout arguments:
    treelayout parameters

    --drawer DRAWER       Circular or Rectangular
    --collapse_level COLLAPSE_LEVEL
                            default collapse level, default is 10
    --ultrametric         ultrametric tree

    Plot arguments:
    Plot parameters

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
    Output parameters

    --interactive         run interactive session
    --port PORT           run interactive session
    --plot PLOT           output as pdf
    -o OUTTREE, --outtree OUTTREE
                            output annotated tree
    --outtsv OUTTSV       output annotated tsv file
