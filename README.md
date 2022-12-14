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

Install MetaTreeProfiler
```
git clone https://github.com/dengzq1234/MetaTreeDrawer
cd MetaTreeDrawer/
# add treeprofiler to path
export PATH=$PATH:$(pwd)
```

### Input files
MetaTreeProfiler takes following file types as input 

| Input    |      Filetype  | 
|----------|-------------   |
| Tree     |      Newick    | 
| Metadata |      TSV       |

### basic usage

Choose your input tree and corresponding metadata, MetaTreeProfiler will map all the metadata into corresponding tree node, and visualize it on the local server browser in interactive interface.

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
### **Simple map metadata into corresponding leaf nodes** 
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
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column_idx 6 --interactive  

# label column index by range, "[star_idx-end_idx]"
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column_idx [1-6] --interactive 
```

Categorical data will be process as counter in each internal node. Users can choose either counter is raw or relative count by using `--counter_stat`
| internal_node properties  |      statistic method  | 
|----------|-------------   |
| `<feature name>`_counter  |      raw, relative    | 

```
# raw count, example internal_node shown as: ```random_type_counter: medium--3||high--2```
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column random_type --counter_stat raw --interactive

# relative count, example internal_node shown as: ```random_type_counter: medium--0.60||high--0.40```
 internal_node example shown as, random_type_counter: medium--3||high--2
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column random_type --counter_stat relative --interactive
```


### Mapping Boolean data
For Boolean dataset, each internal node will count the selected feature(s) of its children nodes as counter as categorical, as shown as `_counter` of suffix feature name(s) of internal node. To label Boolean feature metadata, using following arguments


```
# label boolean data by column name(s) in metadata (for multiple columns, seperate by ","), using --bool_column <header>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column bool_type,bool_type2 --interactive

# label boolean data by column index, using --bool_column_idx <idx>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column_idx 7,8 --interactive  

# label column index by range, "[star_idx-end_idx]"
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column_idx [7-8] --interactive 
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
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column sample1,sample2,sample3,sample4,sample5 --interactive

# label numerical data by column index, using --num_column_idx <idx>
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx 1,2,3,4,5 --interactive  

# label column index by range, "[star_idx-end_idx]"
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx [1-5] --interactive 
```

By default, numerical feature will be calculated all the descriptive statistic, but users can choose specific one to be calculated by using `--num_stat [all, sum, avg, max, min, std] `

--num_stat NUM_STAT   statistic calculation to perform for numerical data in internal nodes, [all, sum, avg, max, min, std] 
```
# by default
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column sample1 --num_stat all --interactive

# only average calculation
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column sample1 --num_stat avg --interactive
```
### Mapping metadata wihtout column names
if metadata doesn't contain column names, please add `--no_colnames` as flag. MetaTreeProfiler will automatically assign feature name by index order

### Taxonomic profiling
If input metadada containcs taxon data, MetaTreeProfiler allows users to process taxonomic annotation with either GTDB or NCBI database.

- `--taxadb`, `NCBI` or `GTDB`, choose the Taxonomic Database for annotation
- `--taxon_column`, choose the column in metadata which representa taxon
- `--taxonomic_profile`, activate taxonomic annotation
- `--taxon_delimiter`, delimiter of taxa columns. default `.`
- `--taxa_field`, field of taxa name after delimiter. default `0`

#### Basic usage on GTDB
Here we demonstrate with `examples/gtdb_example1.nw` and `examples/gtdb_example1.tsv`
```
# in case of gtdb_example1.tsv
python treeprofiler.py --tree examples/gtdb_example1.nw --metadata examples/gtdb_example1.tsv --taxonomic_profile --taxon_column 0 --taxadb GTDB --interactive
```

#### Basic usage on NCBI
For instance of `examples/spongilla_example.nw` and `examples/spongilla_example.tsv`, it contains accession ID such as `83887.comp22273_c0_seq2_m.43352`, hence using `--` 
```
python treeprofiler.py --tree examples/spongilla_example.nw --metadata examples/spongilla_example.tsv --taxonomic_profile --taxon_column name --taxon_delimiter .  --taxa_field 0 --taxadb NCBI --interactive
```


## Visualizing annotated tree with layouts
MetaTreeProfiler provides a several of layout options for visualize features in metadata along with tree, depends on their datatype

### Layouts for categorical data
Users can add the following flag to activate layouts for categorical data
```
--ColorbranchLayout COLORBRANCHLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be plot as Textlayouts
--LabelLayout LABELLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be plot as LabelLayout
--RectangularLayout RECTANGULARLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be plot as RectangularLayout
```

example
```
## target column "random_type" in examples/basic_example1.tsv
# List random_type feature as text in aligned panel using LabelLayout
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column random_type --LabelLayout random_type --interactive

# Label random_type feature on branch with different colors in aligned panel  using --ColorbranchLayout
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column random_type --ColorbranchLayout random_type --interactive

# Label random_type feature with retangular block in aligned panel using --RectangularLayout
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column random_type --ColorbranchLayout random_type --interactive
```
### Layouts for boolean data
Users can add the following flag to activate layouts for Boolean data
```
--BinaryLayout BINARYLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be plot as BinaryLayout, label shown only positive value
--RevBinaryLayout REVBINARYLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be plot as RevBinaryLayout, label shown only negative value
```

```
## target column "bool_type", "bool_type2" in examples/basic_example1.tsv
# List postive bool_type feature in aligned panel using BinaryLayout
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column bool_type --BinaryLayout bool_type --interactive

# List negative bool_type feature in aligned panel using BinaryLayout
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column bool_type --BinaryLayout bool_type --interactive

# multiple columns seperated by ','
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --bool_column bool_type,bool_type2 --BinaryLayout bool_type,bool_type2  --interactive
```

### Layouts for Numerical data
Users can add the following flag to activate layouts for Numerical data
```
--HeatmapLayout HEATMAPLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be read as HeatmapLayout
--BarplotLayout BARPLOTLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be read as BarplotLayouts
```
```
## target column 'sample[1-5]' feature in examples/basic_example1.tsv
# visualize sample1 feature in Barplot
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column sample1 --BarplotLayout sample1 --interactive

# visualize sample1-sample5 in Heatmap
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx [1-5] --HeatmapLayout [1-5] --interactive 
```

### Visualizing annotated internal nodes
If internal nodes are annotated, MetaTreeProfiler is also able to visualize annotated features automatically when layouts are activated

#### Internal nodes of categorical and boolean data
As internal nodes of categorical and boolean data are annotated as counter, hence when activating layouts of categorical or boolean data, it generate pipechart of counter summary at the top of each internal node

#### Internal nodes of numerical data
Internal nodes of numerical data are process descriptive statistic analysis by default, hence when users collapse any branch, BarplotLayout or HeatmapLayout will demonstrate representative value, `avg` by default. representative value can be changed by using `--internal_plot_measure`

example
```
# select max instead of avg as internal node ploting representative
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx [1-5] --HeatmapLayout [1-5] --internal_plot_measure max --interactive
```

### Layouts for Taxonomic data
Activate Taxonomic layout using `--TaxonLayout`


## Conditional query in annotated tree
MetaTreeProfiler allows users to perform conditional processs based on differet circumstances
- `--pruned_by`, prune the annotated tree by conditions, and remove the branches which don't fit.
- `--collapsed_by`, collapse tree branches whose nodes if the conditions
- `--highlighted_by`, select tree nodes which fit the conditions
- `--rank_limit`, prune the taxonomic annotated tree based on rank of classification.

### Query Syntax
#### Basic Query
All the conditional query shared the same syntax, a standard query consists the following 

```
--pruned_by|collapsed_by|highlighted_by "<left_value> <operator> <right_value>"
```
* left value, the property of leaf node or internal node
* operators
    *  `=`
    * `!=`
    * `>` 
    * `>=`
    * `<`
    * `<=`
    * `contains`
* right value, custom value for the condition

Example 
```
# select tree node whose name contains `FALPE` character
...--highlighted_by "name contains FALPE"

# select tree node whose sample1 feature > 0.50
...--highlighted_by "sample1 > 0.50"
```

#### Query in internal nodes
Query in internal nodes' properties is also available, in this case, `left_value` of query will be the internal node property, remember to add the proper suffixes such as `_avg`, `_max`,etc, for the numerical data or `_counter` for categorical and boolean data. 

Example
```
# select tree internal node where sample1_avg feature > 0.50
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx [1-5] --HeatmapLayout [1-5] --highlighted_by "sample1_avg > 0.50" --interactive
```

Special syntax for internal node counter
```
# select tree internal nodes, where `low` relative counter < 0.30 in random_type_counter property
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --text_column random_type --counter_stat relative  --highlighted_by "random_type_counter:low < 0.30" --interactive
```

#### AND and OR conditions
The syntax for the AND condition and OR condition in MetaTreeProfiler is:

AND condition will be under one argument, syntax seperated by `,`, such as 
```
# select tree  node where sample1 feature > 0.50 AND sample2 < 0.2
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx [1-5] --HeatmapLayout [1-5]--highlighted_by "sample1>0.50,sample2<0.2" --interactive
```

OR condition will be used more than one arguments, such as
```
# select tree node where sample1 feature > 0.50 OR sample2 < 0.2
python treeprofiler.py --tree examples/basic_example1.nw --metadata examples/basic_example1.tsv --num_column_idx [1-5] --HeatmapLayout [1-5] --highlighted_by "sample1>0.50" --highlighted_by "sample2<0.2" --interactive
```

### conditional pruning based on taxonomic level
Prune taxonomic annotated tree based on following taxonomic rank level,
`kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `subspecies` 
```
# prune tree by family
... --rank_limit family
```
the above argument equals to 
```
... --pruned_by "rank=family"
```


params:

    input parameters:
        Input parameters

    -t TREE, --tree TREE  Input tree, .nw file, customized tree input
    -d METADATA, --metadata METADATA
                            <metadata.csv> .csv, .tsv. mandatory input
    --annotated_tree      inputtree already annot                 ated by treeprofileer
    --no_colnames         metadata table doesn't contain columns name
    --text_column TEXT_COLUMN
                            <col1,col2> names, column index or index range of columns which need to be read as categorical data
    --num_column NUM_COLUMN
                            <col1,col2> names, column index or index range of columns         which need to be read as numerical data
    --bool_column BOOL_COLUMN
                            <col1,col2> names, column index or index range of columns which need to be read as boolean data
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
                            <col1,col2> names, column index or index range of columns which need to be plot as BinaryLayout
    --RevBinaryLayout REVBINARYLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be plot as RevBinaryLayout
    --ColorbranchLayout COLORBRANCHLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be plot as Textlayouts
    --LabelLayout LABELLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be plot as LabelLayout
    --RectangularLayout RECTANGULARLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be plot as RectangularLayout
    --HeatmapLayout HEATMAPLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be read as HeatmapLayout
    --BarplotLayout BARPLOTLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be read as BarplotLayouts
    --TaxonLayout TAXONLAYOUT
                            <col1,col2> names, column index or index range of columns which need to be read as TaxonLayouts

    Output arguments:
    Output parameters

    --interactive         run interactive session
    --port PORT           run interactive session
    --plot PLOT           output as pdf
    -o OUTTREE, --outtree OUTTREE
                            output annotated tree
    --outtsv OUTTSV       output annotated tsv file
