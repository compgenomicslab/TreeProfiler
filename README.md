# TreeProfiler Tutorial
## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
    - [Quick install](#quick-install-via-pip)
    - [Quick Start](#quick-start-with-examples-dataset)
    - [Manual installation](#manual-installation)
    - [Input files](#input-files)
    - [Basic Usage](#basic-usage)
- [Using TreeProfiler](#using-treeprofiler) 
  - [Parsing Input tree](#parsing-input-tree)
    - [Tree format](#tree-format)
    - [Tree parser](#tree-parser)
  - [Annotate metadata into tree ](#annotate-annotate-metadata-into-tree)
      - [Annotate multiple metadata inputs into tree](#annotate-multiple-metadata-inputs-into-tree)
      - [Annotate metadata into tree internal nodes](#annotate-metadata-into-tree-internal-nodes)
      - [Determine datatype in arguments](#determine-datatype-in-arguments)
      - [Mapping metadata without column names](#mapping-metadata-without-column-names)
      - [Taxonomic annotation](#taxonomic-annotation)
          - [Identify taxon property in metadata](#identify-taxon-property-in-metadata)
          - [Basic usage on GTDB](#basic-usage-on-GTDB)
          - [Basic usage on NCBI](#basic-usage-on-NCBI)
      - [Annotation from eggnog-mapper output](#annotation-from-eggnog-mapper-output)
      - [Annotated tree format](#annotate-tree-format)
  - [Plot annotated tree with layouts](#plot-plot-annotated-tree-with-layouts)
      - [Interactive visualization interface](#interactive-visualization-interface)
      - [Basic options of visualizing layouts](#basic-options-of-visualizing-layouts)
      - [Layouts for categorical data](#layouts-for-categorical-data)
      - [Layouts for boolean data](#layouts-for-boolean-data)
      - [Layouts for numerical data](#layouts-for-numerical-data)
      - [Layouts for multiple text data](#layouts-for-list-data)
      - [Layouts for multiple sequence alignment](#layouts-for-multiple-sequence-alignment)
      - [Layouts for eggnog-mapper pfam annotations](#layouts-for-eggnog-mapper-pfam-annotations)
      - [Layouts for eggnog-mapper smart annotations](#layouts-for-eggnog-mapper-smart-annotations)
      - [Layouts for eggnog-mapper annotations](#layouts-for-eggnog-mapper-annotations)
      - [Visualizing annotated internal nodes](#visualizing-annotated-internal-nodes)
        - [Internal nodes of categorical and boolean data](#internal-nodes-of-categorical-and-boolean-data)
        - [Internal nodes of numerical data](#internal-nodes-of-numerical-data)
      - [Layouts for Taxonomic data](#layouts-for-taxonomic-data)
  - [Conditional query in annotated tree](#conditional-query-in-annotated-tree)
      - [Basic Query](#basic-query)
      - [Query in internal nodes](#query-in-internal-nodes)
      - [AND and OR conditions](#and-and-or-conditions)
      - [Conditional pruning based on taxonomic level](#conditional-limit-based-on-taxonomic-level)
- [Demo1 Explore GTDB and Progenome data](#demo1-explore-gtdb-taxonomic-tree-with-metadata-and-habitat-information-of-progenome3)
- [Demo2 Explore Nif gene family tree with functional annotations data using eggnog-mapper with taxonomic annotation](#demo2-explore-large-nifh-gene-tree-with-functional-and-taxonomic-information)

## Introduction
TreeProfiler is command-line tool for profiling metadata table into phylogenetic tree with descriptive analysis and output visualization

## Installation
### Dependencies
TreeProfiler requires 
  - Python version >= 3.7
  - ETE Toolkit v4
  - biopython
  - selenium
  - scipy
  - matplotlib

### Quick install via pip
```
# Install ETE Toolkit v4
pip install https://github.com/etetoolkit/ete/archive/refs/tags/4.1.0-beta.tar.gz


# Install TreeProfiler dependencies
pip install biopython selenium scipy matplotlib

# Install TreeProfiler tool via pip
pip install TreeProfiler
```


### Quick Start with examples dataset
TreeProfiler provide various example dataset for testing in `examples/` or https://github.com/compgenomicslab/TreeProfiler/tree/main/examples,
each directory consists a demo script `*_demo.sh` for quick starting different functions in TreeProfiler which alreadyh as annotate-plot pipeline of example data. User can fast explore different example tree with different visualizations. Here is the demonstration:

```
# execute demo script of example1
cd examples/basic_example1/
sh ./example1_demo.sh

Annotate example tree with two metadata tables
start parsing...
Time for parse_csv to run:  0.001968860626220703
Time for load_metadata_to_tree to run:  0.0003094673156738281
Time for merge annotations to run:  0.05160331726074219
Time for annotate_taxa to run:  4.76837158203125e-07
Visualize properties categorical data random_type in rectangle_layout, numerical data sample1, sample2 in heatmap_layout and barplot_layout.
Current trees in memory: 0
Added tree example with id 0.
 * Serving Flask app 'ete4.smartview.gui.server' (lazy loading)
 * Environment: production
   WARNING: This is a development server. Do not use it in a production deployment.
   Use a production WSGI server instead.
 * Debug mode: on
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```
As the session starts in local server http://127.0.0.1:5000, annotated tree and selected properties are visualized at the interactive session.
![treeprofiler interface](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/figure1_all.png?raw=true)
Here is detailed introduction of interactive session of visualization([here](#interactive-visualization-interface))

Check other tutorial scripts
```
# display demo script of each example
./examples/basic_example1/example1_demo.sh
./examples/automatic_query/highlight_demo.sh
./examples/automatic_query/collapse_demo.sh
./examples/automatic_query/prune_demo.sh
./examples/basic_example2/example2_demo.sh
./examples/taxonomy_example/ncbi/ncbi_demo.sh
./examples/taxonomy_example/gtdb/gtdb_demo.sh
./examples/pratical_example/progenome3/progenome_demo.sh
./examples/pratical_example/gtdb_r202/gtdbv202full_demo.sh
./examples/pratical_example/gtdb_r202/gtdbv202lite_demo.sh
./examples/pratical_example/emapper/emapper_demo.sh
```


### Manual installation
#### Install ETE v4
Quick way
```
pip install https://github.com/etetoolkit/ete/archive/ete4.zip
```
For local development
To install ETE in a local directory to help with the development, you can:

- Clone this repository (git clone https://github.com/etetoolkit/ete.git)
- Install dependecies
  - If you are using conda: `conda install -c conda-forge cython bottle brotli numpy pyqt`
  - Otherwise, you can install them with `pip install <dependencies>`
  - Build and install ete4 from the repository's root directory: `pip install -e .`

(In Linux there may be some cases where the gcc library must be installed, which can be done with `conda install -c conda-forge gcc_linux-64`)

#### Install TreeProfiler
Install dependencies
```
# install BioPython, selenium, scipy via conda
conda install -c conda-forge biopython selenium scipy matplotlib
# or pip
pip install biopython selenium scipy matplotlib
```

Install TreeProfiler
```
# install TreeProfiler
git clone https://github.com/compgenomicslab/TreeProfiler
cd TreeProfiler/
python setup.py install
```
Or 
```
# install directly
pip install https://github.com/dengzq1234/TreeProfiler/archive/refs/tags/v1.1.0.tar.gz
```

### Input files
TreeProfiler takes following file types as input 

| Input    |      Filetype  | 
|----------|-------------   |
| Tree     |      newick, ete    | 
| Metadata |      tar.gz, tsv       |

- ete format is a novel format developed to solve the situation we encounter in the previous step, annotated tree can be recover easily with all the annotated data without changing the data type. Besides, the ete format optimized the tree file size after mapped with its associated data. Hence it's very handy for programers in their own script. At this moment we can only view the ete format in treeprofiler, but we will make the ete format more universal to other phylogenetic software.
- Metadata input could be single or multiple files, either tar.gz compressed file(s) which contains multiple .tsv or plain .tsv file(s). 

### Basic Usage
TreeProfiler has two main subcommand:
 - annotate
 - plot

The first one `annotate` is used to annotate your input tree and corresponding metadata, TreeProfiler will map all the metadata into corresponding tree node. In this step, annotated tree will be generated in newick and ete format

```
treeprofiler annotate --tree tree.nw --input-type newick --metadata metadata.tsv --outdir ./
```

The second subcommand `plot` is used to visualize tree with associated metadata. By default, treeprofiler will launch an interactive session at localhost for user to explore input tree.

```
treeprofiler plot --tree tree_annotated.nw --input-type newick 
```
or

```
treeprofiler plot --tree tree_annotated.ete --input-type ete 
```



# Using TreeProfiler
In this Tutorial we will use TreeProfiler and demostrate basic usage with data in `examples/`

```
tree examples/
examples/
├── automatic_query
│   ├── basic_example1_metadata1.tsv
│   ├── basic_example1.nw
│   ├── collapse_demo.sh
│   ├── highlight_demo.sh
│   └── prune_demo.sh
├── basic_example1
│   ├── basic_example1_metadata1.tsv
│   ├── basic_example1_metadata2.tsv
│   ├── basic_example1.nw
│   └── example1_demo.sh
├── basic_example2
│   ├── diauxic.array
│   ├── diauxic.nw
│   ├── example2_demo.sh
│   ├── FluA_H3_AA.fas
│   ├── MCC_FluA_H3_Genotype.txt
│   └── MCC_FluA_H3.nw
├── pratical_example
│   ├── emapper
│   │   ├── 7955.ENSDARP00000116736.aln.faa
│   │   ├── 7955.ENSDARP00000116736.nw
│   │   ├── 7955.out.emapper.annotations
│   │   ├── 7955.out.emapper.pfam
│   │   ├── 7955.out.emapper.smart.out
│   │   ├── emapper_demo.sh
│   │   ├── nifH.faa.aln
│   │   ├── nifH.nw
│   │   ├── nifH.out.emapper.annotations
│   │   └── nifH.out.emapper.pfam
│   ├── gtdb_r202
│   │   ├── ar122_metadata_r202_lite.tar.gz
│   │   ├── bac120_metadata_r202_lite.tar.gz
│   │   ├── gtdbv202full_demo.sh
│   │   ├── gtdbv202lite_demo.sh
│   │   ├── gtdbv202.nw
│   │   ├── merge_gtdbtree.py
│   │   └── progenome3.tar.gz
│   └── progenome3
│       ├── progenome3.nw
│       ├── progenome3.tsv
│       └── progenome_demo.sh
└── taxonomy_example
    ├── gtdb
    │   ├── gtdb_demo.sh
    │   ├── gtdb_example1.nw
    │   └── gtdb_example1.tsv
    └── ncbi
        ├── ncbi_demo.sh
        ├── ncbi_example.nw
        └── ncbi_example.tsv

```
## Parsing Input tree
### Tree format
TreeProfiler accpept input tree in `.nw` or `.ete` by putting `--input-type {newick,ete}` flag to identify `[default: ete]`. The difference between `.nw` and `.ete`, 

 - `newick` format is more universal and be able to used in different other phylogenetic software although associated data of tree nodes will be considered as plain text.

 - `ete` format is a novel format developed to solve the situation we encounter in the previous step, annotated tree can be **recover easily with all the annotated data without changing the data type**. Besides, the ete format optimized the tree file size after mapped with its associated data. Hence it's very handy for programers in their own script. At this moment we can only view the ete format in treeprofiler, but we will make the ete format more universal to other phylogenetic software. **Hence using ete format in `plot` subcommand is highly reccomended**

### Tree parser
TreeProfiler provides argument `--internal-parser {name,support}` to specify `newick` tree when it include values in internal node. `[default: name]`

| newick  |      leaves  |  internal_node value |  internal_parser 
|----------|-------------   |-------------   |-------------   |
| (A:0.5, B:0.5)Internal_C:0.5;  |  A, B| Internal_C| `name`| 
| (A:0.5, B:0.5)0.99:0.5;  |   A, B| 0.99| `support`| 


## `annotate`, Annotate metadata into tree 
```
 treeprofiler annotate -h
usage: treeprofiler annotate [-h] -t TREE [--annotated-tree] [--internal-parser {name,support}]
                             [--input-type {newick,ete}] [--prop2type PROP2TYPE]
                             [--rank-limit RANK_LIMIT] [--pruned-by PRUNED_BY]
                             [-d METADATA [METADATA ...]] [--no-colnames] [--aggregate-duplicate]
                             [--text-prop TEXT_PROP [TEXT_PROP ...]]
                             [--multiple-text-prop MULTIPLE_TEXT_PROP [MULTIPLE_TEXT_PROP ...]]
                             [--num-prop NUM_PROP [NUM_PROP ...]]
                             [--bool-prop BOOL_PROP [BOOL_PROP ...]]
                             [--text-prop-idx TEXT_PROP_IDX] [--num-prop-idx NUM_PROP_IDX]
                             [--bool-prop-idx BOOL_PROP_IDX] [--taxadb TAXADB]
                             [--taxon-column TAXON_COLUMN] [--taxon-delimiter TAXON_DELIMITER]
                             [--taxa-field TAXA_FIELD]
                             [--emapper-annotations EMAPPER_ANNOTATIONS]
                             [--emapper-pfam EMAPPER_PFAM] [--emapper-smart EMAPPER_SMART]
                             [--alignment ALIGNMENT] [--taxonomic-profile]
                             [--num-stat {all,sum,avg,max,min,std,none}]
                             [--counter-stat {raw,relative,none}] -o OUTDIR

annotate tree

optional arguments:
  -h, --help            show this help message and exit

SOURCE TREE INPUT:
  Source tree input parameters

  -t TREE, --tree TREE  Input tree, .nw file, customized tree input
  --annotated-tree      input tree already annotated by treeprofiler if you want to skip the
                        annotate part.
  --internal-parser {name,support}
                        To specify how to interpret internal nodes in newick format. [default:
                        name]
  --input-type {newick,ete}
                        Specify input tree format. [newick, ete]. [default: ete]
  --prop2type PROP2TYPE
                        config tsv file where determine the datatype of target properties, if
                        your input tree type is .ete, it's note necessary

Pruning parameters:
  Auto pruning parameters

  --rank-limit RANK_LIMIT
                        TAXONOMIC_LEVEL prune annotate tree by rank limit
  --pruned-by PRUNED_BY
                        target tree pruned by customized conditions, such as --pruned-by "name
                        contains FALPE"

METADATA TABLE parameters:
  Input parameters of METADATA

  -d METADATA [METADATA ...], --metadata METADATA [METADATA ...]
                        <metadata.csv> .csv, .tsv. mandatory input
  --no-colnames         metadata table doesn't contain columns name
  --aggregate-duplicate
                        treeprofiler will aggregate duplicated metadata to a list as a property
                        if metadata contains duplicated row
  --text-prop TEXT_PROP [TEXT_PROP ...]
                        <col1,col2> names, column index or index range of columns which need to
                        be read as categorical data
  --multiple-text-prop MULTIPLE_TEXT_PROP [MULTIPLE_TEXT_PROP ...]
                        <col1,col2> names, column index or index range of columns which need to
                        be read as categorical data which contains more than one value and
                        seperate by ',' such as GO:0000003,GO:0000902,GO:0000904,GO:0003006
  --num-prop NUM_PROP [NUM_PROP ...]
                        <col1,col2> names, column index or index range of columns which need to
                        be read as numerical data
  --bool-prop BOOL_PROP [BOOL_PROP ...]
                        <col1,col2> names, column index or index range of columns which need to
                        be read as boolean data
  --text-prop-idx TEXT_PROP_IDX
                        1,2,3 or [1-5] index of columns which need to be read as categorical data
  --num-prop-idx NUM_PROP_IDX
                        1,2,3 or [1-5] index columns which need to be read as numerical data
  --bool-prop-idx BOOL_PROP_IDX
                        1,2,3 or [1-5] index columns which need to be read as boolean data
  --taxadb TAXADB       <NCBI|GTDB> for taxonomic profiling or fetch taxatree [default: GTDB]
  --taxon-column TAXON_COLUMN
                        <col1> name of columns which need to be read as taxon data
  --taxon-delimiter TAXON_DELIMITER
                        delimiter of taxa columns. [default: None]
  --taxa-field TAXA_FIELD
                        field of taxa name after delimiter. [default: 0]
  --emapper-annotations EMAPPER_ANNOTATIONS
                        attach eggNOG-mapper output out.emapper.annotations
  --emapper-pfam EMAPPER_PFAM
                        attach eggNOG-mapper pfam output out.emapper.pfams
  --emapper-smart EMAPPER_SMART
                        attach eggNOG-mapper smart output out.emapper.smart
  --alignment ALIGNMENT
                        Sequence alignment, .fasta format

Annotation arguments:
  Annotation parameters

  --taxonomic-profile   Activate taxonomic annotation on tree
  --num-stat {all,sum,avg,max,min,std,none}
                        statistic calculation to perform for numerical data in internal nodes,
                        [all, sum, avg, max, min, std, none]. If 'none' was chosen, numerical
                        properties won't be summarized nor annotated in internal nodes. [default:
                        all]
  --counter-stat {raw,relative,none}
                        statistic calculation to perform for categorical data in internal nodes,
                        raw count or in percentage [raw, relative, none]. If 'none' was chosen,
                        categorical and boolean properties won't be summarized nor annotated in
                        internal nodes [default: raw]

OUTPUT options:

  -o OUTDIR, --outdir OUTDIR
                        Directory for annotated outputs.
```
### Annotate multiple metadata inputs into tree
TreeProfiler allows user to annotate more than one metadata inputs to tree. It requires to seperated the input tsv or tar.gz files with space in the arguments `--metadata` or `-d`, such as `--metadata table1.tsv table2.tsv`. Here is the example

Check metadata
```
tree examples/basic_example1/
examples/basic_example1/
├── basic_example1_metadata2.tsv
├── basic_example1_miss.tsv
├── basic_example1_null.tsv
├── basic_example1.nw
├── basic_example1_metadata1.tsv
└── example1_demo.sh

head -5 examples/basic_example1/basic_example1_metadata1.tsv 
#name	sample1	sample2	sample3	sample4	sample5	random_type	*bool_type	bool_type2
Phy003I7ZJ_CHICK	0.05	0.12	0.86	0.01	0.69	medium	1	TRUE
Phy0054BO3_MELGA	0.64	0.67	0.51	0.29	0.14	medium	1	TRUE
Phy00508FR_NIPNI	0.89	0.38	0.97	0.49	0.26	low	1	FALSE
Phy004O1E0_APTFO	0.1	0.09	0.38	0.31	0.41	medium	0	TRUE

head -5 examples/basic_example1/basic_example1_metadata2.tsv
#name	abs_data	list_data
Phy003I7ZJ_CHICK	97	w,t,t
Phy0054BO3_MELGA	16	r,q,s
Phy00508FR_NIPNI	87	z,f,p
Phy004O1E0_APTFO	6	z,t,b
```

Run `annotate` subcommand
```
## annotate tree with more than one metadata tsv, seperated by space
treeprofiler annotate \
--tree examples/basic_example1/basic_example1.nw \
--input-type newick \
--metadata examples/basic_example1/basic_example1_metadata1.tsv examples/basic_example1/basic_example1_metadata2.tsv \
-o examples/basic_example1/
```

### Annotate metadata into tree internal nodes
At the above example, we only mapped metadata to leaf nodes, in this example, we will also profile **internal nodes** annotation and analysis of their children nodes.

TreeProfiler can infer automatically the datatype of each column in your metadata, including 
- `list` (seperate by `,` )
- `string` (categorcial data)
- `numerical`(numerical data, float or integer)
- `booleans` 

Internal node will summurize children nodes information according to their datatypes.

demo tree
```
      ╭╴A
╴root╶┤
      │   ╭╴B
      ╰╴D╶┤
          ╰╴C
```

demo metadata
|  #name | text_property |  multiple_text_property  |   numerical_property  | bool_property| 
|----------|----------|----------|-------------|-------------|
|A|vowel|a,b,c|10|True|
|B|consonant|b,c,d|4|False|
|C|consonant|c,d,e|9|True|

Treeprofiler will infer the datatypes of above metadata and adpot different summary method:
|  - | text_property |  multiple_text_property  |   numerical_property  | bool_property| 
|----------|----------|----------|-------------|-------------|
|datatype|string|list|float|bool|
|method|counter|counter|average,sum,max,min,standard deviation|counter|

After annotation, internal nodes will be summarized. If property was summarize with `counter`, in internal node will be named as ```<property_name>_counter```


Users can choose either counter is raw or relative count by using `--counter-stat`
| internal_node properties  |      statistic method  | 
|----------|-------------   |
| `<feature name>`_counter  |      raw(default), relative    | 

|  internal_node| text_property_counter |  multiple_text_property_counter  | bool_property_counter| 
|----------|----------|-------------|-------------|
|D|consonant--2|b--1\|\|c--2\|\|d--2\|\|e--1|True--1\|\|False--1|
|root|vowel--1\|\|consonant--2|a--2\|\|b--2\|\|c--3\|\|d--2\|\|e--1|True--2\|\|False--1|

After annotation, internal nodes will be summarized.  If property was numerical data, in internal node will be named as 
| internal_node properties  |      statistic method  | 
|----------|-------------   |
| `<feature name>`_avg      |      average    | 
| `<feature name>`_sum      |      sum    | 
| `<feature name>`_max      |      maximum    | 
| `<feature name>`_min      |      minimum    | 
| `<feature name>`_std      |      standard deviation    | 

By default, numerical feature will be calculated all the descriptive statistic, but users can choose specific one to be calculated by using `--num_stat [all, sum, avg, max, min, std] `

In our demo, it would be:
|  internal_node| numerical_property_avg |  numerical_property_sum  | numerical_property_max| numerical_property_max| numerical_property_max| 
|----------|----------|-------------|-------------|-------------|-------------|
|D| 6.5| 13| 9| 4|2.5|
|root| 7.67 | 23| 10| 4| 2.32| 

Excecute example data provided in `examples/`

```sh
treeprofiler annotate \
--tree examples/basic_example1/basic_example1.nw \
--input-type newick \
--metadata examples/basic_example1/basic_example1_metadata1.tsv \
--outdir ./examples/basic_example1/
```

### Determine datatype in arguments
Although TreeProfiler can detect datatype of each column, users still can determine the datatype using the following arguments using

- `--text-prop` and `--text-prop-idx`, to determine columms which need to be read as categorical data, example 

- `--multiple-text-prop`, to determine columns which contains multiple values sperated by `,`, and will be process as list

- `--num-prop` and `--num-prop-idx`, to determine columms which need to be read as numerical data

- `--bool-prop` and `--bool-prop-idx`, to determine columms which need to be read as boolean data

#### Missing value detection
Metadata column which fullfills one of the following criterias will be consider as missing value:

- Entirely non-word characters. Such as `+`, `-`, `~`, `.`, etc.
- The exact strings `none`, `None`, `null`, `Null`, or `NaN`.
- An empty string (zero characters).

Missing value will replaced by string 'NaN' in the corresponding property.

#### Unmapped Tree leaf property detection
If Metadata doesn't cover input tree leaf, tree leaf will be unannotated.  

### Mapping metadata without column names
if metadata doesn't contain column names, please add `--no-colnames` as flag. TreeProfiler will automatically assign feature name by index order

### Taxonomic annotation
If input metadada containcs taxon data, TreeProfiler allows users to process taxonomic annotation with either GTDB or NCBI database.

- `--taxadb`, `NCBI` or `GTDB`, choose the Taxonomic Database for annotation
- `--taxonomic-profile`, activate taxonomic annotation
- `--taxon-column`, choose the column in metadata which representa taxon. default is the first column which should be column of leaf_name
- `--taxon-delimiter`, delimiter of taxa columns. default `''`
- `--taxa-field`, field of taxa name after delimiter. default `0`

 
#### Basic usage on GTDB 
Here we demonstrate with `examples/taxonomy_example/gtdb/gtdb_example1.nw` and `examples/taxonomy_example/gtdb/gtdb_example1.tsv`. Taxonomic accesion IDs are located in the first column which should be the names of leaf. If accesions are located in different columns, using `--taxon-column <column name>` to locate the the column. 

```
# in case of gtdb_example1.tsv
head -3 examples/taxonomy_example/gtdb/gtdb_example1.tsv
name	sample1	sample2	sample3	sample4	sample5	random_type	bool_type	bool_type2
RS_GCF_001560035.1	0.05	0.12	0.86	0.01	0.69	medium	1	True
RS_GCF_001560635.1	0.64	0.67	0.51	0.29	0.14	medium	1	True

# annotate tree with gtdb taxonomic annotation
treeprofiler annotate \
--tree examples/taxonomy_example/gtdb/gtdb_example1.nw \
--input-type newick \
--metadata examples/taxonomy_example/gtdb/gtdb_example1.tsv --taxonomic-profile \
--taxadb GTDB \
--outdir ./examples/taxonomy_example/gtdb/
```

#### Basic usage on NCBI
For instance of `examples/taxonomy_example/ncbi/ncbi_example.nw` and `examples/taxonomy_example/ncbi/ncbi_example.tsv`, it contains accession ID such as `83887.comp22273_c0_seq2_m.43352`, hence using `--taxon-delimiter` and `--taxa-field` to locate the taxonomic accession.

```
treeprofiler annotate \
--tree examples/taxonomy_example/ncbi/ncbi_example.nw \
--metadata examples/taxonomy_example/ncbi/ncbi_example.tsv \
--taxonomic-profile \
--taxon-column name \
--taxadb NCBI \
--taxon-delimiter . \
--taxa-field 0  \
--outdir ./examples/taxonomy_example/ncbi/
```

#### Identify taxon property in metadata
When Taxon properties are embeded in different column or field in metadata, treeprofiler provides `--taxon-column`, `--taxon-delimiter` and `--taxa-field` to identify taxon term in order to process taxonomic annotation sucessfully. Here is summary of different cases with corresponding setting.

| metadata |taxon to be identified |       command line setting  | 
|----------|-------------   | ----|
| `#leafname col1`<br>`9598 wt`     | 9598|     `default` | 
| `#leafname col1`<br>`7739.XP_002609184.1 wt`     |   7739|   `--taxon-column <default> --taxon-delimiter . --taxa-field 0`   | 
| `#leafname ncbi_id`<br>`leaf_A 7739`    | 7739|     `--taxon-column ncbi_id --taxon-delimiter <default> --taxa-field <default> `    | 
| `#leafname ncbi_id`<br>`leaf_A 7739.XP_002609184.1`     |   7739|      `--taxon-column ncbi_id --taxon-delimiter . --taxa-field 0`    | 
|`#leafname col1`<br> `RS_GCF_001560035.1 wt`   | RS_GCF_001560035.1|     `default`   |
| `#leafname gtdb_id`<br>`leaf_A d__Archaea;p__Asgardarchaeota;c__Heimdallarchaeia;o__UBA460;f__Kariarchaeaceae;g__LC-2;s__LC-2 sp001940725`      | s__LC-2 sp001940725|     `--taxon-column gtdb_id --taxon-delimiter ; --taxa-field -1`   |

### Annotation from EggNOG-mapper output
[EggNOG-mapper](http://eggnog-mapper.embl.de/), is a tool for fast functional annotation of novel sequences. It uses precomputed orthologous groups and phylogenies from the eggNOG database (http://eggnog5.embl.de) to transfer functional information from fine-grained orthologs only. 

It generates three kind of ouput file, 

1) Raw standard output, `*.out.emapper.annotations`, that contains functional annotations and prthology predictions, for example:
```
## Mon Feb 27 09:05:50 2023
## emapper-2.1.9
## /data/shared/home/emapper/miniconda3/envs/eggnog-mapper-2.1/bin/emapper.py --cpu 20 --mp_start_method forkserver --data_dir /dev/shm/ -o out --output_dir /emapper_web_jobs/emapper_jobs/user_data/MM_knn6rw6j --temp_dir /emapper_web_jobs/emapper_jobs/user_data/MM_knn6rw6j --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i /emapper_web_jobs/emapper_jobs/user_data/MM_knn6rw6j/queries.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign denovo --num_servers 2 --report_orthologs --decorate_gff yes --excel
##
#query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs
....
## 272 queries scanned
## Total time (seconds): 45.73449420928955
## Rate: 5.95 q/s
```
2) [Pfam](http://pfam.xfam.org/) domain annotations, `*.out.emapper.pfam`, for example:
```
## Mon Feb 27 09:05:52 2023
## emapper-2.1.9
## /data/shared/home/emapper/miniconda3/envs/eggnog-mapper-2.1/bin/emapper.py --cpu 20 --mp_start_method forkserver --data_dir /dev/shm/ -o out --output_dir /emapper_web_jobs/emapper_jobs/user_data/MM_knn6rw6j --temp_dir /emapper_web_jobs/emapper_jobs/user_data/MM_knn6rw6j --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i /emapper_web_jobs/emapper_jobs/user_data/MM_knn6rw6j/queries.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign denovo --num_servers 2 --report_orthologs --decorate_gff yes --excel
##
# query_name	hit	evalue	sum_score	query_length	hmmfrom	hmmto	seqfrom	seqto	query_coverage
...
## 272 queries scanned
## Total time (seconds): 28.74908423423767
## Rate: 9.46 q/s
```

3) [SMART](http://smart.embl-heidelberg.de/) domain annotation, `*.out.emapper.smart.out`, for example:
```
10020.ENSDORP00000023664	MAGE_N	10	63	220000.115599899
10020.ENSDORP00000023664	PTN	44	128	683.160049964146
10020.ENSDORP00000023664	Ephrin_rec_like	73	117	248282.169266432
10020.ENSDORP00000023664	PreSET	87	186	494.036044144428
....
```

TreeProfiler allows users annotate EggNOG-mapper  standard output to target tree with following arguments
 - `--emapper-annotations`, attach eggNOG-mapper output `out.emapper.annotations`.
 - `--emapper-pfam`, attach eggNOG-mapper pfam output `out.emapper.pfams`.
 - `--emapper-smart`, attach eggNOG-mapper smart output `out.emapper.smart`.

 [check EggNOG-mapper annotation example](#demo2-explore-eggnog-mapper-annotations-data-with-taxonomic-annotation)

### **Annotate tree format**
TreeProfiler `annotate` subcommand will generate the following output file

1) `<input_tree>` + *_annotated.nw*, newick format with annotated tree
2) `<input_tree>` + *_annotated.ete*, ete format with annotated tree
3) `<input_tree>` + *_annotated_prop2type.txt*, config file where store the datatype of each annotated properties
4) `<input_tree>` + *_annotated.tsv*,  metadata in tab-separated values format with annotated and summarized internal nodes information. 

In the following `plot` step, users can use either `.nw` or `.ete` by putting `--input-type [newick, ete]` flag to identify. The difference between `.nw` and `.ete` format is 

 - newick file is more universal and be able to used in different other phylogenetic software although associated data of tree nodes will be considered as plain text, so if you use newick format, alongside with the prop2type config file which was generated before by adding `--prop2type <prop2type_file>`

 - ete format is a novel format developed to solve the situation we encounter in the previous step, annotated tree can be **recover easily with all the annotated data without changing the data type**. Besides, the ete format optimized the tree file size after mapped with its associated data. Hence it's very handy for programers in their own script. At this moment we can only view the ete format in treeprofiler, but we will make the ete format more universal to other phylogenetic software. **Hence using ete format in `plot` subcommand is highly reccomended**


## `plot`, Plot annotated tree with layouts
TreeProfiler provides a several of layout options for visualize features in metadata along with tree, depends on their datatype
```
 treeprofiler plot -h
usage: treeprofiler plot [-h] -t TREE [--annotated-tree] [--internal-parser {name,support}]
                         [--input-type {newick,ete}] [--prop2type PROP2TYPE]
                         [--rank-limit RANK_LIMIT] [--pruned-by PRUNED_BY]
                         [--internal-plot-measure INTERNAL_PLOT_MEASURE]
                         [--collapsed-by COLLAPSED_BY] [--highlighted-by HIGHLIGHTED_BY]
                         [--column-width COLUMN_WIDTH] [--barplot-width BARPLOT_WIDTH]
                         [--padding-x PADDING_X] [--padding-y PADDING_Y]
                         [--binary-layout BINARY_LAYOUT [BINARY_LAYOUT ...]]
                         [--revbinary-layout REVBINARY_LAYOUT [REVBINARY_LAYOUT ...]]
                         [--colorbranch-layout COLORBRANCH_LAYOUT [COLORBRANCH_LAYOUT ...]]
                         [--label-layout LABEL_LAYOUT [LABEL_LAYOUT ...]]
                         [--rectangle-layout RECTANGLE_LAYOUT [RECTANGLE_LAYOUT ...]]
                         [--heatmap-layout HEATMAP_LAYOUT [HEATMAP_LAYOUT ...]]
                         [--barplot-layout BARPLOT_LAYOUT [BARPLOT_LAYOUT ...]]
                         [--taxonclade-layout] [--taxonrectangle-layout] [--emapper-layout]
                         [--domain-layout] [--alignment-layout]
                         [--profiling-layout PROFILING_LAYOUT [PROFILING_LAYOUT ...]]
                         [--multi-profiling-layout MULTI_PROFILING_LAYOUT [MULTI_PROFILING_LAYOUT ...]]
                         [--categorical-matrix-layout CATEGORICAL_MATRIX_LAYOUT [CATEGORICAL_MATRIX_LAYOUT ...]]
                         [--numerical-matrix-layout NUMERICAL_MATRIX_LAYOUT [NUMERICAL_MATRIX_LAYOUT ...]]
                         [--verbose] [--port PORT] [--plot PLOT] [--out-colordict]

annotate plot

optional arguments:
  -h, --help            show this help message and exit

SOURCE TREE INPUT:
  Source tree input parameters

  -t TREE, --tree TREE  Input tree, .nw file, customized tree input
  --annotated-tree      input tree already annotated by treeprofiler if you want to skip the
                        annotate part.
  --internal-parser {name,support}
                        To specify how to interpret internal nodes in newick format. [default:
                        name]
  --input-type {newick,ete}
                        Specify input tree format. [newick, ete]. [default: ete]
  --prop2type PROP2TYPE
                        config tsv file where determine the datatype of target properties, if
                        your input tree type is .ete, it's note necessary

Pruning parameters:
  Auto pruning parameters

  --rank-limit RANK_LIMIT
                        TAXONOMIC_LEVEL prune annotate tree by rank limit
  --pruned-by PRUNED_BY
                        target tree pruned by customized conditions, such as --pruned-by "name
                        contains FALPE"

Conditional display arguments:
  Conditional display parameters

  --internal-plot-measure INTERNAL_PLOT_MEASURE
                        statistic measures to be shown in numerical layout for internal nodes,
                        [default: avg]
  --collapsed-by COLLAPSED_BY
                        target tree nodes collapsed by customized conditions
  --highlighted-by HIGHLIGHTED_BY
                        target tree nodes highlighted by customized conditions

Properties' layout arguments:
  Prop layout parameters

  --column-width COLUMN_WIDTH
                        customize column width of each layout.[default: 20]
  --barplot-width BARPLOT_WIDTH
                        customize barplot width of barplot layout.[default: 200]
  --padding-x PADDING_X
                        customize horizontal column padding distance of each layout.[default: 1]
  --padding-y PADDING_Y
                        customize vertical padding distance of each layout.[default: 0]
  --binary-layout BINARY_LAYOUT [BINARY_LAYOUT ...]
                        <prop1,prop2> names of properties which need to be plot as binary-layout
                        which highlights the postives
  --revbinary-layout REVBINARY_LAYOUT [REVBINARY_LAYOUT ...]
                        <prop1,prop2> names of properties which need to be plot as revbinary-
                        layout which highlights the negatives
  --colorbranch-layout COLORBRANCH_LAYOUT [COLORBRANCH_LAYOUT ...]
                        <prop1,prop2> names of properties where branches will be colored based on
                        different values.
  --label-layout LABEL_LAYOUT [LABEL_LAYOUT ...]
                        <prop1,prop2> names of properties where values will be displayed on the
                        aligned panel.
  --rectangle-layout RECTANGLE_LAYOUT [RECTANGLE_LAYOUT ...]
                        <prop1,prop2> names of properties where values will be label as
                        rectangular color block on the aligned panel.
  --heatmap-layout HEATMAP_LAYOUT [HEATMAP_LAYOUT ...]
                        <prop1,prop2> names of numerical properties which need to be read as
                        heatmap-layout
  --barplot-layout BARPLOT_LAYOUT [BARPLOT_LAYOUT ...]
                        <prop1,prop2> names of numerical properties which need to be read as
                        barplot_layouts
  --taxonclade-layout   Activate taxonclade_layout which clades will be colored based on taxonomy
                        of each node.
  --taxonrectangle-layout
                        Activate taxonrectangle-layout which taxonomy of each node will be
                        display as rectangular blocks in aligned panel.
  --emapper-layout      Activate emapper_layout which display all the annotation from EggNOG-
                        mapper.
  --domain-layout       Activate domain_layout which display protein domain annotation in
                        sequence.
  --alignment-layout    Display Multiple Sequence Alignment layout in aligned panel.
  --profiling-layout PROFILING_LAYOUT [PROFILING_LAYOUT ...]
                        <prop1,prop2> names of properties which need to be convert to presence-
                        absence profiling matrix of each value
  --multi-profiling-layout MULTI_PROFILING_LAYOUT [MULTI_PROFILING_LAYOUT ...]
                        <prop1,prop2> names of properties containing values as list which need to
                        be convert to presence-absence profiling matrix
  --categorical-matrix-layout CATEGORICAL_MATRIX_LAYOUT [CATEGORICAL_MATRIX_LAYOUT ...]
                        <prop1,prop2> names which need to be plot as categorical_matrix_layout
                        for categorical values
  --numerical-matrix-layout NUMERICAL_MATRIX_LAYOUT [NUMERICAL_MATRIX_LAYOUT ...]
                        <prop1,prop2> names which need to be plot as numerical_matrix_layout for
                        numerical values

Visualizing output arguments:
  Visualizing output parameters

  --verbose             show detail on prompt when visualizing taget tree.
  --port PORT           run interactive session on custom port.[default: 5000]
  --plot PLOT           output as pdf
  --out-colordict       print color dictionary of each property
```

Here we use `examples/basic_example1/` and `examples/basic_example2/`
```
tree examples/basic_example1/
examples/basic_example1/
├── basic_example1_metadata2.tsv
├── basic_example1_miss.tsv
├── basic_example1_null.tsv
├── basic_example1.nw
├── basic_example1_metadata1.tsv
└── example1_demo.sh

tree examples/basic_example2
examples/basic_example2
├── diauxic.array
├── diauxic.nw
├── example2_demo.sh
├── FluA_H3_AA.fas
├── MCC_FluA_H3_Genotype.txt
└── MCC_FluA_H3.nw


head -5 examples/basic_example1/basic_example1_metadata1.tsv 
#name	sample1	sample2	sample3	sample4	sample5	random_type	*bool_type	bool_type2
Phy003I7ZJ_CHICK	0.05	0.12	0.86	0.01	0.69	medium	1	TRUE
Phy0054BO3_MELGA	0.64	0.67	0.51	0.29	0.14	medium	1	TRUE
Phy00508FR_NIPNI	0.89	0.38	0.97	0.49	0.26	low	1	FALSE
Phy004O1E0_APTFO	0.1	0.09	0.38	0.31	0.41	medium	0	TRUE


head -5 examples/basic_example1/basic_example1_metadata2.tsv
#name	abs_data	list_data
Phy003I7ZJ_CHICK	97	w,t,t
Phy0054BO3_MELGA	16	r,q,s
Phy00508FR_NIPNI	87	z,f,p
Phy004O1E0_APTFO	6	z,t,b

head -3 examples/basic_example2/diauxic.array
#NAMES	col1	col2	col3	col4	col5	col6	col7
YGR138C	-1.23	-0.81	1.79	0.78	-0.42	-0.69	0.58
YPR156C	-1.76	-0.94	1.16	0.36	0.41	-0.35	1.12

head -3 examples/basic_example2/MCC_FluA_H3_Genotype.txt 
#name	PB2	PB1	PA	HA	NP	NA	M	NS
A/Swine/Binh_Duong/03_10/2010	trig	trig	trig	HuH3N2	trig	HuH3N2	trig	trig
A/Swine/Binh_Duong/03_08/2010	trig	trig	trig	HuH3N2	trig	HuH3N2	trig	trig

## *annotate tree with more than one metadata tsv, seperated by `,`
treeprofiler annotate \
--tree examples/basic_example1/basic_example1.nw \
--input-type newick \
--metadata examples/basic_example1/basic_example1_metadata1.tsv examples/basic_example1/basic_example1_metadata2.tsv \
--bool-prop bool_type bool_type2 \
-o examples/basic_example1/

treeprofiler annotate \
--tree examples/basic_example2/diauxic.nw \
--input-type newick \
--metadata examples/basic_example2/diauxic.array \
--outdir examples/basic_example2/

treeprofiler annotate \
--tree examples/basic_example2/MCC_FluA_H3.nw \
--input-type newick \
--metadata examples/basic_example2/MCC_FluA_H3_Genotype.txt \
--outdir examples/basic_example2/
```

*if bool value is 1 or 0, treeprofiler will infer it as numerical data, hence we determine it as boolean value by using `--bool-prop` arguments

### Interactive visualization interface
TreeProfiler uses the new visualization framework implemented in [ETE 4.0](https://github.com/etetoolkit/ete/tree/ete4), which allows for the interactive exploration of huge phylogenies based on a context-based adaptive zooming strategy.

![treeprofiler interface](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/control_panel_page-0001.jpg?raw=true)

Overview of the TreeProfiler visualization interface. (A) The control panel allows users to customize visualization layout and features, and to perform text-based searches. (B) An annotated example tree, from `examples/basic_example1/` after `annotate`, is launched with a command `plot`. Support values (red) and branch distance (grey) are displayed on top of branches. The properties of one of the nodes are shown on the top. The minimap (bottom right) facilitates navigation. (C) The node editor panel provides access to node-specific actions, such as creating subtrees, collapsing, pruning, rooting and more. (D) Visualized properties by order are, categorical data `random_type` in `rectangle-layout`, numerical data `sample1`, `sample2`, `sample3` in `heatmap-layout` and `sample4`, `sample5` in `barplot-layout`, categorical data `random-type` in `profiling-layout` shown as presence-absence matrix. Layouts are shown with the order as input argument order from the command line. Names of properties are shown as titles on the top of each layout. (E) Legends each layout is shown on the top right corner with the same order as the layouts.
  

### Basic options of visualizing layouts
Selected properties of tree will be visualized at the aligned panel alongside with the tree, here is some basic parameters for layouts.
- `--column-width` column width of each property in layout. [default: 20]. 
- `--barplot-width` width of total scale of barplot layout.[default: 200]
- `--padding-x` customize horizontal column padding distance of each
layout.[default: 1]
- `--padding-y` customize vertical padding distance of each layout.[default: 0]


### Layouts for categorical data
Users can add the following flag to activate layouts for categorical data
```
--colorbranch-layout COLORBRANCH_LAYOUT
                        <prop1,prop2> names of properties where branches will be colored based on
                        different values.
--label-layout LABEL_LAYOUT
                      <prop1,prop2> names of properties where values will be displayed on the
                      aligned panel.
--rectangle-layout rectangle_layout
                      <prop1,prop2> names of properties where values will be label as rectangular
                      color block on the aligned panel.
--profiling-layout PROFILING_LAYOUT
                        <prop1,prop2> names of properties which need to be convert to presence-
                        absence profiling matrix of each value.
--categorical-matrix-layout CATEGORICAL_MATRIX_LAYOUT
                        <col1,col2> names, column index which need to be plot as categorical_matrix_layout for categorical columns.
```

example
```

## target column "random_type" in examples/basic_example1/basic_example1_metadata1.tsv
# List random_type feature as text in aligned panel using label_layout
treeprofiler plot --tree examples/basic_example1/basic_example1_annotated.ete --label-layout random_type 

# Label random_type feature on branch with different colors in aligned panel  using --colorbranch-layout
treeprofiler plot --tree examples/basic_example1/basic_example1_annotated.ete  --colorbranch-layout random_type 

# Label random_type feature with retangular block in aligned panel using --rectangle-layout
treeprofiler plot --tree examples/basic_example1/basic_example1_annotated.ete  --rectangle-layout random_type 

# Convert random_type feature into presence-absence profiling matrix using --profiling-layout
treeprofiler plot --tree examples/basic_example1/basic_example1_annotated.ete --profiling-layout random_type

# Label all feature with retangular block in aligned panel using --categorical-matrix-layout
treeprofiler plot --tree examples/basic_example2/MCC_FluA_H3_annotated.ete --categorical-matrix-layout PB2 PB1 PA HA NP NA M NS
```
![label_layout example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_label_layout.jpeg?raw=true)
`label-layout` displays the corresponding value of selected property
of each leaf and categorized with colors. 

![colorbranch_layout example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_colorbranch_layout.jpeg?raw=true)
`colorbranch-layout` categorize values of selected property by coloring the leaf nodes.

![rectangle_layout example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_rectangular_layout.jpeg?raw=true)
`rectangle-layout` categorizes values of selected property by displaying rectangular color block alongside the corresponing leaf node.

![profiling_layout example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_profiling_layout.png?raw=true)
`profiling-layout` convert categorical data of the selected property to presence-absence matrix.

![categorical_matrix_layout example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/cateforical_matrix_layout.png?raw=true)
`categorical-matrix-layout` convert multiple categorical properties into categorical matrix, each value will be represent in different color. In this example we use `MCC_FluA_H3.tree`, time-scaled phylogenetic tree of H3 influenza viruses inferred by BEAST using molecular clock model and `MCC_FluA_H3_Genotype.txt`, Genotype table of the H3 influenza viruses([Yu, Guangchuang et al. (2017)](https://doi.org/10.5061/dryad.v15v0)). 8 gene segments `PB2`,`PB1`,`PA`,`HA`,`NP`,`NA`,`M`,`NS` as properties, and virus strain `trig`, `pdm` and `HuH3N2` are categorized with different colors in the matrix.

### Layouts for boolean data
Users can add the following flag to activate layouts for Boolean data
```
---binary-layout BINARYLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be plot as binary_layout, label shown only positive value
--revbinary-layout REVBINARYLAYOUT
                        <col1,col2> names, column index or index range of columns which need to be plot as revbinary_layout, label shown only negative value
```                      

```
## target column "bool_type", "bool_type2" in examples/basic_example1/basic_example1_metadata1.tsv
# List postive bool_type feature in aligned panel using binary_layout
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
---binary-layout bool_type

# List negative bool_type feature in aligned panel using binary_layout
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--revbinary-layout bool_type2

```
![binary example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_binary_layout.jpeg?raw=true)
`binary_layout` displays postive value as colored block and negative value as grey block

![revbinary example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_revbinary_layout.jpeg?raw=true)
`binary_layout` displays negative value as colored block and postive value as grey block

*Boolean data can be also visualized by categorical layouts, such as 
```
# multiple columns seperated by ','
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--rectangle-layout bool_type
```

### Layouts for Numerical data
Users can add the following flag to activate layouts for Numerical data
```
--heatmap-layout HEATMAP_LAYOUT
                        <prop1,prop2> names of numerical properties which need to be read as heatmap_layout.
--barplot-layout BARPLOT_LAYOUT
                        <prop1,prop2> names of numericalproperties which need to be read as barplot_layouts.
--numerical-matrix-layout NUMERICAL_MATRIX_LAYOUT
                        <prop1,prop2> names which need to be plot as numerical_matrix_layout for numerical values.
```
```
## target column 'sample[1-5]' feature in examples/basic_example1/basic_example1_metadata1.tsv
# visualize sample1 feature in Barplot
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--barplot-layout sample1,sample2,sample3,sample4,sample5

# visualize sample1-sample5 in Heatmap
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--heatmap-layout sample1,sample2,sample3,sample4,sample5

# visualize col1-col7 in diauxic_annotated.nw with numerical profiling

treeprofiler plot \
--tree examples/basic_example2/diauxic_annotated.ete \
--numerical-matrix-layout col1,col2,col3,col4,col5,col6,col7
```
![barplot example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_barplot_layout.jpeg?raw=true)
`barplot-layout` display numerical data to barplot with scale, users can change the length of scale by using argument `--barplot-width [default: 200]`
![heatmap example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_heatmap_layout.png?raw=true)
`heatmap-layout` display numerical data to heatmap, which will automatically scale the minimum and maximum value from white to red.

![numerical profiling example](https://raw.githubusercontent.com/dengzq1234/treeprofiler_gallery/main/plot_numerical_profiling_layout.jpeg)

`numerical-profiling-layout` display multiple numerical data column to numerical data trix, which will automatically scale the minimum and maximum value from blue to red. Comparing to `heatmap-layout`, `numerical-profiling-layout` can afford more data columns with faster memory.  


### Layouts for list data
here we use example in `examples/basic_example1/basic_example1_metadata2.tsv`
```
# column `list_data` contain multiple elements value which can be process as list data in treeprofiler
head examples/basic_example1/basic_example1_metadata2.tsv
#name	abs_data	list_data
Phy003I7ZJ_CHICK	97	w,t,t
Phy0054BO3_MELGA	16	r,q,s
Phy00508FR_NIPNI	87	z,f,p
Phy004O1E0_APTFO	6	z,t,b
Phy004PA1B_ANAPL	72	z,r,p

## annotate tree first
treeprofiler annotate \
--tree examples/basic_example1/basic_example1.nw \
--input-type newick \
--metadata examples/basic_example1/basic_example1_metadata2.tsv \
-o examples/basic_example1/
```

Column `list_data` contain multiple elements value which can be process as list data in treeprofiler. Users can visualize list information using `--multi-profiling-layout`, which will . In this case, we highly reccomend users using `ete` using `--input-type` in order to resume datatype list in annotated tree.

```
# visualize using multi_profiling_layout
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--input-type ete \
--multi-profiling-layout list_data
```


### Layouts for multiple sequence alignment
In order to visualize multiple sequence alignment alongside with the tree, first we need to annotate alignment using `--alignment` in `annotate`. Then activate alignment layout by adding `--alignment-layout` 

```
# annotate
treeprofiler annotate \
--tree examples/basic_example2/MCC_FluA_H3.nw \
--input-type newick \
--alignment  ./examples/basic_example2/FluA_H3_AA.fas \
--outdir examples/basic_example2/

# visualize
treeprofiler plot \
--tree examples/basic_example2/MCC_FluA_H3_annotated.ete \
--alignment-layout
```
![alignment example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_alignment_layout.jpeg?raw=true)
`alignment-layout` displays multiple sequence alignments with a tree. Whole MSA sequences were visualized with a tree in rectangular layout. Sacle of sequence with position roadmark located at the top.

### Layouts for eggnog-mapper pfam annotations
if metadata is pfam annotations from eggnog-mapper, using `--emapper-pfam` to annotate domain information in target tree and **MUST** be with the alignment using `--alignment` to attach corresponding file.

Once tree is annotated, using `--domain-layout` to visualize it.

```
treeprofiler annotate \
--tree examples/pratical_example/emapper/7955.ENSDARP00000116736.nw \
--input-type newick \
--emapper-pfam examples/pratical_example/emapper/7955.out.emapper.pfam \
--alignment examples/pratical_example/emapper/7955.ENSDARP00000116736.aln.faa \
-o examples/pratical_example/emapper/

treeprofiler plot \
--tree examples/pratical_example/emapper/7955.ENSDARP00000116736_annotated.ete \
--domain-layout
```
![pfam example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_domain_layout.jpeg?raw=true)
`domain-layout` displays domain annotation with a tree. It requires sequence infomration `--alignment` in `annotate` step from MSA sequences to locate domain start and end position. 

### Layouts for eggnog-mapper smart annotations
if metadata is smart annotations from eggnog-mapper, using `--emapper-smart` to annotate domain information in target tree and must be with the alignment using `--alignment` to attach corresponding file.

Once tree is annotated, using `--domain-layout` to visualize it.

```
treeprofiler annotate \
--tree examples/pratical_example/emapper/7955.ENSDARP00000116736.nw \
--input-type newick \
--emapper-smart examples/pratical_example/emapper/7955.out.emapper.smart.out \
--alignment examples/pratical_example/emapper/7955.ENSDARP00000116736.aln.faa \
-o examples/pratical_example/emapper/

treeprofiler plot \
--tree examples/pratical_example/emapper/7955.ENSDARP00000116736_annotated.ete \
--domain-layout
```


### Layouts for eggnog-mapper annotations
If metadata is output from eggnog-mapper, using `--emapper-annotations` automatically parse all information as metadata. Program will parse data of all the columns from emapper output. Once tree is annotated, using `--emapper-layout` to visualize tree with all the metadata

```
seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs
```

```
treeprofiler annotate \
--tree examples/pratical_example/emapper/7955.ENSDARP00000116736.nw \
--input-type newick \
--emapper-annotations examples/pratical_example/emapper/7955.out.emapper.annotations \
-o examples/pratical_example/emapper/

treeprofiler plot \
--tree examples/pratical_example/emapper/7955.ENSDARP00000116736_annotated.ete \
--input-type ete \
--emapper-layout
```

[Check eggnogmapper example](#demo2-explore-eggnog-mapper-annotations-data-with-taxonomic-annotation)

### Visualizing annotated internal nodes
If internal nodes are annotated, TreeProfiler is also able to visualize annotated features automatically when layouts are activated

#### Internal nodes of categorical and boolean data
As internal nodes of categorical and boolean data are annotated as counter, for categorical data it generates a stacked bar of counter summary at the top of each internal node. And for boolean data, it generates a heatmap where represent positive(or negative) percentage of total data of each internal node.

**Categorical data**

Before collapsed
![text_uncollapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_text_uncollapsed.jpeg?raw=true)
After collapsed
![text_collapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_text_collapsed.jpeg?raw=true)
In this example, collapsed internal node shows a stacked bar which summarize categorical counter data of children nodes. 1/6 is red, 2/6 is blue and 3/6 is green.


**Boolean data**

Before collapsed
![bool_uncollapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_bool_uncollapsed.jpeg?raw=true)
After collapsed
![bool_collapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_bool_collapsed.jpeg?raw=true)
In this example, collapsed internal node shows a heatmap which represent the gradient level of positive/total ratio.

#### Internal nodes of numerical data
Internal nodes of numerical data are process descriptive statistic analysis by default, hence when users collapse any branch, barplot_layout or heatmap_layout will demonstrate representative value, `avg` by default. representative value can be changed by using `--internal-plot-measure`

**Numerical data**
```
# select max instead of avg as internal node ploting representative
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--heatmap-layout sample1,sample2,sample3,sample4,sample5 \
--internal-plot-measure max 
```
Before collapsed
![heatmap_uncollapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_heatmap_uncollapsed.jpeg?raw=true)

After collapsed
avg as itnernal plot measure
![heatmap_collapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_heatmap_collapsed.jpeg?raw=true)
max as itnernal plot measure
![heatmap_max_collapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_heatmap_collapsed_max.jpeg?raw=true)

### Layouts for Taxonomic data
If target tree was annotated with `--taxonomic-profiler` in previous `annotate` step successfully, now activate Taxonomic layout using `--taxonclade-layout` or `--taxonrectangle-layout` to visualize taxonomic classification. All rank levels will be generated separately and users can switch each of them on/off.

```
## Annotate
# GTDB
treeprofiler annotate \
--tree examples/taxonomy_example/gtdb/gtdb_example1.nw \
--metadata examples/taxonomy_example/gtdb/gtdb_example1.tsv \
--taxon-column name \
--taxonomic-profile \
--taxadb GTDB \
--outdir ./examples/taxonomy_example/gtdb/

# NCBI
treeprofiler annotate \
--tree examples/taxonomy_example/ncbi/spongilla_example.nw \
--metadata examples/taxonomy_example/ncbi/spongilla_example.tsv \
--taxonomic-profile \
--taxon-delimiter .  \
--taxa-field 0 \
--taxadb NCBI \
--outdir ./examples/taxonomy_example/ncbi/


## Visualize 
treeprofiler plot \
--tree ./examples/taxonomy_example/gtdb/gtdb_example1_annotated.nw \ --taxonrectangle-layout

treeprofiler plot \
--tree examples/taxonomy_example/ncbi/spongilla_example_annotated.nw \ --taxonclade-layout
```

![taxarect_collapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_taxarect.jpeg?raw=true)
`taxonrectangle-layout` shows taxonomic classification as rectangular block from root to leaf.

![taxaclade_collapsed](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/plot_taxaclade.jpeg?raw=true)
`taxonclade-layout` color associated clade of each category of each rank. 

## Conditional query in annotated tree
TreeProfiler allows users to perform conditional process based on different circumstances

- Conditional pruning, conditional pruning works both `annotate` and `plot` subcommand
    - `--pruned-by`, prune the annotated tree by conditions, and remove the branches or clades which don't fit the condition.
    - `--rank-limit`, prune the taxonomic annotated tree based on rank of classification.

- Conditional collapsing, conditional collapsing works in `plot` subcommand, allow users to collapsed tree internal nodes to clade under customized conditions
    - `--collapsed-by`, collapse tree branches whose nodes if the conditions, mainly on internal nodes
- Conditional highlight, conditional highlight works in `plot` subcommand, allow users to highlight tree nodes under customized conditions
    - `--highlighted-by`, select tree nodes which fit the conditions
 
### Query Syntax
#### Basic Query
All the conditional query shared the same syntax, a standard query consists the following 

```
--pruned-by|collapsed-by|highlighted-by "<left_value> <operator> <right_value>"
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
## annotate tree 
treeprofiler annotate \
--tree examples/basic_example1/basic_example1.nw \
--input-type newick \
--metadata examples/basic_example1/basic_example1_metadata1.tsv \
--bool-prop bool_type bool_type2 \
--counter-stat relative \
-o examples/basic_example1/ 

# Conditional pruning, prune leaf node whose name contain "FALPE"
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.nw \
--input-type newick \
--pruned-by "name contains FALPE"
```
Left panel is tree before prune, right panel is result after prune
![pruned](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/prune.png?raw=true)
```
# Conditional highlight
# select tree node whose name contains `FALPE` character
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.nw \
--input-type newick \
--highlighted-by "name contains FALPE"
```
![highlighted](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/highlighted.jpeg?raw=true)
```
# select tree node whose sample1 feature > 0.50, here we using ete format which can resume the datatype 
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--input-type ete \
--highlighted-by "sample1 > 0.50" \
--heatmap-layout sample1

# if use tree in newick format, we need to attach the prop2type file which can resume the datatype
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.nw \
--input-type newick \
--prop2type examples/basic_example1/basic_example1_prop2type.txt \
--highlighted-by "sample1 > 0.50" \
--heatmap-layout sample1
```
![highlighted_numeric](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/highlighted_numeric.png?raw=true)

#### Query in internal nodes
Query in internal nodes' properties is also available, in this case, `left_value` of query will be the internal node property, remember to add the proper suffixes such as `_avg`, `_max`,etc, for the numerical data or `_counter` for categorical and boolean data. 

Example
```
# select tree internal node where sample1_avg feature < 0.50
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--input-type ete \
--heatmap-layout sample1 \
--collapsed-by "sample1_avg < 0.50" 
```
![collapsed_numeric](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/collapsed_numeric.png?raw=true)
Syntax for internal node counter data
```
# collapse tree internal nodes, where `high` relative counter > 0.35 in random_type_counter property
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--input-type ete \
--rectangle-layout random_type \
--collapsed-by "random_type_counter:high > 0.35" \
--column-width 70
```
![collapsed_counter](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/collapsed_counter.png?raw=true)
#### AND and OR conditions
The syntax for the AND condition and OR condition in TreeProfiler is:

AND condition will be under one argument, syntax seperated by `,`, such as 
```
# select tree  node where sample1 feature > 0.50 AND sample2 < 0.2
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--input-type ete \
--heatmap-layout sample1 sample2 sample3 sample4 sample5 \
--highlighted-by "sample1>0.50,sample2<0.2" 
```
![highlighted_and](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/highlighted_and.png?raw=true)
OR condition will be used more than one arguments
```
# select tree node where sample1 feature > 0.50 OR sample2 < 0.2
treeprofiler plot \
--tree examples/basic_example1/basic_example1_annotated.ete \
--input-type ete \
--heatmap-layout sample1 sample2 sample3 sample4 sample5 \
--highlighted-by "sample1>0.50" \
--highlighted-by "sample2<0.2" 
```
![highlighted_or](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/highlighted_or.png?raw=true)
### conditional limit based on taxonomic level
Prune taxonomic annotated tree based on following taxonomic rank level,
`kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `subspecies` 
```
# Case in GTDB
# before pruning
treeprofiler plot \
--tree examples/taxonomy_example/gtdb/gtdb_example1_annotated.ete \
--input-type ete \
--taxonclade-layout 
```
before rank limit
![gtdb_before_rank](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/gtdb_taxa.png?raw=true)

```
# prune tree in visualization, rank limit to family level
treeprofiler plot \
--tree examples/taxonomy_example/gtdb/gtdb_example1_annotated.ete \
--input-type ete \
--rank-limit class \
--taxonclade-layout  
```
After rank_limit
![gtdb_class](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/gtdb_taxa_rank_class.png?raw=true)
As you see, class branches of target gtdb tree are all pruned and only left the internal branches which rank as class.  

## Demo1 Explore GTDB taxonomic tree with metadata and habitat information of progenome3
To illustrate the easiness and flexibility of TreeProfiler, we use it to annotate and visualize the version 202 of the GTDB prokaryotic phylogeny, which represents a species tree with 60,000 representative bacterial and archaeal species in [here](https://data.gtdb.ecogenomic.org/releases/release202/). GTDB provides the tree in plain newick format and massive datatable with various associated to such species. Apart from the metadata provided by the GTDB, here we also include annotations of genomes and species clusters to habitats from proGenomes3([Fullam et al. 2023](https://progenomes.embl.de/)). 

Example can be found in directory `./examples/pratical/gtdb_r202/`. We already prepared the gtdb v202 taxonomic tree `gtdbv202.nw` by merging Bacteria and Archaea trees, detailed steps are included in `merge_gtdbtree.py`. Based on the difference of computational capacity, complete steps and pipeline can be found in `gtdbv202full_demo.sh` and `gtdbv202lite_demo.sh`. 

1) A glance of habitat information of progenome3
```
cd examples/pratical/gtdb_r202/

zcat progenome3.tar.gz|head -n 4
progenome3.tsv0000664000175000017500000343216414447266674012351 0ustar  dengdenggtdb_genome_representative	aquatic_habitat	host_associated	soil_habitat
RS_GCF_004210275.1	f	t	f
GB_GCA_014116815.1			
RS_GCF_000730245.1	f	t	f

```

2) Download metadata of archaea and bacteria from gtdb
```
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_metadata_r202.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz
```

3) GTDB metadata annotation

Considering the size of GTDB metadata and phylogeney, here we provide two pipelines for user to choose base on their computational resources. 

- GTDB partial annotation (lightweight), which we will extract only a few of columns from metadata for annotation

```
# Extract genome_size, protein_count, gc_percentage, ncbi_assembly_level, ncbi_genome_category columns from GTDB metadata
tar -xf ar122_metadata_r202.tar.gz -O | cut -f1,14,89,13,46,56 > ar122_metadata_r202_lite.tsv
tar -xf bac120_metadata_r202.tar.gz -O | cut -f1,14,89,13,46,56 > bac120_metadata_r202_lite.tsv

# start annotation
treeprofiler annotate \
--tree gtdbv202.nw \
--input-type newick \
--metadata \
ar122_metadata_r202_lite.tsv bac120_metadata_r202_lite.tsv progenome3.tar.gz \
--taxonomic-profile \
--taxadb GTDB \
-o ./
```  


- GTDB full annotation, which requires **>6G disk space and >15G RAM memory**. 
```
# Annotate metadatas to taxonomic tree(this step may take a few minutes)
treeprofiler annotate \
--tree gtdbv202.nw \
--input-type newick \
--metadata \
ar122_metadata_r202.tar.gz bac120_metadata_r202.tar.gz progenome3.tar.gz \
--taxonomic-profile \
--taxadb GTDB \
-o ./

```

4) Visualizing annotated GTDB tree with GTDB metadata, which are 
- `genome_size`
- `protein_count`
- `gc_percentage` 
- `ncbi_assembly_level` 
- `ncbi_genome_category` 

  and progenome3 habitat information 
- `aquatic_habitat` 
- `host_associated` 
- `soil_habitat`

```
treeprofiler plot \
--tree gtdbv202_annotated.ete \
--input-type ete \
--barplot-layout genome_size protein_count \
--heatmap-layout gc_percentage \
---binary-layout aquatic_habitat host_associated soil_habitat \
--rectangle-layout ncbi_assembly_level ncbi_genome_category \
--taxonclade-layout \
--column-width 70
```

![gtdbv202_general](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/gtdb_v4.png?raw=true)
Here we show the GTDB v202 taxonomy tree (bacteria+archaea, 47894 leaves) in rectangular tree layout, with selected annotated properties which are displayed by order in aligned panel. Numerical data `genome_size` and `protein_count` are visualized as barplot, `gc_percentage` is shown as heatmap. Habitat information of progenome3, `aquatic_habitat`, `host_associated` and `soil_habitat` are shown as binary layout. Two categorical data `ncbi_assembly_level` and `ncbi_genome_category` are visualized as rectangular layout. In order to improve memory effiency, tree has default collapse level (10) hence multiple leaf nodes are collapsed as default, if nodes are collapsed, aligned layouts represented corresponding values of each property of annotated internal nodes. In this level, `taxonclade-layout` of the highest classification `kingdom` is activated, which demonstrate `bacteria` in salmon, `archaea` in blue. 

![gtdbv202_closeup](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/gtdb_v5.png?raw=true)
Once zoom into smaller view in tree, collapse level reduces automatically (or manually) to 1, then leaf nodes are dynamically displayed and rendered. Therefore associated layouts are shown as represending values of annotated leaves. `taxonclade-layout` colored leaf nodes in `specie` rank level.

![progenome3 example](https://raw.githubusercontent.com/dengzq1234/treeprofiler_gallery/main/progenome_example.jpeg)
Annotated tree and layouts can be shown as circular tree layout.

## Demo2 Explore large NifH gene tree with functional and taxonomic information
Here we analyzed the nitrogenase iron protein NifH gene family across bacteria from EggNOG6 with EggNOG-mapper, a tool for functional annotation based on precomputed orthology assignments. TreeProfiler provides options which allows users to directly map EggNOG-mapper outputs including functional annotations and pfam/smart domain predictions. Hence are then able to map these functional annotations to their respective phylogenetic gene trees and them with the evolutionary history, tracing from leaf to root level.

Map emapper annotation, pfam annotation and taxonomic annotation to target tree 
```
treeprofiler annotate \
--tree  examples/pratical_example/emapper/nifH.tree \
--emapper_annotation examples/pratical_example/emapper/nifH.out.emapper.annotations \
--emapper-pfam examples/pratical_example/emapper/nifH.out.emapper.pfam \
--alignment examples/pratical_example/emapper/nifH.faa.aln \
--taxonomic-profile \
--taxadb NCBI \
--taxon-delimiter . \
--taxa-field 0 \
-o examples/pratical_example/emapper/
```

Visualize annotations of emapper, pfam domain and ncbi taxonomy
```
treeprofiler plot \
--tree examples/pratical_example/emapper/nifH_annotated.ete \
--input-type ete \
--emapper-layout \
--domain-layout \
--taxonclade-layout \
--column-width 70
```


![emapper example](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/emapper_nifh_v1.png?raw=true)
visualization of categorical data `seed_orthologs`, `max_annot_lvl`, `COG_category`, `Description`, `Preferred_name`, and numerical data `score`

![emapper example2](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/emapper_nifh_v2.png?raw=true)
visualization of `KEGG_Pathway` profiling


![emapper example3](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/emapper_nifh_v3.png?raw=true)
visualization of `KEGG_ko` profiling

![emapper example4](https://github.com/dengzq1234/treeprofiler_gallery/blob/main/emapper_nifh_v4.png?raw=true)
visualization of `domain` annotation
