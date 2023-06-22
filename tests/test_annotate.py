
import sys
import os
from io import StringIO

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

#from collections import namedtuple
from tempfile import NamedTemporaryFile

from ete4 import Tree
import tree_annotate
import utils

def test_annotate_01():
    # basic annotate categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse('(a);', parser='newick')

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\na\tapple')
        f_annotation.flush()
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree = '(a:1[&&NHX:col1=apple]);'

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree

def test_annotate_02():
    # internal_nodes annotation categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree_no_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);'
    expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root:0[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'


    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    #print(test_tree_annotated.write(properties=[], format=1,format_root_node=True))
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1,format_root_node=True) == expected_tree_with_root

def test_annotate_03():
    # basic annotate numerical data
    # load tree
    test_tree = tree_annotate.ete4_parse('(a);', parser='newick')

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\na\t2')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree = '(a:1[&&NHX:col1=2.0]);'


    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree

def test_annotate_04():
    # internal_nodes annotation numerical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree_no_root = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_avg=3.5:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_sum=7.0])Internal_2:0.5[&&NHX:col1_avg=3.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_sum=9.0]);'
    expected_tree_with_root = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_avg=3.5:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_sum=7.0])Internal_2:0.5[&&NHX:col1_avg=3.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_sum=9.0])Root:0[&&NHX:col1_avg=2.5:col1_max=4.0:col1_min=1.0:col1_std=1.6666666666666667:col1_sum=10.0];'

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1, format_root_node=True) == expected_tree_with_root

def test_annotate_05():
    # test num_stat none and counter_stat none
    # internal_nodes annotation categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\talphabet_type\nA\t1\tvowel\nB\t2\tconsonant\nD\t3\tconsonant\nE\t4\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
    metadata_dict=metadata_dict, node_props=node_props, counter_stat='none', num_stat='none',
    columns=columns, prop2type=prop2type)
    
    expected_tree_no_root = '(A:1[&&NHX:alphabet_type=vowel:col1=1.0],(B:1[&&NHX:alphabet_type=consonant:col1=2.0],(E:1[&&NHX:alphabet_type=vowel:col1=4.0],D:1[&&NHX:alphabet_type=consonant:col1=3.0])Internal_1:0.5)Internal_2:0.5);'
    expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel:col1=1.0],(B:1[&&NHX:alphabet_type=consonant:col1=2.0],(E:1[&&NHX:alphabet_type=vowel:col1=4.0],D:1[&&NHX:alphabet_type=consonant:col1=3.0])Internal_1:0.5)Internal_2:0.5)Root:0;'

    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1,format_root_node=True) == expected_tree_with_root

def test_annotate_06():
    # assign internal node name
    # internal_nodes annotation categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);", parser='newick')

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree_no_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])N4:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])N5:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);'
    expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])N4:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])N5:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root:0[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    #print(test_tree_annotated.write(properties=[], format=1,format_root_node=True))
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1,format_root_node=True) == expected_tree_with_root

def test_annotate_07():
    # internal_nodes annotation boolean data
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);", parser='newick')

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tbool_type\nA\tTrue\nB\tFalse\nD\tTrue\nE\tFalse\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree_no_root = '(A:1[&&NHX:bool_type=True],(B:1[&&NHX:bool_type=False],(E:1[&&NHX:bool_type=False],D:1[&&NHX:bool_type=True])N4:0.5[&&NHX:bool_type_counter=False--1||True--1])N5:0.5[&&NHX:bool_type_counter=False--2||True--1]);'
    expected_tree_with_root = '(A:1[&&NHX:bool_type=True],(B:1[&&NHX:bool_type=False],(E:1[&&NHX:bool_type=False],D:1[&&NHX:bool_type=True])N4:0.5[&&NHX:bool_type_counter=False--1||True--1])N5:0.5[&&NHX:bool_type_counter=False--2||True--1])Root:0[&&NHX:bool_type_counter=False--2||True--2];'
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1,format_root_node=True) == expected_tree_with_root

def test_annotate_08():
    # internal_nodes annotation list data
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);", parser='newick')
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tlist_data\nA\ta,b,c\nB\tc,d\nD\ta,c,d,e\nE\te,d,b\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    expected_tree_no_root = '(A:1[&&NHX:list_data=a|b|c],(B:1[&&NHX:list_data=c|d],(E:1[&&NHX:list_data=e|d|b],D:1[&&NHX:list_data=a|c|d|e])N4:0.5[&&NHX:list_data_counter=a--1||b--1||c--1||d--2||e--2])N5:0.5[&&NHX:list_data_counter=a--1||b--1||c--2||d--3||e--2]);'
    expected_tree_with_root = '(A:1[&&NHX:list_data=a|b|c],(B:1[&&NHX:list_data=c|d],(E:1[&&NHX:list_data=e|d|b],D:1[&&NHX:list_data=a|c|d|e])N4:0.5[&&NHX:list_data_counter=a--1||b--1||c--1||d--2||e--2])N5:0.5[&&NHX:list_data_counter=a--1||b--1||c--2||d--3||e--2])Root:0[&&NHX:list_data_counter=a--2||b--2||c--3||d--3||e--2];'

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1, format_root_node=True) == expected_tree_with_root

def test_annotate_taxnomic_NCBI_01():
    # taxid in the leaf name
    sp_delimiter = ''
    sp_field = ''
    test_tree = tree_annotate.ete4_parse("((9598, 9606), 10090);", parser='newick')
    test_tree_annotated, rank2values = tree_annotate.annotate_taxa(test_tree, db='NCBI', taxid_attr="name", sp_delimiter=sp_delimiter, sp_field=sp_field)
    expected_tree = "((9598:1[&&NHX:common_name=Pan troglodytes:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes:rank=species:sci_name=Pan troglodytes:taxid=9598],9606:1[&&NHX:common_name=Homo sapiens:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens:rank=species:sci_name=Homo sapiens:taxid=9606])Homininae:1[&&NHX:common_name=:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae:rank=subfamily:sci_name=Homininae:taxid=207598],10090:1[&&NHX:common_name=Mus musculus:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus:rank=species:sci_name=Mus musculus:taxid=10090]);"
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree

def test_annotate_taxnomic_NCBI_02():
    # taxid in the leaf name
    test_tree = tree_annotate.ete4_parse("((9598.abc, 9606.nca), 10090.abd);", parser='newick')
    test_tree_annotated, rank2values = tree_annotate.annotate_taxa(test_tree, db='NCBI', taxid_attr="name", sp_delimiter='.', sp_field=0)

    expected_tree = "((9598.abc:1[&&NHX:common_name=Pan troglodytes:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes:rank=species:sci_name=Pan troglodytes:taxid=9598],9606.nca:1[&&NHX:common_name=Homo sapiens:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens:rank=species:sci_name=Homo sapiens:taxid=9606])Homininae:1[&&NHX:common_name=:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae:rank=subfamily:sci_name=Homininae:taxid=207598],10090.abd:1[&&NHX:common_name=Mus musculus:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus:rank=species:sci_name=Mus musculus:taxid=10090]);"
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree

def test_annotate_taxnomic_NCBI_03():
    # ncbi_id in the other column 
    test_tree = tree_annotate.ete4_parse("(A:1,((B:1,D:1):1,E:1):1);", parser='newick')

    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tncbi_id\nA\t7707\nB\t9606\nD\t9598\nE\t10090\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type, taxonomic_profile=True, taxadb='NCBI', taxon_column='ncbi_id', taxon_delimiter='', taxa_field=0)

    expected_tree_no_root = '(A:1[&&NHX:name=A:dist=1.0:support=1.0:taxid=7707:ncbi_id=7707:sci_name=Dendrochirotida:rank=order:lineage=1|131567|2759|33154|33208|6072|33213|33511|7586|133551|7624|7705|7706|7707:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Echinodermata|Eleutherozoa|Echinozoa|Holothuroidea|Dendrochirotacea|Dendrochirotida],((B:1[&&NHX:name=B:dist=1.0:support=1.0:taxid=9606:ncbi_id=9606:sci_name=Homo sapiens:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens],D:1[&&NHX:name=D:dist=1.0:support=1.0:taxid=9598:ncbi_id=9598:sci_name=Pan troglodytes:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes])Homininae:1[&&NHX:name=Homininae:dist=1.0:support=1.0:taxid=207598:sci_name=Homininae:rank=subfamily:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae],E:1[&&NHX:name=E:dist=1.0:support=1.0:taxid=10090:ncbi_id=10090:sci_name=Mus musculus:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus])Euarchontoglires:1[&&NHX:name=Euarchontoglires:dist=1.0:support=1.0:taxid=314146:sci_name=Euarchontoglires:rank=superorder:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires]);'
    expected_tree_with_root = '(A:1[&&NHX:name=A:dist=1.0:support=1.0:taxid=7707:ncbi_id=7707:sci_name=Dendrochirotida:rank=order:lineage=1|131567|2759|33154|33208|6072|33213|33511|7586|133551|7624|7705|7706|7707:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Echinodermata|Eleutherozoa|Echinozoa|Holothuroidea|Dendrochirotacea|Dendrochirotida],((B:1[&&NHX:name=B:dist=1.0:support=1.0:taxid=9606:ncbi_id=9606:sci_name=Homo sapiens:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens],D:1[&&NHX:name=D:dist=1.0:support=1.0:taxid=9598:ncbi_id=9598:sci_name=Pan troglodytes:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes])Homininae:1[&&NHX:name=Homininae:dist=1.0:support=1.0:taxid=207598:sci_name=Homininae:rank=subfamily:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae],E:1[&&NHX:name=E:dist=1.0:support=1.0:taxid=10090:ncbi_id=10090:sci_name=Mus musculus:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus])Euarchontoglires:1[&&NHX:name=Euarchontoglires:dist=1.0:support=1.0:taxid=314146:sci_name=Euarchontoglires:rank=superorder:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires])Deuterostomia:0[&&NHX:name=Deuterostomia:dist=0.0:support=1.0:taxid=33511:sci_name=Deuterostomia:rank=clade:lineage=1|131567|2759|33154|33208|6072|33213|33511:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia];'
    
    show_properties = ['name','dist','support', 'taxid', 'ncbi_id','sci_name','rank','lineage', 'named_lineage']
    assert test_tree_annotated.write(properties=show_properties, format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=show_properties, format=1, format_root_node=True) == expected_tree_with_root

def test_annotate_taxnomic_NCBI_04():
    # ncbi_id in the other column 
    test_tree = tree_annotate.ete4_parse("(A:1,((B:1,D:1):1,E:1):1);", parser='newick')

    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tncbi_id\nA\t7707.abc\nB\t9606.acs\nD\t9598.asd\nE\t10090.sds\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type, taxonomic_profile=True, taxadb='NCBI', taxon_column='ncbi_id', taxon_delimiter='.', taxa_field=0)

    expected_tree_no_root = '(A:1[&&NHX:name=A:dist=1.0:support=1.0:taxid=7707:ncbi_id=7707:sci_name=Dendrochirotida:rank=order:lineage=1|131567|2759|33154|33208|6072|33213|33511|7586|133551|7624|7705|7706|7707:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Echinodermata|Eleutherozoa|Echinozoa|Holothuroidea|Dendrochirotacea|Dendrochirotida],((B:1[&&NHX:name=B:dist=1.0:support=1.0:taxid=9606:ncbi_id=9606:sci_name=Homo sapiens:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens],D:1[&&NHX:name=D:dist=1.0:support=1.0:taxid=9598:ncbi_id=9598:sci_name=Pan troglodytes:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes])Homininae:1[&&NHX:name=Homininae:dist=1.0:support=1.0:taxid=207598:sci_name=Homininae:rank=subfamily:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae],E:1[&&NHX:name=E:dist=1.0:support=1.0:taxid=10090:ncbi_id=10090:sci_name=Mus musculus:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus])Euarchontoglires:1[&&NHX:name=Euarchontoglires:dist=1.0:support=1.0:taxid=314146:sci_name=Euarchontoglires:rank=superorder:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires]);'
    expected_tree_with_root = '(A:1[&&NHX:name=A:dist=1.0:support=1.0:taxid=7707:ncbi_id=7707:sci_name=Dendrochirotida:rank=order:lineage=1|131567|2759|33154|33208|6072|33213|33511|7586|133551|7624|7705|7706|7707:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Echinodermata|Eleutherozoa|Echinozoa|Holothuroidea|Dendrochirotacea|Dendrochirotida],((B:1[&&NHX:name=B:dist=1.0:support=1.0:taxid=9606:ncbi_id=9606:sci_name=Homo sapiens:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens],D:1[&&NHX:name=D:dist=1.0:support=1.0:taxid=9598:ncbi_id=9598:sci_name=Pan troglodytes:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes])Homininae:1[&&NHX:name=Homininae:dist=1.0:support=1.0:taxid=207598:sci_name=Homininae:rank=subfamily:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae],E:1[&&NHX:name=E:dist=1.0:support=1.0:taxid=10090:ncbi_id=10090:sci_name=Mus musculus:rank=species:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus])Euarchontoglires:1[&&NHX:name=Euarchontoglires:dist=1.0:support=1.0:taxid=314146:sci_name=Euarchontoglires:rank=superorder:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires])Deuterostomia:0[&&NHX:name=Deuterostomia:dist=0.0:support=1.0:taxid=33511:sci_name=Deuterostomia:rank=clade:lineage=1|131567|2759|33154|33208|6072|33213|33511:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia];'
    
    show_properties = ['name','dist','support', 'taxid', 'ncbi_id','sci_name','rank','lineage', 'named_lineage']
    assert test_tree_annotated.write(properties=show_properties, format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=show_properties, format=1, format_root_node=True) == expected_tree_with_root

def test_annotate_taxnomic_GTDB_01():
    # taxid in the leaf name
    # (GB_GCA_011358815.1@sample1:1,(RS_GCF_000019605.1@sample2:1,(RS_GCF_003948265.1@sample3:1,GB_GCA_003344655.1@sample4:1)1:0.5)1:0.5);
    test_tree = tree_annotate.ete4_parse("(GB_GCA_011358815.1:1,(RS_GCF_000019605.1:1,(RS_GCF_003948265.1:1,GB_GCA_003344655.1:1):0.5):0.5);", parser='newick')
    test_tree_annotated, rank2values = tree_annotate.annotate_taxa(test_tree, db='GTDB', taxid_attr="name", sp_delimiter='', sp_field=0)
    
    expected_tree_no_root = '(GB_GCA_011358815.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|2176|2177:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|GB_GCA_011358815.1:rank=subspecies:sci_name=s__Korarchaeum cryptofilum:taxid=GB_GCA_011358815.1],(RS_GCF_000019605.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|2176|2178:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_000019605.1:rank=subspecies:sci_name=s__Korarchaeum cryptofilum:taxid=RS_GCF_000019605.1],(RS_GCF_003948265.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|2176|2179:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_003948265.1:rank=subspecies:sci_name=s__Korarchaeum cryptofilum:taxid=RS_GCF_003948265.1],GB_GCA_003344655.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|5196|5197:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum sp003344655|GB_GCA_003344655.1:rank=subspecies:sci_name=s__Korarchaeum sp003344655:taxid=GB_GCA_003344655.1])g__Korarchaeum:0.5[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum:rank=genus:sci_name=g__Korarchaeum:taxid=g__Korarchaeum])g__Korarchaeum:0.5[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum:rank=genus:sci_name=g__Korarchaeum:taxid=g__Korarchaeum]);'
    expected_tree_with_root = '(GB_GCA_011358815.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|2176|2177:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|GB_GCA_011358815.1:rank=subspecies:sci_name=s__Korarchaeum cryptofilum:taxid=GB_GCA_011358815.1],(RS_GCF_000019605.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|2176|2178:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_000019605.1:rank=subspecies:sci_name=s__Korarchaeum cryptofilum:taxid=RS_GCF_000019605.1],(RS_GCF_003948265.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|2176|2179:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_003948265.1:rank=subspecies:sci_name=s__Korarchaeum cryptofilum:taxid=RS_GCF_003948265.1],GB_GCA_003344655.1:1[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175|5196|5197:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum sp003344655|GB_GCA_003344655.1:rank=subspecies:sci_name=s__Korarchaeum sp003344655:taxid=GB_GCA_003344655.1])g__Korarchaeum:0.5[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum:rank=genus:sci_name=g__Korarchaeum:taxid=g__Korarchaeum])g__Korarchaeum:0.5[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum:rank=genus:sci_name=g__Korarchaeum:taxid=g__Korarchaeum])g__Korarchaeum:0[&&NHX:common_name=:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum:rank=genus:sci_name=g__Korarchaeum:taxid=g__Korarchaeum];'
    
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1, format_root_node=True) == expected_tree_with_root

def test_annotate_taxnomic_GTDB_03():
    # taxid in the other column
    # (GB_GCA_011358815.1:1,(RS_GCF_000019605.1:1,(RS_GCF_003948265.1:1,GB_GCA_003344655.1:1):0.5):0.5);
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(C:1,D:1):0.5):0.5);", parser='newick')
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tgtdb_taxid\nA\tGB_GCA_011358815.1\nB\tRS_GCF_000019605.1\nC\tRS_GCF_003948265.1\nD\tGB_GCA_003344655.1\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type, taxonomic_profile=True, taxadb='GTDB', taxon_column='gtdb_taxid', taxon_delimiter='', taxa_field=0)

    show_properties = ['name','dist','support', 'taxid', 'gtdb_taxid','sci_name','rank','lineage', 'named_lineage']

    expected_tree_no_root = '(A:1[&&NHX:name=A:dist=1.0:support=1.0:taxid=GB_GCA_011358815.1:gtdb_taxid=GB_GCA_011358815.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2177:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|GB_GCA_011358815.1],(B:1[&&NHX:name=B:dist=1.0:support=1.0:taxid=RS_GCF_000019605.1:gtdb_taxid=RS_GCF_000019605.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2178:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_000019605.1],(C:1[&&NHX:name=C:dist=1.0:support=1.0:taxid=RS_GCF_003948265.1:gtdb_taxid=RS_GCF_003948265.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2179:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_003948265.1],D:1[&&NHX:name=D:dist=1.0:support=1.0:taxid=GB_GCA_003344655.1:gtdb_taxid=GB_GCA_003344655.1:sci_name=s__Korarchaeum sp003344655:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|5196|5197:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum sp003344655|GB_GCA_003344655.1])g__Korarchaeum:0.5[&&NHX:name=g__Korarchaeum:dist=0.5:support=1.0:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum])g__Korarchaeum:0.5[&&NHX:name=g__Korarchaeum:dist=0.5:support=1.0:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum]);'
    expected_tree_with_root = '(A:1[&&NHX:name=A:dist=1.0:support=1.0:taxid=GB_GCA_011358815.1:gtdb_taxid=GB_GCA_011358815.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2177:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|GB_GCA_011358815.1],(B:1[&&NHX:name=B:dist=1.0:support=1.0:taxid=RS_GCF_000019605.1:gtdb_taxid=RS_GCF_000019605.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2178:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_000019605.1],(C:1[&&NHX:name=C:dist=1.0:support=1.0:taxid=RS_GCF_003948265.1:gtdb_taxid=RS_GCF_003948265.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2179:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_003948265.1],D:1[&&NHX:name=D:dist=1.0:support=1.0:taxid=GB_GCA_003344655.1:gtdb_taxid=GB_GCA_003344655.1:sci_name=s__Korarchaeum sp003344655:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|5196|5197:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum sp003344655|GB_GCA_003344655.1])g__Korarchaeum:0.5[&&NHX:name=g__Korarchaeum:dist=0.5:support=1.0:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum])g__Korarchaeum:0.5[&&NHX:name=g__Korarchaeum:dist=0.5:support=1.0:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum])g__Korarchaeum:0[&&NHX:name=g__Korarchaeum:dist=0.0:support=1.0:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum];'
    
    assert test_tree_annotated.write(properties=show_properties, format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=show_properties, format=1, format_root_node=True) == expected_tree_with_root

def test_annotate_taxnomic_GTDB_04():

    # taxid in the other column
    # (GB_GCA_011358815.1:1,(RS_GCF_000019605.1:1,(RS_GCF_003948265.1:1,GB_GCA_003344655.1:1):0.5):0.5);
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(C:1,D:1):0.5):0.5);", parser='newick')
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tgtdb_taxid\nA\tGB_GCA_011358815.1@sample1\nB\tRS_GCF_000019605.1@sample2\nC\tRS_GCF_003948265.1@sample3\nD\tGB_GCA_003344655.1@sample4\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type, taxonomic_profile=True, taxadb='GTDB', taxon_column='gtdb_taxid', taxon_delimiter='@', taxa_field=0)

    
    show_properties = ['taxid', 'gtdb_taxid','sci_name','rank','lineage', 'named_lineage']
    
    expected_tree_no_root = '(A:1[&&NHX:taxid=GB_GCA_011358815.1:gtdb_taxid=GB_GCA_011358815.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2177:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|GB_GCA_011358815.1],(B:1[&&NHX:taxid=RS_GCF_000019605.1:gtdb_taxid=RS_GCF_000019605.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2178:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_000019605.1],(C:1[&&NHX:taxid=RS_GCF_003948265.1:gtdb_taxid=RS_GCF_003948265.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2179:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_003948265.1],D:1[&&NHX:taxid=GB_GCA_003344655.1:gtdb_taxid=GB_GCA_003344655.1:sci_name=s__Korarchaeum sp003344655:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|5196|5197:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum sp003344655|GB_GCA_003344655.1])g__Korarchaeum:0.5[&&NHX:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum])g__Korarchaeum:0.5[&&NHX:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum]);'
    expected_tree_with_root = '(A:1[&&NHX:taxid=GB_GCA_011358815.1:gtdb_taxid=GB_GCA_011358815.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2177:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|GB_GCA_011358815.1],(B:1[&&NHX:taxid=RS_GCF_000019605.1:gtdb_taxid=RS_GCF_000019605.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2178:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_000019605.1],(C:1[&&NHX:taxid=RS_GCF_003948265.1:gtdb_taxid=RS_GCF_003948265.1:sci_name=s__Korarchaeum cryptofilum:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|2176|2179:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum cryptofilum|RS_GCF_003948265.1],D:1[&&NHX:taxid=GB_GCA_003344655.1:gtdb_taxid=GB_GCA_003344655.1:sci_name=s__Korarchaeum sp003344655:rank=subspecies:lineage=1|2|79|2172|2173|2174|2175|5196|5197:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum|s__Korarchaeum sp003344655|GB_GCA_003344655.1])g__Korarchaeum:0.5[&&NHX:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum])g__Korarchaeum:0.5[&&NHX:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum])g__Korarchaeum:0[&&NHX:taxid=g__Korarchaeum:sci_name=g__Korarchaeum:rank=genus:lineage=1|2|79|2172|2173|2174|2175:named_lineage=root|d__Archaea|p__Thermoproteota|c__Korarchaeia|o__Korarchaeales|f__Korarchaeaceae|g__Korarchaeum];'
    
    assert test_tree_annotated.write(properties=show_properties, format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=show_properties, format=1, format_root_node=True) == expected_tree_with_root

def test_rank_limit():
    test_ncbitree = Tree("(((36324:1[&&NHX:rank=species:sci_name=Psolus fabricii:taxid=36324],42201:1[&&NHX:rank=species:sci_name=Psolus chitonoides:taxid=42201],864673:1[&&NHX:rank=species:sci_name=Psolus antarcticus complex sp. FM-2010:taxid=864673],864674:1[&&NHX:rank=species:sci_name=Psolus dubiosus complex sp. FM-2010:taxid=864674],868587:1[&&NHX:rank=species:sci_name=Psolus phantapus:taxid=868587],2014662:1[&&NHX:rank=species:sci_name=Psolus rufus:taxid=2014662],(1769970:1[&&NHX:rank=species:sci_name=Psolus sp. 11 RWR-2015:taxid=1769970],1769971:1[&&NHX:rank=species:sci_name=Psolus sp. 41 RWR-2015:taxid=1769971],1783523:1[&&NHX:rank=species:sci_name=Psolus sp. RWR-2015:taxid=1783523])2.6225e+06:1[&&NHX:rank=no rank:sci_name=unclassified Psolus:taxid=2622497])36323:1[&&NHX:rank=genus:sci_name=Psolus:taxid=36323],(864665:1[&&NHX:rank=species:sci_name=Psolidium gaini:taxid=864665],864672:1[&&NHX:rank=species:sci_name=Psolidium tenue complex sp. FM-2010:taxid=864672],1902842:1[&&NHX:rank=species:sci_name=Psolidium dorsipes:taxid=1902842],1966443:1[&&NHX:rank=species:sci_name=Psolidium whittakeri:taxid=1966443],1902934:1[&&NHX:rank=species:sci_name=Psolidium sp. AKM-2016:taxid=1902934])864664:1[&&NHX:rank=genus:sci_name=Psolidium:taxid=864664],1902863:1[&&NHX:rank=species:sci_name=Lissothuria nutriens:taxid=1902863])36322:1[&&NHX:rank=family:sci_name=Psolidae:taxid=36322],((28833:1[&&NHX:rank=species:sci_name=Cucumaria miniata:taxid=28833],36326:1[&&NHX:rank=species:sci_name=Cucumaria frondosa:taxid=36326],40245:1[&&NHX:rank=species:sci_name=Cucumaria echinata:taxid=40245],42202:1[&&NHX:rank=species:sci_name=Cucumaria vegae:taxid=42202],42203:1[&&NHX:rank=species:sci_name=Cucumaria pallida:taxid=42203],42204:1[&&NHX:rank=species:sci_name=Cucumaria piperata:taxid=42204],42388:1[&&NHX:rank=species:sci_name=Cucumaria pseudocurata:taxid=42388],122241:1[&&NHX:rank=species:sci_name=Cucumaria salma:taxid=122241],864625:1[&&NHX:rank=species:sci_name=Cucumaria dudexa:taxid=864625],864626:1[&&NHX:rank=species:sci_name=Cucumaria georgiana:taxid=864626],869198:1[&&NHX:rank=species:sci_name=Cucumaria cf. lubrica EAC-2010:taxid=869198],1898586:1[&&NHX:rank=species:sci_name=Cucumaria cf. lubrica KSL-2016:taxid=1898586],(1421346:1[&&NHX:rank=species:sci_name=Cucumaria sp. ZP0016:taxid=1421346],2653558:1[&&NHX:rank=species:sci_name=Cucumaria sp. 17_ECHINO_001_070:taxid=2653558],2653559:1[&&NHX:rank=species:sci_name=Cucumaria sp. 17_ECHINO_001_071:taxid=2653559])2.62272e+06:1[&&NHX:rank=no rank:sci_name=unclassified Cucumaria:taxid=2622715])28832:1[&&NHX:rank=genus:sci_name=Cucumaria:taxid=28832],(42205:1[&&NHX:rank=species:sci_name=Pseudocnus lubricus:taxid=42205],42212:1[&&NHX:rank=species:sci_name=Pseudocnus curatus:taxid=42212],42386:1[&&NHX:rank=species:sci_name=Pseudocnus californicus:taxid=42386],42387:1[&&NHX:rank=species:sci_name=Pseudocnus astigmatus:taxid=42387])42385:1[&&NHX:rank=genus:sci_name=Pseudocnus:taxid=42385],(206683:1[&&NHX:rank=species:sci_name=Aslia lefevrei:taxid=206683],1690249:1[&&NHX:rank=species:sci_name=Aslia forbesi:taxid=1690249],1902836:1[&&NHX:rank=species:sci_name=Aslia pygmaea:taxid=1902836])206682:1[&&NHX:rank=genus:sci_name=Aslia:taxid=206682],(206685:1[&&NHX:rank=species:sci_name=Thyonella gemmata:taxid=206685],2772383:1[&&NHX:rank=species:sci_name=Thyonella sp. USNM IZ 1446416:taxid=2772383])206684:1[&&NHX:rank=genus:sci_name=Thyonella:taxid=206684],396326:1[&&NHX:rank=genus:sci_name=Pentacta:taxid=396326],(861480:1[&&NHX:rank=species:sci_name=Cucumariidae sp. EAC-2010a:taxid=861480],1443624:1[&&NHX:rank=species:sci_name=Cucumariidae sp. RG-2014:taxid=1443624])861479:1[&&NHX:rank=no rank:sci_name=unclassified Cucumariidae:taxid=861479],(864646:1[&&NHX:rank=species:sci_name=Heterocucumis steineni:taxid=864646],1633604:1[&&NHX:rank=species:sci_name=Heterocucumis denticulata:taxid=1633604])864645:1[&&NHX:rank=genus:sci_name=Heterocucumis:taxid=864645],864662:1[&&NHX:rank=genus:sci_name=Psolidiella:taxid=864662],(864668:1[&&NHX:rank=species:sci_name=Staurocucumis liouvillei:taxid=864668],1207538:1[&&NHX:rank=species:sci_name=Staurocucumis turqueti:taxid=1207538],1633605:1[&&NHX:rank=species:sci_name=Staurocucumis krzysztofi:taxid=1633605],1902935:1[&&NHX:rank=species:sci_name=Staurocucumis cf. turqueti AKM-2016:taxid=1902935],1633607:1[&&NHX:rank=species:sci_name=Staurocucumis sp. GP-2015:taxid=1633607])864667:1[&&NHX:rank=genus:sci_name=Staurocucumis:taxid=864667],(1633609:1[&&NHX:rank=species:sci_name=Abyssocucumis abyssorum:taxid=1633609],1777472:1[&&NHX:rank=species:sci_name=Abyssocucumis albatrossi:taxid=1777472])1.63361e+06:1[&&NHX:rank=genus:sci_name=Abyssocucumis:taxid=1633608],(1633611:1[&&NHX:rank=species:sci_name=Cladodactyla crocea:taxid=1633611],1633612:1[&&NHX:rank=species:sci_name=Cladodactyla sicinski:taxid=1633612])1.63361e+06:1[&&NHX:rank=genus:sci_name=Cladodactyla:taxid=1633610],(1633615:1[&&NHX:rank=species:sci_name=Laevocnus laevigatus:taxid=1633615],1633618:1[&&NHX:rank=species:sci_name=Laevocnus perrieri:taxid=1633618],1633619:1[&&NHX:rank=species:sci_name=Laevocnus serratus:taxid=1633619])1.63361e+06:1[&&NHX:rank=genus:sci_name=Laevocnus:taxid=1633613],(1725248:1[&&NHX:rank=species:sci_name=Leptopentacta imbricata:taxid=1725248],2268661:1[&&NHX:rank=species:sci_name=Leptopentacta elongata:taxid=2268661])1.72525e+06:1[&&NHX:rank=genus:sci_name=Leptopentacta:taxid=1725247],(1902849:1[&&NHX:rank=species:sci_name=Colochirus robustus:taxid=1902849],1980634:1[&&NHX:rank=species:sci_name=Colochirus quadrangularis:taxid=1980634],1980635:1[&&NHX:rank=species:sci_name=Colochirus sp_1_GP:taxid=1980635])1.90285e+06:1[&&NHX:rank=genus:sci_name=Colochirus:taxid=1902848],(1633614:1[&&NHX:rank=species:sci_name=Pentactella katrinae:taxid=1633614],1633616:1[&&NHX:rank=species:sci_name=Pentactella leachmani:taxid=1633616],1902875:1[&&NHX:rank=species:sci_name=Pentactella leonina:taxid=1902875],1966442:1[&&NHX:rank=species:sci_name=Pentactella sp. AM-2017:taxid=1966442])1.90287e+06:1[&&NHX:rank=genus:sci_name=Pentactella:taxid=1902874],(1980638:1[&&NHX:rank=species:sci_name=Plesiocolochirus challengeri:taxid=1980638],1980639:1[&&NHX:rank=species:sci_name=Plesiocolochirus ignavus:taxid=1980639],1980640:1[&&NHX:rank=species:sci_name=Plesiocolochirus minaeus:taxid=1980640],1980641:1[&&NHX:rank=species:sci_name=Plesiocolochirus tessellarus:taxid=1980641],1980642:1[&&NHX:rank=species:sci_name=Plesiocolochirus sp_1_GP:taxid=1980642],1980643:1[&&NHX:rank=species:sci_name=Plesiocolochirus sp_2_GP:taxid=1980643])1.98064e+06:1[&&NHX:rank=genus:sci_name=Plesiocolochirus:taxid=1980637],(864629:1[&&NHX:rank=species:sci_name=Echinopsolus charcoti:taxid=864629],864631:1[&&NHX:rank=species:sci_name=Echinopsolus koehleri:taxid=864631],864663:1[&&NHX:rank=species:sci_name=Echinopsolus mollis:taxid=864663])2.77157e+06:1[&&NHX:rank=genus:sci_name=Echinopsolus:taxid=2771573],571821:1[&&NHX:rank=species:sci_name=Athyonidium chilensis:taxid=571821],864642:1[&&NHX:rank=species:sci_name=Cucamba psolidiformis:taxid=864642],864661:1[&&NHX:rank=species:sci_name=Psolicrux coatsi:taxid=864661],864670:1[&&NHX:rank=species:sci_name=Trachythyone bouvetensis:taxid=864670],1622217:1[&&NHX:rank=species:sci_name=Ekmania barthii:taxid=1622217],1690254:1[&&NHX:rank=species:sci_name=Trachasina crucifera:taxid=1690254],2268711:1[&&NHX:rank=species:sci_name=Panningia hyndmani:taxid=2268711],2528974:1[&&NHX:rank=species:sci_name=Pseudocolochirus violaceus:taxid=2528974],2576476:1[&&NHX:rank=species:sci_name=Neocucumis proteus:taxid=2576476],55632:1[&&NHX:rank=species:sci_name=Pseudocnella sykion:taxid=55632],1633606:1[&&NHX:rank=species:sci_name=Psolicucumis nocturna:taxid=1633606],2785017:1[&&NHX:rank=species:sci_name=Pawsonia saxicola:taxid=2785017],2785214:1[&&NHX:rank=species:sci_name=Cercodemas anceps:taxid=2785214],2850318:1[&&NHX:rank=species:sci_name=Paraleptopentacta elongata:taxid=2850318],1902940:1[&&NHX:rank=species:sci_name=Amphicyclus sp. AKM-2016:taxid=1902940])36325:1[&&NHX:rank=family:sci_name=Cucumariidae:taxid=36325],((42391:1[&&NHX:rank=species:sci_name=Pentamera lissoplaca:taxid=42391],861476:1[&&NHX:rank=species:sci_name=Pentamera calcigera:taxid=861476],861516:1[&&NHX:rank=species:sci_name=Pentamera cf. pseudocalcigera EAC-2010:taxid=861516],869196:1[&&NHX:rank=species:sci_name=Pentamera pediparva:taxid=869196],1382469:1[&&NHX:rank=species:sci_name=Pentamera rigida:taxid=1382469],(2172660:1[&&NHX:rank=species:sci_name=Pentamera sp. 3 ML-2018:taxid=2172660],2172661:1[&&NHX:rank=species:sci_name=Pentamera sp. 4 ML-2018:taxid=2172661])2.64079e+06:1[&&NHX:rank=no rank:sci_name=unclassified Pentamera:taxid=2640789])42390:1[&&NHX:rank=genus:sci_name=Pentamera:taxid=42390],(55631:1[&&NHX:rank=species:sci_name=Lipotrapeza vestiens:taxid=55631],1382468:1[&&NHX:rank=species:sci_name=Lipotrapeza eichleri:taxid=1382468])55630:1[&&NHX:rank=genus:sci_name=Lipotrapeza:taxid=55630],(869205:1[&&NHX:rank=species:sci_name=Thyonidium drummondii:taxid=869205],2172724:1[&&NHX:rank=species:sci_name=Thyonidium kurilensis:taxid=2172724],2268686:1[&&NHX:rank=species:sci_name=Thyonidium hyalinum:taxid=2268686],2268729:1[&&NHX:rank=species:sci_name=Thyonidium sp. Echin 6076V:taxid=2268729])869204:1[&&NHX:rank=genus:sci_name=Thyonidium:taxid=869204],(1238286:1[&&NHX:rank=species:sci_name=Phyrella cf. thyonoides FM-2012:taxid=1238286],1238287:1[&&NHX:rank=species:sci_name=Phyrella fragilis:taxid=1238287],1239311:1[&&NHX:rank=species:sci_name=Phyrella mookiei:taxid=1239311])1.23828e+06:1[&&NHX:rank=genus:sci_name=Phyrella:taxid=1238285],(1382473:1[&&NHX:rank=species:sci_name=Massinium magnum:taxid=1382473],2928469:1[&&NHX:rank=species:sci_name=Massinium toyoshiomaruae:taxid=2928469])1.38247e+06:1[&&NHX:rank=genus:sci_name=Massinium:taxid=1382472],(1382475:1[&&NHX:rank=species:sci_name=Phyllophorus brocki:taxid=1382475],1690250:1[&&NHX:rank=species:sci_name=Phyllophorus cebuensis:taxid=1690250],(1382678:1[&&NHX:rank=species:sci_name=Phyllophorus sp. UF_9620:taxid=1382678],1382679:1[&&NHX:rank=species:sci_name=Phyllophorus sp. UF_9621:taxid=1382679],1382680:1[&&NHX:rank=species:sci_name=Phyllophorus sp. UF_9622:taxid=1382680],1382681:1[&&NHX:rank=species:sci_name=Phyllophorus sp. UF_9624:taxid=1382681],1382682:1[&&NHX:rank=species:sci_name=Phyllophorus sp. WAM_Z11501:taxid=1382682],1382683:1[&&NHX:rank=species:sci_name=Phyllophorus sp. WAM_Z21137:taxid=1382683],1382684:1[&&NHX:rank=species:sci_name=Phyllophorus sp. WAM_Z29789:taxid=1382684],1382685:1[&&NHX:rank=species:sci_name=Phyllophorus sp. WAM_Z31837:taxid=1382685])2.62582e+06:1[&&NHX:rank=no rank:sci_name=unclassified Phyllophorus:taxid=2625815])1.38247e+06:1[&&NHX:rank=genus:sci_name=Phyllophorus:taxid=1382474],(1382478:1[&&NHX:rank=species:sci_name=Thyone flindersi:taxid=1382478],1382479:1[&&NHX:rank=species:sci_name=Thyone nigra:taxid=1382479],1382480:1[&&NHX:rank=species:sci_name=Thyone pedata:taxid=1382480],1933098:1[&&NHX:rank=species:sci_name=Thyone fusus:taxid=1933098],2172723:1[&&NHX:rank=species:sci_name=Thyone benti:taxid=2172723],(1902937:1[&&NHX:rank=species:sci_name=Thyone sp. AKM-2016:taxid=1902937],2268728:1[&&NHX:rank=species:sci_name=Thyone sp. Echin 6777V:taxid=2268728])2.6241e+06:1[&&NHX:rank=no rank:sci_name=unclassified Thyone:taxid=2624100],2928468:1[&&NHX:rank=species:sci_name=Thyone toyoshiomaruae:taxid=2928468],2928470:1[&&NHX:rank=species:sci_name=Thyone kyushuensis:taxid=2928470],2928471:1[&&NHX:rank=species:sci_name=Thyone liaoi:taxid=2928471])1.38248e+06:1[&&NHX:rank=genus:sci_name=Thyone:taxid=1382477],(1382476:1[&&NHX:rank=species:sci_name=Phyllophorella kohkutiensis:taxid=1382476],2810320:1[&&NHX:rank=species:sci_name=Phyllophorella liuwutiensis:taxid=2810320])2.77157e+06:1[&&NHX:rank=genus:sci_name=Phyllophorella:taxid=2771572],206687:1[&&NHX:rank=species:sci_name=Neopentadactyla mixta:taxid=206687],1382471:1[&&NHX:rank=species:sci_name=Hemithyone semperi:taxid=1382471],(1382687:1[&&NHX:rank=species:sci_name=Havelockia sp. NMV_F151829:taxid=1382687],1382688:1[&&NHX:rank=species:sci_name=Havelockia sp. NMV_F151830:taxid=1382688])2.62497e+06:1[&&NHX:rank=no rank:sci_name=unclassified Havelockia:taxid=2624973],(1382690:1[&&NHX:rank=species:sci_name=Neothyonidium sp. NMV_F150806:taxid=1382690],1382691:1[&&NHX:rank=species:sci_name=Neothyonidium sp. NMV_F151827:taxid=1382691])2.62759e+06:1[&&NHX:rank=no rank:sci_name=unclassified Neothyonidium:taxid=2627593],(1382693:1[&&NHX:rank=species:sci_name=Stolus sp. NMV_F151822:taxid=1382693],1382694:1[&&NHX:rank=species:sci_name=Stolus sp. UF_9494:taxid=1382694])2.62909e+06:1[&&NHX:rank=no rank:sci_name=unclassified Stolus:taxid=2629094],1690252:1[&&NHX:rank=species:sci_name=Thyonina sp. AB-2015:taxid=1690252])42389:1[&&NHX:rank=family:sci_name=Phyllophoridae:taxid=42389],((7710:1[&&NHX:rank=species:sci_name=Sclerodactyla briareus:taxid=7710],2080421:1[&&NHX:rank=species:sci_name=Sclerodactyla multipes:taxid=2080421])7709:1[&&NHX:rank=genus:sci_name=Sclerodactyla:taxid=7709],(42394:1[&&NHX:rank=species:sci_name=Eupentacta quinquesemita:taxid=42394],1774088:1[&&NHX:rank=species:sci_name=Eupentacta fraudatrix:taxid=1774088],2172693:1[&&NHX:rank=species:sci_name=Eupentacta pseudoquinquesemita:taxid=2172693],(2653564:1[&&NHX:rank=species:sci_name=Eupentacta sp. 17_FA_001:taxid=2653564],2653565:1[&&NHX:rank=species:sci_name=Eupentacta sp. 17_ECHINO_001_087:taxid=2653565],2653566:1[&&NHX:rank=species:sci_name=Eupentacta sp. 17_ECHINO_001_091:taxid=2653566],2653567:1[&&NHX:rank=species:sci_name=Eupentacta sp. 17_ECHINO_001_075:taxid=2653567],2857184:1[&&NHX:rank=species:sci_name=Eupentacta sp. USNM IZ 1503386:taxid=2857184],2857185:1[&&NHX:rank=species:sci_name=Eupentacta sp. USNM IZ 1503387:taxid=2857185],2857186:1[&&NHX:rank=species:sci_name=Eupentacta sp. USNM IZ 1503388:taxid=2857186],2857213:1[&&NHX:rank=species:sci_name=Eupentacta sp. USNM IZ 1523735:taxid=2857213],2857318:1[&&NHX:rank=species:sci_name=Eupentacta sp. USNM IZ 1552526:taxid=2857318])2.65356e+06:1[&&NHX:rank=no rank:sci_name=unclassified Eupentacta:taxid=2653563])42393:1[&&NHX:rank=genus:sci_name=Eupentacta:taxid=42393],(2172754:1[&&NHX:rank=species:sci_name=Pseudothyone levini:taxid=2172754],2268677:1[&&NHX:rank=species:sci_name=Pseudothyone raphanus:taxid=2268677])2.17275e+06:1[&&NHX:rank=genus:sci_name=Pseudothyone:taxid=2172753],206689:1[&&NHX:rank=species:sci_name=Afrocucumis africana:taxid=206689],1902871:1[&&NHX:rank=species:sci_name=Pachythyone rubra:taxid=1902871],1966448:1[&&NHX:rank=species:sci_name=Euthyonidiella huwi:taxid=1966448],2033684:1[&&NHX:rank=species:sci_name=Cladolabes schmeltzii:taxid=2033684])42392:1[&&NHX:rank=family:sci_name=Sclerodactylidae:taxid=42392],(818303:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAB7664:taxid=818303],818304:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAC9852:taxid=818304],818305:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAE3748:taxid=818305],857999:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN112-08:taxid=857999],858000:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN113-08:taxid=858000],858001:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN124-08:taxid=858001],858002:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN125-08:taxid=858002],858003:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN126-08:taxid=858003],858004:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN127-08:taxid=858004],858005:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN138-08:taxid=858005],974485:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAN4944:taxid=974485],2480094:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-07:taxid=2480094],2480095:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-08:taxid=2480095],2480096:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-09:taxid=2480096])734213:1[&&NHX:rank=no rank:sci_name=unclassified Dendrochirotida:taxid=734213],((864639:1[&&NHX:rank=species:sci_name=Crucella hystrix:taxid=864639],864640:1[&&NHX:rank=species:sci_name=Crucella scotiae:taxid=864640],1633603:1[&&NHX:rank=species:sci_name=Crucella susannae:taxid=1633603])864638:1[&&NHX:rank=genus:sci_name=Crucella:taxid=864638],864654:1[&&NHX:rank=species:sci_name=Paracucumis turricata:taxid=864654])864675:1[&&NHX:rank=family:sci_name=Paracucumidae:taxid=864675],((1902853:1[&&NHX:rank=species:sci_name=Echinocucumis hispida:taxid=1902853],2041650:1[&&NHX:rank=species:sci_name=Echinocucumis cf. hispida AM-2017:taxid=2041650])1.90285e+06:1[&&NHX:rank=genus:sci_name=Echinocucumis:taxid=1902852],1902938:1[&&NHX:rank=species:sci_name=Ypsilothuria cf. bitentaculata AKM-2016:taxid=1902938])1.90289e+06:1[&&NHX:rank=family:sci_name=Ypsilothuriidae:taxid=1902891],1902861:1[&&NHX:rank=species:sci_name=Heterothyone alba:taxid=1902861],1902877:1[&&NHX:rank=species:sci_name=Placothuria squamata:taxid=1902877]);", format=1, parser='newick') 
    #print(test_ncbitree.children[0].props)
    pruned_tree = tree_annotate.taxatree_prune(test_ncbitree, rank_limit='family')
    expected_tree_no_root = '(36322:1[&&NHX:rank=family:sci_name=Psolidae:taxid=36322],36325:1[&&NHX:rank=family:sci_name=Cucumariidae:taxid=36325],42389:1[&&NHX:rank=family:sci_name=Phyllophoridae:taxid=42389],42392:1[&&NHX:rank=family:sci_name=Sclerodactylidae:taxid=42392],(818303:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAB7664:taxid=818303],818304:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAC9852:taxid=818304],818305:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAE3748:taxid=818305],857999:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN112-08:taxid=857999],858000:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN113-08:taxid=858000],858001:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN124-08:taxid=858001],858002:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN125-08:taxid=858002],858003:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN126-08:taxid=858003],858004:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN127-08:taxid=858004],858005:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN138-08:taxid=858005],974485:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAN4944:taxid=974485],2480094:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-07:taxid=2480094],2480095:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-08:taxid=2480095],2480096:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-09:taxid=2480096])734213:1[&&NHX:rank=no rank:sci_name=unclassified Dendrochirotida:taxid=734213],864675:1[&&NHX:rank=family:sci_name=Paracucumidae:taxid=864675],1.90289e+06:1[&&NHX:rank=family:sci_name=Ypsilothuriidae:taxid=1902891],1902861:1[&&NHX:rank=species:sci_name=Heterothyone alba:taxid=1902861],1902877:1[&&NHX:rank=species:sci_name=Placothuria squamata:taxid=1902877]);'
    expected_tree_with_root = '(36322:1[&&NHX:rank=family:sci_name=Psolidae:taxid=36322],36325:1[&&NHX:rank=family:sci_name=Cucumariidae:taxid=36325],42389:1[&&NHX:rank=family:sci_name=Phyllophoridae:taxid=42389],42392:1[&&NHX:rank=family:sci_name=Sclerodactylidae:taxid=42392],(818303:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAB7664:taxid=818303],818304:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAC9852:taxid=818304],818305:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAE3748:taxid=818305],857999:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN112-08:taxid=857999],858000:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN113-08:taxid=858000],858001:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN124-08:taxid=858001],858002:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN125-08:taxid=858002],858003:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN126-08:taxid=858003],858004:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN127-08:taxid=858004],858005:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. ECNN138-08:taxid=858005],974485:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. BOLD_AAN4944:taxid=974485],2480094:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-07:taxid=2480094],2480095:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-08:taxid=2480095],2480096:1[&&NHX:rank=species:sci_name=Dendrochirotida sp. WMNH-INV-09:taxid=2480096])734213:1[&&NHX:rank=no rank:sci_name=unclassified Dendrochirotida:taxid=734213],864675:1[&&NHX:rank=family:sci_name=Paracucumidae:taxid=864675],1.90289e+06:1[&&NHX:rank=family:sci_name=Ypsilothuriidae:taxid=1902891],1902861:1[&&NHX:rank=species:sci_name=Heterothyone alba:taxid=1902861],1902877:1[&&NHX:rank=species:sci_name=Placothuria squamata:taxid=1902877]):0;'
    assert pruned_tree.write(properties=[], format=1) == expected_tree_no_root
    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree_with_root

# test pruned_by in order to test if data type is process correctly
def test_pruned_by_00():
    # test "contains" in leaf name
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    expected_tree = '((B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root:1[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    
    condition_inputs = ["name contains A"]
    pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)

    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree

def test_pruned_by_01():
    # test "contains" in leaf node in categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    
    expected_tree = '((B:1[&&NHX:alphabet_type=consonant],(D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root:1[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'
   
    condition_inputs = ["alphabet_type=vowel"]
    pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)
    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree

def test_pruned_by_02():
    # test "contains" in internal node in categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    expected_tree = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root:1[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'
    condition_inputs = ["alphabet_type_counter:consonant < 2"]
    pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)

    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree

def test_pruned_by_03():
    # test operators in leaf node in numerical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    expected_tree = '(((E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_avg=3.5:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_sum=7.0])Internal_2:0.5[&&NHX:col1_avg=3.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_sum=9.0])Root:1[&&NHX:col1_avg=2.5:col1_max=4.0:col1_min=1.0:col1_std=1.6666666666666667:col1_sum=10.0];'
    condition_inputs = ["col1 < 3"]
    pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)
    
    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree


def test_pruned_by_04():
    # test operators in internal node in numerical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", parser='newick')

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    expected_tree = '(A:1[&&NHX:col1=1.0])Root:1[&&NHX:col1_avg=2.5:col1_max=4.0:col1_min=1.0:col1_std=1.6666666666666667:col1_sum=10.0];'
    condition_inputs = ["col1_avg < 3.5"]
    pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)
    
    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree

def test_pruned_by_05():
    # test "contains" in leaf node in list data
    # internal_nodes annotation list data
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);", parser='newick')
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tlist_data\nA\ta,b,c\nB\tc,d\nD\ta,c,d,e\nE\te,d,b\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    
    expected_tree = '((B:1[&&NHX:list_data=c|d],(E:1[&&NHX:list_data=e|d|b])N4:0.5[&&NHX:list_data_counter=a--1||b--1||c--1||d--2||e--2])N5:0.5[&&NHX:list_data_counter=a--1||b--1||c--2||d--3||e--2])Root:1[&&NHX:list_data_counter=a--2||b--2||c--3||d--3||e--2];'
    condition_inputs = ['list_data contains a']
    pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)

    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree

def test_pruned_by_06():
    # test "contains" in internal node in list data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);", parser='newick')
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tlist_data\nA\ta,b,c\nB\tc,d\nD\ta,c,d,e\nE\te,d,b\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    
    expected_tree = '(A:1[&&NHX:list_data=a|b|c])Root:1[&&NHX:list_data_counter=a--2||b--2||c--3||d--3||e--2];'
    condition_inputs = ['list_data_counter:a<2']
    pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)
    
    assert pruned_tree.write(properties=[], format=1, format_root_node=True) == expected_tree

# def test_error():
#     assert 1 == 0, 'Not implemented yet.'