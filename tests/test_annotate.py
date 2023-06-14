# NOTE(JBC): All tests fail. It should be the opposite.

import sys
import os
from io import StringIO
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

#from collections import namedtuple
from tempfile import NamedTemporaryFile

import tree_annotate

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

def test_annotate_taxnomic_NCBI_01():
    # taxid in the leaf name
    test_tree = tree_annotate.ete4_parse("((9598, 9606), 10090);", parser='newick')
    test_tree_annotated, rank2values = tree_annotate.annotate_taxa(test_tree, db='NCBI')
    
    expected_tree = "((9598:1[&&NHX:common_name=Pan troglodytes:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes:rank=species:sci_name=Pan troglodytes:taxid=9598],9606:1[&&NHX:common_name=Homo sapiens:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens:rank=species:sci_name=Homo sapiens:taxid=9606])Homininae:1[&&NHX:common_name=:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae:rank=subfamily:sci_name=Homininae:taxid=207598],10090:1[&&NHX:common_name=Mus musculus:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus:rank=species:sci_name=Mus musculus:taxid=10090]);"
    assert test_tree_annotated.write(properties=[], format=1) == expected_tree

def test_annotate_taxnomic_NCBI_02():
    assert 1 == 0, 'Not implemented yet.'

def test_annotate_taxnomic_GTDB_01():
    # NOTE(JBC): I think you mean this test is not done.
    assert 1 == 0, 'Not implemented yet.'
