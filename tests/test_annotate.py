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

    test_tree_annotated = tree_annotate.run_tree_annotate(test_tree, 
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

    test_tree_annotated = tree_annotate.run_tree_annotate(test_tree, 
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

    test_tree_annotated = tree_annotate.run_tree_annotate(test_tree, 
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

    test_tree_annotated = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    assert test_tree_annotated.write(properties=[], format=1) == expected_tree_no_root
    assert test_tree_annotated.write(properties=[], format=1, format_root_node=True) == expected_tree_with_root

def test_annotate_taxnomic_NCBI():
    assert

def test_annotate_taxnomic_GTDB():
    assert