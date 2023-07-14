
import sys
import os
import tarfile
from io import StringIO, BytesIO
import pytest

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

#from collections import namedtuple
from tempfile import NamedTemporaryFile, TemporaryDirectory

from treeprofiler import tree_annotate
from treeprofiler.src import utils

def test_annotate_01():
    # basic annotate categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse('(a:1);')

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\na\tapple')
        f_annotation.flush()
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    
    expected_tree = '(a:1[&&NHX:col1=apple]);'

    assert test_tree_annotated.write(props=None) == expected_tree

def test_annotate_02():
    # internal_nodes annotation categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    
    expected_tree_no_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);'
    expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'

    assert test_tree_annotated.write(props=None) == expected_tree_no_root
    assert test_tree_annotated.write(props=None,format_root_node=True) == expected_tree_with_root

def test_annotate_03():
    # basic annotate numerical data
    # load tree
    test_tree = tree_annotate.ete4_parse('(a:1);')

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\na\t2')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    expected_tree = '(a:1[&&NHX:col1=2.0]);'

    assert test_tree_annotated.write(props=None) == expected_tree

def test_annotate_04():
    # internal_nodes annotation numerical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    props = ['col1', 'col1_sum','col1_max','col1_min','col1_std','col1_avg']
    expected_tree_no_root = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_sum=7.0:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_avg=3.5])Internal_2:0.5[&&NHX:col1_sum=9.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_avg=3.0]);'    
    expected_tree_with_root = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_sum=7.0:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_avg=3.5])Internal_2:0.5[&&NHX:col1_sum=9.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_avg=3.0])Root[&&NHX:col1_sum=10.0:col1_max=4.0:col1_min=1.0:col1_std=1.6666666666666667:col1_avg=2.5];'

    assert test_tree_annotated.write(props=props) == expected_tree_no_root
    assert test_tree_annotated.write(props=props, format_root_node=True) == expected_tree_with_root

def test_annotate_05():
    # test num_stat none and counter_stat none
    # internal_nodes annotation categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

     # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\talphabet_type\nA\t1\tvowel\nB\t2\tconsonant\nD\t3\tconsonant\nE\t4\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
    metadata_dict=metadata_dict, node_props=node_props, counter_stat='none', num_stat='none',
    columns=columns, prop2type=prop2type)
    
    props = ["alphabet_type", "col1"]
    expected_tree_no_root = '(A:1[&&NHX:alphabet_type=vowel:col1=1.0],(B:1[&&NHX:alphabet_type=consonant:col1=2.0],(E:1[&&NHX:alphabet_type=vowel:col1=4.0],D:1[&&NHX:alphabet_type=consonant:col1=3.0])Internal_1:0.5)Internal_2:0.5);'
    expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel:col1=1.0],(B:1[&&NHX:alphabet_type=consonant:col1=2.0],(E:1[&&NHX:alphabet_type=vowel:col1=4.0],D:1[&&NHX:alphabet_type=consonant:col1=3.0])Internal_1:0.5)Internal_2:0.5)Root;'

    assert test_tree_annotated.write(props=props) == expected_tree_no_root
    assert test_tree_annotated.write(props=props,format_root_node=True) == expected_tree_with_root

def test_annotate_06():
    # assign internal node name
    # internal_nodes annotation categorical data
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree_no_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])N4:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])N5:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);'
    expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])N4:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])N5:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    #print(test_tree_annotated.write(props=None,format_root_node=True))
    assert test_tree_annotated.write(props=None) == expected_tree_no_root
    assert test_tree_annotated.write(props=None,format_root_node=True) == expected_tree_with_root

def test_annotate_07():
    # internal_nodes annotation boolean data
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tbool_type\nA\tTrue\nB\tFalse\nD\tTrue\nE\tFalse\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    expected_tree_no_root = '(A:1[&&NHX:bool_type=True],(B:1[&&NHX:bool_type=False],(E:1[&&NHX:bool_type=False],D:1[&&NHX:bool_type=True])N4:0.5[&&NHX:bool_type_counter=False--1||True--1])N5:0.5[&&NHX:bool_type_counter=False--2||True--1]);'
    expected_tree_with_root = '(A:1[&&NHX:bool_type=True],(B:1[&&NHX:bool_type=False],(E:1[&&NHX:bool_type=False],D:1[&&NHX:bool_type=True])N4:0.5[&&NHX:bool_type_counter=False--1||True--1])N5:0.5[&&NHX:bool_type_counter=False--2||True--1])Root[&&NHX:bool_type_counter=False--2||True--2];'
    
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    assert test_tree_annotated.write(props=None) == expected_tree_no_root
    assert test_tree_annotated.write(props=None, format_root_node=True) == expected_tree_with_root

def test_annotate_08():
    # internal_nodes annotation list data
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);")
    
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tlist_data\nA\ta,b,c\nB\tc,d\nD\ta,c,d,e\nE\te,d,b\n')
        f_annotation.flush()
        
        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    expected_tree_no_root = '(A:1[&&NHX:list_data=a|b|c],(B:1[&&NHX:list_data=c|d],(E:1[&&NHX:list_data=e|d|b],D:1[&&NHX:list_data=a|c|d|e])N4:0.5[&&NHX:list_data_counter=a--1||b--1||c--1||d--2||e--2])N5:0.5[&&NHX:list_data_counter=a--1||b--1||c--2||d--3||e--2]);'
    expected_tree_with_root = '(A:1[&&NHX:list_data=a|b|c],(B:1[&&NHX:list_data=c|d],(E:1[&&NHX:list_data=e|d|b],D:1[&&NHX:list_data=a|c|d|e])N4:0.5[&&NHX:list_data_counter=a--1||b--1||c--1||d--2||e--2])N5:0.5[&&NHX:list_data_counter=a--1||b--1||c--2||d--3||e--2])Root[&&NHX:list_data_counter=a--2||b--2||c--3||d--3||e--2];'

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)

    assert test_tree_annotated.write(props=None) == expected_tree_no_root
    assert test_tree_annotated.write(props=None, format_root_node=True) == expected_tree_with_root
 
def test_annotate_09():
    # specify datatype of each column 
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\tcol2\tcol3\tcol4\nA\tvowel\t1\tTrue\ta,b,c\nB\tconsonant\t2\tFalse\tc,d\nD\tconsonant\t3\tTrue\ta,c,d,e\nE\tvowel\t4\tFalse\te,d,b\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    text_prop = ['col1']
    num_prop = ['col2']
    bool_prop = ['col3']
    multiple_text_prop = ['col4']

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        text_prop=text_prop, multiple_text_prop=multiple_text_prop,
        num_prop=num_prop, bool_prop=bool_prop,
        columns=columns, prop2type=prop2type)

    props = ['col1', 'col2', 'col3', 'col4', 'col1_counter', 'col4_counter', 'col3_counter', 'col2_avg', 'col2_sum', 'col2_max', 'col2_min', 'col2_std']
    expected_tree = '(A:1[&&NHX:col1=vowel:col2=1.0:col3=True:col4=a|b|c],(B:1[&&NHX:col1=consonant:col2=2.0:col3=False:col4=c|d],(E:1[&&NHX:col1=vowel:col2=4.0:col3=False:col4=e|d|b],D:1[&&NHX:col1=consonant:col2=3.0:col3=True:col4=a|c|d|e])Internal_1:0.5[&&NHX:col1_counter=consonant--1||vowel--1:col4_counter=a--1||b--1||c--1||d--2||e--2:col3_counter=False--1||True--1:col2_avg=3.5:col2_sum=7.0:col2_max=4.0:col2_min=3.0:col2_std=0.5])Internal_2:0.5[&&NHX:col1_counter=consonant--2||vowel--1:col4_counter=a--1||b--1||c--2||d--3||e--2:col3_counter=False--2||True--1:col2_avg=3.0:col2_sum=9.0:col2_max=4.0:col2_min=2.0:col2_std=1.0])Root[&&NHX:col1_counter=consonant--2||vowel--2:col4_counter=a--2||b--2||c--3||d--3||e--2:col3_counter=False--2||True--2:col2_avg=2.5:col2_sum=10.0:col2_max=4.0:col2_min=1.0:col2_std=1.6666666666666667];'    
    assert test_tree_annotated.write(props=props, format_root_node=True) == expected_tree

def test_annotate_10():
    # specify datatype of each column index
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\tcol2\tcol3\tcol4\nA\tvowel\t1\tTrue\ta,b,c\nB\tconsonant\t2\tFalse\tc,d\nD\tconsonant\t3\tTrue\ta,c,d,e\nE\tvowel\t4\tFalse\te,d,b\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    text_prop_idx = '1'
    num_prop_idx = '2'
    bool_prop_idx = '3'
    multiple_text_prop = ['col4']

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        text_prop_idx=text_prop_idx, multiple_text_prop=multiple_text_prop,
        num_prop_idx=num_prop_idx, bool_prop_idx=bool_prop_idx,
        columns=columns, prop2type=prop2type)
    props = ['col1', 'col2', 'col3', 'col4', 'col1_counter', 'col4_counter', 'col3_counter', 'col2_avg', 'col2_sum', 'col2_max', 'col2_min', 'col2_std']
    
    expected_tree = '(A:1[&&NHX:col1=vowel:col2=1.0:col3=True:col4=a|b|c],(B:1[&&NHX:col1=consonant:col2=2.0:col3=False:col4=c|d],(E:1[&&NHX:col1=vowel:col2=4.0:col3=False:col4=e|d|b],D:1[&&NHX:col1=consonant:col2=3.0:col3=True:col4=a|c|d|e])Internal_1:0.5[&&NHX:col1_counter=consonant--1||vowel--1:col4_counter=a--1||b--1||c--1||d--2||e--2:col3_counter=False--1||True--1:col2_avg=3.5:col2_sum=7.0:col2_max=4.0:col2_min=3.0:col2_std=0.5])Internal_2:0.5[&&NHX:col1_counter=consonant--2||vowel--1:col4_counter=a--1||b--1||c--2||d--3||e--2:col3_counter=False--2||True--1:col2_avg=3.0:col2_sum=9.0:col2_max=4.0:col2_min=2.0:col2_std=1.0])Root[&&NHX:col1_counter=consonant--2||vowel--2:col4_counter=a--2||b--2||c--3||d--3||e--2:col3_counter=False--2||True--2:col2_avg=2.5:col2_sum=10.0:col2_max=4.0:col2_min=1.0:col2_std=1.6666666666666667];'
    assert test_tree_annotated.write(props=props, format_root_node=True) ==expected_tree 

def test_annotate_11():
    # specify datatype of each column index range
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\tcol2\tcol3\tcol4\tcol5\tcol6\tcol7\nA\tvowel\tvowel\t1\t1\tTrue\tTrue\ta,b,c\nB\tconsonant\tconsonant\t2\t2\tFalse\tFalse\tc,d\nD\tconsonant\tconsonant\t3\t3\tTrue\tTrue\ta,c,d,e\nE\tvowel\tvowel\t4\t4\tFalse\tFalse\te,d,b\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    text_prop_idx = '[1-2]'
    num_prop_idx = '[3-4]'
    bool_prop_idx = '[5-6]'
    multiple_text_prop = ['col7']

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        text_prop_idx=text_prop_idx, multiple_text_prop=multiple_text_prop,
        num_prop_idx=num_prop_idx, bool_prop_idx=bool_prop_idx,
        columns=columns, prop2type=prop2type)
    props = ['col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col1_counter', 'col2_counter', 'col7_counter', 'col5_counter', 'col6_counter', 'col3_avg', 'col3_sum', 'col3_max', 'col3_min', 'col3_std', 'col4_avg', 'col4_sum', 'col4_max', 'col4_min', 'col4_std']
    expected_tree = '(A:1[&&NHX:col1=vowel:col2=vowel:col3=1.0:col4=1.0:col5=True:col6=True:col7=a|b|c],(B:1[&&NHX:col1=consonant:col2=consonant:col3=2.0:col4=2.0:col5=False:col6=False:col7=c|d],(E:1[&&NHX:col1=vowel:col2=vowel:col3=4.0:col4=4.0:col5=False:col6=False:col7=e|d|b],D:1[&&NHX:col1=consonant:col2=consonant:col3=3.0:col4=3.0:col5=True:col6=True:col7=a|c|d|e])Internal_1:0.5[&&NHX:col1_counter=consonant--1||vowel--1:col2_counter=consonant--1||vowel--1:col7_counter=a--1||b--1||c--1||d--2||e--2:col5_counter=False--1||True--1:col6_counter=False--1||True--1:col3_avg=3.5:col3_sum=7.0:col3_max=4.0:col3_min=3.0:col3_std=0.5:col4_avg=3.5:col4_sum=7.0:col4_max=4.0:col4_min=3.0:col4_std=0.5])Internal_2:0.5[&&NHX:col1_counter=consonant--2||vowel--1:col2_counter=consonant--2||vowel--1:col7_counter=a--1||b--1||c--2||d--3||e--2:col5_counter=False--2||True--1:col6_counter=False--2||True--1:col3_avg=3.0:col3_sum=9.0:col3_max=4.0:col3_min=2.0:col3_std=1.0:col4_avg=3.0:col4_sum=9.0:col4_max=4.0:col4_min=2.0:col4_std=1.0])Root[&&NHX:col1_counter=consonant--2||vowel--2:col2_counter=consonant--2||vowel--2:col7_counter=a--2||b--2||c--3||d--3||e--2:col5_counter=False--2||True--2:col6_counter=False--2||True--2:col3_avg=2.5:col3_sum=10.0:col3_max=4.0:col3_min=1.0:col3_std=1.6666666666666667:col4_avg=2.5:col4_sum=10.0:col4_max=4.0:col4_min=1.0:col4_std=1.6666666666666667];'
    assert test_tree_annotated.write(props=props, format_root_node=True) == expected_tree

def test_annotate_12():
    # test missing data and unmapped data they should be see as the same as none
    # r'^(?:\W+|none|None|null|NaN|)$'
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")
    
    # load metadata with missing categorical data
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\talphabet_type\nA\tnone\nB\t-\nD\t\nE\tvowel\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
    
    test_tree_annotated_1, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
    metadata_dict=metadata_dict, node_props=node_props, counter_stat='raw',
    columns=columns, prop2type=prop2type)
    
    # load metadata with unmapped categorical data
    with NamedTemporaryFile(suffix='.tsv') as f_annotation_2:
        f_annotation_2.write(b'#name\talphabet_type\nA\tnone\nD\t\nE\tvowel\n')
        f_annotation_2.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation_2.name])
    
    test_tree_annotated_2, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
    metadata_dict=metadata_dict, node_props=node_props, counter_stat='raw',
    columns=columns, prop2type=prop2type)

    expected_tree = '(A:1[&&NHX:alphabet_type=NaN],(B:1[&&NHX:alphabet_type=NaN],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=NaN])Internal_1:0.5[&&NHX:alphabet_type_counter=NaN--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=NaN--2||vowel--1])Root[&&NHX:alphabet_type_counter=NaN--3||vowel--1];'
    
    assert test_tree_annotated_1.write(props=None,format_root_node=True) == expected_tree
    assert test_tree_annotated_2.write(props=None,format_root_node=True) == expected_tree

def test_annotate_13():
    # test relative on categorical, boolean and list
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\tcol2\tcol3\nA\tvowel\tTrue\ta,b,c\nB\tconsonant\tFalse\tc,d\nD\tconsonant\tTrue\ta,c,d,e\nE\tvowel\tFalse\te,d,b\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, counter_stat='relative',
        columns=columns, prop2type=prop2type)
    props = ['col1', 'col2', 'col3', 'col1_counter', 'col2_counter', 'col3_counter']
    expected_tree = '(A:1[&&NHX:col1=vowel:col2=True:col3=a|b|c],(B:1[&&NHX:col1=consonant:col2=False:col3=c|d],(E:1[&&NHX:col1=vowel:col2=False:col3=e|d|b],D:1[&&NHX:col1=consonant:col2=True:col3=a|c|d|e])Internal_1:0.5[&&NHX:col1_counter=consonant--0.50||vowel--0.50:col2_counter=False--0.50||True--0.50:col3_counter=a--0.14||b--0.14||c--0.14||d--0.29||e--0.29])Internal_2:0.5[&&NHX:col1_counter=consonant--0.67||vowel--0.33:col2_counter=False--0.67||True--0.33:col3_counter=a--0.11||b--0.11||c--0.22||d--0.33||e--0.22])Root[&&NHX:col1_counter=consonant--0.50||vowel--0.50:col2_counter=False--0.50||True--0.50:col3_counter=a--0.17||b--0.17||c--0.25||d--0.25||e--0.17];'
   
    assert test_tree_annotated.write(props=props, format_root_node=True) == expected_tree

def test_annotate_14_a():
    # test different numerical stats
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated_all, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, num_stat='all',
        columns=columns, prop2type=prop2type)
    props = ['col1', 'col1_sum','col1_max','col1_min','col1_std','col1_avg']
    expected_tree_all = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_sum=7.0:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_avg=3.5])Internal_2:0.5[&&NHX:col1_sum=9.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_avg=3.0])Root[&&NHX:col1_sum=10.0:col1_max=4.0:col1_min=1.0:col1_std=1.6666666666666667:col1_avg=2.5];'
    
    assert test_tree_annotated_all.write(props=props, format_root_node=True) == expected_tree_all

def test_annotate_14_b():
    # test different numerical stats
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])


    test_tree_annotated_sum, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, num_stat='sum',
        columns=columns, prop2type=prop2type)

    expected_tree = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_sum=7.0])Internal_2:0.5[&&NHX:col1_sum=9.0])Root[&&NHX:col1_sum=10.0];'

    assert test_tree_annotated_sum.write(props=None, format_root_node=True) == expected_tree

def test_annotate_14_c():
    # test different numerical stats
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated_avg, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, num_stat='avg',
        columns=columns, prop2type=prop2type)

    expected_tree_avg = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_avg=3.5])Internal_2:0.5[&&NHX:col1_avg=3.0])Root[&&NHX:col1_avg=2.5];'
    assert test_tree_annotated_avg.write(props=None, format_root_node=True) == expected_tree_avg

def test_annotate_14_d():
    # test different numerical stats
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated_max, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, num_stat='max',
        columns=columns, prop2type=prop2type)

    expected_tree_max = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_max=4.0])Internal_2:0.5[&&NHX:col1_max=4.0])Root[&&NHX:col1_max=4.0];'
    assert test_tree_annotated_max.write(props=None, format_root_node=True) == expected_tree_max

def test_annotate_14_e():
    # test different numerical stats
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated_min, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, num_stat='min',
        columns=columns, prop2type=prop2type)

    expected_tree_min = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_min=3.0])Internal_2:0.5[&&NHX:col1_min=2.0])Root[&&NHX:col1_min=1.0];'
    
    assert test_tree_annotated_min.write(props=None, format_root_node=True) == expected_tree_min

def test_annotate_14_f():
    # test different numerical stats
    # load tree
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

    # load metadata
    with NamedTemporaryFile(suffix='.tsv') as f_annotation:
        f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
        f_annotation.flush()

        metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

    test_tree_annotated_std, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, num_stat='std',
        columns=columns, prop2type=prop2type)

    expected_tree_std = '(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_std=0.5])Internal_2:0.5[&&NHX:col1_std=1.0])Root[&&NHX:col1_std=1.6666666666666667];'
    
    assert test_tree_annotated_std.write(props=None, format_root_node=True) == expected_tree_std

def test_annotate_tar():
    # test if can read tar.gz file
    # load tree
    test_tree = tree_annotate.ete4_parse('(a);')

    # load metadata
    with TemporaryDirectory() as temp_dir:
        file1_path = temp_dir + '/metadata1.tsv'
        with open(file1_path, 'w') as file1:
            file1.write('#name\tcol1\na\tapple')

        file2_path = temp_dir + '/metadata2.tsv'
        with open(file2_path, 'w') as file2:
            file2.write('#name\tcol2\na\t3')

        with NamedTemporaryFile(suffix='.tar.gz') as temp_tar:
            tar_path = temp_tar.name

            # Create a tarfile and add the files from the temporary directory
            with tarfile.open(tar_path, 'w:gz') as tar:
                tar.add(file1_path, arcname='metadata1.tsv')
                tar.add(file2_path, arcname='metadata2.tsv')

            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([tar_path])
            
    test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
        metadata_dict=metadata_dict, node_props=node_props, 
        columns=columns, prop2type=prop2type)
    props = ['col1', 'col2']
    expected_tree = '(a:1[&&NHX:col1=apple:col2=3.0]);'
    assert test_tree_annotated.write(props=props) == expected_tree

def test_internal_parser_01():
    parser='name'
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", internal_parser=parser)
    expected_tree_paser_1 = "(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;"
    expected_tree_paser_0 = "(A:1,(B:1,(E:1,D:1):0.5[&&NHX:name=Internal_1]):0.5[&&NHX:name=Internal_2]);"

    assert test_tree.write(props=None, parser=1, format_root_node=True) == expected_tree_paser_1
    assert test_tree.write(props=None, parser=0) == expected_tree_paser_0


def test_internal_parser_02():
    parser='support'
    test_tree = tree_annotate.ete4_parse("(A:1,(B:1,(E:1,D:1)1:0.5)1:0.5);", internal_parser=parser)
    expected_tree_paser_0 = "(A:1,(B:1,(E:1,D:1)1:0.5)1:0.5);"
    expected_tree_paser_1 = "(A:1,(B:1,(E:1,D:1):0.5[&&NHX:support=1.0]):0.5[&&NHX:support=1.0]);"

    assert test_tree.write(props=None, parser=1, format_root_node=True) == expected_tree_paser_1
    assert test_tree.write(props=None, parser=0) == expected_tree_paser_0

#pytest.main(['-v'])