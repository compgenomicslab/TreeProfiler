
import sys
import os
from io import StringIO
import unittest

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

#from collections import namedtuple
from tempfile import NamedTemporaryFile

from treeprofiler import tree_annotate
from treeprofiler.src import utils

class TestPrune(unittest.TestCase):
    # test pruned_by in order to test if data type is process correctly
    def test_pruned_by_00(self):
        # test "contains" in leaf name
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
            f_annotation.flush()

            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
        
        expected_tree = '((B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        props = ['alphabet_type', 'alphabet_type_counter']
        condition_inputs = ["name contains A"]
        pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)
        
        self.assertEqual(pruned_tree.write(props=props, parser=parser, format_root_node=True), expected_tree)

    def test_pruned_by_01(self):
        # test "contains" in leaf node in categorical data
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
            f_annotation.flush()

            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        
        expected_tree = '((B:1[&&NHX:alphabet_type=consonant],(D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'
        props = ['alphabet_type', 'alphabet_type_counter']
        condition_inputs = ["alphabet_type=vowel"]
        pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)
        self.assertEqual(pruned_tree.write(props=props, parser=parser, format_root_node=True), expected_tree)

    def test_pruned_by_02(self):
        # test "contains" in internal node in categorical data
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tconsonant\nD\tconsonant\nE\tvowel\n')
            f_annotation.flush()

            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        props = ['alphabet_type', 'alphabet_type_counter']
        expected_tree = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1])Root[&&NHX:alphabet_type_counter=consonant--2||vowel--2];'
        condition_inputs = ["alphabet_type_counter:consonant < 2"]
        pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)

        self.assertEqual(pruned_tree.write(props=props, parser=parser, format_root_node=True), expected_tree)

    def test_pruned_by_03(self):
        # test operators in leaf node in numerical data
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
            f_annotation.flush()

            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        props = ['col1', 'col1_avg', 'col1_sum', 'col1_max', 'col1_min', 'col1_std']
        
        expected_tree = '(((E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_avg=3.5:col1_sum=7.0:col1_max=4.0:col1_min=3.0:col1_std=0.5])Internal_2:0.5[&&NHX:col1_avg=3.0:col1_sum=9.0:col1_max=4.0:col1_min=2.0:col1_std=1.0])Root[&&NHX:col1_avg=2.5:col1_sum=10.0:col1_max=4.0:col1_min=1.0:col1_std=1.6666666666666667];'
        condition_inputs = ["col1 < 3"]
        pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)
        
        self.assertEqual(pruned_tree.write(props=props, parser=parser, format_root_node=True), expected_tree)

    def test_pruned_by_04(self):
        # test operators in internal node in numerical data
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tcol1\nA\t1\nB\t2\nD\t3\nE\t4\n')
            f_annotation.flush()

            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)

        props = ['col1', 'col1_avg', 'col1_sum', 'col1_max', 'col1_min', 'col1_std']

        expected_tree = '(A:1[&&NHX:col1=1.0])Root[&&NHX:col1_avg=2.5:col1_sum=10.0:col1_max=4.0:col1_min=1.0:col1_std=1.6666666666666667];'
        condition_inputs = ["col1_avg < 3.5"]
        pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)

        self.assertEqual(pruned_tree.write(props=props, parser=parser, format_root_node=True), expected_tree)

    def test_pruned_by_05(self):
        # test "contains" in leaf node in list data
        # internal_nodes annotation list data
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5)Root;")

        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tlist_data\nA\ta,b,c\nB\tc,d\nD\ta,c,d,e\nE\te,d,b\n')
            f_annotation.flush()
            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        props = ['list_data', 'list_data_counter']
        expected_tree = '((B:1[&&NHX:list_data=c|d],(E:1[&&NHX:list_data=e|d|b])N4:0.5[&&NHX:list_data_counter=a--1||b--1||c--1||d--2||e--2])N5:0.5[&&NHX:list_data_counter=a--1||b--1||c--2||d--3||e--2])Root[&&NHX:list_data_counter=a--2||b--2||c--3||d--3||e--2];'
        condition_inputs = ['list_data contains a']
        pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)

        self.assertEqual(pruned_tree.write(props=props, parser=parser, format_root_node=True), expected_tree)
    
    def test_pruned_by_06(self):
        # test "contains" in internal node in list data
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1):0.5):0.5);")
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tlist_data\nA\ta,b,c\nB\tc,d\nD\ta,c,d,e\nE\te,d,b\n')
            f_annotation.flush()
            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        props = ['list_data', 'list_data_counter']
        expected_tree = '(A:1[&&NHX:list_data=a|b|c])Root[&&NHX:list_data_counter=a--2||b--2||c--3||d--3||e--2];'
        condition_inputs = ['list_data_counter:a<2']
        pruned_tree = utils.conditional_prune(test_tree_annotated, condition_inputs, prop2type)

        self.assertEqual(pruned_tree.write(props=props, parser=parser, format_root_node=True), expected_tree)

if __name__ == '__main__':
    unittest.main()
#pytest.main(['-v'])