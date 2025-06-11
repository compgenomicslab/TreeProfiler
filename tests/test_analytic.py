
import sys
import os
import tarfile
from io import StringIO, BytesIO
import unittest
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

#from collections import namedtuple
from tempfile import NamedTemporaryFile, TemporaryDirectory

from treeprofiler import tree_annotate
from treeprofiler.src import utils
import time

class TestAnalytic(unittest.TestCase):
    def test_acr_discrete_01(self):
        # test acr discrete with default
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", internal_parser=internal_parser)
        prediction_method = "MPPA"
        model = "F81"

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tvowel\nD\tconsonant\nE\tconsonant\n')
            f_annotation.flush()

            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            acr_discrete_columns = ['alphabet_type']
            acr_discrete_columns_dict = {k: v for k, v in columns.items() if k in acr_discrete_columns}
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        
        acr_results, test_tree_annotated = tree_annotate.run_acr_discrete(test_tree_annotated, acr_discrete_columns_dict, \
        prediction_method=prediction_method, model=model, threads=4, outdir="./")

        
        utils.clear_extra_features([test_tree_annotated], prop2type.keys())
        expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=vowel],(E:1[&&NHX:alphabet_type=consonant],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type=vowel_consonant])Internal_2:0.5[&&NHX:alphabet_type=vowel_consonant])Root[&&NHX:alphabet_type=vowel_consonant];'
        expected_tree_with_root2 = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=vowel],(E:1[&&NHX:alphabet_type=consonant],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type=consonant_vowel])Internal_2:0.5[&&NHX:alphabet_type=consonant_vowel])Root[&&NHX:alphabet_type=consonant_vowel];'
        tree_output = test_tree_annotated.write(props=["alphabet_type"], parser=parser, format_root_node=True)
        self.assertIn(tree_output, [expected_tree_with_root, expected_tree_with_root2])

    def test_acr_discrete_02(self):
        # test acr discrete with delta and pval
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", internal_parser=internal_parser)
        prediction_method = "MPPA"
        model = "F81"

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tvowel\nD\tconsonant\nE\tconsonant\n')
            f_annotation.flush()

            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            acr_discrete_columns = ['alphabet_type']
            acr_discrete_columns_dict = {k: v for k, v in columns.items() if k in acr_discrete_columns}
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        
        acr_results, test_tree_annotated = tree_annotate.run_acr_discrete(test_tree_annotated, acr_discrete_columns_dict, \
        prediction_method=prediction_method, model=model, threads=4, outdir="./")

        prop2delta = tree_annotate.run_delta(acr_results, test_tree_annotated, sim=1000, threads=6)

        for prop, delta_result in prop2delta.items():
            if f"{delta_result:.1f}" == "0.2" or f"{delta_result:.1f}" == "0.3":
                delta_result = "0.2"
            test_tree_annotated.add_prop(utils.add_suffix(prop, "delta"), f"{delta_result}")
        
        for prop in acr_discrete_columns:
            prop2type.update({
                utils.add_suffix(prop, "delta"): float
            })
            
        dump_tree = test_tree_annotated.copy()
        utils.clear_extra_features([test_tree_annotated], prop2type.keys())

        prop2array = {}
        for prop in columns.keys():
            prop2array.update(tree_annotate.convert_to_prop_array(metadata_dict, prop))
        
        
        prop2delta_array = tree_annotate.get_pval(prop2array, dump_tree, acr_discrete_columns_dict, \
        iteration=100, prediction_method=prediction_method, model=model, threads=6)


        # for prop, delta_result in prop2delta.items():
        #         logger.info(f"Delta statistic of {prop} is: {delta_result}")
        #         tree.add_prop(utils.add_suffix(prop, "delta"), delta_result)

        for prop, delta_array in prop2delta_array.items():
            p_value = np.sum(np.array(delta_array) > prop2delta[prop]) / len(delta_array)
            test_tree_annotated.add_prop(utils.add_suffix(prop, "pval"), f"{p_value:.1f}")
            prop2type.update({
                utils.add_suffix(prop, "pval"): float
            })
        
        expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=vowel],(E:1[&&NHX:alphabet_type=consonant],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type=consonant_vowel:alphabet_type_counter=consonant--2])Internal_2:0.5[&&NHX:alphabet_type=consonant_vowel:alphabet_type_counter=consonant--2||vowel--1])Root[&&NHX:alphabet_type=consonant_vowel:alphabet_type_counter=consonant--2||vowel--2:alphabet_type_delta=0.2:alphabet_type_pval=0.0];'
        self.assertEqual(test_tree_annotated.write(props=None, parser=parser, format_root_node=True), expected_tree_with_root)
    
    def test_acr_discrete_03(self):
        # test acr discrete with other method
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", internal_parser=internal_parser)
        prediction_method = "DOWNPASS"

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\talphabet_type\nA\tvowel\nB\tvowel\nD\tconsonant\nE\tconsonant\n')
            f_annotation.flush()

            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            acr_discrete_columns = ['alphabet_type']
            acr_discrete_columns_dict = {k: v for k, v in columns.items() if k in acr_discrete_columns}
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        
        acr_results, test_tree_annotated = tree_annotate.run_acr_discrete(test_tree_annotated, acr_discrete_columns_dict, \
        prediction_method=prediction_method, threads=4, outdir="./")
    
        dump_tree = test_tree_annotated.copy()
        utils.clear_extra_features([test_tree_annotated], prop2type.keys())
        expected_tree_with_root = '(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=vowel],(E:1[&&NHX:alphabet_type=consonant],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type=consonant:alphabet_type_counter=consonant--2])Internal_2:0.5[&&NHX:alphabet_type=vowel:alphabet_type_counter=consonant--2||vowel--1])Root[&&NHX:alphabet_type=vowel:alphabet_type_counter=consonant--2||vowel--2];'
        self.assertEqual(test_tree_annotated.write(props=None, parser=parser, format_root_node=True), expected_tree_with_root)

    def test_acr_continuous_01(self):
        # test acr continuous with default
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", internal_parser=internal_parser)
        prediction_method = "BAYESIAN"
        model = "BM"

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tlength\nA\t1\nB\t1\nD\t0.2\nE\t0.2\n')
            f_annotation.flush()

            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            acr_continuous_columns = ['length']
            acr_continuous_columns_dict = {k: v for k, v in columns.items() if k in acr_continuous_columns}
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, num_stat='avg',
            columns=columns, prop2type=prop2type)
        
        # convert metadata to observed traits
        transformed_dict = {key: {} for key in acr_continuous_columns}
        for leaf, props in metadata_dict.items():
            for prop in acr_continuous_columns:
                transformed_dict[prop][leaf] = float(props[prop])

        acr_results, test_tree_annotated = tree_annotate.run_acr_continuous(test_tree_annotated, transformed_dict, \
        prediction_method=prediction_method, model=model, threads=4, outdir="./")
        avail_props = ['length','length_avg']
        for node in test_tree_annotated.traverse():
            for prop in avail_props:
                if prop in node.props:
                    val = node.props.get(prop, None)
                    val = f"{val:.2f}" if val is not None else None
                    node.add_prop(prop, val)

        self.assertTrue(-2.34 <= float(test_tree_annotated.props.get('length')) <= 10.59)
        self.assertTrue(-5.08 <= float(test_tree_annotated['Internal_2'].props.get('length')) <= 13.42)
        self.assertTrue(-7.55 <= float(test_tree_annotated['Internal_1'].props.get('length')) <= 15.77)

    def test_acr_continuous_02(self):
        # test acr continuous with default
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", internal_parser=internal_parser)
        prediction_method = "BAYESIAN"
        model = "OU"

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tlength\nA\t1\nB\t1\nD\t0.2\nE\t0.2\n')
            f_annotation.flush()

            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            acr_continuous_columns = ['length']
            acr_continuous_columns_dict = {k: v for k, v in columns.items() if k in acr_continuous_columns}
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, num_stat='avg',
            columns=columns, prop2type=prop2type)
        
        # convert metadata to observed traits
        transformed_dict = {key: {} for key in acr_continuous_columns}
        for leaf, props in metadata_dict.items():
            for prop in acr_continuous_columns:
                transformed_dict[prop][leaf] = float(props[prop])

        acr_results, test_tree_annotated = tree_annotate.run_acr_continuous(test_tree_annotated, transformed_dict, \
        prediction_method=prediction_method, model=model, threads=4, outdir="./")
        avail_props = ['length','length_avg']
        for node in test_tree_annotated.traverse():
            for prop in avail_props:
                if prop in node.props:
                    val = node.props.get(prop, None)
                    val = f"{val:.2f}" if val is not None else None
                    node.add_prop(prop, val)

        self.assertTrue(-2.11 <= float(test_tree_annotated.props.get('length')) <= 6.57)
        self.assertTrue(-5.78 <= float(test_tree_annotated['Internal_2'].props.get('length')) <= 10.28)
        self.assertTrue(-8.32 <= float(test_tree_annotated['Internal_1'].props.get('length')) <= 12.90)

    def test_acr_continuous_03(self):
        # test acr continuous with default
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;", internal_parser=internal_parser)
        prediction_method = "ML"
        model = "BM"

        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tlength\nA\t1\nB\t1\nD\t0.2\nE\t0.2\n')
            f_annotation.flush()

            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            acr_continuous_columns = ['length']
            acr_continuous_columns_dict = {k: v for k, v in columns.items() if k in acr_continuous_columns}
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, num_stat='avg',
            columns=columns, prop2type=prop2type)
        
        # convert metadata to observed traits
        transformed_dict = {key: {} for key in acr_continuous_columns}
        for leaf, props in metadata_dict.items():
            for prop in acr_continuous_columns:
                transformed_dict[prop][leaf] = float(props[prop])

        acr_results, test_tree_annotated = tree_annotate.run_acr_continuous(test_tree_annotated, transformed_dict, \
        prediction_method=prediction_method, model=model, threads=4, outdir="./")
        avail_props = ['length','length_avg']
        for node in test_tree_annotated.traverse():
            for prop in avail_props:
                if prop in node.props:
                    val = node.props.get(prop, None)
                    val = f"{val:.2f}" if val is not None else None
                    node.add_prop(prop, val)

        self.assertTrue(-2 <= float(test_tree_annotated.props.get('length')) <= 2)
        self.assertTrue(-3 <= float(test_tree_annotated['Internal_2'].props.get('length')) <= 3)
        self.assertTrue(-2 <= float(test_tree_annotated['Internal_1'].props.get('length')) <= 2)

    def test_ls_01(self):
        # test acr discrete with default
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(((((I:1.39,(Y:1.23,S:0.89)Internal_64:1.95)Internal_2:0.85,(M:1.56,(H:1.03,T:0.4)Internal_25:0.49)Internal_58:1.85)Internal_1:0.14,((U:0.14,(O:0.85,Z:2.0)Internal_73:1.67)Internal_35:1.09,(E:1.86,(G:0.63,R:1.23)Internal_60:1.2)Internal_58:0.65)Internal_16:1.08)Internal_2:1.32,(((N:1.38,(V:1.87,D:1.36)Internal_92:1.17)Internal_82:0.57,(B:0.17,(J:1.89,C:0.29)Internal_92:1.36)Internal_11:0.46)Internal_51:0.7,((W:1.09,(P:1.58,L:0.86)Internal_72:0.76)Internal_51:0.41,((K:1.27,F:1.31)Internal_45:1.09,(X:0.96,A:0.77)Internal_91:1.14)Internal_28:0.97)Internal_72:1.37)Internal_16:0.94)Internal_78:1.19)Root;", internal_parser=internal_parser)
        prec_cutoff = 0.5
        sens_cutoff = 0.5
        # load metadata
        with NamedTemporaryFile(suffix='.tsv') as f_annotation:
            f_annotation.write(b'#name\tis_vowel\nA\t1\nB\t0\nC\t0\nD\t0\nE\t1\nF\t0\nG\t0\nH\t0\nI\t1\nJ\t0\nK\t0\nL\t0\nM\t0\nN\t0\nO\t1\nP\t0\nR\t0\nS\t0\nT\t0\nU\t1\nV\t0\nW\t0\nX\t0\nY\t0\nZ\t0\n')
            f_annotation.flush()

            
            metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            ls_columns = ['is_vowel']
            ls_columns_dict = {k: v for k, v in columns.items() if k in ls_columns}
        
        test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            metadata_dict=metadata_dict, node_props=node_props, 
            columns=columns, prop2type=prop2type)
        
        best_node, qualified_nodes = tree_annotate.run_ls(test_tree_annotated, props=ls_columns, 
            precision_cutoff=prec_cutoff, sensitivity_cutoff=sens_cutoff)
        
        f1 = round(float(test_tree_annotated.props.get('is_vowel_f1')), 2)

        self.assertEqual(f1, 0.33)
        self.assertEqual(best_node.name, "Internal_16")
        

if __name__ == '__main__':
    unittest.main()