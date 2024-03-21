
import sys
import os
import tarfile
from io import StringIO, BytesIO
import unittest
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

#from collections import namedtuple
from tempfile import NamedTemporaryFile, TemporaryDirectory

from ete4 import Tree
from treeprofiler import tree_annotate
from treeprofiler.src import utils

class TestMSA(unittest.TestCase):
    def test_annotate_msa(self):
        # test alignment
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree_msa = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # load metadata
        with NamedTemporaryFile(suffix='.faa') as f_alignment:
            f_alignment.write(b'>A\nMAEIPDETIQQFMALT---HNIAVQYLSEFGDLNEALNSYYASQTDDIKDRREEAH\n>B\nMAEIPDATIQQFMALTNVSHNIAVQY--EFGDLNEALNSYYAYQTDDQKDRREEAH\n>E\nMAEIPDATIQ---ALTNVSHNIAVQYLSEFGDLNEALNSYYASQTDDQPDRREEAH\n>D\nMAEAPDETIQQFMALTNVSHNIAVQYLSEFGDLNEAL--------------REEAH')
            f_alignment.flush()

            #metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])

            test_tree_annotated_msa, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree_msa, 
            alignment=f_alignment.name, emapper_annotations=None
            )
        
            expected_tree_msa = '(A:1[&&NHX:alignment=MAEIPDETIQQFMALT---HNIAVQYLSEFGDLNEALNSYYASQTDDIKDRREEAH],(B:1[&&NHX:alignment=MAEIPDATIQQFMALTNVSHNIAVQY--EFGDLNEALNSYYAYQTDDQKDRREEAH],(E:1[&&NHX:alignment=MAEIPDATIQ---ALTNVSHNIAVQYLSEFGDLNEALNSYYASQTDDQPDRREEAH],D:1[&&NHX:alignment=MAEAPDETIQQFMALTNVSHNIAVQYLSEFGDLNEAL--------------REEAH])Internal_1:0.5[&&NHX:alignment=MAE-PD-TIQQFMALTNVSHNIAVQYLSEFGDLNEALNSYYASQTDDQPDRREEAH])Internal_2:0.5[&&NHX:alignment=MAE-PD-TIQQFMALTNVSHNIAVQYLSEFGDLNEALNSYYA-QTDDQ-DRREEAH])Root[&&NHX:alignment=MAEIPD-TIQQFMALTNVSHNIAVQYLSEFGDLNEALNSYYA-QTDD--DRREEAH];'
            self.assertEqual(test_tree_annotated_msa.write(props=['alignment'], parser=parser, format_root_node=True), expected_tree_msa)

if __name__ == '__main__':
    unittest.main()
#pytest.main(['-v'])