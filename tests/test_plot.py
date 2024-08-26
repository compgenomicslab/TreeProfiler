
import sys
import os
import tarfile
import requests
import time
import json
from collections import defaultdict
import unittest

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

from multiprocessing import Process
from treeprofiler import tree_plot
from treeprofiler.layouts import (
    text_layouts, taxon_layouts, staple_layouts, 
    conditional_layouts, seq_layouts, profile_layouts)
from treeprofiler.src import utils

paried_color = ["red", "darkblue", "lightgreen", "sienna", "lightCoral", "violet", "mediumturquoise",   "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal", "darkyellow"]

def get_parallel(tree, layouts, expected_draw, expected_config):
    p1 = Process(target=tree_session, args=(tree, layouts))
    p2 = Process(target=exam_tree, args=(expected_draw,expected_config))
    p1.start()
    time.sleep(1)
    p2.start()
    time.sleep(2)
    p1.terminate()
    p1.join()

def exam_tree(expected_draw, expected_config):
    url = 'http://127.0.0.1:5002/' 
    response1 = requests.get(url+'trees/0/draw')
    response2 = requests.get(url+'layouts/0')
    # assert str(json.loads(response2.text)) == expected_config
    # assert str(json.loads(response1.text)) == expected_draw

def tree_session(tree, layouts):
    #t = Tree(TREEFILE, format=1)
    tree.explore(layouts=layouts, show_leaf_name=True, 
                show_branch_length=True, show_branch_support=True, port=5002,
                keep_server=True, open_browser=False) 

class TestTreePlots(unittest.TestCase):
    def test_plot_01(self):
        # basic annotate categorical data
        # load tree
        test_tree = utils.ete4_parse('(a);')
        layouts = []
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False}}"
        expected_draw = "[['line', [0, 0.5], [1.0, 0.5], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 1.0, 1.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        get_parallel(test_tree, layouts, expected_draw, expected_config)

    def test_plot_02(self):
        # label_layout
        newick = "(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);"
        test_tree = utils.ete4_parse(newick)
        level = 1
        prop2type = {# start with leaf name
                    'name':str,
                    'dist':float,
                    'support':float,
                    'rank': str,
                    'sci_name': str,
                    'taxid': str,
                    'lineage':str,
                    'named_lineage': str,
                    'evoltype': str,
                    'dup_sp': str,
                    'dup_percent': float,
                    }
        #popup_prop_keys = list(prop2type.keys()) 

        label_layouts, level, color_dict = tree_plot.get_label_layouts(test_tree, ["alphabet_type"], level, prop2type=prop2type)
        layouts = []
        layouts.extend(label_layouts)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Label_alphabet_type': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Label_alphabet_type', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'text_prop': 'alphabet_type', 'column': 1, 'color_dict': {'consonant': '#9a312f', 'vowel': '#9b57d0'}, 'internal_prop': 'alphabet_type_counter', 'width': 70, 'height': None, 'min_fsize': 5, 'max_fsize': 15, 'absence_color': '#EBEBEB', 'padding_x': 1, 'padding_y': 0}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        
        self.assertEqual(str(layouts[0].__dict__), expected_layout)

    def test_plot_03(self):
        # label_layout
        newick = "(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);"
        test_tree = utils.ete4_parse(newick)
        level = 1
        prop2type = {# start with leaf name
                    'name':str,
                    'dist':float,
                    'support':float,
                    'rank': str,
                    'sci_name': str,
                    'taxid': str,
                    'lineage':str,
                    'named_lineage': str,
                    'evoltype': str,
                    'dup_sp': str,
                    'dup_percent': float,
                    }
        #popup_prop_keys = list(prop2type.keys()) 

        colorbranch_layouts, level, color_dict = tree_plot.get_colorbranch_layouts(test_tree, ["alphabet_type"], level, prop2type=prop2type)
        layouts = []
        layouts.extend(colorbranch_layouts)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Colorbranch_alphabet_type': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Colorbranch_alphabet_type', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'text_prop': 'alphabet_type', 'column': 1, 'color_dict': {'consonant': '#9a312f', 'vowel': '#9b57d0'}, 'internal_prop': 'alphabet_type_counter', 'height': None, 'absence_color': '#EBEBEB', 'width': 70, 'padding_x': 1, 'padding_y': 0}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        
        self.assertEqual(str(layouts[0].__dict__), expected_layout)

    def test_plot_04(self):
        # label_layout
        newick = "(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);"
        test_tree = utils.ete4_parse(newick)
        level = 1
        prop2type = {# start with leaf name
                    'name':str,
                    'dist':float,
                    'support':float,
                    'rank': str,
                    'sci_name': str,
                    'taxid': str,
                    'lineage':str,
                    'named_lineage': str,
                    'evoltype': str,
                    'dup_sp': str,
                    'dup_percent': float,
                    }
        #popup_prop_keys = list(prop2type.keys()) 

        rectangular_layouts, level, color_dict = tree_plot.get_rectangle_layouts(test_tree, ["alphabet_type"], level, prop2type=prop2type)
        layouts = []
        layouts.extend(rectangular_layouts)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Rectangular_alphabet_type': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Rectangular_alphabet_type', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'text_prop': 'alphabet_type', 'column': 1, 'color_dict': {'consonant': '#9a312f', 'vowel': '#9b57d0'}, 'absence_color': '#EBEBEB', 'internal_prop': 'alphabet_type_counter', 'width': 70, 'height': None, 'min_fsize': 5, 'max_fsize': 15, 'padding_x': 1, 'padding_y': 0}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        
        self.assertEqual(str(layouts[0].__dict__), expected_layout)

    def test_plot_05(self):
        newick = '(A:1[&&NHX:bool_type=True],(B:1[&&NHX:bool_type=False],(E:1[&&NHX:bool_type=False],D:1[&&NHX:bool_type=True])N4:0.5[&&NHX:bool_type_counter=False--1||True--1])N5:0.5[&&NHX:bool_type_counter=False--2||True--1]);'
        test_tree = utils.ete4_parse(newick)
        level = 1
        prop2type = {# start with leaf name
                    'name':str,
                    'dist':float,
                    'support':float,
                    'rank': str,
                    'sci_name': str,
                    'taxid': str,
                    'lineage':str,
                    'named_lineage': str,
                    'evoltype': str,
                    'dup_sp': str,
                    'dup_percent': float,
                    }
        binary, level, color_dict = tree_plot.get_binary_layouts(test_tree, ["bool_type"], level, prop2type=prop2type)
        layouts = []
        layouts.extend(binary)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Binary_bool_type': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        #expected_layout = "{'name': 'Binary_bool_type', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'bool_prop': 'bool_type', 'column': 1, 'color': '#9f3fbf', 'negative_color': '#EBEBEB', 'prop_colour_dict': {'False': 'red', 'True': 'darkblue'}, 'internal_prop': 'bool_type_counter', 'reverse': False, 'radius': 25, 'padding_x': 1, 'padding_y': 0, 'width': 70, 'height': None, 'min_fsize': 5, 'max_fsize': 10}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        #assert str(layouts[0].__dict__) == expected_layout

    def test_plot_06(self):
        newick = "(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_avg=3.5:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_sum=7.0])Internal_2:0.5[&&NHX:col1_avg=3.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_sum=9.0]);"
        test_tree = utils.ete4_parse(newick)
        level = 1
        prop2type = {# start with leaf name
                    'name':str,
                    'dist':float,
                    'support':float,
                    'rank': str,
                    'sci_name': str,
                    'taxid': str,
                    'lineage':str,
                    'named_lineage': str,
                    'evoltype': str,
                    'dup_sp': str,
                    'dup_percent': float,
                    }
        heatmap_layouts, level = tree_plot.get_heatmap_layouts(test_tree, ["col1"], level)
        layouts = []
        layouts.extend(heatmap_layouts)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Heatmap_col1': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Heatmap_col1_min-max', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'heatmap_prop': 'col1', 'internal_prop': 'col1_avg', 'column': 1, 'value_color': {1.0: '#fff5f0', 2.0: '#fca689', 3.0: '#dd2a25', 4.0: '#67000d', 3.5: '#af1117'}, 'value_range': [1.0, 4.0], 'color_range': {1: '#fff5f0', 2: '#ffece4', 3: '#fee4d8', 4: '#fdd7c6', 5: '#fdc7b2', 6: '#fcb79c', 7: '#fca689', 8: '#fc9474', 9: '#fc8464', 10: '#fb7252', 11: '#f96044', 12: '#f34c37', 13: '#ed392b', 14: '#dd2a25', 15: '#cf1c1f', 16: '#be151a', 17: '#af1117', 18: '#9a0c14', 19: '#800610', 20: '#67000d'}, 'absence_color': '#EBEBEB', 'maxval': 4.0, 'minval': 1.0, 'width': 70, 'height': None, 'padding_x': 1, 'padding_y': 0}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        self.assertEqual(str(layouts[0].__dict__), expected_layout)

    def test_plot_07(self):
        newick = "(A:1[&&NHX:col1=1.0],(B:1[&&NHX:col1=2.0],(E:1[&&NHX:col1=4.0],D:1[&&NHX:col1=3.0])Internal_1:0.5[&&NHX:col1_avg=3.5:col1_max=4.0:col1_min=3.0:col1_std=0.5:col1_sum=7.0])Internal_2:0.5[&&NHX:col1_avg=3.0:col1_max=4.0:col1_min=2.0:col1_std=1.0:col1_sum=9.0]);"
        test_tree = utils.ete4_parse(newick)
        level = 1
        prop2type = {# start with leaf name
                    'name':str,
                    'dist':float,
                    'support':float,
                    'rank': str,
                    'sci_name': str,
                    'taxid': str,
                    'lineage':str,
                    'named_lineage': str,
                    'evoltype': str,
                    'dup_sp': str,
                    'dup_percent': float,
                    }
        barplot_layouts, level, color_dict = tree_plot.get_barplot_layouts(test_tree, ["col1"], level, prop2type)
        layouts = []
        layouts.extend(barplot_layouts)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Barplot_col1': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Barplot_col1', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'width': 70, 'position': 'aligned', 'column': 1, 'scale': True, 'padding_x': 10, 'padding_y': 0, 'internal_rep': 'avg', 'prop': 'col1', 'size_prop': 'col1', 'color_prop': None, 'size_range': [0, 4.0], 'color': '#9b57d0', 'colors': None, 'color_gradient': None}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        self.assertEqual(str(layouts[0].__dict__), expected_layout)

    def test_plot_08(self):
        newick = '(A:1[&&NHX:alignment=MAEIPDETIQQFMALT---HNIAVQYLSEFGDLNEALNSYYASQTDDIKDRREEAH],(B:1[&&NHX:alignment=MAEIPDATIQQFMALTNVSHNIAVQY--EFGDLNEALNSYYAYQTDDQKDRREEAH],(E:1[&&NHX:alignment=MAEIPDATIQ---ALTNVSHNIAVQYLSEFGDLNEALNSYYASQTDDQPDRREEAH],D:1[&&NHX:alignment=MAEAPDETIQQFMALTNVSHNIAVQYLSEFGDLNEAL--------------REEAH])Internal_1:0.5[&&NHX:alignment=MAE-PD-TIQQFMALTNVSHNIAVQYLSEFGDLNEALNSYYASQTDDQPDRREEAH])Internal_2:0.5[&&NHX:alignment=MAE-PD-TIQQFMALTNVSHNIAVQYLSEFGDLNEALNSYYA-QTDDQ-DRREEAH])Root:0[&&NHX:alignment=MAEIPD-TIQQFMALTNVSHNIAVQYLSEFGDLNEALNSYYA-QTDD--DRREEAH];'
        test_tree = utils.ete4_parse(newick)
        level = 1
        lengh = len(max(utils.children_prop_array(test_tree, 'alignment'),key=len))
        aln_layout = seq_layouts.LayoutAlignment(name='Alignment_layout', 
                            alignment_prop='alignment', column=level, scale_range=lengh,
                            summarize_inner_nodes=True)
        layouts = []
        layouts.append(aln_layout)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Alignment_layout': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Alignment_layout', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'alignment_prop': 'alignment', 'width': 700, 'height': 15, 'column': 1, 'format': 'seq', 'scale_range': (0, 56), 'summarize_inner_nodes': True}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        self.assertEqual(str(layouts[0].__dict__), expected_layout)

    def test_plot_09(self):
        # label_layout
        newick = "(A:1[&&NHX:alphabet_type=vowel],(B:1[&&NHX:alphabet_type=consonant],(E:1[&&NHX:alphabet_type=vowel],D:1[&&NHX:alphabet_type=consonant])Internal_1:0.5[&&NHX:alphabet_type_counter=consonant--1||vowel--1])Internal_2:0.5[&&NHX:alphabet_type_counter=consonant--2||vowel--1]);"
        test_tree = utils.ete4_parse(newick)
        level = 1
        layouts = []
        profiling_prop = 'alphabet_type'
        matrix, all_values = tree_plot.single2profile(test_tree, profiling_prop)
        profile_layout = profile_layouts.LayoutProfile(name=f'Profiling_{profiling_prop}', mode='multi',alignment=matrix, seq_format='profiles', profiles=all_values, column=level, summarize_inner_nodes=True)
        layouts.append(profile_layout)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Profiling_alphabet_type': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Profiling_alphabet_type', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'mode': 'multi', 'width': 40, 'height': 20, 'column': 1, 'seq_format': 'profiles', 'profiles': ['consonant', 'vowel'], 'length': 2, 'scale_range': (0, 2), 'value_range': [], 'value_color': {}, 'summarize_inner_nodes': True}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        layout_dict = layouts[0].__dict__
        del layout_dict['alignment']
        self.assertEqual(str(layout_dict), expected_layout)

    def test_plot_10(self):
        # seqgroup is wrong
        newick = '(A:1[&&NHX:list_data=a|b|c],(B:1[&&NHX:list_data=c|d],(E:1[&&NHX:list_data=e|d|b],D:1[&&NHX:list_data=a|c|d|e])N4:0.5[&&NHX:list_data_counter=a--1||b--1||c--1||d--2||e--2])N5:0.5[&&NHX:list_data_counter=a--1||b--1||c--2||d--3||e--2]);'
        test_tree = utils.ete4_parse(newick)
        level = 1
        layouts = []
        profiling_prop = 'list_data'
        matrix, all_values = tree_plot.multiple2profile(test_tree, profiling_prop)
        profile_layout = profile_layouts.LayoutProfile(name=f'Profiling_{profiling_prop}', mode='multi',
        alignment=matrix, seq_format='profiles', profiles=all_values, column=level, summarize_inner_nodes=False)
        level += 1
        layouts.append(profile_layout)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'Profiling_list_data': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'Profiling_list_data', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'mode': 'multi', 'width': 80, 'height': 20, 'column': 1, 'seq_format': 'profiles', 'profiles': ['a|b|c', 'a|c|d|e', 'c|d', 'e|d|b'], 'length': 4, 'scale_range': (0, 4), 'value_range': [], 'value_color': {}, 'summarize_inner_nodes': False}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        layout_dict = layouts[0].__dict__
        del layout_dict['alignment']
        self.assertEqual(str(layout_dict), expected_layout)

    def test_plot_11(self):
        # categorical to matrix
        newick = '(A:1[&&NHX:col1=vowel:col2=vowel:col3=1.0:col4=1.0:col5=True:col6=True:col7=a|b|c],(B:1[&&NHX:col1=consonant:col2=consonant:col3=2.0:col4=2.0:col5=False:col6=False:col7=c|d],(E:1[&&NHX:col1=vowel:col2=vowel:col3=4.0:col4=4.0:col5=False:col6=False:col7=e|d|b],D:1[&&NHX:col1=consonant:col2=consonant:col3=3.0:col4=3.0:col5=True:col6=True:col7=a|c|d|e])Internal_1:0.5[&&NHX:col1_counter=consonant--1||vowel--1:col2_counter=consonant--1||vowel--1:col3_avg=3.5:col3_max=4.0:col3_min=3.0:col3_std=0.5:col3_sum=7.0:col4_avg=3.5:col4_max=4.0:col4_min=3.0:col4_std=0.5:col4_sum=7.0:col5_counter=False--1||True--1:col6_counter=False--1||True--1:col7_counter=a--1||b--1||c--1||d--2||e--2])Internal_2:0.5[&&NHX:col1_counter=consonant--2||vowel--1:col2_counter=consonant--2||vowel--1:col3_avg=3.0:col3_max=4.0:col3_min=2.0:col3_std=1.0:col3_sum=9.0:col4_avg=3.0:col4_max=4.0:col4_min=2.0:col4_std=1.0:col4_sum=9.0:col5_counter=False--2||True--1:col6_counter=False--2||True--1:col7_counter=a--1||b--1||c--2||d--3||e--2])Root:0[&&NHX:col1_counter=consonant--2||vowel--2:col2_counter=consonant--2||vowel--2:col3_avg=2.5:col3_max=4.0:col3_min=1.0:col3_std=1.6666666666666667:col3_sum=10.0:col4_avg=2.5:col4_max=4.0:col4_min=1.0:col4_std=1.6666666666666667:col4_sum=10.0:col5_counter=False--2||True--2:col6_counter=False--2||True--2:col7_counter=a--2||b--2||c--3||d--3||e--2];'
        test_tree = utils.ete4_parse(newick)
        level = 1
        layouts = []
        profiling_props = ['col1', 'col2'] 
        matrix, value2color = tree_plot.categorical2matrix(test_tree, profiling_props)
        # profile_layout = profile_layouts.LayoutProfile(name='categorical_matrix_layout', 
        # mode='single', alignment=matrix, seq_format='categories', profiles=profiling_props, 
        # value_color=value2color, column=level)
        profile_layout = profile_layouts.LayoutPropsMatrixOld(name='categorical_matrix_layout',
                matrix=matrix, matrix_type='categorical', matrix_props=profiling_props,
                value_color=value2color, column=level)

        level += 1
        layouts.append(profile_layout)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'categorical_matrix_layout': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'categorical_matrix_layout', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'matrix': {'A': ['vowel', 'vowel'], 'B': ['consonant', 'consonant'], 'E': ['vowel', 'vowel'], 'D': ['consonant', 'consonant']}, 'matrix_type': 'categorical', 'matrix_props': ['col1', 'col2'], 'is_list': False, 'width': 40, 'height': 20, 'column': 1, 'length': 2, 'scale_range': (0, 2), 'value_range': [], 'value_color': {'consonant': '#9a312f', 'vowel': '#9b57d0'}, 'summarize_inner_nodes': False}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        layout_dict = layouts[0].__dict__

        self.assertEqual(str(layout_dict), expected_layout)

    def test_plot_12(self):
        # categorical to matrix
        newick = '(A:1[&&NHX:col1=vowel:col2=vowel:col3=1.0:col4=1.0:col5=True:col6=True:col7=a|b|c],(B:1[&&NHX:col1=consonant:col2=consonant:col3=2.0:col4=2.0:col5=False:col6=False:col7=c|d],(E:1[&&NHX:col1=vowel:col2=vowel:col3=4.0:col4=4.0:col5=False:col6=False:col7=e|d|b],D:1[&&NHX:col1=consonant:col2=consonant:col3=3.0:col4=3.0:col5=True:col6=True:col7=a|c|d|e])Internal_1:0.5[&&NHX:col1_counter=consonant--1||vowel--1:col2_counter=consonant--1||vowel--1:col3_avg=3.5:col3_max=4.0:col3_min=3.0:col3_std=0.5:col3_sum=7.0:col4_avg=3.5:col4_max=4.0:col4_min=3.0:col4_std=0.5:col4_sum=7.0:col5_counter=False--1||True--1:col6_counter=False--1||True--1:col7_counter=a--1||b--1||c--1||d--2||e--2])Internal_2:0.5[&&NHX:col1_counter=consonant--2||vowel--1:col2_counter=consonant--2||vowel--1:col3_avg=3.0:col3_max=4.0:col3_min=2.0:col3_std=1.0:col3_sum=9.0:col4_avg=3.0:col4_max=4.0:col4_min=2.0:col4_std=1.0:col4_sum=9.0:col5_counter=False--2||True--1:col6_counter=False--2||True--1:col7_counter=a--1||b--1||c--2||d--3||e--2])Root:0[&&NHX:col1_counter=consonant--2||vowel--2:col2_counter=consonant--2||vowel--2:col3_avg=2.5:col3_max=4.0:col3_min=1.0:col3_std=1.6666666666666667:col3_sum=10.0:col4_avg=2.5:col4_max=4.0:col4_min=1.0:col4_std=1.6666666666666667:col4_sum=10.0:col5_counter=False--2||True--2:col6_counter=False--2||True--2:col7_counter=a--2||b--2||c--3||d--3||e--2];'
        test_tree = utils.ete4_parse(newick)
        level = 1
        layouts = []
        profiling_props = ['col3', 'col4'] 
        matrix, minval, maxval, _, _, _, _ = tree_plot.numerical2matrix(test_tree, profiling_props)
        # profile_layout = profile_layouts.LayoutProfile(name='numerical_matrix_layout', mode='numerical', 
        # alignment=matrix, seq_format='gradients', profiles=profiling_props, 
        # value_range=[minval, maxval], column=level)   
        profile_layout = profile_layouts.LayoutPropsMatrixOld(name='numerical_matrix_layout',
                matrix=matrix, matrix_type='categorical', matrix_props=profiling_props,
                column=level)
        level += 1
        layouts.append(profile_layout)
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'numerical_matrix_layout': True}}"
        expected_draw = "[['line', [0, 2.0], [2.0, 2.0], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 4.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'numerical_matrix_layout', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'matrix': {'Root': [None, None], 'A': [1.0, 1.0], 'Internal_2': [None, None], 'B': [2.0, 2.0], 'Internal_1': [None, None], 'E': [4.0, 4.0], 'D': [3.0, 3.0]}, 'matrix_type': 'categorical', 'matrix_props': ['col3', 'col4'], 'is_list': False, 'width': 40, 'height': 20, 'column': 1, 'length': 2, 'scale_range': (0, 2), 'value_range': [], 'value_color': {}, 'summarize_inner_nodes': False}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        layout_dict = layouts[0].__dict__
        self.assertEqual(str(layout_dict), expected_layout)

    def test_plot_13(self):
        # taxonomic layouts
        newick = "((9598.abc:1[&&NHX:common_name=Pan troglodytes:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9596|9598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Pan|Pan troglodytes:rank=species:sci_name=Pan troglodytes:taxid=9598],9606.nca:1[&&NHX:common_name=Homo sapiens:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae|Homo|Homo sapiens:rank=species:sci_name=Homo sapiens:taxid=9606])Homininae:1[&&NHX:common_name=:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|9443|376913|314293|9526|314295|9604|207598:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae|Homininae:rank=subfamily:sci_name=Homininae:taxid=207598],10090.abd:1[&&NHX:common_name=Mus musculus:lineage=1|131567|2759|33154|33208|6072|33213|33511|7711|89593|7742|7776|117570|117571|8287|1338369|32523|32524|40674|32525|9347|1437010|314146|314147|9989|1963758|337687|10066|39107|10088|862507|10090:named_lineage=root|cellular organisms|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Deuterostomia|Chordata|Craniata|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Dipnotetrapodomorpha|Tetrapoda|Amniota|Mammalia|Theria|Eutheria|Boreoeutheria|Euarchontoglires|Glires|Rodentia|Myomorpha|Muroidea|Muridae|Murinae|Mus|Mus|Mus musculus:rank=species:sci_name=Mus musculus:taxid=10090]);"
        test_tree = utils.ete4_parse(newick)
        level = 1
        layouts = []
        taxon_color_dict = {}
        taxa_layouts = []
        rank2values = {}

        # generate a rank2values dict for pre taxonomic annotated tree
        if not rank2values:
            rank2values = defaultdict(list)
            for n in test_tree.traverse():
                if n.props.get('rank') and n.props.get('rank') != 'Unknown':
                    rank2values[n.props.get('rank')].append(n.props.get('sci_name',''))
        else:
            pass
        
        # assign color for each value of each rank
        for rank, value in sorted(rank2values.items()):
            color_dict = {} 
            nvals = len(value)
            for i in range(0, nvals):
                if nvals <= 14:
                    color_dict[value[i]] = paried_color[i]
                else:
                    color_dict[value[i]] = random_color(h=None)
            
            
            taxa_layout = taxon_layouts.TaxaClade(name='TaxaClade_'+rank, level=level, rank = rank, color_dict=color_dict)
            taxa_layouts.append(taxa_layout)

            
            taxa_layout = taxon_layouts.TaxaRectangular(name = "TaxaRect_"+rank, rank=rank ,color_dict=color_dict, column=level)
            taxa_layouts.append(taxa_layout)
                #level += 1
            taxon_color_dict[rank] = color_dict
            
        #taxa_layouts.append(taxon_layouts.TaxaRectangular(name = "Last Common Ancester", color_dict=taxon_color_dict, column=level))
        taxa_layouts.append(taxon_layouts.LayoutSciName(name = 'Taxa Scientific name', color_dict=taxon_color_dict))
        taxa_layouts.append(taxon_layouts.LayoutEvolEvents(name='Taxa Evolutionary events', prop="evoltype",
            speciation_color="blue", 
            duplication_color="red", node_size = 2,
            legend=True))
        layouts = layouts + taxa_layouts
        level += 1
        
        expected_config = "{'default': {'Branch length': True, 'Branch support': True, 'Leaf name': True, 'Number of leaves': False, 'TaxaClade_species': True, 'TaxaRect_species': True, 'TaxaClade_subfamily': True, 'TaxaRect_subfamily': True, 'Taxa Scientific name': True, 'Taxa Evolutionary events': True}}"
        expected_draw = "[['line', [0, 1.5], [2.0, 1.5], '', [], {'stroke': '#000000', 'stroke-width': 0.5, 'fill': '#e5e5e5', 'fill-opacity': 0.3, 'type': 'solid'}], ['nodebox', [0, 0, 2.0, 3.0], '', {'name': '', 'dist': 0.0, 'support': 1.0}, [], [], {'fill': 'transparent'}]]"
        expected_layout = "{'name': 'TaxaClade_species', 'active': True, 'aligned_faces': True, 'description': '', 'legend': True, 'always_render': False, 'ts': None, 'ns': None, 'activate': False, 'column': 1, 'rank': 'species', 'color_dict': {'Mus musculus': 'red', 'Pan troglodytes': 'darkblue', 'Homo sapiens': 'lightgreen'}}"
        get_parallel(test_tree, layouts, expected_draw, expected_config)
        layout_dict = layouts[0].__dict__
        self.assertEqual(str(layout_dict), expected_layout)

if __name__ == '__main__':
    unittest.main()
#pytest.main(['-v'])