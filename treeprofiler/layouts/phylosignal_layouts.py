from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from treeprofiler.layouts.general_layouts import get_piechartface, get_stackedbarface

class LayoutACRDiscrete(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, legend=True, width=70, padding_x=1, padding_y=0):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.legend = legend
        self.height = None
        self.absence_color = "#EBEBEB"
        self.width = width
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=5, max_fsize=15, padding_x=self.padding_x, width=self.width, rotation=315)
        #tree_style.aligned_panel_header.add_face(text, column=self.column)
        if self.legend:
            if self.color_dict:
                self.color_dict["Undecided Ancestral Character State"] = "black"
                tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        if node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                if self.color_dict:
                    # node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                    # padding_x=self.padding_x),column=0, position="branch_right")
                    node.sm_style["hz_line_color"] = self.color_dict.get(prop_text,"")
                    node.sm_style["hz_line_width"] = 2
                    node.sm_style["vt_line_color"] = self.color_dict.get(prop_text,"")
                    node.sm_style["vt_line_width"] = 2
                    node.sm_style["outline_color"] = self.color_dict.get(prop_text,"black")
                    node.sm_style["size"] = 3
                    node.sm_style["fgcolor"] = self.color_dict.get(prop_text,"black")

        # # Delta statistic
        # if node.props.get(self.text_prop+"_delta"):
        #     prop_text = "%.2f" % float(node.props.get(self.text_prop+"_delta"))
        #     if prop_text:

        #         node.add_face(TextFace(prop_text, color = "red", 
        #         padding_x=self.padding_x),column=0, position="branch_right")

class LayoutACRContinuous(TreeLayout):
    def __init__(self, name, column, color_dict, score_prop, value_range=None, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.score_prop = score_prop
        self.color_dict = color_dict
        self.legend = legend
        self.absence_color = "#EBEBEB"
        self.value_range = value_range
        self.line_width = 5

    def set_tree_style(self, tree, tree_style):
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.name,
                                    variable='continuous',
                                    value_range=self.value_range,
                                    color_range=['#F62D2D','#FFFFFF', "#1034A6"],
                                    )

    def set_node_style(self, node):
        prop_score = node.props.get(self.score_prop)
        if prop_score:
            if self.color_dict:
                # node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                # padding_x=self.padding_x),column=0, position="branch_right")
                node.sm_style["hz_line_color"] = self.color_dict.get(prop_score,"")
                node.sm_style["hz_line_width"] = self.line_width
                node.sm_style["vt_line_color"] = self.color_dict.get(prop_score,"")
                node.sm_style["vt_line_width"] = self.line_width
                node.sm_style["outline_color"] = self.color_dict.get(prop_score,"black")
