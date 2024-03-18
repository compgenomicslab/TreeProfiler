from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from treeprofiler.layouts.general_layouts import get_piechartface, get_stackedbarface

class LayoutACRDiscrete(TreeLayout):
    def __init__(self, name, column, color_dict, acr_prop, legend=True, width=70, padding_x=1, padding_y=0):
        super().__init__(name)
        self.aligned_faces = True
        self.acr_prop = acr_prop
        self.delta_prop = acr_prop+"_delta"
        self.pval_prop = acr_prop+"_pval"
        self.column = column
        self.color_dict = color_dict
        self.legend = legend
        self.height = None
        self.absence_color = "black"
        self.width = width
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.acr_prop, min_fsize=5, max_fsize=15, padding_x=self.padding_x, width=self.width, rotation=315)
        #tree_style.aligned_panel_header.add_face(text, column=self.column)
        if self.legend:
            if self.color_dict:
                #self.color_dict["Undecided Ancestral Character State"] = self.absence_color 
                tree_style.add_legend(title=self.name,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )
    def format_p_value(self, pval):
        if pval < 0.001:
            return "P < 0.001"
        elif pval < 0.01:
            return f"P = {pval:.3f}"
        else:
            # Handle rounding issue for P values greater than .01
            if 0.049 <= pval < 0.050 or 0.0049 <= pval < 0.005:
                return f"P = {pval:.3f}"
            else:
                return f"P = {pval:.2f}"

    def set_node_style(self, node):
        if node.props.get(self.acr_prop):
            prop_text = node.props.get(self.acr_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                if self.color_dict:
                    node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                    padding_x=self.padding_x),column=0, position="branch_right")
                    node.sm_style["hz_line_color"] = self.color_dict.get(prop_text,"")
                    node.sm_style["hz_line_width"] = 3
                    node.sm_style["vt_line_color"] = self.color_dict.get(prop_text,"")
                    node.sm_style["vt_line_width"] = 3
                    node.sm_style["outline_color"] = self.color_dict.get(prop_text, self.absence_color)
                    node.sm_style["size"] = 4
                    node.sm_style["fgcolor"] = self.color_dict.get(prop_text, self.absence_color)

        # # Delta statistic
        if node.props.get(self.delta_prop):
            prop_text = "%.2f" % float(node.props.get(self.delta_prop))
            if prop_text:
                output = u"\u0394" + f"-{self.acr_prop}: " + prop_text
                node.add_face(TextFace(output, color = "red", 
                padding_x=self.padding_x*5), column=0, position="branch_right")
            # p_value
            if node.props.get(self.pval_prop) is not None:
                pval = float(node.props.get(self.pval_prop))

                prop_text = self.format_p_value(pval)
                if prop_text:
                    node.add_face(TextFace(prop_text, color = "red", 
                    padding_x=self.padding_x*5), column=0, position="branch_right")

class LayoutACRContinuous(TreeLayout):
    def __init__(self, name, column, color_dict, score_prop, value_range=None, color_range=None, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.score_prop = score_prop
        self.color_dict = color_dict
        self.legend = legend
        self.absence_color = "black"
        self.value_range = value_range
        self.color_range = color_range
        self.line_width = 5

    def set_tree_style(self, tree, tree_style):
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.name,
                                    variable='continuous',
                                    value_range=self.value_range,
                                    color_range=self.color_range,
                                    )

    def set_node_style(self, node):
        prop_score = node.props.get(self.score_prop)
        if prop_score:
            prop_score = float(prop_score)
            if self.color_dict:
                # node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                # padding_x=self.padding_x),column=0, position="branch_right")
                node.sm_style["hz_line_color"] = self.color_dict.get(prop_score,"")
                node.sm_style["hz_line_width"] = self.line_width
                node.sm_style["vt_line_color"] = self.color_dict.get(prop_score,"")
                node.sm_style["vt_line_width"] = self.line_width
                node.sm_style["outline_color"] = self.color_dict.get(prop_score, self.absence_color)

class LayoutLineageSpecific(TreeLayout):
    def __init__(self, name, ls_prop, color, legend=True, active=True):
        super().__init__(name)
        self.ls_prop = ls_prop
        self.color = color
        self.legend = legend
        self.active = active
    def set_tree_style(self, tree, tree_style):
        if self.legend:
            
            tree_style.add_legend(title=self.name,
                                variable='discrete',
                                colormap={
                                    self.ls_prop: self.color,
                                }
                                )

    def set_node_style(self, node):
        if node.props.get(self.ls_prop):
            node.sm_style["bgcolor"] = self.color # highligh clade
            node.sm_style["outline_color"] = self.color