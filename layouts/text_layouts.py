from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from layouts.general_layouts import get_piechartface

#paried_color = ["red", "darkblue", "darkgreen", "darkyellow", "violet", "mediumturquoise", "sienna", "lightCoral", "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal"]

class LayoutText(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, width=70, min_fsize=5, max_fsize=15, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'
        self.legend = legend
        self.width = width
        self.min_fsize = min_fsize 
        self.max_fsize = max_fsize

    # def set_tree_style(self, tree, tree_style):
    #     super().set_tree_style(tree, tree_style)
    #     text = TextFace(self.name, max_fsize=11, padding_x=2)
    #     tree_style.aligned_panel_header.add_face(text, column=self.column)
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=2, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                if self.color_dict:
                    prop_face = TextFace(prop_text, color=self.color_dict.get(prop_text, 'black'),min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=2, width=self.width )
                else:
                    prop_face = TextFace(prop_text, color='black', min_fsize=self.min_fsiz, max_fsize=self.max_fsize, padding_x=2, width=self.width )
            node.add_face(prop_face, column=self.column, position="aligned")
            
        elif node.is_leaf() and node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=True)

class LayoutColorbranch(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'
        self.legend = legend

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                if self.color_dict:
                    node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                    padding_x=2),column=0, position="branch_right")
                    node.sm_style["hz_line_color"] = self.color_dict.get(prop_text,"")
                    node.sm_style["hz_line_width"] = 2
            
        elif node.is_leaf() and node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=True)

class LayoutRect(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, width=70, height=None, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'
        self.legend = legend
        self.width = width
        self.height = height
        self.min_fsize = 5
        self.max_fsize = 15
        self.padding_x = 1
        self.padding_y = 0
    # def set_tree_style(self, tree, tree_style):
    #     super().set_tree_style(tree, tree_style)
    #     text = TextFace(self.name, max_fsize=11, padding_x=1)
    #     tree_style.aligned_panel_header.add_face(text, column=self.column)
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=2, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )
                                    
    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass

                tooltip = ""
                if node.name:
                    tooltip += f'<b>{node.name}</b><br>'
                if self.text_prop:
                    tooltip += f'<br>{self.text_prop}: {prop_text}<br>'
                
                if self.color_dict:
                    color = self.color_dict.get(str(prop_text),"")
                    prop_face = RectFace(width=self.width, height=self.height, color=color, \
                        padding_x=self.padding_x , padding_y=self.padding_y, tooltip=tooltip)
                    node.add_face(prop_face, column=self.column, position="aligned")
            
        elif node.is_leaf() and node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=True)