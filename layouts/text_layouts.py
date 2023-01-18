from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from layouts.general_layouts import get_piechartface

#paried_color = ["red", "darkblue", "darkgreen", "darkyellow", "violet", "mediumturquoise", "sienna", "lightCoral", "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal"]

class LayoutText(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'

    # def set_tree_style(self, tree, tree_style):
    #     super().set_tree_style(tree, tree_style)
    #     text = TextFace(self.name, max_fsize=11, padding_x=2)
    #     tree_style.aligned_panel_header.add_face(text, column=self.column)

    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                # human_orth = " ".join(human_orth.split('|'))
                # human_orth_face = RectFace(width=50,height=50, color=self.color)
                # node.add_face(human_orth_face, column=self.column, position="aligned")
                if self.color_dict:
                    prop_face = TextFace(prop_text, color=self.color_dict[prop_text])
                else:
                    prop_face = TextFace(prop_text, color='blue')
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
    def __init__(self, name, column, color_dict, text_prop):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'

    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if self.color_dict:
                    node.sm_style["hz_line_color"] = self.color_dict[prop_text]
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
    def __init__(self, name, column, color_dict, text_prop):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'

    # def set_tree_style(self, tree, tree_style):
    #     super().set_tree_style(tree, tree_style)
    #     text = TextFace(self.name, max_fsize=11, padding_x=1)
    #     tree_style.aligned_panel_header.add_face(text, column=self.column)

    def set_node_style(self, node):
        if node.is_leaf() and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                tooltip = ""
                if node.name:
                    tooltip += f'<b>{node.name}</b><br>'
                if self.text_prop:
                    tooltip += f'<br>{self.text_prop}: {node.props.get(self.text_prop)}<br>'
                
                if self.color_dict:
                    color = self.color_dict[str(prop_text)]
                    prop_face = RectFace(width=50,height=50, color=color, \
                        padding_x=1, padding_y=1, tooltip=tooltip)
                    node.add_face(prop_face, column=self.column, position="aligned")
            
        elif node.is_leaf() and node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            piechart_face = get_piechartface(node, self.internal_prop, self.color_dict)
            node.add_face(piechart_face, column = self.column, position = "branch_top")
            node.add_face(piechart_face, column = self.column, position = "aligned", collapsed_only=True)