from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace
from ete4.smartview  import RectFace, CircleFace, SeqMotifFace, TextFace, OutlineFace
from treeprofiler.layouts.general_layouts import get_piechartface, get_stackedbarface

"""
label_layout, colorbranch_layout, rectangular_layout   
"""
#paried_color = ["red", "darkblue", "darkgreen", "darkyellow", "violet", "mediumturquoise", "sienna", "lightCoral", "lightSkyBlue", "indigo", "tan", "coral", "olivedrab", "teal"]

class LayoutText(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, width=70, min_fsize=5, max_fsize=15, padding_x=1, padding_y=0, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'
        self.legend = legend
        self.width = width
        self.height = None
        self.min_fsize = min_fsize 
        self.max_fsize = max_fsize
        self.absence_color = "#EBEBEB"
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        if node.is_leaf and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                if self.color_dict:
                    prop_face = TextFace(prop_text, color=self.color_dict.get(prop_text, 'black'),min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width )
                else:
                    prop_face = TextFace(prop_text, color='black', min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width )
            node.add_face(prop_face, column=self.column, position="aligned")
            
        elif node.is_leaf and node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=True)
        else:
            #prop_face = CircleFace(radius=self.radius, color='grey', padding_x=self.padding_x, padding_y=self.padding_y)
            prop_face = RectFace(width=self.width, height=self.height, color=self.absence_color, \
                    padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
            node.add_face(prop_face, column=self.column, position="aligned", collapsed_only=True)

class LayoutColorbranch(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, legend=True, width=70, padding_x=1, padding_y=0):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.internal_prop = text_prop+'_counter'
        self.legend = legend
        self.height = None
        self.absence_color = "#EBEBEB"
        self.width = width
        self.padding_x = padding_x
        self.padding_y = padding_y

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=5, max_fsize=15, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)
        if self.legend:
            if self.color_dict:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )

    def set_node_style(self, node):
        if node.is_leaf and node.props.get(self.text_prop):
            prop_text = node.props.get(self.text_prop)
            if prop_text:
                if type(prop_text) == list:
                    prop_text = ",".join(prop_text)
                else:
                    pass
                if self.color_dict:
                    
                    node.add_face(TextFace(node.name, color = self.color_dict.get(prop_text,""), 
                    padding_x=self.padding_x),column=0, position="branch_right")
                    node.sm_style["hz_line_color"] = self.color_dict.get(prop_text,"")
                    node.sm_style["hz_line_width"] = 2
                    node.add_face(RectFace(width=self.width, height=None, color=self.absence_color, \
                        padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None),column=self.column, position="aligned")
            
        elif node.is_leaf and node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=True)
        
        else:
            prop_face = RectFace(width=self.width, height=self.height, color=self.absence_color, \
                    padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
            node.add_face(prop_face, column=self.column, position="aligned", collapsed_only=True)

class LayoutRect(TreeLayout):
    def __init__(self, name, column, color_dict, text_prop, width=70, height=None, padding_x=1, padding_y=0, legend=True):
        super().__init__(name)
        self.aligned_faces = True
        self.text_prop = text_prop
        self.column = column
        self.color_dict = color_dict
        self.absence_color = "#EBEBEB"
        self.internal_prop = text_prop+'_counter'
        self.legend = legend
        self.width = width
        self.height = height
        self.min_fsize = 5
        self.max_fsize = 15
        self.padding_x = padding_x
        self.padding_y = padding_y
    
    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
        text = TextFace(self.text_prop, min_fsize=self.min_fsize, max_fsize=self.max_fsize, padding_x=self.padding_x, width=self.width, rotation=315)
        tree_style.aligned_panel_header.add_face(text, column=self.column)

        if self.legend:
            if self.color_dict:
                self.color_dict['NA'] = self.absence_color
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap=self.color_dict
                                    )
            else:
                tree_style.add_legend(title=self.text_prop,
                                    variable='discrete',
                                    colormap={'NA':self.absence_color}
                                    )
    def set_node_style(self, node):
        if node.is_leaf:
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
            else:
                #prop_face = CircleFace(radius=self.radius, color='grey', padding_x=self.padding_x, padding_y=self.padding_y)
                prop_face = RectFace(width=self.width, height=self.height, text="NA", color=self.absence_color, \
                        padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
                node.add_face(prop_face, column=self.column, position="aligned")
        
        elif node.is_leaf and node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=False)

        elif node.props.get(self.internal_prop):
            stackedbar_face = get_stackedbarface(node, self.internal_prop, self.color_dict, width=self.width, padding_x=self.padding_x, padding_y=self.padding_y)
            node.add_face(stackedbar_face, column = self.column, position = "aligned", collapsed_only=True)

        else:
            
            prop_face = RectFace(width=self.width, height=self.height, color=self.absence_color, \
                    padding_x=self.padding_x , padding_y=self.padding_y, tooltip=None)
            node.add_face(prop_face, column=self.column, position="aligned", collapsed_only=True)