from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, ScaleFace



__all__ = [ "LayoutBarplot" ]


class LayoutPlot(TreeLayout):
    def __init__(self, name=None, prop=None, width=200, size_prop=None, color_prop=None, 
            position="aligned", column=0, padding_x=10, internal_rep='avg'):
        super().__init__(name, aligned_faces=True if position == "aligned" else False)
    
        self.width = width
        self.position = position
        self.column = column
        
        self.padding_x = padding_x

        self.internal_rep = internal_rep
        self.prop = prop
        # if not (size_prop or color_prop):
            # raise InvalidUsage("Either size_prop or color_prop required")

        self.size_prop = size_prop
        self.color_prop = color_prop

    def set_tree_style(self, tree, tree_style):
        def update_vals(metric, node):
            p, minval, maxval, uniqvals = vals[metric]
            prop = node.props.get(p)
            try:
                prop = float(prop)
                if type(prop) in [int, float]:
                    vals[metric][1] = min(minval, prop)
                    vals[metric][2] = max(maxval, prop)
                elif prop is None:
                    return
                else:
                    uniqvals.add(prop)
            except:
                pass

        vals = { 
            "size": [ self.size_prop, 0, 0, set() ],    # min, max, unique
            "color":  [ self.color_prop,  0, 0, set() ] # min, max, unique
            }

        for node in tree.traverse():
            if self.size_prop:
                update_vals("size", node)

            if self.color_prop:
                update_vals("color", node)
                
        if self.size_prop:
            self.size_range = vals["size"][1:3]
            print(vals["size"][1:3])

    def get_size(self, node, prop):
        if self.size_range != [0, 0]:
            minval, maxval = self.size_range
        else:
            minval, maxval = 0,1

        return float(node.props.get(prop, 0)) / float(maxval) * self.width


class LayoutBarplot(LayoutPlot):
    def __init__(self, name=None, prop=None, width=200, size_prop=None,
            color_prop=None, position="aligned", column=0, padding_x=10, internal_rep='avg'):

        name = name or f'Barplot_{size_prop}_{color_prop}'
        super().__init__(name=name, prop=prop, width=width, size_prop=size_prop,
                color_prop=color_prop, position=position, column=column,
                padding_x=padding_x, internal_rep=internal_rep)

    def set_tree_style(self, tree, tree_style):
        super().set_tree_style(tree, tree_style)
            
        if self.width:
            face = ScaleFace(width=self.width, scale_range=self.size_range, 
                    formatter='%.2f',
                    padding_x=self.padding_x, padding_y=10)
            tree_style.aligned_panel_header.add_face(face, column=self.column)

    def set_node_style(self, node):
        internal_prop = self.prop + '_' + self.internal_rep
        if node.is_leaf() and node.props.get(self.prop):
            width = self.get_size(node, self.size_prop)
            color = self.color_prop
            face = RectFace(width, None, color=color, padding_x=self.padding_x)
            node.add_face(face, position=self.position, column=self.column,
                    collapsed_only=False)
                    
        elif node.is_leaf() and node.props.get(internal_prop):
            width = self.get_size(node, internal_prop)
            color = self.color_prop
            face = RectFace(width, None, color=color, padding_x=self.padding_x)
            node.add_face(face, position=self.position, column=self.column,
                    collapsed_only=False)

        elif node.props.get(internal_prop):
            width = self.get_size(node, internal_prop)
            color = self.color_prop
            face = RectFace(width, None, color=color, padding_x=self.padding_x)
            node.add_face(face, position=self.position, column=self.column,
                    collapsed_only=True)