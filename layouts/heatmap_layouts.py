from ete4.smartview import TreeStyle, NodeStyle, TreeLayout
from ete4.smartview  import RectFace, ScaleFace
import colorsys


def color_gradient(hue, intensity, granularity):
    min_lightness = 0.35 
    max_lightness = 0.9
    base_value = intensity

    # each gradient must contain 100 lightly descendant colors
    colors = []   
    rgb2hex = lambda rgb: '#%02x%02x%02x' % rgb
    l_factor = (max_lightness-min_lightness) / float(granularity)
    l = min_lightness
    while l <= max_lightness:
        l += l_factor
        rgb =  rgb2hex(tuple(map(lambda x: int(x*255), colorsys.hls_to_rgb(hue, l, base_value))))
        colors.append(rgb)
        
    colors.append("#ffffff")
    return list(reversed(colors))

def heatmap_layout(prop, level, internal_rep='avg'):
    internal_prop = prop+'_'+internal_rep
    def layout_fn(node):
        #if not node.is_root() and node.props.get(prop):
        if node.is_leaf() and node.props.get(prop):
            # heatmap
            redgradient = color_gradient(0.95, 0.6, 10)
            relative_abundance = float(node.props.get(prop))
            color_idx = int(relative_abundance*10)
            color = redgradient[color_idx]

            identF = RectFace(width=50,height=50,text="%.1f" % (relative_abundance*100), color=color, 
            padding_x=1, padding_y=1)
            #face_name = TextFace(node.props.get('name'), color="red")
            #face_name = TextFace("%.1f" % (relative_abundance*100), color=color)
            node.add_face(identF, column = level,  position = 'branch_right')
            #node.add_face(identF, column = level,  position = 'branch_right', collapsed_only=True)
        
        elif node.props.get(internal_prop):
            # heatmap
            redgradient = color_gradient(0.95, 0.6, 10)
            relative_abundance = float(node.props.get(internal_prop))
            color_idx = int(relative_abundance*10)
            color = redgradient[color_idx]

            identF = RectFace(width=50,height=50,text="%.1f" % (relative_abundance*100), color=color, 
            padding_x=1, padding_y=1)
            #face_name = TextFace(node.props.get('name'), color="red")
            #face_name = TextFace("%.1f" % (relative_abundance*100), color=color)
            node.add_face(identF, column = level+2,  position = 'branch_right', collapsed_only=True)

    return layout_fn
    return