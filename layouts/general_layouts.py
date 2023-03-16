from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace, LegendFace, RectFace

import numpy as np
from distutils.util import strtobool
import matplotlib as mpl

from utils import to_code, call, counter_call
from utils import check_nan

def get_piechartface(node, prop, colour_dict=None, radius=40, tooltip=None):
    pair_delimiter = "--"
    item_seperator = "||"
    piechart_data = []
    counter_props = node.props.get(prop).split(item_seperator)
    for counter_prop in counter_props:
        k, v = counter_prop.split(pair_delimiter)
        piechart_data.append([k,float(v),colour_dict.get(k,None),None])
        
    if piechart_data:
        # tooltip = ""
        # if node.name:
        #     tooltip += f'<b>{node.name}</b><br>'
        # if prop:
        #     tooltip += f'<br>{prop}: {piechart_data}<br>' # {counter_props}
        piechart_face = PieChartFace(radius=radius, data=piechart_data, padding_x=5, tooltip=tooltip)
        
        return piechart_face
    else:
        return None

def color_gradient(c1, c2, mix=0):
    """ Fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1) """
    # https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def get_heatmapface(node, prop, color, tooltip=None, width=70, height=50, padding_x=2, padding_y=2):
    counter_props = node.props.get(prop).split('||')
    
    total = 0
    positive = 0
    for counter_prop in counter_props:
        k, v = counter_prop.split('--')
        if not check_nan(k):
            if strtobool(k):
                positive = float(v)
        total += float(v)
    total = int(total)
    ratio = positive / total
    if ratio < 0.05 and ratio != 0: # show minimum color for too low
        ratio = 0.05
    c1 = 'white'
    c2 = color
    gradient_color = color_gradient(c1, c2, mix=ratio)
    text = f"{positive} / {total}"
    # gradientFace = RectFace(width=100,height=50,text="%.1f" % (ratio*100), color=gradient_color, 
    #         padding_x=1, padding_y=1)

    if not tooltip:
        if node.name:
            tooltip = f'<b>{node.name}</b><br>'
        else:
            tooltip = ''
        if prop:
            tooltip += '<br>{}: {} / {} <br>'.format(prop, positive, total)

    gradientFace = RectFace(width=width, height=height, 
                            #text=text, 
                            color=gradient_color, 
                            padding_x=padding_x, padding_y=padding_y, tooltip=tooltip)
    return gradientFace