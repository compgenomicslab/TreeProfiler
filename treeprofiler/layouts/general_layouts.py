from __future__ import annotations
import Bio
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
import numpy as np
from distutils.util import strtobool
import matplotlib as mpl
import matplotlib.colors as mcolors
from itertools import chain
import math

from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace, LegendFace, RectFace
from ete4.smartview.renderer.draw_helpers import *

from treeprofiler.src.utils import to_code, call, counter_call, check_nan
from treeprofiler.src import utils

Box = namedtuple('Box', 'x y dx dy')  # corner and size of a 2D shape

def get_piechartface(node, prop, color_dict=None, radius=20, tooltip=None):
    pair_delimiter = "--"
    item_seperator = "||"
    piechart_data = []
    counter_props = node.props.get(prop).split(item_seperator)
    for counter_prop in counter_props:
        k, v = counter_prop.split(pair_delimiter)
        piechart_data.append([k,float(v),color_dict.get(k,None),None])
        
    if piechart_data:
        piechart_face = PieChartFace(radius=radius, data=piechart_data, padding_x=5, tooltip=tooltip)
        
        return piechart_face
    else:
        return None

def get_aggregated_heatmapface(node, prop, min_color="#EBEBEB", max_color="#971919", tooltip=None,
                               width=70, height=None, padding_x=1, padding_y=0, count_missing=True, max_count=0):
    counter_props = node.props.get(prop).split('||')
    total = 0
    positive = 0
    for counter_prop in counter_props:
        k, v = counter_prop.split('--')
        if count_missing:
            if not check_nan(k):
                if strtobool(k):
                    positive += float(v)
            total += float(v)  # Consider missing data in total
        else:
            if not check_nan(k):
                total += float(v)  # Doesn't consider missing data in total
                if strtobool(k):
                    positive += float(v)
    
    total = int(total)
    # ratio = positive / total if total != 0 else 0
    # if ratio < 0.05 and ratio != 0:  # Show minimum color for too low
    #     ratio = 0.05
    
    # Adjust the maximum color based on 'total' to simulate darkening
    adjusted_max_color = utils.make_color_darker_scaled(max_color, positive, max_count, base=10, scale_factor=10)
    #adjusted_max_color = make_color_darker(max_color, darkening_factor=0.01)  # Example factor
    #gradient_color = color_gradient(min_color, adjusted_max_color, mix=ratio)

    if not tooltip:
        tooltip = f'<b>{node.name}</b><br>' if node.name else ''
        if prop:
            tooltip += f'<br>{prop}: {positive} / {total} <br>'
    if positive == 0:
        aggregateFace = RectFace(width=width, text=int(positive), height=height, color=min_color, padding_x=padding_x, padding_y=padding_y, tooltip=tooltip)
    else:
        aggregateFace = RectFace(width=width, text=int(positive), height=height, color=adjusted_max_color, padding_x=padding_x, padding_y=padding_y, tooltip=tooltip)
    return aggregateFace

def get_heatmapface(node, prop, min_color="#EBEBEB", max_color="#971919", tooltip=None, width=70, height=None, padding_x=1, padding_y=0, count_missing=True, reverse=False):
    counter_props = node.props.get(prop).split('||')
    total = 0
    positive = 0
    for counter_prop in counter_props:
        k, v = counter_prop.split('--')
        if count_missing:
            if not check_nan(k):
                if strtobool(k):
                    positive = float(v)
            total += float(v) # here consider missing data in total
        else:
            if not check_nan(k):
                total += float(v) # here doesn't consider missing data in total
                if strtobool(k):
                    positive = float(v)
            
    total = int(total)
    if total != 0:
        ratio = positive / total
    else:
        ratio = 0
    
    if reverse:
        ratio = 1 - ratio

    if ratio < 0.05 and ratio != 0: # show minimum color for too low
        ratio = 0.05

    c1 = min_color
    c2 = max_color
    gradient_color = utils.color_gradient(c1, c2, mix=ratio)
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


SeqRecord = Bio.SeqRecord.SeqRecord
def get_consensus_seq(filename: Path | str, threshold=0.7) -> SeqRecord:
    #https://stackoverflow.com/questions/73702044/how-to-get-a-consensus-of-multiple-sequence-alignments-using-biopython
    common_alignment = MultipleSeqAlignment(
        chain(*AlignIO.parse(filename, "fasta"))
    )
    summary = SummaryInfo(common_alignment)
    consensus = summary.dumb_consensus(threshold, "-")
    return consensus

def get_stackedbarface(node, prop, color_dict=None, width=70, height=None, padding_x=1, padding_y=0, tooltip=None):
    pair_delimiter = "--"
    item_seperator = "||"
    stackedbar_data = []
    absence_color = "#EBEBEB"
    counter_props = node.props.get(prop).split(item_seperator)
    tooltip = ""
    total = 0
    
    for counter_prop in counter_props:    
        k, v = counter_prop.split(pair_delimiter)
        if v:
            total += float(v)
        stackedbar_data.append([k,float(v),color_dict.get(k,absence_color),None])
        
    if stackedbar_data:
        tooltip = ""
        if node.name:
            tooltip += f'<b>{node.name}</b><br>'
        
        if counter_props:
            for counter_prop in counter_props:
                k, v = counter_prop.split(pair_delimiter)
                tooltip += f'<b>{k}</b>:  {v}/{int(total)}<br>'

        stackedbar_face = StackedBarFace(width=width, height=None, data=stackedbar_data, padding_x=padding_x, padding_y=padding_y, tooltip=tooltip)
        
        return stackedbar_face
    else:
        return None



class StackedBarFace(RectFace):
    def __init__(self, width, height, data=None, name="", opacity=0.7, 
                min_fsize=6, max_fsize=15, ftype='sans-serif',
                padding_x=1, padding_y=0, tooltip=None):
        
        RectFace.__init__(self, width=width, height=height, name=name, color=None, 
            min_fsize=min_fsize, max_fsize=max_fsize, 
            padding_x=padding_x, padding_y=padding_y, tooltip=tooltip)
        
        self.width = width
        self.height = height

        # data = [ [name, value, color, tooltip], ... ]
        # self.data = [
        #     ['first', 10, 'red', None],
        #     ['second', 40, 'blue', None],
        #     ['green', 50, 'green', None]
        #     ]
        self.data = data
        

    def __name__(self):
        return "StackedBarFace"

    def draw(self, drawer):

        # Draw RectFace if only one datum

        if len(self.data) == 1:
            self.color = self.data[0][2]
            yield from RectFace.draw(self, drawer)

        else:
            total_value = sum(d[1] for d in self.data)
            start_x, start_y, dx, dy = self._box        
            
            for i in range(len(self.data)):
                i_value = self.data[i][1]
                color = self.data[i][2]
                
                if i > 0:
                    start_x += new_dx # start with where last segment ends
                
                new_dx = i_value/total_value * dx # width of segment
                
                self._box = Box(start_x,start_y,new_dx,dy)
                style = { 'fill': color }
                yield draw_rect(self._box,
                self.name,
                style=style,
                tooltip=self.tooltip)

# def color_gradient(c1, c2, mix=0):
#     """ Fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1) """
#     # https://stackoverflow.com/questions/25668828/how-to-create-colour-gradient-in-python
#     c1 = np.array(mpl.colors.to_rgb(c1))
#     c2 = np.array(mpl.colors.to_rgb(c2))
#     return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

# def make_color_darker_log(hex_color, total, base=10):
#     """Darkens the hex color based on a logarithmic scale of the total."""
#     # Calculate darkening factor using a logarithmic scale
#     darkening_factor = math.log(1 + total, base) / 50  # Adjust base and divisor as needed
#     return make_color_darker(hex_color, darkening_factor)
    
# def make_color_darker(hex_color, darkening_factor):
#     """Darkens the hex color by a factor. Simplified version for illustration."""
#     # Simple darkening logic for demonstration
#     c = mcolors.hex2color(hex_color)  # Convert hex to RGB
#     darker_c = [max(0, x - darkening_factor) for x in c]  # Darken color
#     return mcolors.to_hex(darker_c)

# def make_color_darker_scaled(hex_color, positive, maximum, base=10, scale_factor=10, min_darkness=0.6):
#     """
#     Darkens the hex color based on the positive count, maximum count, and a scaling factor.
    
#     :param hex_color: The original color in hex format.
#     :param positive: The current count.
#     :param maximum: The maximum count achievable, corresponding to the darkest color.
#     :param base: The base for the logarithmic calculation, affecting darkening speed.
#     :param scale_factor: Factor indicating how much darker the color can get at the maximum count.
#     :param min_darkness: The minimum darkness level allowed.
#     :return: The darkened hex color.
#     """
#     if positive > maximum:
#         raise ValueError("Positive count cannot exceed the maximum specified.")
    
#     # Calculate the normalized position of 'positive' between 0 and 'maximum'
#     normalized_position = positive / maximum if maximum != 0 else 0
    
#     # Calculate the logarithmic scale position
#     log_position = math.log(1 + normalized_position * (scale_factor - 1), base) / math.log(scale_factor, base)
    
#     # Ensure the log_position respects the min_darkness threshold
#     if log_position >= min_darkness:
#         log_position = min_darkness

#     # Convert hex to RGB
#     rgb = mcolors.hex2color(hex_color)
    
#     # Apply the darkening based on log_position
#     darkened_rgb = [(1 - log_position) * channel for channel in rgb]
    
#     return mcolors.to_hex(darkened_rgb)
