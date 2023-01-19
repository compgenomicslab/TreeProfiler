from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace, LegendFace

def get_piechartface(node, prop, colour_dict=None, radius=20, tooltip=None):
    piechart_data = []
    counter_props = node.props.get(prop).split('||')
    for counter_prop in counter_props:
        k, v = counter_prop.split('--')
        piechart_data.append([k,float(v),colour_dict[k],None])

    if piechart_data:
        tooltip = ""
        if node.name:
            tooltip += f'<b>{node.name}</b><br>'
        if prop:
            tooltip += f'<br>{prop}: {piechart_data}<br>' # {counter_props}
            
        piechart_face = PieChartFace(radius=radius, data=piechart_data, padding_x=5, tooltip=tooltip)
        return piechart_face
    else:
        return None