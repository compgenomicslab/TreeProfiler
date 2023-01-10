from ete4.smartview import TreeStyle, NodeStyle, TreeLayout, PieChartFace

def get_piechartface(node, prop, colour_dict=None, radius=2000,):
    piechart_data = []
    counter_props = node.props.get(prop).split('||')
    for counter_prop in counter_props:
        k, v = counter_prop.split('--')
        piechart_data.append([k,float(v),colour_dict[k],None])

    if piechart_data:
        piechart_face = PieChartFace(radius=radius, data=piechart_data, padding_x=5)
        return piechart_face
    else:
        return None