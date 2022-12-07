from ete4 import Tree

import sys,os
import csv

t = Tree(sys.argv[1])

# output 
with open('../claudia_tree/annotated_concatenate_pg3_faa_ft.tsv','w') as f:
    fnames = ['name', 'GC', 'GCA', 'GTDB-tk', 'aquatic_habitat', 'host_associated', 'size', 'soil_habitat']
    writer = csv.DictWriter(f, fieldnames=fnames,delimiter='\t')
    writer.writeheader()
    for node in t.traverse('postorder'):
        node2annotations = {}
        if node.is_leaf():
            node.name = node.name +'.1'
            old_GCA = node.props.get('GCA')
            node.add_prop('GCA', 'GB_'+old_GCA+'.1') 
            node2annotations = node.props
            del node2annotations["dist"]
            del node2annotations["support"]
            writer.writerow(node2annotations)

t.write(outfile='../claudia_tree/concatenate_pg3_faa_ft.nw', properties=None)