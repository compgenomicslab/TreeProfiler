#!/usr/bin/env python
from ete4 import GTDBTaxa, NCBITaxa
from ete4 import Tree, PhyloTree
import os, sys
import csv

"""
    input: 
        1) metadata table from mOTUs/gtdb-tk/general tsv file
        2) user's tree or pre taxnomic annotated tree from GTDB/NCBI
    output:
        1) annotated tree
        2) annotated data
"""

input_file = sys.argv[1]

# GTDB ids to species/order/genus/family/order/class/phylum/kingdom
def parse_csv(input_file):
    metadata = []
    treenodes = []
    header = []
    with open(input_file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        #next(reader, None)  # skip the headers
        count = 0
        for line in reader:
            if count == 0:
                headers = line
                count += 1
            else:
                metadata.append(line)
                treenodes.append(line[0])
    return headers, metadata, treenodes

def extract_lineage(treenodes):
    gtdb = GTDBTaxa()
    gtdb_lineages = {}
    for d in gtdb.get_name_lineage(treenodes):
        #print(d.values())
        #print(';'.join(list(d.values())[0][1:]))
        lineage = [i for i in list(d.values())[0] if cut_subspecies(i)]
        gtdb_lineages[list(d.keys())[0]] = ';'.join(lineage)

    return gtdb_lineages

def cut_subspecies(taxa):
    levels = [  
                "d__",
                "p__",
                "c__",
                "o__",
                "f__",
                "g__",
                "s__",
            ]
    try:
        prefix = taxa[:3]
        if prefix in levels:
            return True
        else:
            return False
    except IndexError:
        return False


def add_taxon(input_file):
    headers, metadata, treenodes = parse_csv(input_file)
    gtdb_lineages = extract_lineage(treenodes)
    for m in metadata:
        m.append(gtdb_lineages[m[0]])
    return metadata


# headers.append('Taxon')
# # headers = ['ID', 'fraction_uncultivated_cultivated','abundance', 'lineage']
# sys.stdout.write('\t'.join(headers)+'\n')
# for m in metadata:
#     m.append(gtdb_lineages[m[0]])
#     sys.stdout.write('\t'.join(m)+'\n')

print(add_taxon(input_file))
