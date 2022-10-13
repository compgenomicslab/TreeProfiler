#!/usr/bin/env python
from ete4 import GTDBTaxa, NCBITaxa
from ete4 import Tree, PhyloTree

import os, sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict

"""
    input: 
        1) metadata table from mOTUs/gtdb-tk/general tsv file
        2) user's tree or pre taxnomic annotated tree from GTDB/NCBI
    output:
        1) annotated tree
        2) annotated data
"""

input_file = sys.argv[1]
taxa_annotation=True
taxa_column= -1
taxa_separator=';'
rank_limit='species'



short2levels = {
            "d__":"kingdom",
            "p__":"phylum",
            "c__":"class",
            "o__":"order",
            "f__":"family",
            "g__":"genus",
            "s__":"species",
        }

levels2short = {
            "kingdom":"d__",
            "phylum":"p__",
            "class":"c__",
            "order":"o__",
            "family":"f__",
            "genus":"g__",
            "species":"s__",
        }

# GTDB ids to species/order/genus/family/order/class/phylum/kingdom
def parse_df(input_file, taxa_annotation=True, taxa_column=-1, taxa_separator=';'):
    df = pd.read_csv(input_file, sep='\t')
    origin_headers = list(df.columns)
    if taxa_annotation:
        taxon = df.iloc[:, taxa_column] 
        taxa_df = taxon.str.split(taxa_separator, expand=True)
        taxa_levels = list(levels2short.keys())
        for idx, name in enumerate(taxa_levels):
            df[name] = taxa_df[idx]
        return origin_headers, df
    else:
        add_taxon(metadata)

    

# def parse_csv(input_file):
#     metadata = []
#     treenodes = []
#     header = []
#     with open(input_file, newline='') as csvfile:
#         reader = csv.reader(csvfile, delimiter='\t')
#         #next(reader, None)  # skip the headers
#         count = 0
#         for line in reader:
#             if count == 0:
#                 headers = line
#                 count += 1
#             else:
#                 metadata.append(line)
#                 treenodes.append(line[0])
    
#     columns = defaultdict(list) # each value in each column is appended to a list
#     with open(input_file) as f:
#         reader = csv.DictReader(f, delimiter="\t") # read rows into a dictionary format
#         for row in reader: # read a row as {column1: value1, column2: value2,...}
#             for (k,v) in row.items(): # go over each column name and value 
#                 columns[k].append(v) # append the value into the appropriate list
#                                     # based on column name k

#     return metadata, columns

# def parse_csv(input_file):
#     metatable = []
#     tsv_file = open(input_file)
#     read_tsv = csv.DictReader(tsv_file, delimiter="\t")

#     for row in read_tsv:
#         metatable.append(row)
#     tsv_file.close()

#     columns = defaultdict(list) # each value in each column is appended to a list
#     with open(input_file) as f:
#         reader = csv.DictReader(f, delimiter="\t") # read rows into a dictionary format
#         for row in reader: # read a row as {column1: value1, column2: value2,...}
#             for (k,v) in row.items(): # go over each column name and value 
#                 columns[k].append(v) # append the value into the appropriate list
#                                     # based on column name k
#     return metatable, columns

def extract_lineage(treenodes, rank_limit='subspecies'):
    gtdb = GTDBTaxa()
    gtdb_lineages = {}
    for d in gtdb.get_name_lineage([treenodes]):
        lineage = [i for i in list(d.values())[0] ] # until species

        gtdb_lineages[list(d.keys())[0]] = ';'.join(lineage)

    return gtdb_lineages

# def cut_taxa(taxa, rank_limit='subspecies'):
#     if taxa.startswith('root'):
#         return False
#     else:
#         try:
#             prefix = taxa[:3]
#             if short2levels[prefix] == rank_limit:
#                 return False
#             elif prefix in short2levels.keys():
#                 return True
#             else:
#                 return False

#         except IndexError:
#             return False


def add_taxon(metadata, name_column=0):
    #headers, metadata, treenodes = parse_csv(input_file)
    treenodes = extract(metadata, name_column)
    gtdb_lineages = extract_lineage(treenodes)
    for m in metadata:
        m.append(gtdb_lineages[m[0]])
    return metadata

##### if it's taxa #############
def convert_by_taxa(origin_headers, metadata_df, level='species'):
    taxa = list(levels2short.keys())
    new_metadata_df = pd.DataFrame()
    numeric_headers = list(df.select_dtypes(include=np.number).columns)
    text_headers = list(filter(lambda v: v not in taxa, df.select_dtypes(exclude=np.number).columns))[1:] # filter not the taxa list

    collapsed_headers = 'collapse_leaves'
    g1 = df.groupby(level)[level].count()
    new_metadata_df[level] = g1.index
    new_metadata_df[collapsed_headers] = g1.values

    ##### properties
    grouped_df = df.groupby(level)[numeric_headers]
    describe_df = grouped_df.describe()

    ##### numeric ################
    for header in numeric_headers:
        #print(describe_df[header].mean())
        header_sum = header + '_sum'
        header_mean = header + '_mean'
        new_metadata_df[header_sum] = grouped_df.sum()[header].values
        new_metadata_df[header_mean] = describe_df[header]['mean'].values

    ################ taxa ##########################
    taxa_header = origin_headers[taxa_column]
    new_metadata_df[taxa_header] = new_metadata_df[level].apply(lambda x: list(extract_lineage(x, rank_limit=level).values())[0])

    return new_metadata_df
    
########## utils ##################################
 
def extract(lst, pos):
    return [item[pos] for item in lst]

###################################################
# metadata, columns = parse_csv(input_file)
#origin_headers, df = parse_df(input_file, taxa_column=1, taxa_separator='|')

# headers = list(columns.keys())

# if 'Taxon' not in headers:
#     metadata = add_taxon(metadata)
#     columns['Taxon'] = extract(metadata,-1)
#     headers.append('Taxon')


############Write it out##############

# sys.stdout.write('\t'.join(headers)+'\n')
# for m in metadata:
#     #m.append(gtdb_lineages[m[0]])
#     sys.stdout.write('\t'.join(m)+'\n')


###########get taxa tree ######################
#tree = convert_by_taxa(headers, metadata, 'species', 6)
#print(tree)

origin_headers, df = parse_df(input_file, taxa_column=taxa_column, taxa_separator=taxa_separator)
output_df  = convert_by_taxa(origin_headers, df, rank_limit)
output_df.to_csv('demo/test_species', sep='\t')