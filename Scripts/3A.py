# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:42:07 2020

This code performs association rule mining using the python package apriori for
different threshold values given.

This code requires the following packages: pandas, networkx and efficient_apriori.
All these packages can be downloaded via pip.

Input
------
(a) MEDUSA database organism reads counts where the columns are
sample IDs and the rows are organisms.
(b) Organisms found in our Supplementary File 1 (as a pickle file) 

Output
------
Saves the food web file in the working directory.

@author: Aarthi Ravikrishnan
"""

import pickle
import os
import pandas as pd
import networkx as nx
from efficient_apriori import apriori


# Initialise with all the thresholds
pathname = os.chdir('.')
input_data_file = '../Data/CollapsedSpeciesCount.count'
cfile = pd.read_csv(input_data_file, '\t', index_col='SPID')
list_name = cfile.columns
different_thresholds = [0.01] #Append entries here in case other thresholds are to be tried 
# Loading the data and removing the empty rows
# Row sum
cfile1 = (cfile.loc[:, list_name].sum(axis=0) == 0)
row_zero_indices = cfile1.index[cfile1 == True].tolist()

# Column_sum
cfile2 = (cfile.loc[:, list_name].sum(axis=1) == 0)
col_zero_indices = cfile2.index[cfile2 == True].tolist()

# All zero indices
zero_indices = row_zero_indices + col_zero_indices
# Main counts dataframe
cfile = cfile.drop(zero_indices)

# Normalise the read counts by total count in samples and multiply by 10^8
cfile_new = cfile.copy()
total_reads = pd.DataFrame()
total_reads['sum'] = cfile_new.loc[:, list_name].sum(axis=0)
cfile_new = cfile_new.T
cfile_normalised = cfile_new.div(total_reads['sum'], axis='rows')
cfile_normalised = cfile_normalised.apply(lambda x: x*pow(10, 8))

# Discretise to get presence absence
#cfile_normalised = cfile_normalised.T
max_min_val_df = pd.DataFrame()
max_min_val_df['maximum'] = cfile_normalised.quantile(0.95) 
max_min_val_df['minimum'] = cfile_normalised.min(axis=0)
max_min_val_df['difference'] = max_min_val_df['maximum'] - \
    max_min_val_df['minimum']
    

for threshold in different_thresholds:
    print(threshold)
    max_min_val_df['thresoldingvalueonepart'] = max_min_val_df['difference'].apply(lambda x: x*threshold)
    max_min_val_df['thresoldingvalue'] = max_min_val_df['thresoldingvalueonepart'].add(max_min_val_df['minimum'])
    discrete_profiles = cfile_normalised.ge(max_min_val_df['thresoldingvalue'], axis='columns').astype(int)
    discrete_profiles_trans = discrete_profiles.copy()

    # Aprori
    organisms = list(discrete_profiles_trans.columns)
    list_discrete_profiles = []
    # xidx is the number of organisms, yidx is the sample
    [xidx, yidx] = discrete_profiles_trans.shape

    for i in range(0, xidx):
        list_discrete_profiles_tmp = []
        for j in range(0, yidx):
            if discrete_profiles_trans.values[i, j] == 1:
                list_discrete_profiles_tmp.append(organisms[j])
                list_discrete_profiles_tmp.sort()
        list_discrete_profiles.append(list_discrete_profiles_tmp)

    [itemsets, association_rules] = apriori(
        list_discrete_profiles, min_support=0.05, min_confidence=0.95, max_length=2)
    rules = list(association_rules)
    rules_rhs = filter(lambda rule: len(rule.lhs) ==
                       1 and len(rule.rhs) == 1, rules)
    all_rules = []
    for ru in rules:
        all_rules.append([ru.lhs[0], ru.rhs[0]])

    # Constructing graph to remove bidirectional edges
    G = nx.DiGraph()
    for orgs in all_rules:
        G.add_edges_from([(orgs[0], orgs[1])])
    pred = G.predecessors
    succ = G.successors
    bidirectional_edges = []
    for (u,v) in G.edges():
        if u in G[v] :
            bidirectional_edges.append((u,v))
        
    G.remove_edges_from(bidirectional_edges) 
    with open("../Data/all_organisms_non_zero.pickle", 'rb') as f:
        master_organism_list = pickle.load(f)

    # Get our list of organisms
    nodes_in_our_graph = G.nodes()
    nodes_to_retain = set(nodes_in_our_graph).intersection(
        set(master_organism_list))
    nodes_to_remove = nodes_in_our_graph - nodes_to_retain
    G.remove_nodes_from(list(nodes_to_remove))
    H = nx.transitive_reduction(G)

    # Association list based on threshold
    file_name_to_store_1 = 'food_web_' + \
        + str(threshold) + '.txt'
    file_name_to_store = os.path.join(file_name_to_store_1)
    with open(file_name_to_store, 'w') as f:
        f.write('Source' + '\t' + 'Target' + '\n')
        for eg in H.edges():
            f.write(eg[1] + '\t' + eg[0])
            f.write('\n')
    print('Done')
    print('------------------------')
