#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Top-level script to make list of genes amenable for Cytoscape/Reactome analysis
# Author:   John Joseph Valletta
# Date:     12/09/2017
# Comment:  Big change: port all code to Python3, GPy-1.7.7, paramz-0.7.4, GPclust-0.1.0 
# Data:     Microarray Illumina Mouse WG6 v2 (45,281 probe sets representing 30,854 genes)
# From:     Jingwen Lin and Jean Langhorne, Crick Institute 
#*************************************************************************************************#
#*************************************************************************************************#

#*************************************************************************************************#
# Preamble
#*************************************************************************************************#
# Libraries
import os
import pandas as pd # version 0.20.1
from config import setup_folders, ORGANS, STRAINS

#*************************************************************************************************#
# Setup/Create folders and thus all paths (stored as a structured dictionary)
#*************************************************************************************************#
path = setup_folders()

#*************************************************************************************************#
# *_list = list of gene symbols
#*************************************************************************************************#
for organ in ORGANS:
    for strain in STRAINS:
        data = pd.read_csv(os.path.join(path['Clust']['GeneList'], organ + strain + '.csv')) 
        for clust in data.keys():
            df = pd.DataFrame({'name': data[clust].dropna().str.upper()})
            filePath = os.path.join(path['PathwayAnalysis']['Reactome'], clust + '_list.txt')
            df.to_csv(filePath, sep='\t', index=False)