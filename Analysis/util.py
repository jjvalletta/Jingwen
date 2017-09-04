#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Utility functions for mice microarray data
# Author:   John Joseph Valletta
# Date:     31/08/2017
# Comment:  Big change: port all code to Python3, GPy-1.7.7, paramz-0.7.4, GPclust-0.1.0 
# Data:     Microarray Illumina Mouse WG6 v2 (45,281 probe sets representing 30,854 genes)
# From:     Jingwen Lin and Jean Langhorne, Crick Institute 
#*************************************************************************************************#
#*************************************************************************************************#

#*************************************************************************************************#
# Preamble
#*************************************************************************************************#
# Libraries
import pandas as pd # version 0.20.1
import os
from config import TOP_RANKED

#*************************************************************************************************#
# Functions
#*************************************************************************************************# 
def normalise(x):
    """
    Normalise a list or numpy array to be between 0 and 1 
    
    =========
    x - a list or numpy array
    

    =========
    normalised list/array [0, 1]
    
    """
    return (x - min(x))/(max(x) - min(x))

#*************************************************************************************************#
def top_ranked_genes(organ, strain, path, topRanked=TOP_RANKED, IDType='ProbeID'):
    """
    Find top ranked genes/probeID

    Arguments
    =========
    organ - Blood/Spleen
    strain - AS/CB
    path - a structured dictionary of paths returned by config
    topRanked - no. of top ranked genes to consider
    returnType - whether to return a set of 'ProbeID' or 'Symbol'
    
    Returns
    =========
    Set of top ranked probeID/genesymbol

    """    
    # Find top ranked genes 
    data = pd.read_csv(os.path.join(path['GPFit']['Metrics'], organ + strain + ".csv"))
    data = data.sort_values(by='rank', ascending=False) # sort by rank    

    return set(data.iloc[:topRanked][IDType])