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
from sklearn.decomposition import PCA # version 0.18.1
import pandas as pd # version 0.20.1
import numpy as np # version 1.12.1
import matplotlib.pyplot as plt # version 2.0.2
import os
from functools import reduce
from collections import OrderedDict # important for consistent legend
# "Local" libraries
import input_output as io
from config import TOP_RANKED, COL, LWD, STRAINS, ORGANS

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
    IDType - whether to return a set of 'ProbeID' or 'Symbol'
    
    Returns
    =========
    Set of top ranked probeID/genesymbol

    """    
    # Find top ranked UNIQUE genes
    # For multiple probes that point to the same gene symbol retain the 
    # probe that had the highest ranking
    data = pd.read_csv(os.path.join(path['GPFit']['Metrics'], organ + strain + ".csv"))   
    data = data.sort_values(by='rank', ascending=False).drop_duplicates(subset=['Symbol'], keep='first')

    return set(data.iloc[:topRanked][IDType])

#*************************************************************************************************#
def _pca_plot(obj, embedding, markerList=None, colList=None, legDict=None, bLegend=False):
    """
    Low-level PCA plot routine. 
    
    Note: the ORDER of markerList and colList HAS to match the order of the 
    matrix used when creating the PCA object
         
    Arguments
    =========
    obj - sklearn.decomposition.pca.PCA object
    embedding - PCA(n_components=2).fit(y)
    markerList - list of marker strings ['o', '^', etc..] ORDER important
    colList - list of colours ['r', 'g', 'g', ...] ORDER important
    legDict - dictionary of legend name and colour e.g ['Day1': 'black, 'Day2': 'blue']
    bLegend - a Boolean whether to show/hide legend
  
    Returns
    =========
    (hWt, hPCA) - figure handles to weights and PCA figure
    
    """ 
    # Create markers and cols if None
    if markerList is None:
        markerList = ['o']*obj.n_samples_
    if colList is None:
        colList = ['grey']*obj.n_samples_
    
    # Extract explained variance and transformed variables 
    PC = [np.round(obj.explained_variance_ratio_[i]*100) for i in range(obj.n_components)]
    
    # Plot weights
    hWt, hAx = plt.subplots(2, sharex=True, figsize=(12, 12))
    for i in range(len(PC)):
        hAx[i].plot(np.sort(obj.components_[i, :]**2)[::-1]*100)
        hAx[i].set_title('PC{}: {:.0f}%'.format(i+1, PC[i]))
        hAx[i].set_ylabel('PCA Weight ($\sum$ weights = 1)')
        hAx[i].set_xlabel('No. of genes')
        hAx[i].axhline(y=0, color=COL['Zero'], linestyle='--', linewidth=LWD['M'])


    # Create PCA figure
    # Note: marker cannot be a list so have to loop over unique markers (argh!)
    hPCA = plt.figure(figsize=(8, 8)); hAx = plt.gca()
    markers = np.array(markerList)
    cols = np.array(colList)
    for m in np.unique(markers):
        hAx.scatter(embedding[markers==m, 0], embedding[markers==m, 1], linewidth=LWD['M'], 
                    edgecolor='black', c=cols[markers==m], s=300, marker=m)

    # Set x/y limits
    xMin, yMin = np.min(embedding, axis=0)
    xMax, yMax = np.max(embedding, axis=0)
    hAx.set_xlim(xMin-(0.1*abs(xMin)), xMax+(0.1*abs(xMax))) 
    hAx.set_ylim(yMin-(0.1*abs(yMin)), yMax+(0.1*abs(yMax)))  
    hAx.set_xlabel('PC1: {:.0f}%'.format(PC[0]), fontsize=24)
    hAx.set_ylabel('PC2: {:.0f}%'.format(PC[1]), fontsize=24)
    plt.tick_params(labelsize=16) # change tixk font 
    plt.tight_layout()
    plt.axes().set_aspect('equal', 'datalim')
    
    # Add legend if wanted
    if bLegend and legDict is not None:
        faff = [] # 'Hack' to show legend 
        for key, value in legDict.items():
            faff.append(hAx.scatter([], [], edgecolor='black', c=value[0], marker=value[1]))
        hAx.legend(faff, legDict.keys(), fontsize=20, ncol=1, markerscale=2, 
                   handlelength=0.2, borderpad=-0.4, loc=0, frameon=False)
        
    return hWt, hPCA

#*************************************************************************************************#
def pca_plot(organ, strain, path, topRanked=None, bLegend=False):
    """
    Performs and plots PCA for organ(s) and strain.
         
    Arguments
    =========
    organ - a list with one or two organs e.g ['Spleen'] or ['Blood', 'Spleen']
    strain - single strain AS/CB (string)
    path - a structured dictionary of paths returned by config
    topRanked - how many genes to consider (int) or None to consider everything
    bLegend - a Boolean whether to show/hide legend
  
    Returns
    =========
    (hWt, hPCA) - figure handles to weights and PCA figure
    
    """ 
    # Assert inputs are ok
    set(organ).issubset(ORGANS), '{} not a valid organ! Use {}'.format(organ, ORGANS)
    assert strain in STRAINS, '{} not a valid strain! Use {}'.format(strain, STRAINS)
    
    # Set marker dictionary
    marker = {'Bl': 'o', 'Sp': '^'}
        
    # Read data
    df = [] # list of data frames
    for org in organ:
        temp = pd.read_csv(os.path.join(path['Data'], org + strain + '.csv'))
        temp = temp.set_index('ProbeID')
        temp = temp.drop(labels='Symbol', axis=1)
        temp = temp.add_prefix(org[:2] + '_')
        if topRanked is not None:
            probeID = top_ranked_genes(org, strain, path)
            temp = temp.loc[probeID]
        df.append(temp)
        
    # Merge into single data frame
    data = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), df)
    del(df)
    y = data.filter(regex='\.').as_matrix().T # samples/days=rows; features/genes=columns
    
    # Perform PCA
    obj = PCA(n_components=2).fit(y)
    embedding = obj.transform(y)
        
    # Create marker + colour list; preserve ORDER
    # And legend dictionary
    legDict = OrderedDict()
    markerList = []
    colList = []
    for col in data.columns:
        # Marker+col list
        org = col.split('_')[0]
        day = int(col.split('_')[1].split('.')[0])
        markerList.append(marker[org])
        colList.append(COL['PCA'][int(day/2)])
        # Legend 
        x = col.split('.')[0] # ignore replicate
        if x not in legDict:
            legDict[x] = [COL['PCA'][int(day/2)], marker[org]]        
            
    # Return figure handles
    return _pca_plot(obj, embedding, markerList, colList, legDict, bLegend)

#*************************************************************************************************#
def cluster_membership_boxplot(organ, strain, path):
    """
    Boxplot of posterior probability that a gene pertains to that cluster

    Arguments
    =========   
    organ - Blood/Spleen
    strain - AS/CB
    path - dictionary with all results paths
    
    Returns
    =========
    None - Figure is saved to 'Misc' folder
    """
    # Read model fit
    fit = io.load_pickle(os.path.join(path['Clust']['Model'], organ + strain + ".pickle"))    
    fit.reorder()
    labels = np.argmax(fit.phi, axis=1) # 0, 1, 2, etc.
    # Read cluster names e.g 0 --> AS_Bl_01
    data = pd.read_csv(os.path.join(path['Clust']['Centres'], organ + strain + ".csv"), sep=",")      
    
    # Create a list of cluster membership
    membership = [] 
    clustName = []
    for label in np.unique(labels):
        membership.append(fit.phi[labels==label, label])
        clustName.append(data['Cluster'][label])

    # Create boxplot
    hFig = plt.figure(figsize=(16, 5))    
    hAx = hFig.gca()    
    hAx.boxplot(membership, labels=clustName, patch_artist=True, notch=True,
                boxprops={"facecolor": COL[organ]},
                medianprops={"color": "black", "linewidth": 3},
                whiskerprops={"color": "black"})
    hAx.set_ylabel("Posterior probability of assigning gene to cluster")
    hAx.set_ylim(0, 1)
    
    # Save figure to file
    filePath = os.path.join(path['Misc'], "ClusterMembership" + organ + strain + '.pdf')    
    io.save_pdf(filePath, hFig)
    