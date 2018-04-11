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
# "Local" libraries
import input_output as io
from config import TOP_RANKED, COL, LWD

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

#*************************************************************************************************#
def pca_plot(organ, strain, path, bLegend=False):
    """
    Performs and plots PCA for organ and strain. The plot is saved to path['Misc']
         
    Arguments
    =========
    organ - Blood/Spleen
    strain - AS/CB
    path - a structured dictionary of paths returned by config
    bLegend - a Boolean whether to show/hide legend
  
    Returns
    =========
    None - PCA plot and PCA weights figure is saved as a .pdf to path['Misc']
    
    """ 
    # Read data
    # Note: PCA for Log2FC or Raw is the same as we're only subtracting a constant
    data = pd.read_csv(os.path.join(path['Data'], organ + strain + '.csv'))
    #probeID = top_ranked_genes(organ, strain, path) # if we only want to consider TOP_RANKED
    #data = data[data['ProbeID'].isin(probeID)] # subset data
    x = np.floor(data.filter(regex='\.').columns.values.astype('float')).astype('int')
    y = data.filter(regex='\.').as_matrix().T # samples/days=rows; features/genes=columns
    
    # Colour dictonary and col list
    colDict = {}
    for i, value in enumerate(np.unique(x)):
        colDict[value] = COL['PCA'][i]
    myCol = [colDict[value] for value in x]
    
    # Perform PCA        
    pcaOut = PCA(n_components=2).fit(y)
    PC = [np.round(pcaOut.explained_variance_ratio_[i]*100) for i in range(pcaOut.n_components)]
    embedding = pcaOut.transform(y)
    
    # Plot weights
    hFig, hAx = plt.subplots(2, sharex=True, figsize=(12, 12))
    for i in range(len(PC)):
        hAx[i].plot(np.sort(pcaOut.components_[i, :]**2)[::-1]*100)
        hAx[i].set_title('PC{}: {:.0f}%'.format(i, PC[i]))
        hAx[i].set_ylabel('PCA Weight ($\sum$ weights = 100)')
        hAx[i].set_xlabel('No. of genes')
    io.save_pdf(os.path.join(path['Misc'], 'PCAWeights' + organ + strain + '.pdf'), hFig)
    
    # Create PCA figure
    hFig = plt.figure(figsize=(8, 8)); hAx = plt.gca()
    hAx.scatter(embedding[:, 0], embedding[:, 1], linewidth=LWD['M'], edgecolor='black', 
                c=myCol, s=300)
    
    # Set x/y limits
    xMin, yMin = np.min(embedding, axis=0)
    xMax, yMax = np.max(embedding, axis=0)
    hAx.set_xlim(xMin-(0.1*abs(xMin)), xMax+(0.1*abs(xMax))) 
    hAx.set_ylim(yMin-(0.1*abs(yMin)), yMax+(0.1*abs(yMax)))  
    hAx.set_xlabel('PC1: {:.0f}%'.format(PC[0]), fontsize=24)
    hAx.set_ylabel('PC2: {:.0f}%'.format(PC[1]), fontsize=24)
    plt.tick_params(labelsize=16) # change tixk font 
    plt.tight_layout()
    #plt.axes().set_aspect('equal', 'datalim')
    
    # Add legend if wanted
    if bLegend:
        col = []
        lab = []
        for key, value in colDict.items():
            lab.append('Day {}'.format(key))
            col.append(value)
        # 'Hack' to show legend
        faff = [hAx.scatter([], [], edgecolor='black', c=c) for c in col]        
        hAx.legend(faff, lab, fontsize=20, ncol=1, markerscale=2, 
                   handlelength=0.2, borderpad=-0.4, loc=0, frameon=False)
    
    # Write to disk
    io.save_pdf(os.path.join(path['Misc'], "PCA" + organ + strain + ".pdf"), hFig) 

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
    hFig = plt.figure(figsize=(16/1.2, 9/2))    
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
    