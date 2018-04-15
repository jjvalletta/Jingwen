#*************************************************************************************************#
#*************************************************************************************************#
# Title:    A/B comparison functions for mice microarray data
# Author:   John Joseph Valletta
# Date:     02/09/2017
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
import numpy as np # version 1.12.1
import os
import warnings
import matplotlib.pyplot as plt # version 2.0.2
from matplotlib_venn import venn2 # venn diagram
import seaborn as sns # version 0.7.1
# "Local" libraries
import input_output as io
from config import COL, XLIM, YLIM, LWD, ALPHA

#*************************************************************************************************#
# Functions
#*************************************************************************************************# 
def genes(organs, strains, path, geneDF, desc=''):
    """
    Visually compare arbitrary genes across different conditions organ/strain A/B 
    
    =========
    organs - a list of two organs e.g ['Blood', 'Spleen']
    strains - a list of two strains e.g ['AS', 'CB']
    path - a structured dictionary of paths returned by config
    geneDF - a pandas series/data frame with column either ProbeID or Symbol
    desc - optional description string e.g cluster name
    
    Returns
    =========
    hFig - figure handle
    
    """
    # Simple data checking
    if len(organs) != 2 or len(strains) !=2 or len(geneDF) <= 0:
        warnings.warn("Arguments not compatible with function requirements!", UserWarning)
        hFig = None
    else:
        # Create figure (a4 = 11.69 x 8.27 inches)
        hFig, hAx = plt.subplots(1, 2, sharex=True, sharey=False, figsize=(11.69, 8.27/1.5))
        
        # Extract gene list to be plotted
        col = geneDF.columns[0] # pick first column by default to be col id
        geneList = list(geneDF.dropna()[col])
        for i, (organ, strain) in enumerate(zip(organs, strains)):
            # Load smooth expression data
            data = pd.read_csv(os.path.join(path['GPFit']['SmoothExprs'], organ + strain + '.csv'))
            smoothExprs = data[data[col].isin(geneList)]
            x = smoothExprs.filter(regex='\.').columns.values.astype('float')
            y = smoothExprs.filter(regex='\.').as_matrix()
            labels = list(smoothExprs[col])
            
            # Plot data
            hAx[i].set_prop_cycle('color', COL['Cycle']) # set the cycling colours
            hLines = hAx[i].plot(x, y.T)
            hAx[i].axhline(y=0, color=COL['Zero'], linestyle='--', linewidth=LWD['M'])
            
            # Set limits
            hAx[i].set_xlim(XLIM)
            hAx[i].set_ylim(YLIM)
            hAx[i].set_title('{}{} {} (#{} probes/genes)'.format(organ, strain, desc, y.shape[0]))
            hAx[i].set_ylabel('$\log_2$(FC) wrt naive')
            hAx[i].set_xlabel('Days')
            #hAx[i].set_yticks(np.arange(-6, 7, 2))
              
            # Show legend if less than 20 genes
            if len(hLines)<=20 and i==1:
                hAx[i].legend(hLines, labels, ncol=2, loc=0)
                
    return hFig

#*************************************************************************************************#
def clusters(organs, strains, path, clustIndices):
    """
    Visually compare arbitrary clusters across different conditions organ/strain A/B 
    
    =========
    organs - a list of two organs e.g ['Blood', 'Spleen']
    strains - a list of two strains e.g ['AS', 'CB']
    path - a structured dictionary of paths returned by config
    clustIndices - a list of cluster indices e.g [[1, 4, 5], [6, 7]]
    
    Returns
    =========
    hFig - figure handle
    
    """
    # Simple data checking
    if len(organs) != 2 or len(strains) !=2 or len(clustIndices) != 2:
        warnings.warn("Arguments not compatible with function requirements!", UserWarning)
        hFig = None
    else:
        # Find common probes i.e to plot in grey
        shared, left, right = shared_genes(organs, strains, path, clustIndices)
        
        # Create figure (a4 = 11.69 x 8.27 inches)
        hFig, hAx = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(16, 8.27/1.8)) 
        title = ['{}_{}{}'.format(strains[i], organs[i][:2], clustIndices[i]) for i in range(2)]
        
        # Plot shared and only genes
        cluster_subplot(organs[0], strains[0], path, shared, left, title[0], hAx[0])
        cluster_subplot(organs[1], strains[1], path, shared, right, title[1], hAx[2]) # beware 2
        
        # Plot venn diagram
        plt.sca(hAx[1])
        hVen = venn2([left.shape[0], right.shape[0], shared.shape[0]], set_colors=(COL[organs[0]], 
                     COL[organs[1]]), set_labels=("", ""), alpha=0.9)
        for label in hVen.subset_labels:
            if hasattr(label, 'set_fontsize'):
                label.set_fontsize(20)
        if shared.shape[0] > 0:
            hVen.get_patch_by_id('11').set_color(COL['Shade']) # make intersection grey
        
        # Create folder where to put results
        folderName = '{}_{}'.format(title[0], title[1])
        thisDir = io.create_folder(path['Comparison'], folderName)
        
        # Save plot and gene lists
        io.save_pdf(os.path.join(thisDir, 'timeplot.pdf'), hFig)
        shared.to_csv(os.path.join(thisDir, 'shared.csv'), index=False)
        left.to_csv(os.path.join(thisDir, 'left.csv'), index=False)
        right.to_csv(os.path.join(thisDir, 'right.csv'), index=False)
        
#*************************************************************************************************#
def shared_genes(organs, strains, path, clustIndices):
    """
    Find shared genes/probes across arbitrary clusters
    
    =========
    organs - a list of two organs e.g ['Blood', 'Spleen']
    strains - a list of two strains e.g ['AS', 'CB']
    path - a structured dictionary of paths returned by config
    clustIndices - a list of cluster indices e.g [[1, 4, 5], [6, 7]]
    
    Returns
    =========
    hFig - figure handle
    
    """
    # Extract probe lists
    geneDF = [] # data frame of combined genes/probes
    for organ, strain, clustIndex in zip(organs, strains, clustIndices):
        # Read probe and gene list
        probe = pd.read_csv(os.path.join(path['Clust']['ProbeList'], organ + strain + '.csv'))
        symbol = pd.read_csv(os.path.join(path['Clust']['GeneList'], organ + strain + '.csv'))
        name = ['{}_{}_{:02d}'.format(strain, organ[:2], i) for i in clustIndex] # cluster names
        # Create data frame from subset of 
        df = pd.DataFrame({'ProbeID': probe[name].stack(), 'Symbol': symbol[name].stack()})
        geneDF.append(df)
    
    # Find shared/left/right genes
    shared = pd.merge(geneDF[0], geneDF[1], how='inner', on=['ProbeID', 'Symbol'])
    left = pd.merge(geneDF[0], geneDF[1], how='outer', on=['ProbeID', 'Symbol'], indicator=True)\
                    .query('_merge == "left_only"').drop(['_merge'], axis=1)
    right = pd.merge(geneDF[0], geneDF[1], how='outer', on=['ProbeID', 'Symbol'], indicator=True)\
                    .query('_merge == "right_only"').drop(['_merge'], axis=1)
    
    return shared, left, right

#*************************************************************************************************#
def cluster_subplot(organ, strain, path, shared, only, title, hAx): 
    """
    Helper function to plot shared and left/right genes
    
    =========
    organ - Blood/Spleen
    strain - AS/CB
    path - a structured dictionary of paths returned by config
    shared - pandas data frame of shared ProbeID/Symbols
    only - pandas data frame of left/right ProbeID/Symbols
    title - plot title
    hAx - handle to figure axis
    
    Returns
    =========
    None - axis handle is modified
    
    """
    # Load smooth expression data
    data = pd.read_csv(os.path.join(path['GPFit']['SmoothExprs'], organ + strain + '.csv'))
    
    # Plot "only" genes first
    smoothExprs = data[data['ProbeID'].isin(only['ProbeID'])]
    x = smoothExprs.filter(regex='\.').columns.values.astype('float')
    y = smoothExprs.filter(regex='\.').as_matrix()
    hAx.plot(x, y.T, color=COL[organ], linewidth=LWD['S'], alpha=ALPHA)
    
    # Plot "shared" genes
    smoothExprs = data[data['ProbeID'].isin(shared['ProbeID'])]
    x = smoothExprs.filter(regex='\.').columns.values.astype('float')
    y = smoothExprs.filter(regex='\.').as_matrix()
    hAx.plot(x, y.T, color=COL['Shade'], linewidth=LWD['M'], alpha=ALPHA)

    # Set limits
    hAx.axhline(y=0, color=COL['Zero'], linestyle='--', linewidth=LWD['M'])
    hAx.set_xlim(XLIM)
    hAx.set_ylim(YLIM)
    hAx.set_title(title)
    hAx.set_ylabel('$\log_2$(FC) wrt naive')
    hAx.set_xlabel('Days')
        
#*************************************************************************************************#
def corr2_coeff(A, B):
    """
    Computing the correlation coefficient between two multi-dimensional arrays:
    http://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays

    Arguments
    =========
    A - a 2D numpy array (e.g N1 x D1)
    B - a 2D numpy array (e.g N2 x D2)    
    
    Returns
    =========
    Correlation matrix (N1 x N2)

    """
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);

    # Finally get corr coeff
    return np.dot(A_mA,B_mB.T)/np.sqrt(np.dot(ssA[:, None],ssB[None]))
    
#*************************************************************************************************#
def cluster_mean_pearson(organA, strainA, organB, strainB, path):
    """
    Correlate the MEAN of each cluster of organA/strainA with every cluster 
    of organB/strainB and plots a heatmap
    
    This could be improved by correlating the individual genes in each cluster 
    rather than the cluster mean and then report the average, but I think using 
    the mean profile of the cluster is reasonable for visualisation purposes  
    

    Arguments
    =========
    organ - Blood/Spleen
    strain - AS/CB
    path - dictionary with all results paths
    
    Returns
    =========
    None - Figure is saved to path['Misc']

    """
    # Load clust centres A
    data = pd.read_csv(os.path.join(path['Clust']['Centres'], 
                                    organA + strainA + ".csv"), sep=",")  
    data = data.set_index('Cluster') # name rows as cluster name
    clustCentreA = data.values
    clustNameA = data.index
    clustNameA.name = organA # change 'Cluster' to 'Spleen' or 'Blood'
    
    # Load clust centres B
    data = pd.read_csv(os.path.join(path['Clust']['Centres'], 
                                    organB + strainB + ".csv"), sep=",")  
    data = data.set_index('Cluster') # name rows as cluster name
    clustCentreB = data.values
    clustNameB = data.index
    clustNameB.name = organB # change 'Cluster' to 'Spleen' or 'Blood'
    
    # Compute correlation coefficient between pairwise cluster centres
    corrMatrix = corr2_coeff(clustCentreA, clustCentreB)
    #corrMatrix[corrMatrix<0.6] = 0
    #corrMatrix[corrMatrix>=0.6] = 1
    df = pd.DataFrame(data=corrMatrix, index=clustNameA, columns=clustNameB) 
    
    # Plot heatmap        
    # method = http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    # metric = http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
    #obj = sns.clustermap(df, method='ward', metric='euclidean', cmap='gray_r')#, vmin=-1, vmax=1)
    obj = sns.clustermap(df, method='ward', metric='euclidean', vmin=-1, vmax=1, 
                         annot=True, fmt=".2f", figsize=(12, 12))
    # Rotate ticks
    plt.setp(obj.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(obj.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # Save to file
    fileName = 'Heatmap_' + organA + strainA + '_' + organB + strainB + '.pdf' 
    filePath = os.path.join(path['Misc'], fileName)    
    obj.savefig(filePath) 
    plt.close(obj.fig) # close figure          