#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Top-level script to cluster mice micro-array data (computationally expensive)
# Author:   John Joseph Valletta
# Date:     30/08/2017
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
# "Local" libraries
from config import setup_folders, ORGANS, STRAINS, SEED
import gaussian_process as gp
import util

#*************************************************************************************************#
# Setup/Create folders and thus all paths (stored as a structured dictionary)
#*************************************************************************************************#
print('Setup/Create folders and thus all paths (stored as a structured dictionary)...\n')
path = setup_folders()

#*************************************************************************************************#
# Fit Gaussian Process model to every probe and call differential expression
#*************************************************************************************************#
print('Fit Gaussian Process model to every probe and call differential expression...\n')
for organ in ORGANS:
    for strain in STRAINS:
        # Load data
        data = pd.read_csv(os.path.join(path['Data'], organ + strain + ".csv"))
        yTrain = data.filter(regex='\\.')
        xTrain = np.floor(yTrain.columns.values.astype('float'))[:, None]
        yTrain = yTrain.as_matrix() # convert to numpy array
        xTestInterp = np.linspace(start=np.min(xTrain), stop=np.max(xTrain), num=49)[:, None] # for interpolation
        xTestDEG = np.unique(xTrain)[1:, None] # to evaluate DEGs - ignore day 0 
        
        # Compute gene-wise models
        print("Computing probe-wise models: %s%s...\n" % (organ, strain))        
        NProbes = yTrain.shape[0]  
        
        # Initialise arrays where to store parameters, metrics and DEGs
        params = np.empty((NProbes, 3), dtype='float') # rbfVar, rbfScale, noiseVar
        metrics = np.empty((NProbes, 5), dtype='float') # maxLogFC, SNR, score, netLogFC, rank
        smoothExprs = np.empty((yTrain.shape[0], xTestInterp.shape[0]))
        bDEGUp = np.empty((yTrain.shape[0], xTestDEG.shape[0]), dtype=bool)
        bDEGDown = np.empty((yTrain.shape[0], xTestDEG.shape[0]), dtype=bool)
        
        # Loop through all probes and fit a Gaussian Process Regression Model 
        for i in range(NProbes):
            if ((i+1)%1000==0):
                print("%i/%i\n" % ((i+1), NProbes))
            GPModel = gp.gene_model(xTrain, yTrain[i, :][:, None])        
            
            # Get posterior mean and posterior variance to create a smooth expression profile
            mu, var = GPModel.predict(xTestInterp)
            smoothExprs[i, :] = mu.T 
            
            # Get hyperparameters and metrics  
            # REMEMBER: http://python.net/~goodger/projects/pycon/2007/idiomatic/handout.html#other-languages-have-variables          
            params[i, :] = gp.get_gene_params(GPModel.copy()) # rbfVar, rbfScale, noiseVar
            metrics[i, 0:4] = gp.compute_gene_metrics(GPModel.copy()) # maxLogFC, SNR, score, netLogFC
            
            # Get posterior mean and posterior variance to check for DEG
            mu, var = GPModel.predict(xTestDEG)
            # Is logFC = 0 is within the 2sigma intervals if no and logFC>1.5, then gene is DEG
            lower = mu - 2*np.sqrt(var)
            upper = mu + 2*np.sqrt(var)
            bDEGUp[i, :] = ((lower > 0) & (upper > 0) & (mu > 1.5)).T
            bDEGDown[i, :] = ((lower < 0) & (upper < 0) & (mu < -1.5)).T
        
        # Compute rank for time-series differential expression as a whole (i.e not individual days)
        metrics[:, 4] = util.normalise(metrics[:, 0]) + util.normalise(metrics[:, 1]) + util.normalise(metrics[:, 2])
        
        # Write all results to disk (a helper function to shorten code in this loop)
        gp.write_gene_model_results(data, organ, strain, path, params, metrics, smoothExprs, 
                                    bDEGUp, bDEGDown, xTestInterp, xTestDEG)
        
#*************************************************************************************************#
# Cluster time-profiles using a Mixtures of Hierarchical Gaussian Processes model
#*************************************************************************************************#
print('Cluster time-profiles using a Mixtures of Hierarchical Gaussian Processes model...\n')
for organ in ORGANS:
    for strain in STRAINS:
        probesToCluster = util.top_ranked_genes(organ, strain, path)
        gp.cluster_with_MOHGP(probesToCluster, organ, strain, 2, 1e-100, path, seed=SEED)