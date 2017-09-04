#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Gaussian process models for mice microarray data
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
import GPy # Gaussian Process module version 1.7.7
import GPclust # version 0.1.0
import numpy as np # version 1.12.1
import pandas as pd # version 0.20.1
import os
import matplotlib.pyplot as plt # version '2.0.2'
# "Local" libraries
import config
import input_output as io

#*************************************************************************************************#
# Functions
#*************************************************************************************************# 
def gene_model(x, y):
    """
    Fit a Gaussian Process Regression Model to a single gene

    Arguments
    =========
    x - 1 x D numpy array (e.g days = (0, 2, 4, ...))
    y - 1 x D numpy array gene expression data (e.g log 2 fold change)
    
    Returns
    =========
    GPModel - a Gaussian Process Regression model (GPy - package; SheffieldML)

    """
    # Define covariance function (mean function = 0) and create GPR Model
    # Prior to optimising hyperparameters set variance = set it to variance of data
    covFun = GPy.kern.RBF(input_dim=1, variance=np.var(y), lengthscale=config.LENGTHSCALE)
    GPModel = GPy.models.GPRegression(X=x, Y=y, kernel=covFun)
    GPModel.Gaussian_noise.variance = config.NOISE_VARIANCE
    # Bound parameters
    GPModel.Gaussian_noise.variance.constrain_bounded(config.LOWER_BOUND_NOISE, 
                                                      config.UPPER_BOUND_NOISE, warning=False)
    GPModel.rbf.lengthscale.constrain_bounded(config.LOWER_BOUND_LENGTHSCALE, 
                                              config.UPPER_BOUND_LENGTHSCALE, warning=False)
    # Optimise hyperparameters
    GPModel.optimize()
    
    return GPModel

#*************************************************************************************************#
def get_gene_params(GPModel):        
    """
    Get gene_model(x, y) hyperparameters = [rbfVar, rbfScale, noiseVar]

    Arguments
    =========
    GPModel - a gene-wise GP model returned by gene_model(x, y)
    
    Returns
    =========
    [rbfVar, rbfScale, noiseVar]

    """        
    
    return GPModel.rbf.variance, GPModel.rbf.lengthscale, GPModel.Gaussian_noise.variance  

#*************************************************************************************************#
def compute_gene_metrics(GPModel):
    """
    Compute four metrics for a single gene model:
    1) max absolute log fold change (maxLogFC)
    2) signal-to-noise-ratio (SNR) 
    3) score (as per Kalaitzis et al.(2011))
    4) net log fold change (netLogFC) i.e area under curve

    Arguments
    =========
    GPRModel - a gene-wise GP model returned by gene_model(x, y)
    
    Returns
    =========
    [maxLogFC, SNR, score, netLogFC]

    """
    # Max log fold change
    xTest = np.linspace(start=min(GPModel.X), stop=max(GPModel.X), num=100)[:, None]   
    yTest = GPModel.predict(xTest)[0]
    maxLogFC = np.max(np.abs(yTest))
    netLogFC = np.trapz(yTest, xTest, axis=0)       
    
    # SNR (signal variance/noise variance)
    SNR = float(GPModel.rbf.variance/GPModel.Gaussian_noise.variance)
    
    # Score as per Kalaitzis et al. (2011)
    logLikeNotQuiet = GPModel.log_likelihood()
    # Fit assuming gene is quiet    
    # Fix hyperparameters assuming data is just noise
    # lengthscale --> inf, rbf.variance=0, noise.variance=var(signal)  
    GPModel.rbf.lengthscale.fix(value=1000, warning=False)
    GPModel.rbf.variance.fix(value=0, warning=False)
    GPModel.Gaussian_noise.variance.fix(value=np.var(GPModel.Y), warning=False)
    logLikeQuiet = GPModel.log_likelihood()
    score = np.log(np.exp(logLikeNotQuiet)/np.exp(logLikeQuiet))
    
    return maxLogFC, SNR, score, netLogFC

#*************************************************************************************************#     
def write_gene_model_results(data, organ, strain, path, params, metrics, smoothExprs, bDEGUp, 
                             bDEGDown, xTestInterp, xTestDEG):  
    """
    Write all the results from fitting gene_model() and calling differential expression to disk
    This is a helper function to shorten code in the main script

    Arguments
    =========
    data - clean pandas data frame of the microarray data
    organ - Blood/Spleen 
    strain - AS/CB
    path - a structured dictionary of paths returned by config
    params - numpy array of [rbfVar, rbfScale, noiseVar]
    metrics - numpy array of [maxLogFC, SNR, score, netLogFC, rank]
    smoothExprs - numpy array of interpolated time-series
    bDEGUp/Down - numpy array of booleans DEG yes/no
    
    Returns
    =========
    None - data saved to disk

    """     
    # Write GP params to disk
    df = pd.DataFrame(data=params, columns=['rbfVar', 'rbfScale', 'noiseVar'])
    df['ProbeID'] = data['ProbeID']; df['Symbol'] = data['Symbol']
    df.to_csv(os.path.join(path['GPFit']['Params'], organ + strain + '.csv'), index=False)
    
    # Write GP time-series metrics to disk
    df = pd.DataFrame(data=metrics, columns=['maxLogFC', 'SNR', 'score', 'netLogFC', 'rank'])
    df['ProbeID'] = data['ProbeID']; df['Symbol'] = data['Symbol']
    df.to_csv(os.path.join(path['GPFit']['Metrics'], organ + strain + '.csv'), index=False)
    
    # Write smooth expression to disk
    df = pd.DataFrame(data=smoothExprs, columns=list(map(str, xTestInterp.flatten())))
    df['ProbeID'] = data['ProbeID']; df['Symbol'] = data['Symbol']
    df.to_csv(os.path.join(path['GPFit']['SmoothExprs'], organ + strain + '.csv'), index=False)
    
    # Write DEG results to disk
    colName = ['Day%d' % j for j in xTestDEG]  
    # Up regulated
    df = pd.DataFrame(data=bDEGUp, columns=colName)
    df['ProbeID'] = data['ProbeID']; df['Symbol'] = data['Symbol']
    df.to_csv(os.path.join(path['GPFit']['DEG'], 'b' + organ + strain + 'Up.csv'), index=False)        
    # Down regulated    
    df = pd.DataFrame(data=bDEGDown, columns=colName)
    df['ProbeID'] = data['ProbeID']; df['Symbol'] = data['Symbol']
    df.to_csv(os.path.join(path['GPFit']['DEG'], 'b' + organ + strain + 'Down.csv'), index=False)  
    
    # Get gene list per day and store those
    geneUp = []
    geneDown = []        
    for j in range(bDEGUp.shape[1]):
        geneUp.append(data.loc[bDEGUp[:, j], 'Symbol'])  
        geneDown.append(data.loc[bDEGDown[:, j], 'Symbol']) 
        
    # Write gene lists to csv
    io.write_list_to_csv(os.path.join(path['GPFit']['DEG'], organ + strain + "Up.csv"), colName, geneUp)
    io.write_list_to_csv(os.path.join(path['GPFit']['DEG'], organ + strain + "Down.csv"), colName, geneDown)
    
#*************************************************************************************************#  
def cluster_with_MOHGP(probesToCluster, organ, strain, K, alpha, path, seed=0):
    """
    Cluster time-series profiles with Mixtures of Hierarchical Gaussian Processes (MOHGP)

    Arguments
    =========
    probesToCluster - a set of unique probeIDs to cluster
    organ - Blood/Spleen
    strain - AS/CB
    K - init no. of clusters
    alpha - concentration parameter/strength parameter of the Dirichlet Process Prior
    path - a structured dictionary of paths returned by config
    seed - to reproduce results due to multiple local optima
    
    Returns
    =========
    None - a Mixture of Hierarchical Gaussian Process model is fitted and saved to disk

    """
    # To reproduce results
    np.random.seed(seed) 
    
    # Load gene expression data and assign training data
    # Remember column names are of the type Day.NReplicate e.g 3.4 = day 3, replicate 4
    data = pd.read_csv(os.path.join(path['Data'], organ + strain + ".csv"))
    yTrain = data[data['ProbeID'].isin(probesToCluster)] # subset data
    # NOTE: ordering of ProbeID and probesToCluster is different hence store probeID/geneSymbol
    geneID = yTrain[['ProbeID', 'Symbol']]
    yTrain = yTrain.drop(['ProbeID', 'Symbol'], axis=1)
    xTrain = np.floor(yTrain.columns.values.astype('float'))[:, None]
    yTrain = yTrain.as_matrix()
        
    # Mixture of Hierarchical Gaussian Process model
    # Define the covariance functions for the hierarchical GP structure
    # The model of any cluster of genes has a hierarchical structure, with the unknown 
    # cluster-specific mean drawn from a GP, and then each gene in that cluster being drawn from 
    # a GP with said unknown mean function.
    # Covariance function for the latent function that describes each cluster. 
    covFunClust = GPy.kern.RBF(input_dim=1, variance=np.var(yTrain.ravel()), lengthscale=config.LENGTHSCALE) 
    # Covariance function that describes how eacg time-course (gene) deviates from the cluster
    covFunGene = GPy.kern.RBF(input_dim=1, variance=np.var(yTrain.ravel())/10, lengthscale=config.LENGTHSCALE) + \
                 GPy.kern.White(1, variance=config.NOISE_VARIANCE)
    # Set-up the clustering problem NB: For large alpha P resembles Po (i.e the base distribution)
    fit = GPclust.MOHGP(X=xTrain, kernF=covFunClust, kernY=covFunGene, Y=yTrain, K=K, prior_Z='DP', alpha=alpha)   
    # Constrain lengthscales (to avoid very short lengthscales as per Topa et al. (2012) on arXiv)
    fit.rbf.lengthscale.constrain_bounded(config.LOWER_BOUND_LENGTHSCALE, 
                                          config.UPPER_BOUND_LENGTHSCALE , warning=False)
    fit.sum.rbf.lengthscale.constrain_bounded(config.LOWER_BOUND_LENGTHSCALE, 
                                              config.UPPER_BOUND_LENGTHSCALE , warning=False)    
    fit.hyperparam_opt_interval = 1000 # how often to optimize the hyperparameters
    
    # Optimise hyperparameters
    fit.optimize()
    fit.systematic_splits(verbose=False)
    
    # Name and reorder cluster columns i.e cluster 1 = the biggest one etc..   
    fit.name = organ + strain
    fit.reorder()
    
    # Write results to disk
    write_MOHGP_results(organ, strain, path, fit, geneID)
    
    
#*************************************************************************************************#     
def write_MOHGP_results(organ, strain, path, fit, geneID):
    """
    A helper function akin to 'write_gene_model_results' to save MOHGP results to disk

    Arguments
    =========
    organ - Blood/Spleen
    strain - AS/CB
    path - a structured dictionary of paths returned by config
    fit - a Mixture of Hierarchical Gaussian Process model
    geneID - a pandas data frame containing ordered probeID/geneSymbols
            NOTE: this is different from probesToCluster order
        
    Returns
    =========
    None - data saved to disk

    """
    # Extract the cluster assigned to each probe
    clustNum = np.argmax(fit.phi, axis=1) + 1 # cluster number
    clustName = [strain + '_' + organ[:2] + '_%02d' % i for i in clustNum] # cluster name
    geneID['Cluster'] = clustName # add to data frame
    
    # Extract the gene and probe list
    geneList = []; probeList = []; header = []
    for name in np.unique(clustName):
        bWant = geneID['Cluster'] == name
        geneList.append(list(geneID.loc[bWant, 'Symbol']))
        probeList.append(list(geneID.loc[bWant, 'ProbeID']))
        header.append(name)
   
    # Save to disk
    io.write_list_to_csv(os.path.join(path['Clust']['GeneList'], organ + strain + '.csv'), 
                         header, geneList)
    io.write_list_to_csv(os.path.join(path['Clust']['ProbeList'], organ + strain + '.csv'), 
                         header, probeList) # Probe list 
    
    # Save model and standard plot    
    io.save_pickle(os.path.join(path['Clust']['Model'], organ + strain + ".pickle"), fit)
    io.save_pdf(os.path.join(path['Clust']['Plot'], organ + strain + ".pdf"), standard_plot(fit))
    
    # Compute cluster predictions for xTest where xTest is taken from SmoothExprs
    data = pd.read_csv(os.path.join(path['GPFit']['SmoothExprs'], organ + strain + ".csv"))  
    xTest = data.drop(['ProbeID', 'Symbol'], axis=1).columns.values.astype('float64')[:, None]
    mu, var = fit.predict_components(xTest) # Compute posterior mean and posterior variance
    # Write to disk (mu row ordering is biggest to smallest cluster)
    df = pd.DataFrame(data=np.array(mu), columns=list(map(str, xTest.flatten())))
    df['Cluster'] = header # header = cluster name
    df.to_csv(os.path.join(path['Clust']['Centres'], organ + strain + '.csv'), index=False)
    clustCentre = df # for readability
    
    # Merge smooth expression data frame with gene ID
    smoothExprs = pd.merge(geneID, data, how='left', on=['ProbeID', 'Symbol'])
    
    # Produce alternate plot
    hFig = alternate_plot(smoothExprs, clustCentre, config.COL[organ])    
    io.save_pdf(os.path.join(path['Clust']['Plot'], organ + strain + '2.pdf'), hFig)
   
#*************************************************************************************************# 
def standard_plot(fit):
    """
    Use the GPy plotting functionality to produce the standard plot for time-series by cluster

    Arguments
    =========
    fit - a Mixture of Hierarchical Gaussian Process model
    
    Returns
    =========
    hFig - a handle to the figure

    """
    hFig = plt.figure(figsize=(11.69, 8.27))    
    fit.plot(on_subplots=True, colour=True, in_a_row=False, newfig=False, 
             min_in_cluster=0.99, joined=False, ylim=config.YLIM)
    
    return hFig

#*************************************************************************************************# 
def alternate_plot(smoothExprs, clustCentre, colour):
    """
    Manually replicating 'standard_plot()' with a few tweaks

    Arguments
    =========
    smoothExprs - a pandas data frame containing data to plot and cluster name
    clustCentre - a pandas data frame containing clust centres computed from MOHGP fit
    colour - colour of cluster centre
    
    Returns
    =========
    hFig - a handle to the figure

    """
    # Initialise plot values (how many subplots etc.)     
    Ntotal = len(np.unique(smoothExprs['Cluster']))
    Nx = np.floor(np.sqrt(Ntotal))
    Ny = int(np.ceil(Ntotal/Nx))
    Nx = int(Nx) 
    hFig = plt.figure(figsize=(11.69, 8.27)) # w,h tuple in inches
    
    # For every cluster plot smoothed fi
    for i, name in enumerate(np.unique(smoothExprs['Cluster'])):
        thisExprs = smoothExprs[smoothExprs['Cluster']==name].filter(regex='\.') # keep numeric col
        xTest = thisExprs.columns.values.astype('float')
        
        # Plot all probes of that label
        hAx = hFig.add_subplot(Nx, Ny, i+1)
        hAx.plot(xTest, thisExprs.as_matrix().T, linewidth=config.LWD['S'], color=config.COL['Shade'])
        
        # Plot mean of that cluster
        thisCentre = clustCentre[clustCentre['Cluster']==name].filter(regex='\.')
        hAx.plot(xTest, thisCentre.as_matrix().squeeze(), linewidth=config.LWD['L'], color=colour)
        
        # Set limits and title/ Limits should be adaptive but for now it's fine        
        hAx.axhline(y=0, color=config.COL['Zero'], linestyle='--', linewidth=config.LWD['M'])
        hAx.set_title("%s (%d probes)" % (name, thisExprs.shape[0]), fontsize=9)
        hAx.set_xlim(config.XLIM)
        hAx.set_ylim(config.YLIM)
        #hAx.set_ylabel("$\log_2$(FC) wrt naive", fontsize=8)
        #hAx.set_xlabel("Days", fontsize=8)
        hAx.grid(axis="y")
     
    return hFig  

#*************************************************************************************************#     