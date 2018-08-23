#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Top-level script to analyse results from clustering mice micro-array data
# Author:   John Joseph Valletta
# Date:     13/04/2018
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
import matplotlib.pyplot as plt # version 2.0.2
#from matplotlib_venn import venn2 # venn diagram
import pandas as pd # version 0.20.1
import time
from matplotlib.backends.backend_pdf import PdfPages # to store Pdf
import seaborn as sns
import numpy as np
from functools import reduce
# "Local" libraries
from config import setup_folders, ORGANS, STRAINS, COL, TOP_RANKED
import util
import input_output as io
import enrichr
import compare

#*************************************************************************************************#
# Setup/Create folders and thus all paths (stored as a structured dictionary)
#*************************************************************************************************#
path = setup_folders()

#*************************************************************************************************#
# PCA plot considering all measured probes and blood and spleen together
#*************************************************************************************************#
# Consider all tx
strain = 'AS'
for organ in ORGANS:
    hWt, hPCA = util.pca_plot([organ], strain, path, bLegend=True)
    io.save_pdf(os.path.join(path['Misc'], 'PCA' + organ + strain + 'All.pdf'), hPCA)
    io.save_pdf(os.path.join(path['Misc'], 'PCAWts' + organ + strain + 'All.pdf'), hWt)
# Spleen + Blood together
hWt, hPCA = util.pca_plot(ORGANS, strain, path, bLegend=True)
io.save_pdf(os.path.join(path['Misc'], 'PCA' + "".join(ORGANS) + strain + 'All.pdf'), hPCA)
io.save_pdf(os.path.join(path['Misc'], 'PCAWts' + "".join(ORGANS) + strain + 'All.pdf'), hWt)
    
# Consider ONLY top ranked
for organ in ORGANS:
    hWt, hPCA = util.pca_plot([organ], strain, path, topRanked=TOP_RANKED, bLegend=True)
    io.save_pdf(os.path.join(path['Misc'], 'PCA' + organ + strain + 'Top.pdf'), hPCA)
    io.save_pdf(os.path.join(path['Misc'], 'PCAWts' + organ + strain + 'Top.pdf'), hWt)
# Spleen + Blood together
hWt, hPCA = util.pca_plot(ORGANS, strain, path, topRanked=TOP_RANKED, bLegend=True)
io.save_pdf(os.path.join(path['Misc'], 'PCA' + "".join(ORGANS) + strain + 'Top.pdf'), hPCA)
io.save_pdf(os.path.join(path['Misc'], 'PCAWts' + "".join(ORGANS) + strain + 'Top.pdf'), hWt)

#*************************************************************************************************#
# Heatmap of top 200 ranked genes (pick shared blood/spleen)
#*************************************************************************************************#
# Find top 200 ranked genes
blood = util.top_ranked_genes('Blood', 'AS', path, 1650, 'ProbeID')
spleen = util.top_ranked_genes('Spleen', 'AS', path, 1650, 'ProbeID')
top = blood & spleen # this gives me 200 genes

# Get data
strain = 'AS'
df = [] # list of data frames
for organ in ORGANS:
    temp = pd.read_csv(os.path.join(path['Data'], organ + strain + '.csv'))
    temp = temp.set_index('ProbeID')
    temp = temp.loc[top]
    temp = temp.set_index('Symbol')
    temp = temp.add_prefix(organ[:2] + '_')
    df.append(temp)
data = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), df)
del(df)

obj = sns.clustermap(data, method='ward', metric='euclidean', cmap="RdBu_r", 
                     col_cluster=False, figsize=(12, 10))
plt.setp(obj.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.setp(obj.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=4)
obj.savefig(os.path.join(path['Misc'], 'HetmapTop200.pdf'))
plt.close(obj.fig) # close figure
        
#*************************************************************************************************#
# Ranking of DE genes [unique genes rather than probes]
#*************************************************************************************************#
hFig, hAxs= plt.subplots(2, 2, sharex=True, sharey=False, figsize=(11.69, 8.27))
metrics = ['score', 'SNR', 'maxLogFC', 'rank']
for i, hAx in enumerate(hFig.axes):
    for organ in ORGANS:
        for strain in ["AS"]:
            data = pd.read_csv(os.path.join(path['GPFit']['Metrics'], organ + strain + ".csv"), sep=",")  
            data = data.set_index('ProbeID')
            uniqueProbes = util.top_ranked_genes(organ, strain, path, len(np.unique(data['Symbol'])), 'ProbeID')
            data = data.loc[uniqueProbes]
            metric = data.sort_values(by=metrics[i], ascending=True)[metrics[i]].as_matrix()
            hAx.plot(metric, label=organ, linewidth=2, color=COL[organ])            
            if metrics[i] == 'rank':            
                hAx.axvline(x=len(metric) - TOP_RANKED, c=COL['Zero'], linewidth=2, linestyle='--')
                hAx.legend(loc=2)
    hAx.set_xlabel("No. of Genes")
    hAx.set_ylabel(metrics[i])
    hAx.set_xlim((0, len(metric)))
    plt.tight_layout()
io.save_pdf(os.path.join(path['Misc'], 'DEG.pdf'), hFig)

#*************************************************************************************************#
# Boxplot of posterior probablity that gene is in cluster X
#*************************************************************************************************#
for organ in ORGANS:
    for strain in STRAINS:
        util.cluster_membership_boxplot(organ, strain, path)
        
#*************************************************************************************************#
# EnrichR gene set enrichment analysis all clusters
#*************************************************************************************************#
for organ in ORGANS:
    for strain in STRAINS:
        genes = pd.read_csv(os.path.join(path['Clust']['GeneList'], organ + strain + '.csv'))      
        for cluster in genes.keys():
            time.sleep(10) # to avoid Connection reset by peer
            print(cluster)
            enrichr.enrichment(genes[cluster].dropna().tolist(), cluster, 
                               path['PathwayAnalysis']['Enrichr'])

#*************************************************************************************************#
# Compare cluster mean by Pearson correlation coefficient
#*************************************************************************************************#
compare.cluster_mean_pearson('Blood', 'AS', 'Spleen', 'AS', path)
            
#*************************************************************************************************#
# Compare time-profiles of each gene in every cluster across each tissue
# i.e Every gene in Bl_AS_01 how do they look in the Spleen and vice-versa
#*************************************************************************************************#
strain = 'AS'
for organA, organB in zip(ORGANS, ORGANS[::-1]):
    data = pd.read_csv(os.path.join(path['Clust']['ProbeList'], organA + strain + '.csv'))
    filePath = os.path.join(path['Comparison'], 
                            'Compare' + organA + strain + 'vs' + organB + strain + '.pdf')
    with PdfPages(filePath) as pdf:
        for cluster in data.keys():
            df = pd.DataFrame({'ProbeID': data[cluster].dropna()})
            hFig = compare.genes([organA, organB], [strain, strain], path, df, desc=cluster)
            pdf.savefig(hFig) # save figure
            plt.close(hFig) 

#*************************************************************************************************#
#
#*************************************************************************************************#

#
##*************************************************************************************************#
## Venn diagram (Blood vs Spleen)
##*************************************************************************************************#
#for IDType in ['ProbeID', 'Symbol']:    
#    blood = util.top_ranked_genes('Blood', 'AS', path, IDType=IDType)
#    spleen = util.top_ranked_genes('Spleen', 'AS', path, IDType=IDType)
#    hFig = plt.figure()
#    hVen = venn2([spleen, blood], set_colors=(COL['Spleen'], COL['Blood']), 
#                 set_labels=('Spleen', 'Blood'), alpha=0.9)
#    plt.title(IDType)    
#    io.save_pdf(os.path.join(path['Misc'], "VennDiagAS" + IDType + ".pdf"), hFig)
#
##*************************************************************************************************#
## Venn diagram (Blood/Spleen vs Immune Genes)
##*************************************************************************************************#
#immune = pd.read_csv(IMMUNE_PATH)
#immune = set(immune['GeneSymbol'].str.upper())
#
#for organ in ORGANS:
#    genes = util.top_ranked_genes(organ, 'AS', path, IDType='Symbol')
#    #genes = pd.read_csv(os.path.join(path['Data'], organ + 'AS.csv'))
#    #genes = set(genes['Symbol'])
#    genes = set(map(str.upper, genes)) # capitalize to compare like-with-like
#    hFig = plt.figure()
#    venn2([genes, immune], set_labels=('All ' + organ, 'Immune'), alpha=0.9)
#    io.save_pdf(os.path.join(path['Misc'], 'VennDiag' + organ + 'vsImmuneAS.pdf'), hFig)