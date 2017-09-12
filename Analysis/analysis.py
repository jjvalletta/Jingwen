#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Top-level script to analyse results from clustering mice micro-array data
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
import matplotlib.pyplot as plt # version 2.0.2
from matplotlib_venn import venn2 # venn diagram
import pandas as pd # version 0.20.1
import time
from matplotlib.backends.backend_pdf import PdfPages # to store Pdf
# "Local" libraries
from config import setup_folders, ORGANS, STRAINS, COL, IMMUNE_PATH
import util
import input_output as io
import enrichr
import compare

#*************************************************************************************************#
# Setup/Create folders and thus all paths (stored as a structured dictionary)
#*************************************************************************************************#
path = setup_folders()

#*************************************************************************************************#
# PCA plot considering all measured probes
#*************************************************************************************************#
for organ in ORGANS:
    for strain in STRAINS:
        util.pca_plot(organ, strain, path, True)
    
#*************************************************************************************************#
# Venn diagram (Blood vs Spleen)
#*************************************************************************************************#
blood = util.top_ranked_genes('Blood', 'AS', path, IDType='Symbol')
spleen = util.top_ranked_genes('Spleen', 'AS', path, IDType='Symbol')
hFig = plt.figure()
hVen = venn2([spleen, blood], set_colors=(COL['Spleen'], COL['Blood']), 
             set_labels=('Spleen', 'Blood'), alpha=0.9)
plt.title('Genes')    
io.save_pdf(os.path.join(path['Misc'], "VennDiagAS.pdf"), hFig)

#*************************************************************************************************#
# Venn diagram (Blood/Spleen vs Immune Genes)
#*************************************************************************************************#
immune = pd.read_csv(IMMUNE_PATH)
immune = set(immune['GeneSymbol'].str.upper())

for organ in ORGANS:
    genes = util.top_ranked_genes(organ, 'AS', path, IDType='Symbol')
    #genes = pd.read_csv(os.path.join(path['Data'], organ + 'AS.csv'))
    #genes = set(genes['Symbol'])
    genes = set(map(str.upper, genes)) # capitalize to compare like-with-like
    hFig = plt.figure()
    venn2([genes, immune], set_labels=('All ' + organ, 'Immune'), alpha=0.9)
    io.save_pdf(os.path.join(path['Misc'], 'VennDiag' + organ + 'vsImmuneAS.pdf'), hFig)

#*************************************************************************************************#
# Pathway analysis all clusters
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
# Compare time-profile of every cluster
#*************************************************************************************************#
strainA = 'AS'; strainB = 'AS'
for organA, organB in zip(ORGANS, ORGANS[::-1]):
    data = pd.read_csv(os.path.join(path['Clust']['ProbeList'], organA + strain + '.csv'))
    filePath = os.path.join(path['Comparison'], 
                            'Compare' + organA + strainA + 'vs' + organB + strainB + '.pdf')
    with PdfPages(filePath) as pdf:
        for cluster in data.keys():
            df = pd.DataFrame({'ProbeID': data[cluster].dropna()})
            hFig = compare.genes([organA, organB], [strainA, strainB], path, df, desc=cluster)
            pdf.savefig(hFig) # save figure
            plt.close(hFig) 

#*************************************************************************************************#
#
#*************************************************************************************************#