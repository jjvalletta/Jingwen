#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Config file containing any constants required for the analysis of mice data
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
import os
from palettable.tableau import Tableau_20 # version 3.1.0
from palettable.colorbrewer.qualitative import Set1_7 # PCA colours
# "Local" libraries
import input_output as io

# Constants
# Paths
DATA_PATH = os.path.join('..', 'Data', 'Log2FC') # location of clean data Blood/Spleen AS/CB data
BASE_PATH = os.path.join('..', 'Results', 'Analysis030917') # results base path
IMMUNE_PATH = os.path.join('..', 'Data', 'ImmuneGenes', 'immune_genes.csv') # list of immune genes

# All experimental conditions
ORGANS = ['Blood', 'Spleen'] # organs
STRAINS = ['AS', 'CB'] # Plasmodium chabaudi strain

# To reproduce results
SEED = 120
TOP_RANKED = 4528 # top 10% of probes

# Gaussian process initial conditions and bounds 
LENGTHSCALE = 4 # 1/3 of the length of time series i.e 12 days
NOISE_VARIANCE = 0.01 # arbitrary measurement noise
# Bound parameter values as per Topa et al. (2012) on arXiv
# This is chosen to focus on frquencies below the Nyquist frequency
LOWER_BOUND_LENGTHSCALE = 1.6 # lengthscale bounded to 0.8199 x deltaT = 0.82*2days = 1.6 days
UPPER_BOUND_LENGTHSCALE = 24 # a lengthscale of 12 is pretty much flat already
LOWER_BOUND_NOISE = 0.0001
UPPER_BOUND_NOISE = 2

# Plotting control (I want all plots to look the same)
XLIM = (0, 12)  
YLIM = (-7, 7)
COL = {'Blood': '#e41a1c', 'Spleen': '#377eb8', 
       'Shade': '#525252', 'Zero': 'black',
       'Cycle': Tableau_20.hex_colors,
       'PCA': Set1_7.hex_colors}
COL['PCA'][5] = COL['PCA'][0] # don't like yellow use pink instead
COL['PCA'][0] = '#000000' # Naive/Control == Black
LWD = {'L': 4, 'M': 2, 'S': 0.5} # linewidth; large, medium, small
ALPHA = 0.5 # set transparency

# Pathway analysis
SPECIES = 'mouse' 
ENRICHR = {'Libraries': ['GO_Biological_Process_2015', 'GO_Cellular_Component_2015', 
                         'GO_Molecular_Function_2015', 'WikiPathways_2016', 'KEGG_2016', 
                         'Panther_2016', 'Reactome_2016', 'PPI_Hub_Proteins', 'BioCarta_2016', 
                         'OMIM_Disease', 'Mouse_Gene_Atlas', 'Transcription_Factor_PPIs', 
                         'GeneSigDB', 'MSigDB_Computational'], 
            'UploadURL': 'http://amp.pharm.mssm.edu/Enrichr/addList', 
            'ExportURL': 'http://amp.pharm.mssm.edu/Enrichr/export', 
            'Query': '?userListId=%s&filename=%s&backgroundType=%s'}

#*************************************************************************************************#
# Functions
#*************************************************************************************************#
def setup_folders(basePath=BASE_PATH, dataPath=DATA_PATH):
    """
    Create all folders where analysis results will be stored
    
    Arguments
    =========
    basePath - typically of the form "/AnalysisDDMMYY"
    dataPath - location of clean data Blood/Spleen AS/CB .csv
    
    Returns
    =========
    path - a structured dictionary of paths

    """
    # Create empty dictionary
    path = {}
    # RawData    
    path['Data'] = dataPath
    # GPFit
    path['GPFit'] = {}
    path['GPFit']['SmoothExprs'] = io.create_folder(os.path.join(basePath, 'GPFit'), 'SmoothExprs')
    path['GPFit']['Metrics'] = io.create_folder(os.path.join(basePath, 'GPFit'), 'Metrics')
    path['GPFit']['Params'] = io.create_folder(os.path.join(basePath, 'GPFit'), 'Params')
    path['GPFit']['DEG'] = io.create_folder(os.path.join(basePath, 'GPFit'), 'DEG')
    # Clust
    path['Clust'] = {}
    path['Clust']['Model'] = io.create_folder(os.path.join(basePath, 'Clust'), 'Model')
    path['Clust']['GeneList'] = io.create_folder(os.path.join(basePath, 'Clust'), 'GeneList')
    path['Clust']['ProbeList'] = io.create_folder(os.path.join(basePath, 'Clust'), 'ProbeList')
    path['Clust']['Plot'] = io.create_folder(os.path.join(basePath, 'Clust'), 'Plot')
    path['Clust']['Centres'] = io.create_folder(os.path.join(basePath, 'Clust'), 'Centres')
    # Comparison
    path['Comparison'] = io.create_folder(basePath, 'Comparison')
    # PathwayAnalysis
    path['PathwayAnalysis'] = {}
    path['PathwayAnalysis']['Enrichr'] = io.create_folder(os.path.join(basePath, 'PathwayAnalysis'), 'Enrichr')
    path['PathwayAnalysis']['Reactome'] = io.create_folder(os.path.join(basePath, 'PathwayAnalysis'), 'Reactome')
    # Misc
    path['Misc'] = io.create_folder(basePath, 'Misc')
    
    return path