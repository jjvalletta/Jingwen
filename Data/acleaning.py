#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Cleaning mice micro-array data, specifically:
#           a) Remove meta data columns apart from ProbeID and Symbol
#           b) Fix issue with gene names auto renamed by Excel (i.e Mar-01 --> MARCH1)
#           c) Split data by P. chabaudi strain (AS, CB)
# Author:   John Joseph Valletta
# Date:     28/08/2017
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
import re # regular expression: version 2.2.1
import os

# Constants
DATA_PATH = {'Blood': os.path.join('Original', 'BL-all.entities.txt'), 
             'Spleen': os.path.join('Original', 'SP-all.entities.txt')}
ANNOTATION_PATH = os.path.join('Original', "MouseWG-6_V2_0_R3_11278593_A.txt")
ORGANS = ['Blood', 'Spleen'] # tissue samples
STRAINS = ['AS', 'CB'] # Plasmodium Chabaudi strain
# Reg expr for sample identifier (e.g Bl_AS_N0.3) --> Organ=Blood/Strain=AS/N/Day=0/Replicate=3
ID_REGEXP = "(?P<organ>\w{2,2})_(?P<strain>\w{2,2})_(?P<naive>\w{1,1})(?P<day>\d{1,2}).(?P<NReplicate>\d{1,1})" 

#*************************************************************************************************#
# Functions
#*************************************************************************************************#
def cleanup(organ):
    """
    Cleanup Blood/Spleen data files i.e remove columns apart from ProbeID, Symbol and data
         
    Arguments
    =========
    organ - Blood or Spleen
  
    Returns
    =========
    data - clean pandas data frame; remove meta data columns apart from ProbeID and Symbol

    """
    # Read raw data file
    data = pd.read_table(DATA_PATH[organ], sep="\t", low_memory=False)
    
    # Initiate dictionary of columns to be renamed
    if organ == 'Blood':
        renameDict = {'Bl_CB__D8.1': 'Bl_CB_D8.1'} # for some reason it had two underscores
    elif organ == 'Spleen':
        renameDict = {}

    # Find GeneSymbol and ProbeID columns and add to dictionary (weird structure of .txt file)
    iSymbol = np.where(data.iloc[0]=='Symbol')[0][0]  
    iProbe = np.where(data.iloc[0]=='ProbeID')[0][0]       
    renameDict[data.keys()[iSymbol]] = 'Symbol'
    renameDict[data.keys()[iProbe]] = 'ProbeID' 

    # Rename columns
    data.rename(columns=renameDict, inplace=True) 

    # Delete first row containing sample names (8762536052 B(normalized) etc.)
    data.drop(data.index[0], 0, inplace=True) # axis = 0 rows and 1 cols

    # Delete the unnamed columns
    delColumns = list(data.filter(regex='Unnamed*')) # list of columns to delete
    data.drop(delColumns, axis=1, inplace=True)
    
    # Change dtype to float/int/string
    typeDict = {'ProbeID': 'int', 'Symbol': 'str'}
    for col in data.keys():
        if col not in ['ProbeID', 'Symbol']:
            typeDict[col] = 'float'
    data = data.astype(typeDict)

    return data

#*************************************************************************************************#
def fix_excel_names(data, annoPath=ANNOTATION_PATH):
    """
    Fix issue with gene names auto renamed by Excel (i.e Mar-01 --> MARCH1)
         
    Arguments
    =========
    data - clean pandas data frame returned by cleanup()
  
    Returns
    =========
    data - same as data but with gene names fixed

    """
    # Load annotation file:
    # Interesting columns are:
    # Entrez_Gene_ID --> EntrezID (e.g 212772)
    # Symbol --> GeneSymbol (e.g Thrsp)
    # Probe_Id --> IlluminaID (e.g ILMN_204164)
    # Array_Address_Id --> ProbeID (e.g 2600193)
    anno = pd.read_table(annoPath, header=8, sep="\t") # read gene annotation file
    anno = anno[['Array_Address_Id', 'Symbol']] # keep required cols
    anno.rename(columns={'Array_Address_Id': 'ProbeID'}, inplace=True) # rename col
    anno.dropna(inplace=True) # drop NAs
    anno = anno.astype({'ProbeID': 'int'}) # change to int

    # Rename Symbol to ExcelSymbol i.e possibly wrong/old
    data = data.rename(columns={'Symbol': 'ExcelSymbol'})

    # Merge data sets
    data = pd.merge(data, anno, how='left', on='ProbeID')

    return data

#*************************************************************************************************#
def compute_log2FC(data, strain):
    """
    Split the data by strain (AS/CB) and compute log 2 fold change (logFC)
    
    Arguments
    =========
    data - pandas dataframe returned by cleanup(); fix_excel_names();
    strain - "AS" or "CB"
    
    Returns
    =========
    dfRaw - pandas dataframe raw data with columns renamed to days
    dfLog2FC - pandas dataframe log 2 fold change with columns renamed to days
    
    """
    # Initialise variables
    renameDict = {} # renaming column names e.g Sp_CB_D8.2 --> 8.2
    dropCols = ['ExcelSymbol'] # list of columns to drop
     
    # Loop through all columns and decide whether to drop or rename
    cols = set(data.keys()) - set(['ExcelSymbol', 'Symbol', 'ProbeID']) 
    for col in cols:
        token = re.match(ID_REGEXP, col, flags=0) # reg expression : )
        if token.group('strain') in strain:
            # Check if it's a naive sample
            # Differentiate between naive day0 and day12
            if token.group('naive') in 'N' and token.group('day') not in '0' :
                day = 0
                NReplicate = int(token.group('NReplicate')) + 10
            else:
                day = token.group('day')
                NReplicate = token.group('NReplicate')
            # Rename columns
            renameDict[col] = str(day) + '.' + str(NReplicate)
        else:
            dropCols.append(col) # drop columns
             
    # Create dfRaw by dropping/renaming columns
    dfRaw = data.drop(dropCols, axis=1, inplace=False)
    dfRaw.rename(columns=renameDict, inplace=True)
    
    # Create log 2 fold change data frame
    meanNaive = dfRaw.filter(regex='^0\\.').mean(axis=1)
    dfLog2FC = dfRaw.filter(regex='\\.').subtract(meanNaive, axis=0)
    dfLog2FC['ProbeID'] = dfRaw['ProbeID']
    dfLog2FC['Symbol'] = dfRaw['Symbol']

    return dfRaw, dfLog2FC

#*************************************************************************************************#
# Main()
#*************************************************************************************************# 
# Loop through raw blood and spleen data file and clean up
for organ in ORGANS:
    print('a) Remove meta data columns apart from ProbeID and Symbol...\n')
    data = cleanup(organ)

    print('b) Fix issue with gene names auto renamed by Excel (i.e Mar-01 --> MARCH1)...\n')
    data = fix_excel_names(data)
    data.to_csv(os.path.join('Processed', organ + '.csv'), index=False) # write to file

    print('c) Split data by P. chabaudi strain (AS, CB)...\n')
    for strain in STRAINS:
        dfRaw, dfLog2FC = compute_log2FC(data, strain) # Note: day 0 + day 12 naive --> day 0
        dfRaw.to_csv(os.path.join('Processed', organ + strain + '.csv'), index=False)
        dfLog2FC.to_csv(os.path.join('Log2FC', organ + strain + '.csv'), index=False)