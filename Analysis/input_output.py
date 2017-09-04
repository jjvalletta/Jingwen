#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Input/Output utility functions for mice microarray data
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
import warnings
import csv # version 1.0
from itertools import zip_longest
from matplotlib.backends.backend_pdf import PdfPages # to store Pdf
import matplotlib.pyplot as plt # version 2.0.2
import pickle # to save any object

#*************************************************************************************************#
# Functions
#*************************************************************************************************# 
def create_folder(basePath, folderName):
    """
    Creates a folder (folderName) in basePath

    Arguments
    =========
    basePath   - path where the folder will be created
    folderName - name of folder

    Returns
    =========
    None (folder is created)

    """    
    # Create folder if it exists throw a warning
    thisDir = os.path.join(basePath, folderName)
    if not os.path.exists(thisDir):
        os.makedirs(thisDir)
    else:
        warnings.warn("Folder has already been created!", UserWarning)
    # Note:    
    # makedirs() creates all the intermediate directories if they don't exist (just like mkdir -p in bash).
    # mkdir() only creates a single sub-directory; throws an exception if intermediate directories aren't specified.
    
    return thisDir

#*************************************************************************************************#     
def write_to_csv(filePath, header, dataMatrix):
    """
    Write a .csv file in the following format
    
    ---------------
    | header      |
    ---------------
    | dataMatrix  |  
    ---------------

    Arguments
    =========
    filePath    - complete file path where .csv is to be saved
    header      - column names corresponding to data e.g day=0, day=0.25, etc..
    dataMatrix  - needs to match size of header

    Returns
    =========
    None (.csv is created)

    """
    with open(name=filePath, mode='w') as fID:
        # Write header        
        writer = csv.writer(fID, delimiter=',')             
        writer.writerow(list(header))
        # Write data line by line
        for iLine in range(dataMatrix.shape[0]): 
            writer.writerow(list(dataMatrix[iLine, :]))
            
#*************************************************************************************************#     
def write_list_to_csv(filePath, header, data):
    """
    Write a list to csv. Specifically the following:
        
    If data = [[2,3,4,8],[5,6]] and header = ["x1", "x2"]
    Then written csv looks like this:
    x1, x2    
    2, 5
    3, 6
    4,
    8,

    Arguments
    =========
    filePath    - complete file path where .csv is to be saved
    header      - name of columns (1 x D)
    data        - list of data

    Returns
    =========
    None (.csv is created)

    """
    with open(filePath, "w") as fID:
        writer = csv.writer(fID)
        writer.writerow(header)    
        for iLine in zip_longest(*data):
            writer.writerow(iLine)

#*************************************************************************************************# 
def save_pdf(filePath, hFig):
    """
    Saves a figure as a .PDF file  

    Arguments
    =========
    filePath - path where pdf is to be saved
    hFig - figure handle

    Returns
    =========
    None (pickle file is created)

    """     
    with PdfPages(filePath) as pdf:     
        pdf.savefig(hFig) # save figure
        plt.close(hFig) # close figure
        
#*************************************************************************************************#  
def save_pickle(filePath, anyObject):
    """
    Saves anyObject to filePath as a pickle  

    Arguments
    =========
    filePath   - path where to save the pickled object
    anyObject - the object we want to pickle and save

    Returns
    =========
    None (pickle file is created)

    """    
    with open(filePath, 'wb') as fID:
        pickle.dump(anyObject, fID)

#*************************************************************************************************#  
def load_pickle(filePath):
    """
    Loads a pickled object found in filePath  

    Arguments
    =========
    filePath   - path where pickled object is

    Returns
    =========
    anyObject - the loaded pickled object

    """    
    with open(filePath, 'rb') as fID:
        anyObject = pickle.load(fID)
    
    return anyObject

#*************************************************************************************************# 