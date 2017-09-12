#*************************************************************************************************#
#*************************************************************************************************#
# Title:    Enrichr gene set enrichment analysis for mice microarray data
# Author:   John Joseph Valletta
# Date:     01/09/2017
# Comment:  Big change: port all code to Python3, GPy-1.7.7, paramz-0.7.4, GPclust-0.1.0 
# Data:     Microarray Illumina Mouse WG6 v2 (45,281 probe sets representing 30,854 genes)
# From:     Jingwen Lin and Jean Langhorne, Crick Institute 
#*************************************************************************************************#
#*************************************************************************************************#

#*************************************************************************************************#
# Preamble
#*************************************************************************************************#
# Libraries
import requests # URL get requests version 2.14.2
import json # deserialize instance containing JSON --> Python object version 2.0.9 
import os
import mygene # http://mygene.info/ and https://pypi.python.org/pypi/mygene and http://mygene-py.readthedocs.org/en/latest/
# "Local" libraries
import input_output as io
from config import ENRICHR, SPECIES

#*************************************************************************************************# 
# Functions
#*************************************************************************************************# 
def enrichment(geneList, description='test', basePath=os.getcwd(), library=ENRICHR['Libraries']):
    """
    Gene set enrichment using Enrichr (http://amp.pharm.mssm.edu/Enrichr/) and export results
    
    Note: **This only supports human and mouse** 

    Arguments
    =========
    geneList - list of symbols (e.g geneList = ['PHF14', 'RBM3', 'MSL1', 'PHF21A']; case insensitive)
    description - description of gene set, which is also name of created folder
    basePath - path were folder is created and results exported to
    library - a list of libraries to perform enrichment on, which is also the name of exported file 
            (http://amp.pharm.mssm.edu/Enrichr/#stats)
    
    Returns
    =========
    None - A folder is created with results in it e.g basePath/description/KEGG_2016.txt

    """    
    # Upload gene list to enrichr and get a jobID
    geneString = '\n'.join(geneList) # convert list to string as needed by the url request
    payload = {'list': (None, geneString), 'description': (None, description)} # query parameters
    response = requests.post(ENRICHR['UploadURL'], files=payload)
    if not response.ok:
      raise Exception('Error analyzing gene list: %s' % description)
    jobID = json.loads(response.text) # enrichr jobID
    
    # Create folder where to put results
    path = io.create_folder(basePath, description)  
    
    # Compute enrichment and export results for every library
    userListID = str(jobID['userListId'])
    for lib in library:
        # Retrieve results
        url = ENRICHR['ExportURL'] + ENRICHR['Query'] % (userListID , lib, lib)
        response = requests.get(url, stream=True)
        # Write to text file        
        fullPath = os.path.join(path, lib + '.txt')         
        with open(fullPath, 'wb') as fID:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    fID.write(chunk)

#*************************************************************************************************#
def gene_description(geneList, description='test', basePath=os.getcwd()):
    """
    Uses mygene to get generic information about individual genes. Strictly speaking it's a utility function
    that should be in a separate module, but for now, it makes sense to put it here, even though it has nothing
    to do with enrichr per se. 
    
    Arguments
    =========
    geneList - list of gene symbols (e.g geneList = ['PHF14', 'RBM3', 'MSL1', 'PHF21A'])
    description - description of gene set, which is also name of folder
    basePath - path were folder is created and results exported to
    
    Returns
    =========
    None - A text file is created with results in it e.g basePath/AS_Sp_01.txt

    """    
    # Run query 
    mg = mygene.MyGeneInfo() # instantiate MyGeneInfo class
    out =  mg.querymany(qterms=geneList, scopes="symbol", fields=['name', 'generif'], species=SPECIES)   
    
    # Open file and start writing gene description in it
    with open(name=os.path.join(basePath, description + '.txt'), mode='w') as fID:
        # Loop through all genes
        for iGene in range(len(out)): 
            # Header with gene name            
            fID.write("\n#====================================================================================#\n")
            fID.write("# %s\n" % out[iGene]['query'])
            fID.write("#====================================================================================#\n")
            if 'name' in out[iGene]:
                fID.write("Name: %s\n\n" % out[iGene]['name'])
            else:
                fID.write("Name: Not found\n\n")
            # Generic information
            if 'generif' in out[iGene]:
                fID.write("Generif:\n\n")
                for iTerm in range(len(out[iGene]['generif'])):
                    fID.write("Pubmed: %i\n" % out[iGene]['generif'][iTerm]['pubmed'])
                    fID.write("%s\n\n" % out[iGene]['generif'][iTerm]['text'])
            else:
                fID.write("Generif: Not found\n\n")
                    
#*************************************************************************************************#                