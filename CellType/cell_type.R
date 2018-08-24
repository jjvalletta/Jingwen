#***************************************************************************************#
# Title:    Estimate cell type proportions for each sample
# Author:   John Joseph Valletta
# Date:     24/08/2018
# Ref: https://www.nature.com/articles/srep40508
# Ref: https://github.com/chenziyi/ImmuCC/blob/master/Microarray_Deconvolution.R
# Output: Blood/Spleen AS/CB .txt
#***************************************************************************************#

#***************************************************************************************#
# Clear workspace
#***************************************************************************************#
rm(list = setdiff(ls(), lsf.str())) 

#***************************************************************************************#
# Libraries
#***************************************************************************************#
source('cibersort.R')

#***************************************************************************************#
# Constants
#***************************************************************************************#
ORGANS <- c('Blood', 'Spleen') # tissue samples
STRAINS <- c('AS', 'CB') # Plasmodium Chabaudi strain
DATA_PATH <- file.path("..", "Data", "DeconReady")
SIG_MATRIX <- file.path(DATA_PATH, "SigMatrix.txt")

# No. of permutations to compute  an empirical P value for the deconvolution using 
# Monte Carlo sampling. To test the null hypothesis that no cell types in the signature matrix 
# are present in a given gene expression profile mixture 
NPerm <- 201 

#***************************************************************************************#
# Main function
#***************************************************************************************#
for (organ in ORGANS)
{
    for (strain in STRAINS)
    {
        mixFile <- file.path(DATA_PATH, paste0(organ, strain, '.txt'))
        dummy <- CIBERSORT(sig_matrix=SIG_MATRIX, mixture_file=mixFile, 
                           perm=NPerm, desc=paste0(organ, strain, '.txt'))
    }
}

#***************************************************************************************#
# Save sessionInfo() to disk
#***************************************************************************************#
writeLines(capture.output(sessionInfo()), "session_info.txt")