#***************************************************************************************#
# Title:    Plot cell type results
# Author:   John Joseph Valletta
# Date:     23/08/2018
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
library(tidyverse)
library(stringr) # str_replace_all()

#***************************************************************************************#
# Constants
#***************************************************************************************#
ORGANS <- c('Blood', 'Spleen') # tissue samples
STRAINS <- c('AS', 'CB') # Plasmodium Chabaudi strain

#***************************************************************************************#
# Main function
#***************************************************************************************#
for (organ in ORGANS)
{
    for (strain in STRAINS)
    {
        # Read CIBERSORT results
        res <- read.table(paste0(organ, strain, '.txt'), header=T, sep="\t")
        res <- gather(res, "CellType", "Proportion", 2:26)
        
        # Get rid of the X's so that I can numerically sort the levels
        res$Mixture <- str_replace_all(res$Mixture, 'X', '')
        newLevels <- sort(as.numeric(unique(res$Mixture)))
        res$Mixture <- factor(res$Mixture, levels=as.character(newLevels))
        
        g <- ggplot(res, aes(x=Mixture, y=Proportion, fill=CellType, alpha=1-P.value)) + 
            geom_bar(position='fill', stat='identity')
        
        ggsave(filename=paste0(organ, strain, '.pdf'), plot=g, width=297, height=210, units='mm')
    }
}