#***************************************************************************************#
# Title:    Convert MSigDB's .gmt files to list of genes
# Author:   John Joseph Valletta
# Date:     29/06/2017
# Source: http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C7
#***************************************************************************************#

#***************************************************************************************#
# Clear workspace
#***************************************************************************************#
rm(list = setdiff(ls(), lsf.str())) 

#***************************************************************************************#
# Libraries
#***************************************************************************************#
library(qusage) # read.gmt

#***************************************************************************************#
# Main function
#***************************************************************************************#
df <- read.gmt(file.path("ImmuneGenes", "c7.all.v6.0.symbols.gmt"))
geneSymbol <- unlist(df)
geneSymbol <- unique(geneSymbol)
write.csv(data.frame(GeneSymbol=geneSymbol), 
          file.path("ImmuneGenes", "immune_genes.csv"), row.names=FALSE)
