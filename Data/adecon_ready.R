#***************************************************************************************#
# Title:    Make microarray data amenable for deconvolution i.e set it up in the right format
# Author:   John Joseph Valletta
# Date:     24/08/2018
# Ref: https://www.nature.com/articles/srep40508
# Ref: https://github.com/chenziyi/ImmuCC/blob/master/Microarray_Deconvolution.R
# Output: Rows are samples/Columns are gene symbols ONLY
#***************************************************************************************#

#***************************************************************************************#
# Clear workspace
#***************************************************************************************#
rm(list = setdiff(ls(), lsf.str())) 

#***************************************************************************************#
# Constants
#***************************************************************************************#
DATA_PATH <- "Processed" # Original data cleaned up by acleaning.py - no changes to data just removing unwanted cols etc.  
RANK_PATH <- file.path("..", "Results", "Analysis110418", "GPFit", "Metrics")
OUT_PATH <- "DeconReady"
SIG_PATH <- "DeconReady/srep40508-s1.csv" # signature matrix path downloaded from 
ORGANS <- c('Blood', 'Spleen') # tissue samples
STRAINS <- c('AS', 'CB') # Plasmodium Chabaudi strain

#***************************************************************************************#
# Main function
#***************************************************************************************#
for (organ in ORGANS)
{
    for (strain in STRAINS)
    {
        # Read data set
        dta <- read.csv(file.path(DATA_PATH, paste0(organ, strain, '.csv')), header=T, stringsAsFactors=FALSE)
        
        # Read ranking table
        rk <- read.csv(file.path(RANK_PATH, paste0(organ, strain, '.csv')), header=T, stringsAsFactors=FALSE)
        rk <- rk[order(rk$rank, decreasing=T), ] # sort by rank
        rk <- rk[!duplicated(rk$Symbol), ] # remove duplicates i.e only keep top ranked probe
        
        # Create data frame that only contains top ranked probes
        df <- dta[dta$ProbeID %in% rk$ProbeID, ]
        
        # Some gene symbols have been updated since the release of:
        # http://emea.support.illumina.com/array/array_kits/mousewg-6_v2_expression_beadchip_kit/downloads.html
        # MouseWG-6_V2_0_R3_11278593_A.txt
        # e.g Il8RB is now Cxcr2 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=CXCR2)
        # In order to maximise the number of genes matching those in the signature matrix 
        # we will update these gene symbols manually here
        # TODO: update these gene symbols everywhere (i.e in cluster definitions too)
        df$Symbol[df$Symbol %in% 'Il8rb'] <- 'Cxcr2'
        df$Symbol[df$Symbol %in% 'LOC100044531'] <- 'Pirb'
        df$Symbol[df$Symbol %in% 'Lin28'] <- 'Lin28a'
        df$Symbol[df$Symbol %in% 'LOC100045268'] <- 'Gzma'
        df$Symbol[df$Symbol %in% '2010317E24Rik'] <- 'Sapcd2'
        df$Symbol[df$Symbol %in% 'Tm7sf4'] <- 'Dcstamp'
        df$Symbol[df$Symbol %in% 'Indo'] <- 'Ido1'
        df$Symbol[df$Symbol %in% 'LOC100048807'] <- 'Gp1ba'
        df$Symbol[df$Symbol %in% 'LOC100044182'] <- 'Cd200r4' 
        df$Symbol[df$Symbol %in% 'Tmem146'] <- 'Catsperd'
        df$Symbol[df$Symbol %in% 'LOC100045607'] <- 'Cd5' 
        df$Symbol[df$Symbol %in% 'LOC100048556'] <- 'Ccl12'
        df$Symbol[df$Symbol %in% 'Edg8'] <- 'S1pr5'
        df$Symbol[df$Symbol %in% 'LOC100044160'] <- 'Sema3d'
        df$Symbol[df$Symbol %in% 'LOC100038897'] <- 'Klra8'
        df$Symbol[df$Symbol %in% 'LOC100046930'] <- 'Sh2d1a'

        # Make symbol row name - so that later I can anti-log the data.frame
        df$ProbeID <- NULL
        rownames(df) <- df$Symbol
        df$Symbol <- NULL
        
        # Data should be in non-log space [see CIBERSORT docs]
        df <- 2^df
        df <- data.frame(Symbol=rownames(df), df)
        rownames(df) <- NULL
        
        # Write to disk
        write.table(df, file.path(OUT_PATH, paste0(organ, strain, '.txt')), sep="\t",
                    row.names=FALSE, quote=FALSE)
    }
}

# Write sig matrix as .txt file
df <- read.csv(SIG_PATH, header=T, row.names=1, check.names=F)
df <- data.frame(Symbol=rownames(df), df)
rownames(df) <- NULL
write.table(df, file.path(OUT_PATH, 'SigMatrix.txt'), sep="\t",
            row.names=FALSE, quote=FALSE)

#***************************************************************************************#
# Genes present in signature matrix but not in Jinwgwen's data
#***************************************************************************************#
# Cd300ld - no hit (have Cd300lg/lf/a etc.. but not ld)
# Ai839979 - no hit
# Rab44 - no hit (have Rasd1/2 but not Rasd3 --> Rab44)
# Lct - no hit
# Rasal3 - no hit (have Rasal2)
# Ighg2c - no hit (have Ighg though)
# Zbtbd6 - no hit
# Ighg2b - no hit (have Ighg though)
# Mzb1 - no hit
# Spns3 - no hit (have Spns1)
# Gm9706 - no hit
# Gimap3 - no hit (lots of Gimap but not 3)
# Pyhin1 - no hit 
# Fam71b - no hit (lots of Fam but not 71b)
# Ighg1 - no hit (have Ighg though)
# Herc6 - no hit (lots of Herc but not 6)
# Gm20199 - no hit
# Nmrk1 - no hit (only got Nrk)
# H2-Ea-ps - no hit
# Gm11346 - no hit
# Gm15708 - no hit
# Ear1 - no hit (Ear2 present)
# A430088P11Rik - no hit
# Trac - no hit
# Tpsb2 - no hit
# Endou - no hit
# Fam26f - no hit
# Trbj2-1 - no hit
# Mcpt-ps1 - no hit
# Lacc1 - no hit
# Ighg3 - no hit (have Ighg though)
# Fam129c - no hit (have Fam129a/b though)
# Tfec - no hit
# ProbeID: 6620722 mismatch between Jingwen's file [Tpsg1] and MouseWG annotation file [Cacna1h]
# Tcrg-C4 - no hit
# Cd163l1 - no hit (Cd163 present)
# Gm4951 - no hit
# Prss30 - no hit
# Hpgds - no hit
# Zg16 - no hit
# Nxpe2 - no hit 
# Gbp2b - no hit (but have Gbp2)
# 2500002B13Rik - no hit
# Gm14047 - no hit
# Trdc - no hit (but have Tcrd-V1)
# Gbp8 - no hit (lots of Gbp but not 8)
# Myl10 - no hit (Myl but not 10)
# Fam167a - no hit
# 4930428O21Rik - no hit
# Themis - no hit
# Gm11110 - no hit
# Mylk3 - no hit
# Trav12-3 - no hit
# Spef2 - no hit
# Trbj2-2 - no hit
# Pglyrp2 - no hit
# Gm5547 - no hit
# Trbv14 - no hit
# Khdc1b - no hit (have Khdc1c)
# 1810073O08Rik - no hit
# Trbv29 - no hit
# Art2a-ps - no hit (have Art2a though)
# Gm8369 - no hit
# Ddx43 - no hit
# Mucl1 - no hit
# Trbv5 - no hit
# Khdc1a - no hit (have Khdc1c though)
# A430108G06Rik - no hit
# 1700061F12Rik - no hit
# 5430437J10Rik - no hit