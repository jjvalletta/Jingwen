# Mice Microarray Data Analysis

Source code for analysing gene expression time-series data for mice experiments performed by Jingwen Lin/Jean Langhorne at 
Francis Crick Institute, London

* Authors: John Joseph Valletta and Mario Recker
* Data: Microarray Illumina Mouse WG6 v2 (45,281 probe sets representing 30,854 genes)
* Contacts: Jingwen Lin and Jean Langhorne, Francis Crick Institute, London	

# Data
Contains all data sets used in the analysis.

`acleaning.py`: Cleans up `Original` data to produce `Processed` and `Log2FC`

## Original
Data as received from Jingwen or downloaded from original source.

* `BL/SP-all.entities.txt`: Microarray data normalised to the median across all the chips and log2 transformed supplied by Jingwen Lin
* `BL/SP-2ANOVA.csv`: Two-way ANOVA across time-points from GeneSpring supplied by Jingwen Lin
* `MouseWG-6_V2_0_R3_11278593_A.txt`: Annotation file downloaded from Illumina [here](http://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mousewg-6/mousewg-6_v2_0_r3_11278593_a_txt.zip) to map probe IDs to gene symbols. **Note**:
	+ `Entrez_Gene_ID` --> `EntrezID` (e.g 212772)
	+ `Symbol` --> `GeneSymbol` (e.g Thrsp)
	+ `Probe_Id` --> `IlluminaID` (e.g )
	+ `Array_Address_Id` --> `ProbeID` (e.g 2600193)

## Processed
`Original` data cleaned up by `acleaning.py`

* `Blood/Spleen.csv`: Reduced versions of `BL/SP-all.entities.txt` as follows:
	1. Remove meta data columns apart from ProbeID, Symbol and data
	2. Fix issue with gene names see [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1044-7) and [here](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-80)
    3. Note: `ExcelSymbol` = original gene symbol (potentially incorrect), `Symbol` = corrected gene symbol

* `Blood/Spleen AS/CB.csv`: Same as `Blood/Spleen.csv` but split by strain, `ExcelSymbol` removed and columns renamed as Day.NReplicate i.e 3.2 -> day 3, replicate 2 

## Log2FC
* `Blood/Spleen AS/CB.csv`: Same as `Processed` data but with log 2 fold change computed (day 0 and day 12 naive mice pooled together) 

# Analysis

## Dependencies

* `GPy-1.7.7`, `paramz-0.7.4`, `GPclust-0.1.0`: Core modules for Gaussian Process Modelling ([Sheffield Machine Learning Group](https://github.com/SheffieldML)) 

## Top-level scripts

* `script.py`: Top-level script to cluster gene expression data (computationally expensive)

## Custom modules

All the scripts used to analyse the data

* `config.py`: Configuration file containing all constants used
* `input_output.py`: Input/output functions (creating folders, saving to pdf, pickle etc)
* `gaussian_process.py`: Wrapper functions for Gaussian Process Models
* `util.py`: Various utility functions