# susieR WDL 

## Overview: 
This WDL is designed to run on terra with the use of data tables.The inputs are the following:
1. `GenotypeDosages` file: This can be generated from a vcf using `bcftools dose`
2. `GenotypeDosageIndex`: Tabix `.tbi` file for the dosage file
3. `TensorQTLPermutations`: Permutation P-values output from tensorQTL
4. `PhenotypeBed`: File corresponding to the gene to be fine mapped
5. `CisDistance`: integer that determines the window to be added to either side of the TSS for fine-mapping
6. `PhenotypeID`: String containing the ID of gene or protein to be fine-mapped
7. `QTLCovariates`: Table containing the covariates used in QTL calling, this is the same file that in given to tensorQTL 
8. `SampleList`: List of sampleIDs that are being used in finemapping



## Tasks:
1. `PrepInputs`: Using `PhenotypeID` this task subsets all of the input files to data just pertaining to the selected gene.
2. `susieR`: 



