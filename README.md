# susieR WDL 

## Overview: 
This WDL is designed to run on terra with the use of data tables.The inputs are the following:
1. `GenotypeDosages` file: This can be generated from a vcf using `bcftools dose`\n
2. `GenotypeDosageIndex`: Tabix `.tbi` file for the dosage file\n
3. `TensorQTLPermutations`: Permutation P-values output from tensorQTL\n
4. `PhenotypeBed`: File corresponding to the gene to be fine mapped\n
5. `CisDistance`: integer that determines the window to be added to either side of the TSS for fine-mapping\n
6. `PhenotypeID`: String containing the ID of gene or protein to be fine-mapped\n
7. `QTLCovariates`: Table containing the covariates used in QTL calling, this is the same file that in given to tensorQTL\n 
8. `SampleList`: List of sampleIDs that are being used in finemapping, needs a header but can be any column name\n 



### Optional Inputs:
1. `MAF`: Float that specifies a MAF cutoff for variants that should be used in the genotype file. If this is supplied then AncestryMetadata will be needed, the workflow calculates the MAF in each population and tests wether the variant meets the MAF cutoff in any population. Note that doing this will remove individuals from MAF calculation that dont belong to any specific population.  
2. `AncestryMetadata`: file that contains the ancestry metadata for individuals being analyzed. Currently requires a column called `ancestry_pred_oth` to assign individuals to each ancestry.
3. `VariantList`: Single column file that contains variants formatted as `chr_pos_ref_alt`, if this is provided then it subsets variants by id in the genotype matrix. Providing this would supercede MAF option

## Tasks:
1. `PrepInputs`: Using `PhenotypeID` this task subsets all of the input files to data just pertaining to the selected gene.\n
2. `susieR`: fine-maps individual phenotypeID, this tasks takes in all of the subsetted files except for the bed file. Creating a `summarizedExperiment` object with just a single entry fails so this cant take in a subsetted bed file\n

## Outputs:
1. SusieParquet - contains fine-mapping information for all variants that are in credible sets\n
2. SusielbfParquet - contains log-bayes factor for each credible set (also present in SusieParquet)\n
3. FullSusieParquet - contains all fine-mapping information for every variant tested in a cis window\n



