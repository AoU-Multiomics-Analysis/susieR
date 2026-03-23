# susieR WDL

## Overview

This repository provides WDL workflows and R scripts for running [SusieR](https://github.com/stephenslab/susieR) fine-mapping on Terra. The workflows are designed to work with data tables on Terra and perform Bayesian fine-mapping of QTL signals identified by tensorQTL.

## Repository Structure

```
.
├── scripts/                    # Standalone R scripts
│   ├── susie.R                 # Main fine-mapping script
│   ├── ComputeR2Susie.R        # Cross-validation R² computation
│   ├── PrepSusieCVPCs.R        # Prepare cross-validation PCs/covariates
│   └── merge_susie.R           # Merge sharded susie output files
├── susieR/                     # WDL workflows and supporting files
│   ├── susieR.wdl              # Full pipeline: input prep + fine-mapping
│   ├── susieRonly.wdl          # Fine-mapping only (no input prep)
│   ├── prepInputsSusieR.wdl    # Input preparation only
│   ├── ComputeR2Susie.wdl      # Cross-validation R² workflow
│   ├── Dockerfile              # Docker image definition
│   └── utils/                  # R utility functions used by the scripts
└── .github/workflows/
    └── docker-image.yml        # CI/CD: build and push Docker image
```

---

## Scripts

### `scripts/susie.R`

The main fine-mapping script. It loads genotype dosages, expression/phenotype data, covariates, and a list of phenotypes to fine-map, then runs SusieR on each phenotype in the provided cis window. Results are written to three Parquet files:

- `<out_prefix>.parquet` — variants belonging to credible sets (corresponds to the `SusieParquet` WDL output)
- `<out_prefix>.lbf_variable.parquet` — log-Bayes factors per credible set (corresponds to `SusielbfParquet`)
- `<out_prefix>.full_susie.parquet` — fine-mapping statistics for all tested variants (corresponds to `FullSusieParquet`)

**Key inputs (command-line arguments):**

| Argument | Description |
|---|---|
| `--genotype_matrix` | Tabix-indexed genotype dosage file |
| `--expression_matrix` | BED-format phenotype (expression) file |
| `--phenotype_list` | TensorQTL permutation output (defines phenotypes to fine-map) |
| `--covariates` | Covariate matrix used in QTL calling |
| `--sample_meta` | List of sample IDs to include |
| `--out_prefix` | Output file prefix |
| `--cisdistance` | Cis window (bp) added to each side of the TSS |
| `--MAF` | (Optional) Minor allele frequency cutoff |
| `--AncestryMetadata` | (Optional) Ancestry assignment file for per-population MAF filtering |
| `--VariantList` | (Optional) Single-column file of variants (`chr_pos_ref_alt`) to restrict analysis |

---

### `scripts/ComputeR2Susie.R`

Runs cross-validation (CV) to evaluate fine-mapping predictive performance (R²). For a given gene/phenotype, it runs SusieR on each CV fold, generates predicted vs. observed expression values, and outputs a TSV of predictions. Two runs are performed per fold: one using all variants and one restricted to variants passing a hardcoded 1% allele-frequency filter applied via `filterMAF()` (optionally further restricted by `VariantList`).

**Note:** Input expression values should be CPM-normalized but **not** rank-normalized; rank normalization is applied within each fold to avoid feature leakage. The `PrepSusieCVPCs.R` script should be used to generate the required `CVMetadata` file.

**Output:** `<out_prefix>_SusiePredictions.tsv` — per-fold predicted and observed expression values with allele-frequency threshold and fold annotations.

---

### `scripts/PrepSusieCVPCs.R`

Prepares the cross-validation metadata and principal components needed by `ComputeR2Susie.R`. It:

1. Assigns samples to CV folds using ancestry information.
2. CPM-normalizes raw count data using edgeR.
3. Computes expression PCs within each fold to be used as covariates during CV fine-mapping.

**Note:** The input to this script should be raw counts (not rank-normalized), as normalization is handled internally per fold to prevent data leakage.

---

### `scripts/merge_susie.R`

Merges multiple sharded SusieR output Parquet files (e.g., from scatter-gather workflows) into a single combined file. Accepts a text file listing paths to all shards and outputs a merged Parquet and a gzipped TSV.

**Command-line arguments:**

| Argument | Description |
|---|---|
| `--FilePaths` | Path to a single-column text file listing all input Parquet shard paths |
| `--OutputPrefix` | Prefix for output files |

**Outputs:**
- `<OutputPrefix>_SusieMerged.parquet`
- `<OutputPrefix>_SusieMerged.tsv.gz`

---

## WDL Workflows

### `susieR/susieR.wdl` — Full Pipeline

Runs both input preparation and fine-mapping in a single workflow. First calls `PrepInputs` to subset all input files to the region around the target phenotype, then calls `susieR` to perform fine-mapping on those subsetted files.

**Inputs:**

| Input | Type | Description |
|---|---|---|
| `GenotypeDosages` | File | Tabix-indexed genotype dosage file (can be generated from VCF using `bcftools dose`) |
| `GenotypeDosageIndex` | File | Tabix `.tbi` index for the dosage file |
| `TensorQTLPermutations` | File | Permutation p-values output from tensorQTL |
| `PhenotypeBed` | File | BED file for the gene/phenotype to be fine-mapped |
| `CisDistance` | Int | Window size (bp) added to each side of the TSS |
| `PhenotypeID` | String | ID of the gene or protein to fine-map |
| `QTLCovariates` | File | Covariate table used in QTL calling (same file given to tensorQTL) |
| `SampleList` | File | Sample IDs used in fine-mapping (requires a header) |
| `susie_rscript` | File | Path to the `susie.R` script |
| `memory` | Int | Memory (GB) for the fine-mapping task |
| `NumPrempt` | Int | Number of preemptible retries |
| `OutputPrefix` | String | Prefix for merged output files |
| `MAF` | Float | Minor allele frequency cutoff (required by the WDL; pass `0` to include all variants) |

**Outputs:**

| Output | Description |
|---|---|
| `SusieParquet` | Variants in credible sets |
| `SusielbfParquet` | Log-Bayes factors per credible set |
| `FullSusieParquet` | Fine-mapping statistics for all tested variants |
| `SubsetBed` | Subsetted BED file for the phenotype region |
| `SubsetDosages` | Subsetted genotype dosage file |
| `SubsetDosagesIndex` | Index for the subsetted dosage file |

---

### `susieR/susieRonly.wdl` — Fine-Mapping Only

Runs only the fine-mapping step (no input preparation). Accepts pre-subsetted or full input files. Useful when inputs are already prepared (e.g., from `prepInputsSusieR.wdl`) or when re-running fine-mapping with different parameters.

**Optional inputs (in addition to the required inputs shared with `susieR.wdl`):**

| Input | Type | Description |
|---|---|---|
| `MAF` | Float? | Minor allele frequency cutoff |
| `VariantList` | File? | Single-column file of variants (`chr_pos_ref_alt`) to restrict analysis |
| `AncestryFile` | File? | Ancestry metadata for per-population MAF filtering |
| `AdditionalGenotypesBed` | File? | Additional genotype BED file |

**Outputs:** Same as `susieR.wdl` (SusieParquet, SusielbfParquet, FullSusieParquet).

---

### `susieR/prepInputsSusieR.wdl` — Input Preparation Only

Subsets genotype dosages, the phenotype BED file, and TensorQTL permutation results to just the region surrounding the target phenotype (`PhenotypeID`). This is the first step of `susieR.wdl` exposed as a standalone workflow, useful for preparing inputs once before running fine-mapping multiple times.

**Outputs:**

| Output | Description |
|---|---|
| `SubsetBed` | Subsetted BED file for the phenotype region |
| `SubsetDosages` | Subsetted genotype dosage file |
| `SubsetDosagesIndex` | Tabix index for the subsetted dosage file |

---

### `susieR/ComputeR2Susie.wdl` — Cross-Validation R²

Runs the `ComputeR2Susie.R` script in a WDL workflow to evaluate fine-mapping predictive R² via cross-validation. Requires a `CVMetadata` file generated by `PrepSusieCVPCs.R`.

**Additional inputs (beyond standard fine-mapping inputs):**

| Input | Type | Description |
|---|---|---|
| `CVMetadata` | File | Cross-validation fold and PC metadata generated by `PrepSusieCVPCs.R` |
| `VariantList` | File | Variants to include in the 1% AF threshold analysis |

**Output:** `<PhenotypeID>_SusiePredictions.tsv` — per-fold predicted and observed expression values.

---

## GitHub Actions Workflow

### `.github/workflows/docker-image.yml` — Docker Image CI

Automatically builds and pushes the Docker image (defined in `susieR/Dockerfile`) to the GitHub Container Registry (`ghcr.io`) on every push or pull request to the `main` branch.

- **Trigger:** Push or pull request to `main`
- **Registry:** `ghcr.io/aou-multiomics-analysis/susier`
- **Tags:** Automatically generated from branch/tag metadata using `docker/metadata-action`

The Docker image contains all R dependencies and utility scripts required by the WDL workflows.

---

## Common Inputs (Shared Across Workflows)

| Input | Description |
|---|---|
| `GenotypeDosages` | Genotype dosage file; can be generated from a VCF using `bcftools dose` |
| `GenotypeDosageIndex` | Tabix `.tbi` index for the dosage file |
| `TensorQTLPermutations` | Permutation p-values output from tensorQTL |
| `PhenotypeBed` | BED file for the gene/phenotype to be fine-mapped |
| `CisDistance` | Window (bp) added to each side of the TSS for fine-mapping |
| `PhenotypeID` | ID of the gene or protein to fine-map |
| `QTLCovariates` | Covariate table used in QTL calling (same file given to tensorQTL) |
| `SampleList` | List of sample IDs used in fine-mapping (requires a header) |

### Optional Inputs

| Input | Description |
|---|---|
| `MAF` | Float; MAF cutoff for variants. Requires `AncestryMetadata`. MAF is calculated per population and a variant must pass the cutoff in at least one population. Note: individuals not assigned to any population are excluded from the MAF calculation. |
| `AncestryMetadata` | Ancestry metadata file; requires a column `ancestry_pred_oth` for population assignment |
| `VariantList` | Single-column file of variants formatted as `chr_pos_ref_alt`; restricts analysis to listed variants. Takes precedence over `MAF` filtering. |



