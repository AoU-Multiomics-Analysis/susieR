# R Scripts

The workflow entrypoint scripts live in `R/scripts/`. Shared helper functions live in `R/utils/` and are copied into the Docker image at `/opt/r/lib`.

## `R/scripts/susie.R`

The main fine-mapping script. It loads genotype dosages, expression or phenotype data, covariates, and a list of phenotypes to fine-map, then runs susieR on each phenotype in the provided cis window.

Results are written to three Parquet files:

| Output | Description |
|---|---|
| `<out_prefix>.parquet` | Variants belonging to credible sets. Corresponds to the `SusieParquet` WDL output. |
| `<out_prefix>.lbf_variable.parquet` | Log-Bayes factors per credible set. Corresponds to `SusielbfParquet`. |
| `<out_prefix>.full_susie.parquet` | Fine-mapping statistics for all tested variants. Corresponds to `FullSusieParquet`. |

Key command-line arguments:

| Argument | Description |
|---|---|
| `--genotype_matrix` | Tabix-indexed genotype dosage file. |
| `--expression_matrix` | BED-format phenotype or expression file. |
| `--phenotype_list` | TensorQTL permutation output that defines phenotypes to fine-map. |
| `--covariates` | Covariate matrix used in QTL calling. |
| `--sample_meta` | List of sample IDs to include. |
| `--out_prefix` | Output file prefix. |
| `--cisdistance` | Cis window in bp added to each side of the TSS. |
| `--MAF` | Optional minor allele frequency cutoff. |
| `--AncestryMetadata` | Optional ancestry assignment file for per-population MAF filtering. |
| `--VariantList` | Optional single-column file of variants formatted as `chr_pos_ref_alt` to restrict analysis. |

## `R/scripts/ComputeR2Susie.R`

Runs cross-validation to evaluate fine-mapping predictive performance. For a given gene or phenotype, it runs susieR on each CV fold, generates predicted vs. observed expression values, and outputs a TSV of predictions.

Two runs are performed per fold: one using all variants and one restricted to variants passing a hardcoded 1% allele-frequency filter applied via `filterMAF()`. The 1% analysis can also be restricted by `VariantList`.

Input expression values should be CPM-normalized but not rank-normalized. Rank normalization is applied within each fold to avoid feature leakage. Use `PrepSusieCVPCs.R` to generate the required `CVMetadata` file.

Output: `<out_prefix>_SusiePredictions.tsv`, with per-fold predicted and observed expression values, allele-frequency threshold, and fold annotations.

## `R/scripts/PrepSusieCVPCs.R`

Prepares the cross-validation metadata and principal components needed by `ComputeR2Susie.R`.

The script:

1. Assigns samples to CV folds using ancestry information.
2. CPM-normalizes raw count data using edgeR.
3. Computes expression PCs within each fold to be used as covariates during CV fine-mapping.

The input should be raw counts rather than rank-normalized values, because normalization is handled internally per fold to prevent data leakage.

## `R/scripts/MergeSusie.R`

Merges multiple sharded susieR output Parquet files, such as outputs from scatter-gather workflows, into a single combined file. It accepts a text file listing paths to all shards and outputs a merged Parquet and a gzipped TSV.

| Argument | Description |
|---|---|
| `--FilePaths` | Path to a single-column text file listing all input Parquet shard paths. |
| `--OutputPrefix` | Prefix for output files. |

Outputs:

| Output | Description |
|---|---|
| `<OutputPrefix>_SusieMerged.parquet` | Merged Parquet output. |
| `<OutputPrefix>_SusieMerged.tsv.gz` | Merged gzipped TSV output. |

## `R/scripts/AnnotateSusie.R`

Annotates a merged Susie TSV with external annotation resources used by the post-analysis workflows.
