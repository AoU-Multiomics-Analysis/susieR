# R Scripts

The workflow entrypoint scripts live in `R/scripts/`. Shared helper functions live in `R/utils/` and are copied into the Docker image at `/opt/r/lib`.

## `R/scripts/susie.R`

The main fine-mapping script. It loads genotype dosages, expression or phenotype data, covariates, and a list of phenotypes to fine-map, then runs susieR on each phenotype in the provided cis window. The feature coordinate is the midpoint of the phenotype BED interval, so 1 bp site/TSS intervals and already-expanded phenotype windows are centered consistently.

Results are written to three Parquet files and two diagnostic TSV files:

| Output | Description |
|---|---|
| `<out_prefix>.parquet` | Variants belonging to credible sets. Corresponds to the `SusieParquet` WDL output. |
| `<out_prefix>.lbf_variable.parquet` | Log-Bayes factors per credible set. Corresponds to `SusielbfParquet`. |
| `<out_prefix>.full_susie.parquet` | Fine-mapping statistics for all tested variants. Corresponds to `FullSusieParquet`. |
| `<out_prefix>.variant_position_summary.tsv` | Per-phenotype summary of tested variant counts, upstream/downstream counts, min/max variant distance from the feature, and requested cis-window size. Corresponds to `VariantPositionSummary`. |
| `<out_prefix>.covariate_explained_variants.tsv` | Per-phenotype list of variants almost fully explained by covariates. Corresponds to `CovariateExplainedVariants`. Contains headers even when no variants are flagged or no credible sets are found. |

Key command-line arguments:

| Argument | Description |
|---|---|
| `--genotype_matrix` | Tabix-indexed genotype dosage file. |
| `--expression_matrix` | BED-format phenotype or expression file. Columns 2-3 are converted to a midpoint feature coordinate before cis-window expansion. |
| `--phenotype_list` | TensorQTL permutation output that defines phenotypes to fine-map. |
| `--covariates` | Covariate matrix used in QTL calling. |
| `--sample_meta` | List of sample IDs to include. |
| `--out_prefix` | Output file prefix. |
| `--cisdistance` | Cis window in bp added to each side of the phenotype BED interval midpoint. |
| `--reuse_genotype_matrix` | Reuse one residualized genotype matrix when selected phenotype windows merge into one region. |
| `--select_top_phenotype_per_cluster` | Fine-map only the strongest FDR-passing phenotype per parsed LeafCutter `clu_*` cluster. |
| `--top_phenotype_pvalue_column` | Preferred TensorQTL p-value column for choosing the representative phenotype per cluster. Defaults to `qval`. |
| `--MAF` | Optional minor allele frequency cutoff. |
| `--AncestryMetadata` | Optional ancestry assignment file for per-population MAF filtering. |
| `--VariantList` | Optional single-column file of variants formatted as `chr_pos_ref_alt` to restrict analysis. |

Ordinary covariate adjustment uses a factored calculation that avoids a dense
sample-by-sample residual matrix. The same adjustment is used for phenotypes
and genotypes. The covariate reader uses a faster path for ordinary TSV input
and uses the original parser for unsupported or ambiguous formats. Marginal
z-scores are calculated in bulk after the SuSiE fit. Unusual or singular input
uses the native calculation. The fit still uses 10 effects.

For compatibility, a residual SD at most `1e-8` times the original column
scale triggers the original dense adjustment for both genotypes and phenotype.
These fits also use native marginal z-scores. The scale is the larger of the
original SD and absolute mean; an all-zero input does not trigger this guard.
This threshold selects the calculation method. It does not filter variants.
If a later phenotype triggers this guard after genotype reuse has started,
the runner rebuilds the full merged window with the original calculation and
reuses that result. This rare path retains the earlier fast cache as well as
the rebuilt result and dense matrix, so it can use more memory than the original
method. It preserves the previous numerical behavior for almost fully
covariate-explained inputs.

The diagnostic TSV contains `molecular_trait_id`, `variant`, `input_scale`,
`residual_sd`, and `residual_sd_ratio`. It lists a variant when its input scale
is positive and its residual SD is at most `1e-8` times that scale. The input is
the dosage after the existing row-mean imputation and sample-mean subtraction,
for the final aligned samples. The scale is the larger of its sample SD and
absolute sample mean. The residual SD is measured after the selected covariate
adjustment, including the original dense calculation when the guard uses it.
`residual_sd_ratio` is the residual SD divided by the input scale. The list can
include variants outside credible sets. Each saved fit also contains this
diagnostic table as `covariate_explained_variants`.

The primary credible-set Parquet has a logical
`has_covariate_explained_variant` column. It is `TRUE` on every row of a
credible set if any original member is in that phenotype's diagnostic list.
The internal credible-set summary has the same column. Fits saved before this
diagnostic field was added use `FALSE`. These annotations do not filter variants
or change the existing fit statistics. The full-SuSiE Parquet columns and the
checkpointed runner's table columns stay the same.

One CPU remains the workflow default. The local synthetic benchmark used
10,000 samples and 20,000 variants. Its results apply to that local test and
do not establish cloud performance. The complete workflow has not been tested
on Terra for this change.

Run the regression tests from the repository root:

```sh
Rscript tests/run_tests.R
```

For a local test library, set `BENCH_LIBRARY=/path/to/library`. By default,
tests use helpers in `R/utils`. Set `SUSIER_FUNCTIONS_PATH` to test another
helper directory. GitHub Actions runs the tests against `/opt/r/lib` in the
built image, with the repository test files mounted read-only.

## `R/scripts/run_checkpointed_window_susie.R`

This runner processes one prepared window. It trusts `window_phenotypes.tsv`.
The shared manifest contains `outcome_key`, the original `phenotype_id`, and
`p_value`. The runner selects rows with a finite P value and uses `outcome_key`
to match the combined phenotype BED. It uses every usable supplied variant. It
does not calculate a phenotype-centered interval.

Use `--window-id`, `--window-dosage`, `--window-phenotypes`, and
`--phenotype-data` for the prepared input bundle. Use `--covariate-files` and
`--covariate-modalities` for matching comma-separated arrays. Use
`--checkpoint-root` for a checkpoint prefix. In production, give a writable
`gs://` prefix. Use `--output-dir` for local final WDL outputs.

The runner writes a complete fitted `susie` object to an RDS file for each fit.
`window_fit_index.tsv` locates every fitted RDS. It also records the terminal
state for every ordered phenotype.

The checkpoint root contains one directory for the analysis ID and one
directory for the window ID. It contains `window_manifest.json`,
`window_fit_index.tsv`, final result tables, and one `phenotypes/<key>/`
directory for each phenotype. Each phenotype directory contains
`fit_manifest.json`, its RDS, and its result tables.

The runner uploads GCS payloads before `fit_manifest.json`. It uploads
`window_manifest.json` last. Normal resume reads only the window cursor and the
latest fit. It does not list the full checkpoint directory. `SKIPPED` and
`NONCONVERGED` are terminal states. They do not block later phenotypes.

The runner loads raw window dosage once. It prepares one residualized genotype
for each exact retained-sample mask and covariate design. It reuses that cache
for matching phenotypes. Each fit manifest records raw, overlap, and retained
sample counts; input and retained variant counts; removed covariates; design
identity; cache identity; runtime source hashes; and image identities.

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

Optional allelic fold-change data can be supplied with `--AllelicFoldChangeData`. The table must include `pid`, `sid`, `log2_aFC`, `log2_aFC_lower`, and `log2_aFC_upper`; `sid` values such as `chr1:169874249_C_T` are normalized to the SuSiE `chr1_169874249_C_T` variant format. The output keeps the three `log2_aFC` columns and drops aFC key/location columns.

## `R/scripts/ComputeAncestrySkew.R`

Computes ancestry skew for variants in an annotated TSV or gzipped TSV. The script filters variants by PIP, recalculates the GVS max subpopulation from minor allele frequency, and runs Fisher tests comparing the max subpopulation to the remaining cohort. Fisher tests add a pseudocount of one to each valid count cell so odds ratios remain finite when one cell is zero.

The script also reports a second ancestry skew calculation with admixed samples removed from the cohort background. By default, the admixed subpopulation is `oth`, but this can be changed with a comma-separated value such as `oth,amr`.

Required annotation columns:

| Column pattern | Description |
|---|---|
| `variant` | Variant identifier. |
| `pip` | Posterior inclusion probability used for filtering. |
| `gvs_all_ac`, `gvs_all_an` | Cohort-level alternate allele count and allele number. |
| `gvs_<subpop>_af` | Subpopulation allele frequency columns, such as `gvs_afr_af`, `gvs_amr_af`, and `gvs_eur_af`. |
| `gvs_<subpop>_ac`, `gvs_<subpop>_an` | Matching subpopulation AC/AN columns for every `gvs_<subpop>_af` column. |

Key command-line arguments:

| Argument | Description |
|---|---|
| `--AnnotationData` | Annotated variant table. Plain TSV and gzip-compressed TSV inputs are supported. |
| `--OutputPrefix` | Output file prefix. |
| `--PipThreshold` | Minimum PIP threshold used to select variants. Defaults to `0.9`. |
| `--AdmixedSubpops` | Comma-separated GVS subpopulation labels to remove for the no-admixed skew calculation. Defaults to `oth`. |
| `--KeepInputColumns` | Keep all input annotation columns and append only the compact ancestry-skew result columns. |

Output: `<OutputPrefix>.AncestrySkew.tsv.gz`.

Standalone output columns are `variant`, `pip`, `gvs_max_subpop`, recomputed MAF-based `gvs_max_af`, `gvs_odds_ratio`, `gvs_p_value`, `gvs_no_admixed_odds_ratio`, and `gvs_no_admixed_p_value`.
