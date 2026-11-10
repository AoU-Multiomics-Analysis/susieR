# SuSiE 2.0 PGEN workflow

[SusieR2.wdl](../workflows/SusieR2.wdl) is a WDL 1.0 workflow for Terra's managed Cromwell service. It plans batches, applies shared PLINK QC, creates the batch PGEN subsets, fits each selected outcome, and combines the batch tables. Each original input triplet is localized once by its preparation task. Each SuSiE task localizes only its own smaller batch PGEN/PVAR/PSAM files. The default batch size is 50 genes. Set `SusieR2.BatchSize` to another positive integer, such as 100, to change this limit. All outcomes for one gene stay in one batch.

The complete workflow has **not completed validation on Terra**. Local WDL checks and command tests do not establish Terra compatibility. No cloud jobs were submitted for these checks.

## Inputs

Copy [the input example](../examples/inputs/SusieR2.inputs.json) and replace each placeholder cloud path. Supply a headered PGEN TSV using [the manifest example](../examples/inputs/SusieR2.pgen_manifest.tsv). Its four columns, in order, are `chromosome`, `pgen`, `pvar`, and `psam`. Add one row for each chromosome that contains a selected phenotype. A shared PSAM path can appear in several rows. Genotype paths must be absolute cloud or local paths; they are not resolved relative to the manifest.

| Input | Type | Use |
|---|---|---|
| `PgenManifest` | File | Headered TSV with chromosome and matching PGEN, PVAR, and PSAM locations. Replaces the previous `ChromosomePgens` array input. |
| `PhenotypeBed` | File | Phenotype BED with its original header, first four feature columns, and sample columns. |
| `QtlResults` | File | QTL results for the existing phenotype selection rules (`qval < 0.05`). |
| `Covariates` | File | Covariate table for the existing data loader. |
| `SampleList` | Optional File | Explicit sample list for the existing data loader. If omitted, use the genotype/phenotype sample-ID intersection in PSAM order. |
| `LegacyMissingCode` | String, required | Use `NA` or `-1`. Specify how the legacy genotype code represents missing values. The workflow does not infer this policy. |
| `GeneMap` | Optional File | TSV with `phenotype_id` and `gene_id`. Without this input, each selected phenotype is its own gene and outcome. |
| `Ancestry` | Optional File | Ancestry table for ancestry MAF filtering. |
| `VariantList` | Optional File | Variant selection list, with the existing loader's selection rules. |
| `AdditionalGenotypes` | Optional File | Additional genotype BED for the existing data loader. |
| `MAC` | Integer, default `0` | Minimum minor allele count for PGEN variants. `0` disables this filter. |
| `MinCallRate` | Float, default `0.95` | Minimum fraction of input samples with called genotypes. For diploid autosomes this is `AN / (2 × sample count)`. `0` disables this filter. |

Initial PGEN support is for diploid, biallelic variants on autosomes 1–22. The planner accepts chromosome aliases such as `1` and `chr1`. A gene cannot span chromosomes. The planner rejects missing chromosome inputs, duplicate chromosome IDs, ambiguous mappings, and duplicate phenotype IDs.

The phenotype BED must use tabs and have a header. Its first four columns are chromosome, start, end, and phenotype ID. The chromosome header can be `chr`, `chrom`, or `chromosome`, with an optional `#`. Coordinate headers are `start` and `end`; `chromStart` and `chromEnd` are also accepted. The reader skips blank lines and comment, `track`, or `browser` lines before the header.

Plain text, gzip, and BGZF files use the same parser, regardless of the filename suffix. Quoted IDs and surrounding field spaces are decoded consistently during planning, splitting, and PGEN batch loading. Original sample values and selected BED row text are preserved. Each data row must have the same number of tab-delimited fields as the header. Embedded tabs and line breaks within fields are not supported.

The planner checks the header, metadata, row count, and selected-row counts again during splitting. Errors report the failing line or metadata row. An ID mismatch alone does not establish that the file changed.

The planner validates the manifest without opening its genotype locations and preserves its row order. It writes unique selected PGEN/PVAR/PSAM triplets and integer indices that link each batch to its source and subset. A triplet shared by several chromosomes or batches has one preparation task. Different companion files make distinct inputs. Unused manifest entries do not start preparation. WDL converts the source locations into typed Files. Preparation receives those Files and an explicit array of worklist Files. Its batch subset outputs remain typed Files through localization by the fit task. The planner's JSON index arrays are output metadata, not script argument wrappers.

When `SampleList` is omitted, sample selection uses IDs rather than phenotype values. Genotype-only and phenotype-only samples are excluded, while phenotype missing values retain their existing handling. Covariates must cover every retained sample; they do not silently reduce the cohort. A supplied list keeps its existing order and selection behavior. The PGEN reader and ancestry MAF calculation retain the complete genotype cohort before the phenotype subset. Additional-genotype sample restrictions still apply afterward. Provenance records the initial selection mode, PSAM/list order rule, sample count, and cohort hash before any additional-genotype restriction.

PGEN-only fits use `pgenlibr::ReadList(meanimpute = TRUE)` with the exact analysis
samples. The fitting helper skips repeated R imputation. A native all-missing
variant error removes that variant and retries the affected read. Other native
errors propagate. Existing empty-window handling is unchanged.

Each batch selects a prepared-genotype cache mode before reading genotypes.
The cache stores mean-imputed, centered, covariate-adjusted columns; phenotype
residualization and SuSiE fitting remain separate for each feature. Selection
uses physical PGEN indices and the actual gene/outcome execution order. It does
not depend on the assay name. Cache contents do not enter saved fitted models.

`GenotypeCacheMode = "auto"` first selects whole-batch caching when the combined
variants fit the budget and reading them once avoids at least 50% of repeated
variant reads. Otherwise it selects rolling caching when the largest window
fits and retaining the previous window avoids at least 50% of reads. Otherwise
it uses no caching. Whole mode prepares each distinct batch variant once.
Rolling mode releases variants outside the next window, retains prepared
overlap, and reads and prepares only new variants. No-cache mode uses the
existing per-feature path.

The automatic cache limit is 20% of `BatchMemoryGB`. Set
`GenotypeCacheMemoryGB` to a positive decimal-GB limit to replace it. The estimate
includes numeric cache columns (`8 * samples * variants` bytes) and bounded
preparation buffers. The fit's working memory is separate; these rules are not a
peak-memory guarantee. Explicit modes `whole`, `rolling`, and `none` support
controlled comparisons. Forced caching must fit its memory limit. Automatic
mode uses no caching with active ancestry MAF, a variant list, or additional
genotypes; forced caching rejects these preparation paths.

All cache modes use the individual-data SuSiE fitter. Whole and rolling caches
reuse prepared genotype columns; they do not build a shared genotype
cross-product or select a different fitting method. Each outcome retains its
own phenotype residualization and SuSiE fit.

Cached fits preserve the existing numerical safeguard for a genotype or
phenotype almost fully explained by covariates. Such a feature uses the existing
full-window calculation. Other features can continue using the prepared cache.
The cache is released when the batch returns or fails. Logs and
`GenotypeCacheReports` record the selected mode, reason, budget, estimated
storage, read/preparation counts, reuse, eviction, and guarded recomputations.

For cached fits, cache access and preparation are included in the
`data_preparation` timer. Loading, preparation, and fitting events remain in the
outcome logs and tables. A guarded recomputation adds its stages to the earlier
cache-check time. No-cache timing scopes retain their existing behavior.

Active ancestry MAF filtering needs the original full-cohort calls before
imputation. That input path retains the raw read, filters the variants, and
then reads variants with missing calls through native imputation on the fit
samples. Inputs with additional genotypes retain their existing merged-matrix
preparation. Temporary subset readers close after each call. `ReadList()` uses
stored dosages, so native imputation expects the input dosages to match the
integer hard calls. PGEN-only filtering uses `NA` missing calls before native
imputation, so the legacy `-1` sentinel cannot cancel observed ALT counts.

## Settings

| Input | Default |
|---|---:|
| `BatchSize` | 50 genes |
| `CisDistance` | 1000000 bases |
| `BaseSeed` | 1 |
| `MemoryDebug` | false |
| `GenotypeCacheMode` | auto |
| `GenotypeCacheMemoryGB` | 0 (20% of task memory) |
| `MAF` | 0 |
| `MAC` | 0 (off) |
| `MinCallRate` | 0.95 |
| `SelectTopPhenotypePerCluster` | false |
| `TopPhenotypePvalueColumn` | `qval` |
| `DockerImage` | `ghcr.io/aou-multiomics-analysis/susier/pgen:main` |
| `Cpu` | 1 |
| `Threads` | 0 (detect the task CPU allowance) |
| `BatchMemoryGB` / `BatchDiskGB` | 16 / 100 |
| `PlinkMemoryMB` | 1024 |
| `PlanMemoryGB` / `PlanDiskGB` | 8 / 50 |
| `QcMemoryGB` / `QcDiskGB` | 8 / 100 |
| `AggregateMemoryGB` / `AggregateDiskGB` | 8 / 50 |
| `PreemptibleAttempts` / `MaxRetries` | 3 / 1 |

Select an image digest through `DockerImage` when a fixed image is required. `QcMemoryGB` and `QcDiskGB` now apply to preparation. Its disk must hold the original triplet, the shared QC output when enabled, and all batch subsets for that input. `BatchDiskGB` covers one localized subset and its gene outputs. Keep `PlinkMemoryMB` below the preparation memory limit. Preparation uses one PLINK thread. The number of subset triplets equals the existing planned batch count. Larger `BatchSize` values reduce file count and increase each batch's input size. Preparation logs its batch count, variant counts and total subset bytes.

`Cpu` requests resources from Terra. `Threads = 0` detects the CPU capacity available to the batch process, including CPU affinity and container quotas. A positive `Threads` value sets a limit, capped at the detected capacity. For the local batch CLI, use `--threads`. The worker applies this count to supported BLAS/OpenMP operations, data.table, and Arrow. Gene fits remain sequential. Standalone batch extraction also uses this count for PLINK; prepared workers do not run PLINK. The log and provenance report the detected limits, selected count, BLAS library, and effective thread counts. A serial BLAS library can remain at one thread. The planner, preparation, and aggregate tasks keep their existing thread settings.

CPU detection does not request a larger VM or bypass a task CPU limit. Set `Cpu` higher when the task needs more allocated CPU capacity. Multithreaded numerical libraries can change the last floating-point bits. Set `Threads = 1` for the previous thread setting.

Detailed memory reports are off by default. Set `SusieR2.MemoryDebug` to `true` to enable object-size reports and two extra garbage collections per fit. For a local PGEN batch, add `--memory-debug`. The legacy `susie.R` command accepts `--memory_debug true`. Direct R callers can use `options(susier.memory_debug = TRUE)`. The default keeps progress and error logs, variant diagnostics, and the three garbage collections after large matrices are removed. This setting does not change model calculations. A batch restores the caller's R option when it returns or fails.

The planner sorts genes by chromosome and position. It derives each feature position from `floor((BED start + BED end) / 2)`. Cis windows use inclusive coordinates, with the query start clamped to 1. PLINK uses `--extract range` for these 1-based regions. Disjoint windows retain their gaps. Extraction retains the complete genotype sample cohort. Ancestry MAF filtering therefore runs before the phenotype sample subset. The fit keeps the existing phenotype missing-value behavior.

`MinCallRate` defaults to `0.95`: a variant must have called genotypes in at least 95% of the input PGEN samples. A value exactly at the threshold passes. `MAC` defaults to `0`; set it to a positive integer, such as `"SusieR2.MAC": 5`, to require a minor allele count of at least 5 as well. These filters use the complete input cohort, including nonfounders and samples without phenotypes, before imputation or phenotype sample selection.

Preparation uses PLINK's native `--mac`, `--geno` and `--nonfounders` filters. `--geno` receives `1 - MinCallRate`, or `0.05` by default. Its missingness rule uses hard calls. Native `--mac` uses dosages when present; for hardcall-only inputs these are the same allele counts. The QC table keeps descriptive hardcall `AC`, `AN` and `MAC` values and adds `FILTER_MAC` for the native allele count used by the filter. PVAR INFO annotations do not supply these counts. The retained native PGEN records preserve phase and dosage data.

Each unique selected triplet has one shared QC run, followed by one extraction per assigned batch. QC uses three PLINK commands: count reporting, native filter selection, and writing the retained records. Temporary row IDs prevent duplicate or missing raw variant IDs from selecting the wrong records. Retained PVAR records and PSAM contents are preserved. Fresh counts are in `GenotypeQcCounts`; original INFO annotations remain unchanged. Unsupported records outside the analysed windows are kept with QC marked as not applied. The existing reader still rejects unsupported variants if an analysed window selects them.

Set both `MAC = 0` and `MinCallRate = 0` to skip shared QC. Preparation still extracts the planned subsets from the original files. `MAC` must be a nonnegative integer; `MinCallRate` must be from 0 to 1. A filter that removes all supported variants stops with a clear error. SuSiE workers open their prepared subsets directly and do not run PLINK or repeat QC. The existing ancestry MAF or variant-list filter runs during fitting; it cannot restore a variant removed by shared QC. `AdditionalGenotypes` keeps its existing handling.

Gene provenance binds the shared QC report, source file records, exact worklist, shard hashes, cis windows and ordered sample cohort. A mismatch is rejected before phenotype loading or fitting. These reports are structured outputs, not JSON argument configurations. For standalone preparation, call `prepare_pgen_susie.R` with named genotype arguments and a task-local `--worklists-file` containing readable local worklist paths. Use a new or empty output directory; existing products are rejected to prevent stale shards or QC reports. Give its subset triplet and `--prepared-batch-manifest` to the batch CLI, with `--mac 0`. When preparation applied shared QC, also supply its `--genotype-qc-report`; the worker requires this companion file to verify the recorded thresholds. Direct batch calls without the preparation manifest retain their previous extraction behavior.

Seeds depend on the base seed, original gene ID, and original outcome ID. They do not depend on batch size, input order, or output paths. Changing batch size preserves the seed for each outcome. Each RDS contains a named list of that gene's outcome fits.

The fitting stage reuses identical matrix-vector products, expected residual-error values, and raw column standard deviations within one fit. It also uses `matrixStats::colVars()` for SuSiE's constant-column warning check. The model settings, optimizer, convergence checks, and iteration logs are unchanged. Each fit uses private copies of the audited SuSiE functions; the installed namespace is unchanged. Cached input references are cleared after success or failure. Unrecognized package sources, unsupported matrix types, and unsupported model modes use the original SuSiE path. Source guards match the pinned image source; local source checks do not constitute an image or Terra execution test.

Post-fit marginal scores reuse the same raw input means, standard deviations, and finite-value validation when the complete input matrix matches. Existing numerical fallback rules remain in place. These statistics stay local to the fit and are excluded from saved models. Each PGEN outcome's extracted tables are shared by its credible-set count and final output formatting; scientific output tables keep their existing schemas and ordering.

The outcome index adds four elapsed-time columns, in seconds: `data_loading_seconds`, `data_preparation_seconds`, `fine_mapping_seconds`, and `outcome_total_seconds`. Loading includes the gene's PGEN window and MAF or variant-list filtering. Preparation includes phenotype normalization, sample alignment, imputation, standardization, and covariate adjustment. Fine mapping includes the SuSiE fit and marginal z-score calculation. Total outcome time also includes result extraction. An error records partial time for the active stage; a stage that did not start has `NA`.

Each gene log reports its stage totals and total elapsed time across outcomes and scientific output writing. The gene total stops before final artifact hashing and provenance writing. Timing stays outside the saved fit objects. The batch log reports shared PLINK extraction, shared input loading, and input verification separately. These shared costs and Terra file localization are excluded from per-gene stage times. A preempted attempt can end before it writes a complete timing record.

The [IKZF1 cis benchmark](../benchmarks/synthetic_scale/reports/2026-11-06-real-ikzf1/README.md) compared complete results for 8,895 samples and 23,376 fitted variants. All saved scientific artifacts matched byte for byte. Mean fit and marginal-score time decreased from 98.0 to 53.4 seconds; the gene step, including input-verification overhead, decreased from 121.9 to 78.2 seconds.

## Localization and commands

All files that tasks open use WDL `File`, `File?`, or `Array[File]` inputs. Genotype files remain typed in the workflow's chromosome structs and in the batch task's separate inputs. The workflow has no file-writing expressions. The planner writes worklists and phenotype BEDs inside its task. It collects these files with explicit `glob()` outputs so Terra uploads them before the planner VM exits. Reading a list of filenames with `read_lines()` is not sufficient for this upload step. The command checks file counts, matching batch IDs, and file order against the planner lists. The workflow pairs both file arrays and chromosome indexes with `zip()`, which rejects unequal counts before batch localization. These files remain typed outputs and inputs for the fit tasks.

Use separate task `File` inputs when a command applies `sub()` to escape paths. In Cromwell 85, `sub(struct.member, ...)` can retain the cloud URI even after the file downloads. The batch task uses `sub(pgen, ...)`, `sub(pvar, ...)`, and `sub(psam, ...)` so the local paths reach the command. The unresolved-URI and readable-file checks still run before R or PLINK.

Each task uses named CLI arguments. Optional arguments are added only when their Files are present. The aggregate task writes newline-delimited lists of the already localized companion Files inside its command. Cloud URIs in a table are metadata and cannot cause automatic localization. Required inputs must be supplied as explicit WDL Files. The commands and R scripts check that local input files are readable before computation. They report an unresolved cloud URI as a localization error.

The command placeholders use WDL 1.0 `write_lines()` to create these lists after input localization. Their generated paths remain typed Files until command rendering, so Cromwell can map them into the container. No file-writing function runs at workflow scope. CI validates the complete PGEN descriptor and input metadata with pinned WOMtool 85, in addition to miniwdl and the command-localization tests.

## Outputs and retries

The workflow returns these combined tables:

- `SusieParquet`: selected variant results.
- `AggregateFullSuSiEParquet`: one combined `full_susie.parquet` with all rows already present in the per-gene full outputs, including variants outside retained credible sets, `X_column_scale_factors`, and `phenotype_residual_SD`.
- `CredibleSetsParquet`: credible-set results.
- `OutcomeIndex`: every expected outcome, including outcomes with no credible sets and fits that did not converge.
- `GenotypeCacheReports`: one JSON decision and runtime-counter report per batch.
- `CovariateExplainedVariants`: covariate-explained variant diagnostics.
- `ExpectedOutcomes` and `BatchIndex`: planner tables for checking batch ownership.
- `PreparedPgenManifest`: a headered manifest with the selected chromosome locations for later reuse. Set both QC thresholds to 0 when using already prepared files to avoid another QC pass.
- `PreparedPgens`, `PreparedPvars`, and `PreparedPsams`: one effective triplet per unique selected input, with original files used when shared QC is disabled.
- `BatchPgens`, `BatchPvars`, `BatchPsams`, and `BatchPreparationReports`: one subset triplet and preparation report per planned batch, in global batch order.
- `GenotypeQcCounts` and `GenotypeQcReports`: shared QC results, empty when shared QC is disabled. `GenotypeQcLogs` and `GenotypeQcErrors` contain preparation task stdout and stderr, including when only extraction runs.
- `ChromosomeGeneManifests`: one TSV per processed chromosome, named `chr<N>.gene_outputs.tsv`.

Each gene has a directory at `batch/genes/<gene-key>/`. Its key is `gene-` plus the SHA256 hash of its original ID. Separate output arrays contain each gene's RDS, three Parquet files, TSV files, log, and provenance JSON. The workflow also returns batch stdout, stderr, and PLINK logs.

Credible-set tables and both variant Parquets include the logical field
`feature_converged`. It records whether that row's molecular-trait fit converged.
All credible sets from a fit that did not converge have `feature_converged = false`.
These rows remain in the outputs. Filter on `feature_converged` when an analysis
requires converged fits. This flag does not describe credible-set purity.

Aggregation takes convergence from the validated outcome index, so it can add
the flag to older Parquets without another fit. A supplied flag must agree with
the index. Formatting an old saved fit with no convergence field records a
missing flag until an outcome index supplies that value. Empty outputs retain a
logical column. Shared checkpointed-workflow schemas are unchanged.

The `pip` column in primary and full variant tables is the native fitted
SuSiE variant PIP across effects. It is copied from the saved fit's `pip`
vector. The `alpha1` through `alpha10` columns in full tables contain the
separate effect probabilities. A variant's `cs_index` does not select or
change its PIP. Full tables also retain native PIPs for variants outside
retained credible sets.

**Output correction, 2026-11-09:** Earlier exporters put the assigned effect's
alpha value in `pip`. Full-table variants outside retained credible sets
received missing PIPs. Those exports need to be regenerated from the original
saved SuSiE RDS fits with the corrected formatter, then aggregated again.
Saved fits already contain native PIPs, so re-export does not require a new
model fit. Primary summary tables alone do not contain enough information to
recover the native PIPs. Existing downloaded files are not changed by this
code correction.

Each chromosome manifest has one row per gene. Its columns are `chromosome`, `gene_id`, `gene_key`, `batch_id`, `rds_uri`, `summary_parquet_uri`, `lbf_parquet_uri`, `full_parquet_uri`, `covariate_explained_variants_uri`, `variant_position_summary_uri`, `outcome_index_uri`, `log_uri`, and `provenance_uri`. Chromosome values use `1` through `22`. Original gene IDs are preserved. A gene with several outcomes has one row and one RDS containing those outcome fits. Genes with no credible sets or nonconverged outcomes remain in the manifests; their linked outcome index gives the per-outcome status.

The manifests use the completed gene output locations supplied by Cromwell. On Terra these are cloud URIs; local CLI tests can use local output paths. The manifest locations remain output-only string metadata. The full files are supplied separately as typed `Array[File]` inputs, so Cromwell localizes them before the aggregate command reads them. RDS, LBF, and log locations in the metadata lists are not opened. The aggregate task checks all nine manifest roles and writes the manifests alongside the combined tables.

For a direct call to `aggregate_pgen_susie.R`, supply all five optional metadata-list arguments together: `--gene-rds-uris-file`, `--gene-parquet-uris-file`, `--gene-tsv-uris-file`, `--gene-log-uris-file`, and `--gene-provenance-uris-file`. Each argument names a local newline-delimited list of completed output locations. These list contents are metadata, never file inputs to open. The WDL creates all five lists during command rendering. Direct CLI calls without them keep their previous outputs. Use an image containing the updated aggregation script when running the updated WDL.

The main workflow also supplies `--full-susie-files`, a task-local list of
localized full-result files. The aggregator writes one Parquet file through
bounded Arrow batches. It preserves the native columns and adds outcome
ownership from the validated index. It checks every expected gene file,
including empty files, and rejects missing, duplicate, foreign, or incompatible
inputs. The final file appears only after all inputs pass validation.

Aggregation does not recover rows absent from the source full files. The current
formatter can write an empty full file for an outcome with no retained credible
set. Such outcomes remain in `OutcomeIndex`. Rows outside credible sets in
nonempty full files are retained, including their missing CS assignments.

To combine completed outputs without running any model tasks, use
`workflows/AggregateFullSusie.wdl`. Its inputs are the completed combined
`outcome_index`, an explicit `Array[File] full_susie_files`, and the updated PGEN
image. Its single task returns `FullSuSiEParquet`. A chromosome gene-output
manifest provides the full-file URI values; supply those values as typed File
inputs. A manifest alone does not localize files.

For a local or task-level call using already localized files:

```bash
Rscript aggregate_pgen_susie.R \
  --full-only \
  --outcome-index completed/outcome_index.tsv \
  --full-susie-files localized_full_files.txt \
  --output-dir aggregate
```

The list needs one readable local path per line and one full file per expected
gene. URI strings inside that list are rejected. The full-only mode does not
open genotypes or invoke PLINK or SuSiE. Ordinary direct aggregation without
`--full-susie-files` keeps its previous outputs.

Full files can be much larger than summaries. Set aggregation disk space to
cover the localized full inputs, the combined Parquet, and temporary buffers.
Streaming limits R memory use but does not reduce localized disk requirements.
The new workflows have not been run on Terra; syntax, localization, and small
installed-code checks are separate from a complete Terra execution.

The full and credible-set summary Parquets record `phenotype_residual_SD` on
each variant row. This is the sample SD of the final residual phenotype vector
passed to SuSiE, after the pipeline's normalization, sample selection, and
covariate adjustment. The saved fit also contains this value. The individual-data
path uses this definition in every genotype cache mode. This metadata does not
change model settings, effects, PIPs, or credible sets.

For effects in phenotype residual SDs per additional ALT copy, divide both
`posterior_mean` and `posterior_sd` by `phenotype_residual_SD`. For effects per SD
of the covariate-adjusted genotype, multiply both by `X_column_scale_factors`.
Apply both operations for SD units on both variables. These conversions do not
require another model fit. A zero phenotype SD cannot be used as a divisor.
The fitted SuSiE error SD and the genotype diagnostic `residual_sd` are different
quantities.

Older saved fits and Parquets do not contain phenotype residual SD. Their rows
receive numeric missing values when formatted or aggregated. A mixed set of old
and new files retains the known scales and leaves the old scales missing. The
aggregator does not estimate a missing scale from the fitted error variance.

For `AdditionalGenotypes`, each gene uses the union of feature rows for its selected phenotype IDs. Every outcome for that gene receives the same union, as in the legacy per-gene task. Features from other genes are excluded. Source row order, sample order, and missing codes are preserved. A gene with no matching rows still uses the additional file's sample columns for the legacy sample intersection.

An explicitly associated additional variant keeps its supplied chromosome, including when it differs from the gene's chromosome. In the variant summary, `chromosome` describes the variant. In the outcome index, `chromosome` describes the gene. Aggregation checks each value against the corresponding identity.

Each gene's provenance records scientific files by argument role, with explicit presence and SHA256 values. This includes optional ancestry, variant-list, and additional-genotype files. Preparation records the source triplet once; each worker verifies only its own subset and companion inputs. The record includes MAF, effective variant-filter mode, cis distance, phenotype selection, missing-code policy, seed settings, package versions, and SHA256 values for the pipeline scripts and functions. Variant lists take precedence over ancestry MAF. The image reference comes from the typed WDL image input; local runs record `local` unless `SUSIER_RUNTIME_IMAGE` is set. Upstream source SHAs come from the image environment when available. An image tag records the requested reference; use a digest for a fixed image identity.

The aggregator reads explicit batch summary, credible-set, index, and diagnostic Files. It does not read gene RDS files. It rejects missing or duplicate outcomes, failed statuses, foreign outcomes, unreadable inputs, and incompatible table schemas. Fits that did not converge remain visible in the index.

Cromwell retries operate on complete batch tasks. There is no durable checkpoint for a completed gene within a preempted batch. A restarted batch repeats its genes. Use smaller batches when the cost of repeating a batch is too high.

## Checks

Run these checks with Python and miniwdl installed:

```sh
miniwdl check workflows/SusieR2.wdl
python3 tests/test_pgen_wdl.py .
```

The tests reject workflow-scope file writers. They also render and run every task command after simulated cloud-to-local File replacement. A recording CLI checks exact scalar arguments and file-list contents, including optional inputs, spaces, apostrophes, dollar expressions, and backticks. The WDL CI job also uses the command renderer in the pinned WOMtool 85 JAR to check the batch command with cloud File values and a local-path mapper. This catches the struct-member path error that a test with pre-localized inputs cannot detect. These checks use no Docker build or Terra submission. The separate image workflow builds and tests the installed R and PLINK code in GitHub Actions.

Normal image builds check CPU detection in a container limited to one CPU. They also run thread-setting tests, the standalone MAC check, the shared QC count tests, and one small installed-code smoke test: 48 samples, shared preparation, one cis-window extraction, one SuSiE fit, and aggregation. They check count and call-rate boundaries, duplicate IDs, retained phase/dosage, source hashes, extracted genotypes, inferred sample order, model outputs, timing records, thread settings, and reusable manifests. The full regression and legacy parity suites remain available. To run them in GitHub Actions, start the **SuSiE PGEN image CI** workflow manually and enable `full_checks`. To run them in an existing image, use `bash tests/test_pgen_container.sh --full`. Automatic image builds still require a change to source installed in the image.

Every image smoke run also checks exported primary and full-table PIPs against
a real multi-effect fitted model. This check includes non-CS variants,
cached extraction, and a Parquet write/read round trip.
