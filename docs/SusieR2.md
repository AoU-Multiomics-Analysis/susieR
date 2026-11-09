# SuSiE 2.0 PGEN workflow

[SusieR2.wdl](../workflows/SusieR2.wdl) is a WDL 1.0 workflow for Terra's managed Cromwell service. It plans batches, extracts the cis windows with PLINK 2, fits each selected outcome, and combines the batch tables. The default batch size is 50 genes. Set `SusieR2.BatchSize` to another positive integer, such as 100, to change this limit. All outcomes for one gene stay in one batch.

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

Initial PGEN support is for diploid, biallelic variants on autosomes 1–22. The planner accepts chromosome aliases such as `1` and `chr1`. A gene cannot span chromosomes. The planner rejects missing chromosome inputs, duplicate chromosome IDs, ambiguous mappings, and duplicate phenotype IDs.

The phenotype BED must use tabs and have a header. Its first four columns are chromosome, start, end, and phenotype ID. The chromosome header can be `chr`, `chrom`, or `chromosome`, with an optional `#`. Coordinate headers are `start` and `end`; `chromStart` and `chromEnd` are also accepted. The reader skips blank lines and comment, `track`, or `browser` lines before the header.

Plain text, gzip, and BGZF files use the same parser, regardless of the filename suffix. Quoted IDs and surrounding field spaces are decoded consistently during planning, splitting, and PGEN batch loading. Original sample values and selected BED row text are preserved. Each data row must have the same number of tab-delimited fields as the header. Embedded tabs and line breaks within fields are not supported.

The planner checks the header, metadata, row count, and selected-row counts again during splitting. Errors report the failing line or metadata row. An ID mismatch alone does not establish that the file changed.

The planner validates the manifest without opening its genotype locations and preserves its row order. It writes a small normalized TSV. WDL reads that table and converts each PGEN, PVAR, and PSAM location into a typed File inside a `ChromosomePgen` struct. Each batch receives the three members as separate `File` inputs. Cromwell can then replace each cloud URI with its local path before the command escapes it for the shell. Each batch's zero-based chromosome index refers to the manifest row order. The planner's `batch_chromosome_indices.json` is output metadata; it is not a script argument wrapper.

When `SampleList` is omitted, sample selection uses IDs rather than phenotype values. Genotype-only and phenotype-only samples are excluded, while phenotype missing values retain their existing handling. Covariates must cover every retained sample; they do not silently reduce the cohort. A supplied list keeps its existing order and selection behavior. The PGEN reader and ancestry MAF calculation retain the complete genotype cohort before the phenotype subset. Additional-genotype sample restrictions still apply afterward. Provenance records the initial selection mode, PSAM/list order rule, sample count, and cohort hash before any additional-genotype restriction.

## Settings

| Input | Default |
|---|---:|
| `BatchSize` | 50 genes |
| `CisDistance` | 1000000 bases |
| `BaseSeed` | 1 |
| `MemoryDebug` | false |
| `MAF` | 0 |
| `MAC` | 0 (off) |
| `SelectTopPhenotypePerCluster` | false |
| `TopPhenotypePvalueColumn` | `qval` |
| `DockerImage` | `ghcr.io/aou-multiomics-analysis/susier/pgen:main` |
| `Cpu` | 1 |
| `BatchMemoryGB` / `BatchDiskGB` | 16 / 100 |
| `PlinkMemoryMB` | 1024 |
| `PlanMemoryGB` / `PlanDiskGB` | 8 / 50 |
| `AggregateMemoryGB` / `AggregateDiskGB` | 8 / 50 |
| `PreemptibleAttempts` / `MaxRetries` | 3 / 1 |

Select an image digest through `DockerImage` when a fixed image is required. Keep `PlinkMemoryMB` below the batch memory limit. Disk size must cover the localized chromosome files, extracted batch files, and gene outputs. CPU, memory, disk, and retry inputs are Terra runtime settings. PLINK extraction uses one thread.

Detailed memory reports are off by default. Set `SusieR2.MemoryDebug` to `true` to enable object-size reports and two extra garbage collections per fit. For a local PGEN batch, add `--memory-debug`. The legacy `susie.R` command accepts `--memory_debug true`. Direct R callers can use `options(susier.memory_debug = TRUE)`. The default keeps progress and error logs, variant diagnostics, and the three garbage collections after large matrices are removed. This setting does not change model calculations. A batch restores the caller's R option when it returns or fails.

The planner sorts genes by chromosome and position. It derives each feature position from `floor((BED start + BED end) / 2)`. Cis windows use inclusive coordinates, with the query start clamped to 1. PLINK uses `--extract range` for these 1-based regions. Disjoint windows retain their gaps. Extraction retains the complete genotype sample cohort. Ancestry MAF filtering therefore runs before the phenotype sample subset. The fit keeps the existing phenotype missing-value behavior. The initial image and batch command use one BLAS thread for reproducible arithmetic.

Set `SusieR2.MAC` to a positive integer to enable the minor allele count filter. For example, `"SusieR2.MAC": 5` keeps PGEN variants with a count of at least 5. The batch adds `--mac 5 --nonfounders` to PLINK during cis-window extraction. Counts use all samples in the input PGEN, including nonfounders, before phenotype sample selection. The filter uses PLINK's allele dosage sums; these are allele counts for hard-call genotypes. For the supported biallelic variants, PLINK's default nonmajor count is the minor allele count. See the [PLINK allele count filter](https://www.cog-genomics.org/plink/2.0/filter#mac).

`MAC` must be a nonnegative integer. With `MAC = 0`, the batch adds no PLINK MAC filter. The MAC filter runs before the existing ancestry MAF or variant-list filter, so a variant list cannot restore a variant removed by MAC. It applies to PGEN variants; variants supplied through `AdditionalGenotypes` keep their existing handling. Each gene's provenance records the threshold and sample scope. The direct batch CLI uses `--mac 5` and also defaults to `0`.

Seeds depend on the base seed, original gene ID, and original outcome ID. They do not depend on batch size, input order, or output paths. Changing batch size preserves the seed for each outcome. Each RDS contains a named list of that gene's outcome fits.

The fitting stage reuses identical matrix-vector products, expected residual-error values, and raw column standard deviations within one fit. It also uses `matrixStats::colVars()` for SuSiE's constant-column warning check. The model settings, optimizer, convergence checks, and iteration logs are unchanged. Each fit uses private copies of the audited SuSiE functions; the installed namespace is unchanged. Cached input references are cleared after success or failure. Unrecognized package sources, unsupported matrix types, and unsupported model modes use the original SuSiE path. Source guards match the pinned image source; local source checks do not constitute an image or Terra execution test.

Post-fit marginal scores reuse the same raw input means, standard deviations, and finite-value validation when the complete input matrix matches. Existing numerical fallback rules remain in place. These statistics stay local to the fit and are excluded from saved models. Each PGEN outcome's extracted tables are shared by its credible-set count and final output formatting; output schemas and ordering are unchanged.

The [IKZF1 cis benchmark](../benchmarks/synthetic_scale/reports/2026-11-06-real-ikzf1/README.md) compared complete results for 8,895 samples and 23,376 fitted variants. All saved scientific artifacts matched byte for byte. Mean fit and marginal-score time decreased from 98.0 to 53.4 seconds; the gene step, including input-verification overhead, decreased from 121.9 to 78.2 seconds.

## Localization and commands

All files that tasks open use WDL `File`, `File?`, or `Array[File]` inputs. Genotype files remain typed in the workflow's chromosome structs and in the batch task's separate inputs. The workflow has no file-writing expressions. The planner writes worklists and phenotype BEDs inside its task. It collects these files with explicit `glob()` outputs so Terra uploads them before the planner VM exits. Reading a list of filenames with `read_lines()` is not sufficient for this upload step. The command checks file counts, matching batch IDs, and file order against the planner lists. The workflow pairs both file arrays and chromosome indexes with `zip()`, which rejects unequal counts before batch localization. These files remain typed outputs and inputs for the fit tasks.

Use separate task `File` inputs when a command applies `sub()` to escape paths. In Cromwell 85, `sub(struct.member, ...)` can retain the cloud URI even after the file downloads. The batch task uses `sub(pgen, ...)`, `sub(pvar, ...)`, and `sub(psam, ...)` so the local paths reach the command. The unresolved-URI and readable-file checks still run before R or PLINK.

Each task uses named CLI arguments. Optional arguments are added only when their Files are present. The aggregate task writes newline-delimited lists of the already localized companion Files inside its command. Cloud URIs in a table are metadata and cannot cause automatic localization. Required inputs must be supplied as explicit WDL Files. The commands and R scripts check that local input files are readable before computation. They report an unresolved cloud URI as a localization error.

The command placeholders use WDL 1.0 `write_lines()` to create these lists after input localization. Their generated paths remain typed Files until command rendering, so Cromwell can map them into the container. No file-writing function runs at workflow scope. CI validates the complete PGEN descriptor and input metadata with pinned WOMtool 85, in addition to miniwdl and the command-localization tests.

## Outputs and retries

The workflow returns these combined tables:

- `SusieParquet`: selected variant results.
- `CredibleSetsParquet`: credible-set results.
- `OutcomeIndex`: every expected outcome, including outcomes with no credible sets and fits that did not converge.
- `CovariateExplainedVariants`: covariate-explained variant diagnostics.
- `ExpectedOutcomes` and `BatchIndex`: planner tables for checking batch ownership.
- `ChromosomeGeneManifests`: one TSV per processed chromosome, named `chr<N>.gene_outputs.tsv`.

Each gene has a directory at `batch/genes/<gene-key>/`. Its key is `gene-` plus the SHA256 hash of its original ID. Separate output arrays contain each gene's RDS, three Parquet files, TSV files, log, and provenance JSON. The workflow also returns batch stdout, stderr, and PLINK logs.

Each chromosome manifest has one row per gene. Its columns are `chromosome`, `gene_id`, `gene_key`, `batch_id`, `rds_uri`, `summary_parquet_uri`, `lbf_parquet_uri`, `full_parquet_uri`, `covariate_explained_variants_uri`, `variant_position_summary_uri`, `outcome_index_uri`, `log_uri`, and `provenance_uri`. Chromosome values use `1` through `22`. Original gene IDs are preserved. A gene with several outcomes has one row and one RDS containing those outcome fits. Genes with no credible sets or nonconverged outcomes remain in the manifests; their linked outcome index gives the per-outcome status.

The manifests use the completed gene output locations supplied by Cromwell. On Terra these are cloud URIs; local CLI tests can use local output paths. The locations are output-only string metadata: the aggregate task never opens or localizes those gene files. It reads the existing small batch files and its validated outcome index, checks all nine artifact roles, and writes the manifests in the same task. No extra VM, gene-file copy, output-directory guess, or new workflow input is needed. Files that the aggregate task opens remain typed WDL Files.

For a direct call to `aggregate_pgen_susie.R`, supply all five optional metadata-list arguments together: `--gene-rds-uris-file`, `--gene-parquet-uris-file`, `--gene-tsv-uris-file`, `--gene-log-uris-file`, and `--gene-provenance-uris-file`. Each argument names a local newline-delimited list of completed output locations. These list contents are metadata, never file inputs to open. The WDL creates all five lists during command rendering. Direct CLI calls without them keep their previous outputs. Use an image containing the updated aggregation script when running the updated WDL.

For `AdditionalGenotypes`, each gene uses the union of feature rows for its selected phenotype IDs. Every outcome for that gene receives the same union, as in the legacy per-gene task. Features from other genes are excluded. Source row order, sample order, and missing codes are preserved. A gene with no matching rows still uses the additional file's sample columns for the legacy sample intersection.

An explicitly associated additional variant keeps its supplied chromosome, including when it differs from the gene's chromosome. In the variant summary, `chromosome` describes the variant. In the outcome index, `chromosome` describes the gene. Aggregation checks each value against the corresponding identity.

Each gene's provenance records scientific files by argument role, with explicit presence and SHA256 values. This includes optional ancestry, variant-list, and additional-genotype files. Source chromosome files and extracted batch files are hashed once per batch. The record includes MAF, effective variant-filter mode, cis distance, phenotype selection, missing-code policy, seed settings, package versions, and SHA256 values for the pipeline scripts and functions. Variant lists take precedence over ancestry MAF. The image reference comes from the typed WDL image input; local runs record `local` unless `SUSIER_RUNTIME_IMAGE` is set. Upstream source SHAs come from the image environment when available. An image tag records the requested reference; use a digest for a fixed image identity.

The aggregator reads explicit batch summary, credible-set, index, and diagnostic Files. It does not read gene RDS files. It rejects missing or duplicate outcomes, failed statuses, foreign outcomes, unreadable inputs, and incompatible table schemas. Fits that did not converge remain visible in the index.

Cromwell retries operate on complete batch tasks. There is no durable checkpoint for a completed gene within a preempted batch. A restarted batch repeats its genes. Use smaller batches when the cost of repeating a batch is too high.

## Checks

Run these checks with Python and miniwdl installed:

```sh
miniwdl check workflows/SusieR2.wdl
python3 tests/test_pgen_wdl.py .
```

The tests reject workflow-scope file writers. They also render and run every task command after simulated cloud-to-local File replacement. A recording CLI checks exact scalar arguments and file-list contents, including optional inputs, spaces, apostrophes, dollar expressions, and backticks. The WDL CI job also uses the command renderer in the pinned WOMtool 85 JAR to check the batch command with cloud File values and a local-path mapper. This catches the struct-member path error that a test with pre-localized inputs cannot detect. These checks use no Docker build or Terra submission. The separate image workflow builds and tests the installed R and PLINK code in GitHub Actions.
