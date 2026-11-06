# SuSiE 2.0 PGEN workflow

[SusieR2.wdl](../workflows/SusieR2.wdl) is a WDL 1.0 workflow for Terra's managed Cromwell service. It plans batches, extracts the cis windows with PLINK 2, fits each selected outcome, and combines the batch tables. The default batch size is 50 genes. Set `SusieR2.BatchSize` to another positive integer, such as 100, to change this limit. All outcomes for one gene stay in one batch.

The complete workflow has **not been tested on Terra**. Local WDL checks and command tests do not establish Terra compatibility. No cloud jobs were submitted for these checks.

## Inputs

Copy [the input example](../examples/inputs/SusieR2.inputs.json) and replace each placeholder cloud path. The example includes chromosomes 1 and 2. Add an entry for each chromosome that contains a selected phenotype.

| Input | Type | Use |
|---|---|---|
| `ChromosomePgens` | Array of `ChromosomePgen` structs | Each entry has a chromosome string and three typed Files: `pgen`, `pvar`, and `psam`. Use one entry per chromosome. |
| `PhenotypeBed` | File | Phenotype BED with its original header, first four feature columns, and sample columns. |
| `QtlResults` | File | QTL results for the existing phenotype selection rules (`qval < 0.05`). |
| `Covariates` | File | Covariate table for the existing data loader. |
| `SampleList` | File | Sample list for the existing data loader. |
| `LegacyMissingCode` | String, required | Use `NA` or `-1`. Specify how the legacy genotype code represents missing values. The workflow does not infer this policy. |
| `GeneMap` | Optional File | TSV with `phenotype_id` and `gene_id`. Without this input, each selected phenotype is its own gene and outcome. |
| `Ancestry` | Optional File | Ancestry table for ancestry MAF filtering. |
| `VariantList` | Optional File | Variant selection list, with the existing loader's selection rules. |
| `AdditionalGenotypes` | Optional File | Additional genotype BED for the existing data loader. |

Initial PGEN support is for diploid, biallelic variants on autosomes 1–22. The planner accepts chromosome aliases such as `1` and `chr1`. A gene cannot span chromosomes. The planner rejects missing chromosome inputs, duplicate chromosome IDs, ambiguous mappings, and duplicate phenotype IDs.

`ChromosomePgens` is an ordered array because WDL 1.0 does not have a `keys(Map)` function. A workflow scatter extracts only the chromosome strings. The planner receives those strings and no PGEN Files. It writes each batch's zero-based index into the input array. The fit scatter uses that index to select the typed PGEN, PVAR, and PSAM Files. The planner's `batch_chromosome_indices.json` is output metadata; it is not a script argument wrapper.

## Settings

| Input | Default |
|---|---:|
| `BatchSize` | 50 genes |
| `CisDistance` | 1000000 bases |
| `BaseSeed` | 1 |
| `MAF` | 0 |
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

The planner sorts genes by chromosome and position. It derives each feature position from `floor((BED start + BED end) / 2)`. Cis windows use inclusive coordinates, with the query start clamped to 1. PLINK uses `--extract range` for these 1-based regions. Disjoint windows retain their gaps. Extraction retains the complete genotype sample cohort. Ancestry MAF filtering therefore runs before the phenotype sample subset. The fit keeps the existing phenotype missing-value behavior. The initial image and batch command use one BLAS thread for reproducible arithmetic.

Seeds depend on the base seed, original gene ID, and original outcome ID. They do not depend on batch size, input order, or output paths. Changing batch size preserves the seed for each outcome. Each RDS contains a named list of that gene's outcome fits.

## Localization and commands

All files that tasks open use WDL `File`, `File?`, or `Array[File]` inputs. Files remain typed inside chromosome structs until command rendering. The workflow has no file-writing expressions. The planner writes worklists and phenotype BEDs inside its task. These files become typed outputs for the fit tasks.

Each task uses named CLI arguments. Optional arguments are added only when their Files are present. The aggregate task writes newline-delimited lists of the already localized companion Files inside its command. Cloud URIs in a table are metadata and cannot cause automatic localization. Required inputs must be supplied as explicit WDL Files. The commands and R scripts check that local input files are readable before computation. They report an unresolved cloud URI as a localization error.

## Outputs and retries

The workflow returns these combined tables:

- `SusieParquet`: selected variant results.
- `CredibleSetsParquet`: credible-set results.
- `OutcomeIndex`: every expected outcome, including outcomes with no credible sets and fits that did not converge.
- `CovariateExplainedVariants`: covariate-explained variant diagnostics.
- `ExpectedOutcomes` and `BatchIndex`: planner tables for checking batch ownership.

Each gene has a directory at `batch/genes/<gene-key>/`. Its key is `gene-` plus the SHA256 hash of its original ID. Separate output arrays contain each gene's RDS, three Parquet files, TSV files, log, and provenance JSON. The workflow also returns batch stdout, stderr, and PLINK logs.

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

The tests reject workflow-scope file writers. They also render and run every task command after simulated cloud-to-local File replacement. A recording CLI checks exact scalar arguments and file-list contents, including optional inputs, spaces, apostrophes, dollar expressions, and backticks. These checks use no Docker build or Terra submission. The separate image workflow builds and tests the installed R and PLINK code in GitHub Actions.
