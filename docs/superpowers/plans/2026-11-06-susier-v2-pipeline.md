# SuSiE 2.0 pipeline implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run PGEN SuSiE analyses in configurable batches of 50 sorted genes, with PLINK multi-region extraction, separate gene artifacts and complete summary aggregation.

**Architecture:** A metadata planner selects the same phenotypes as the current loader, groups complete genes, and streams small batch phenotype BEDs. Each batch task extracts its interval union into a local PGEN and runs the existing statistical path once per outcome. The aggregator reads explicitly supplied summary Files and checks complete outcome ownership.

**Tech Stack:** R/tidyverse, pgenlibr 0.7.0, susieR 0.16.6, PLINK 2, Arrow, WDL 1.0 and GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-05-susier-v2-pgen-design.md` (approved).

## Global constraints

- Default to 50 genes per batch; batch size is configurable. Sort within chromosome; keep all outcomes for a gene together.
- Preserve the existing statistical settings, compatibility fallback, and covariate-explained variant/credible-set diagnostics.
- Keep each gene's RDS and supporting files separate. RDS remains a named list of outcome fits.
- Use the complete genotype cohort for ancestry MAF before phenotype sample subsetting. Do not add complete-case phenotype filtering.
- Require an explicit `NA` or `-1` legacy missing-code policy. Do not infer mixed original missing codes.
- Use WDL 1.0 for Terra. Keep inputs typed File, File? or Array[File], including structs, until command rendering.
- No workflow-scope write_lines/write_tsv/write_map/write_json. Use named CLI arguments and task-local newline lists for file arrays; no JSON CLI argument wrapper.
- Check localization/readability before computation. No local Docker build and no Terra job submission.
- Use a pinned micromamba image and pinned conda-forge/bioconda dependencies. Image CI rebuild triggers cover contained source changes, not unrelated files.
- Work on `codex/susier-v2-pgen`; preserve unrelated files and PR #53. Agents do not commit or push; controller handles Git.

## Review focus

- Wrong chromosome/file association or off-by-one cis bounds must fail or retain the exact expected variants, including overlapping and disjoint intervals.
- Shuffling genes/outcomes or changing batch size must preserve each gene/outcome's seed and numerical fit.
- A batch may contain a gene with multiple outcomes, zero-CS fits or nonconverged fits; none may disappear from the status index.
- Aggregation must reject missing/duplicate expected outcomes, failed statuses, unreadable files or incompatible schemas; never treat missing summaries as empty success.
- WDL must pass localized paths, including optional files and aggregation arrays, without serializing cloud File inputs before command rendering.

---

### Task 1: Per-gene PGEN adapter and shared output formatting

**Files:** Create `R/utils/PgenGeneFunctions.R`, `R/utils/SusieOutputFunctions.R`, `tests/test_pgen_gene.R`, `tests/run_pgen_gene_tests.R`. Modify `R/scripts/susie.R`, `tests/test_covariate_explained_outputs.R` only as needed to share the current formatter.

**Interfaces:**
- Consume reviewed `openPgenReader()`, `readPgenWindow()` and `closePgenReader()` in `PgenFunctions.R`.
- `formatSusieOutputs(results, region_df)` returns `in_cs_variant_df`, `lbf_df`, `variant_df`, `cs_df`, `covariate_explained_variants` using current transforms and schemas. `results` is a named list, keyed by unique original outcome IDs; region_df has phenotype_id (these keys) and region. Preserve old-fit and empty behavior.
- `writeSusieOutputs(results, tables, output_prefix)` writes the existing named-list `_susie.rds`, three Parquets and `.covariate_explained_variants.tsv`; returns named paths. No process-level quit.
- `finemapPgenPhenotype(phenotype_id, se, reader, covariates, cis_distance, missing_code, ancestry_data = NULL, MAF = 0, variant_list = NULL, additional_genotypes = NULL, variant_report_file = NULL)` reads only that phenotype's cis window, applies filterMAF over the complete raw genotype cohort, and calls unchanged finemapPhenotype on the filtered double matrix. Clamp query start to 1 only; preserve requested coordinates in reports. Reject subset readers when ancestry MAF needs the complete cohort. Preserve the current NA-phenotype behavior.
- `geneArtifactKey(gene_id)` returns `gene-` plus the complete SHA256 of the UTF-8 ID. `geneOutcomeSeed(base_seed, gene_id, outcome_id)` hashes a versioned length-delimited tuple and converts the first seven hexadecimal digits to an R integer. Seed identity excludes batch, list order and output path.
- `runPgenGene(gene_id, outcomes, reader, output_root, base_seed, missing_code, provenance = list())`: outcomes is a named list keyed by outcome_id, each with phenotype_id, se, covariates, cis_distance, optional ancestry_data/MAF/variant_list/additional_genotypes and region. Validate unique IDs. Preserve/restore caller RNG kind/state around each fit. Write `output_root/genes/<gene_key>/<gene_key>` artifacts, plus `.variant_position_summary.tsv`, `.outcome_index.tsv`, `.provenance.json` and `.log`. Always write diagnostic headers. On failed fit write failure status/log then stop, without a successful manifest. Return only compact `summary`, `credible_sets`, `diagnostics`, `index` and named `files`, not full fitted objects.
- Index columns in this order: gene_id, gene_key, outcome_id, phenotype_id, seed (integer), status (`CONVERGED`, `NONCONVERGED`, `FAILED`), converged (logical), credible_set_count (integer), reason. Batch runner adds batch_id and chromosome. Successful no-CS fits remain in index/RDS. Provenance records artifact basenames/hashes and supplied structured input identities; do not hash an entire batch PGEN once per gene.

- [ ] Add and run failing tests for the absent adapter/formatter and identity functions.
- [ ] Implement shared formatting without changing numerical values, the checkpointed QR kernel or current CLI behavior.
- [ ] Implement adapter and gene artifact/RNG handling. Existing filtering and dense fallback stay in finemapPhenotype.
- [ ] Test real PGEN versus matching raw double-matrix reference and real Tabix fixture when available, both missing policies, MAF with extra PSAM samples, VariantList order/precedence, sample order, genotype/response fallback, multiple outcomes, reverse order, zero-CS schemas and failure status. Compare all native fit fields with fixed seeds and all decoded old output columns. Exercise the existing entrypoint's writer call too.
- [ ] Run `BENCH_LIBRARY=/private/tmp/susier-scale-benchmark.Src9T0/lib Rscript tests/run_pgen_gene_tests.R .`, existing regression tests and lint. Self-review and report exact RED/GREEN evidence.

### Task 2: Sorted batch planner, interval writer and strict aggregator

**Files:** Create `R/utils/PgenBatchFunctions.R`, `R/scripts/plan_pgen_susie.R`, `R/scripts/aggregate_pgen_susie.R`, `tests/test_pgen_batches.R`, `tests/run_pgen_batch_tests.R`.

**Interfaces:**
- `planPgenBatches(phenotype_bed, qtl_results, chromosomes, cis_distance, output_dir, batch_size = 50L, gene_map = NULL, select_top_phenotype_per_cluster = FALSE, top_phenotype_pvalue_column = "qval")`. Read first four BED columns only for planning; reuse current QTL selection semantics (qval < .05, optional strongest cluster). Default gene_id/outcome_id to phenotype_id; optional mapping TSV has phenotype_id, gene_id. Reject duplicate/ambiguous IDs, cross-chromosome genes and missing selected mappings. Derive feature position as floor((BED start + end)/2), as current loader does. Sort genes by numeric autosome and minimum selected position, then gene ID; group up to batch_size distinct genes within chromosome. Reject selected chromosomes absent from supplied chromosome names; initial reader supports autosomes 1–22.
- Batch worklist exact columns: batch_id, chromosome, gene_id, outcome_id, phenotype_id, phenotype_pos, requested_window_start, requested_window_end, window_start, window_end, region. Inclusive query bounds use max(1, requested start). region preserves current chromosome/requested bounds formatting. Keep every selected outcome for each gene.
- Write `plan/batches/batch_000001.tsv` etc and stream matching full phenotype rows into corresponding `.bed` files with original header, sample order and values. Avoid loading the genome-wide sample matrix. Emit `expected_outcomes.tsv` (all worklist rows), `batch_index.tsv`, and ordered newline files `batch_worklists.txt`, `batch_phenotype_beds.txt`, `batch_chromosomes.txt`; each path line is a task-output path relative to task cwd. The WDL will declare those generated worklists/BEDs as Array[File] outputs via read_lines, preserving their order.
- `writePgenBatchRegions(worklist, output_file)` writes PLINK bed1-compatible chromosome/start/end/unique-safe-label rows, no header. Include all actual outcome windows; retain gaps between disjoint windows. No genotype/MAF filtering here.
- `aggregatePgenBatches(expected_outcomes, summary_files, credible_set_files, index_files, diagnostic_files, output_dir)` reads only explicit local paths, rejects unresolved URIs, and validates each expected (batch_id,gene_id,outcome_id,phenotype_id) exactly once, matching chromosome and no FAILED status. Validate supplied batch file sets, schemas, ownership of all rows and diagnostic booleans. Keep NONCONVERGED status visible. Handle empty tables with valid schemas. Write `susie.parquet`, `credible_sets.parquet`, `outcome_index.tsv`, `covariate_explained_variants.tsv`. Do not read per-gene RDS/full/LBF files or follow paths hidden in a TSV.
- Planner CLI: --phenotype-bed, --qtl-results, --chromosomes-file (newline strings), --cis-distance, --output-dir, optional --batch-size, --gene-map, --select-top-phenotype-per-cluster, --top-phenotype-pvalue-column. Aggregator CLI: --expected-outcomes, --summary-files, --credible-set-files, --index-files, --diagnostic-files (newline lists of already localized files), --output-dir. Use optparse; preserve named arguments and clear log messages.

- [ ] Add and run failing planner/aggregator tests.
- [ ] Implement planner, streamed BED splitting and interval writer. Preserve original phenotype row data and complete windows.
- [ ] Implement strict aggregator and CLIs. Validate duplicate/missing outcomes even when all summary tables are empty.
- [ ] Test 101 genes -> 50/50/1, configurable size including 100, chromosome boundaries, ties, multi-outcome genes counting once, gene mapping, cluster selection, clipped/disjoint/overlapping windows, missing mappings, and exact streamed rows. Test aggregation success/no-CS plus missing/duplicate/failed/foreign-outcome/unreadable/schema errors. Include paths with spaces and cloud-URI failures.
- [ ] Run `BENCH_LIBRARY=/private/tmp/susier-scale-benchmark.Src9T0/lib Rscript tests/run_pgen_batch_tests.R .`, parse/lint and report RED/GREEN evidence. Do not change Task 1 files.

### Task 3: Batch CLI, PLINK execution, WDL and installed image smoke tests

**Files:** Create `R/scripts/run_pgen_susie_batch.R`, `workflows/SusieR2.wdl`, `containers/SusieR2/Dockerfile`, `containers/SusieR2/environment.yml`, `.github/workflows/susier-v2-image.yml`, `tests/test_pgen_pipeline.R`, `tests/test_pgen_wdl.py`, `tests/test_pgen_container.sh`, `docs/SusieR2.md`, `examples/inputs/SusieR2.inputs.json`. Update .dockstore.yml and static WDL CI registration.

**Interfaces:**
- Batch CLI: --pgen, --pvar, --psam, --worklist, --phenotype-bed, --qtl-results, --covariates, --sample-list, --cis-distance, --missing-code (required), --output-dir; optional --ancestry, --variant-list, --additional-genotypes, --maf (0), --base-seed (1), --plink2. Also preserve optional cluster-selection arguments. Check every input path before PLINK. Extract regions from worklist, invoke PLINK via safely quoted process arguments: --pgen <File> --pvar <File> --psam <File> --extract bed1 <regions> --make-pgen --out <batch prefix>. No sample filtering or allele changes. Load small batch molecular data once using current LoadData options, create SE and explicit work-item outcomes, open full-cohort batch reader once, call runPgenGene sequentially, discard large results between genes. Write compact batch outputs `batch.susie.parquet`, `batch.credible_sets.parquet`, `batch.outcome_index.tsv`, `batch.covariate_explained_variants.tsv`; index adds batch_id/chromosome. Log PLINK invocation/version, gene start/end, dimensions, fit status and batch completion.
- WDL `struct ChromosomePgen { String chromosome; File pgen; File pvar; File psam }`; workflow input `Array[ChromosomePgen]+ ChromosomePgens`, typed common phenotype/QTL/covariate/sample Files, required LegacyMissingCode, optional GeneMap/Ancestry/VariantList/AdditionalGenotypes, batch size 50, cis distance 1e6, seed 1, MAF 0, image String and configurable CPU/memory/disk/preemptible attributes. A workflow scatter lifts chromosome strings for the planner, no PGEN payload. Planner outputs ordered Array[File] worklists/BEDs, chromosome labels and zero-based input-array indices, plus File expected_outcomes. WDL 1.0 does not support keys(Map); the task-produced chromosome-index JSON is actual output metadata read as Array[Int], not an argument wrapper or a serialized input File. Scatter by index selects typed chromosome structs. Batch outputs include compact summaries/index/diagnostics and separate gene artifact File arrays. Aggregator gets explicit arrays of Files; create path lists only inside task commands. Workflow returns aggregate tables, complete index, separate flattened per-gene RDS/Parquet/TSV/log/provenance arrays. Preserve File types through commands; no workflow file writing.
- Pin dependency versions and eQTLUtils source SHA (26242562a4e244334fd9691d03bc1ef4d2d6c1d9). Verify package availability from primary sources; install pgenlibr/susieR source if needed. Do not build locally. CI builds only in GitHub Actions, runs real synthetic PLINK -> batch -> aggregate installed-code tests and legacy parity, and publishes only on main pushes. Trigger builds only when contained source changes (explicit paths), with manual dispatch allowed.

- [ ] Add static WDL/cloud-localization and real CLI smoke fixtures before implementation. Ensure tests fail for missing workflow/CLI.
- [ ] Implement CLI and WDL using preceding interfaces. Quote user strings safely; reject malformed lists and local/cloud path confusion. Use task-local safely rendered argument files only for scalar text or newline array lists, never serialize File inputs before localization.
- [ ] Run local real PLINK subset round trips and pipeline tests without Docker. Verify overlapping/disjoint intervals, missingness, REF/ALT/sample/variant order, single versus batched fits, gene artifacts, summary parity and masks. Test cloud-to-local command rendering for every new task, including optional paths/arrays; static scan must reject workflow-scope writers.
- [ ] Run miniwdl syntax checks, static regression tests, R lint and existing tests. Add installed-image Actions smoke tests and usage/input documentation. State clearly that Terra execution is untested.
- [ ] Review complete changes independently, resolve findings, commit and push the feature branch, and open a draft PR for review with concrete test evidence. Monitor relevant GitHub smoke checks; do not submit cloud analysis jobs.
