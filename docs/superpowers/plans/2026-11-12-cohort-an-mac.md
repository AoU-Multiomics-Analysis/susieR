# Cohort AN and MAC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add observed fit-sample AN and MAC to compact variant Parquets.

**Architecture:** Compute count vectors before imputation, retain physical row
alignment through fits and caches, and extend the existing compact MAF export.
Native PGEN reads use an opt-in raw-dosage sidecar without changing model matrices.

**Tech Stack:** R, tidyverse, matrixStats, pgenlibr, arrow; existing GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-12-cohort-an-mac-design.md`

## Global Constraints

- Names: `cohort_an`, `cohort_mac`; numeric vectors and Parquet columns.
- Final fit samples; observed diploid dosages before imputation or adjustment.
- Preserve model matrices, PIPs, moments, CSs, RNG and existing MAF values.
- Preserve full/LBF/CS/public extraction schemas and checkpointed defaults.
- No WDL change, local Docker build, dependency installation or Terra submission.

## Review Focus

- Missing calls and high ALT frequency: observed AN and the minor REF count.
- Fractional dosage: no hardcall rounding and no imputed-count inflation.
- Repeated labels: metadata follows physical rows, including cache eviction.
- Old and empty outputs: independently missing fields remain numeric.
- Native-first final sample subset: counts use the same subset as the fit.

### Task 1: Capture count metadata

**Files:** `R/utils/SusieFunctions.R`, `R/utils/PgenFunctions.R`,
`R/utils/PgenGeneFunctions.R`, `R/utils/PgenCacheFunctions.R`, dosage/native/cache
and fit-path tests.

**Interfaces:** `cohortAlleleCountsFromDosages(mat)` returns a list with numeric
`cohort_an` and `cohort_mac` vectors. `readPgenVariants(...,
include_cohort_counts=FALSE)` preserves the existing default return; opt-in returns
the two vectors beside retained physical `variant_indices`. Fits and cache windows
carry the same two named metadata fields.

- [ ] Add hand-checked missing/ALT-major/fractional/all-missing/subset tests; run RED.
- [ ] Add native-reader/cache/fit-path propagation tests; run available RED cases.
- [ ] Implement raw count capture and positional propagation; run GREEN.
- [ ] Run existing dosage, cache and statistical/output checks; use CI for missing native dependencies.

### Task 2: Export and aggregate counts

**Files:** `R/utils/SusieOutputFunctions.R`, `R/utils/InitFunctions.R`,
`R/utils/PgenBatchFunctions.R`, compact-output and aggregation tests, relevant
schema expectations, `docs/SusieR2.md`.

**Interfaces:** Fit metadata from Task 1 becomes columns after `cohort_maf` and
before `feature_converged`. Each missing old field gets numeric NA independently.

- [ ] Add tests for physical row alignment, cached extraction, old fits, zero/empty
  counts, Parquet round trips, mixed batches and preserved scientific outputs; run RED.
- [ ] Extend compact export, empty schema and aggregate fallback; run GREEN.
- [ ] Document sample scope, missingness, fractional MAC and output scope.

### Task 3: Review and open PR

- [ ] Run repository validation, R lint, available regression suites and diff checks.
- [ ] Request independent review; address material findings with regression tests.
- [ ] Commit and push only this feature; open and attach the PR.
- [ ] Check GitHub Actions image/installed-code smoke tests and report their state.
