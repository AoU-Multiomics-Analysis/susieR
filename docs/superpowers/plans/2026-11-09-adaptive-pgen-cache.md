# Adaptive PGEN Cache Implementation Plan

> **For agentic workers:** Use test-driven development, parallel ownership, and independent review.

**Goal:** Automatically select prepared-genotype whole, rolling, or no caching per batch.

**Architecture:** A separate cache module plans from metadata, stores prepared
columns in bounded blocks, and returns exact feature windows. The PGEN adapter
fits these matrices or invokes the existing numerical guard path. The batch
runner owns cache lifetime and reporting. WDL passes typed settings and returns
the report without changing file localization.

**Tech stack:** R, pgenlibr 0.7.0, SuSiE, WDL 1.0, existing GitHub image CI.

**Spec:** `docs/superpowers/specs/2026-11-09-adaptive-pgen-cache-design.md`

## Constraints and review focus

- Reuse the feature checkout; preserve previous preparation/imputation edits and unrelated files.
- Whole/rolling/none use the same fixed cohort and covariates, exact PGEN indices, feature windows, and seeds.
- No new imputation fallback or scientific model settings. Preserve current dense numerical safeguards.
- Unsupported preparation options keep the no-cache path in auto mode.
- Budget estimates cover cached columns and bounded preparation buffers; do not claim a peak-RAM guarantee.
- No genotype cache object may enter the saved fit or scientific output tables.
- WDL 1.0, typed file localization, named CLI args, and task-only file writing.
- No local Docker build or cloud submission; GitHub image CI owns image tests.

## Task 1: Core planning, indexed reading, and storage

- [x] Add failing real PGEN tests in `tests/test_pgen_cache.R`.
- [x] Implement `planPgenGenotypeCache`, `createPgenGenotypeCache`, `getPgenGenotypeCacheWindow`, `clearPgenGenotypeCache`, and `pgenGenotypeCacheSummary` in `R/utils/PgenCacheFunctions.R`.
- [x] Expose `readPgenVariants` in `R/utils/PgenFunctions.R`, returning matrix and retained physical indices; retain the existing window reader contract.
- [x] Verify policy thresholds, budgets, overlapping/disjoint/backtracking windows, exclusions, physical duplicates, counts, and cleanup.

## Task 2: Fits and batch ownership

- [x] Write complete-fit cache parity tests and observe failure.
- [x] Add optional cache inputs to `finemapPgenPhenotype` and `runPgenGene`.
- [x] Use prepared columns with per-feature phenotype residualization; preserve genotype/response guard recomputation.
- [x] Add CLI settings and per-batch plan/cache lifetime/reporting to `run_pgen_susie_batch.R`.
- [x] Verify ordinary scientific parity, exact guarded fits, RNG, output tables, and no-cache behavior.

## Task 3: WDL, image, documentation, and review

- [x] Add typed WDL settings, task named arguments, and report output, with regression tests first.
- [x] Install the cache module in the image and restrict automatic image triggers to installed source changes.
- [x] Wire cache checks into the existing GitHub image smoke path.
- [x] Run relevant R, CLI, WDL/Cromwell, image configuration, smoke, and lint checks.
- [x] Request independent review and resolve material findings.

## Execution record

- User explicitly requested implementation of the discussed automatic policy and prepared-genotype cache.
- Core reader/cache and WDL/image integration are delegated with separate file ownership. Root owns fit and batch integration.

- Verified the PR source against current main (`007b7957`) and preserved CPU/thread controls and batch, gene, and outcome timing logs.
- The integrated cache suite passes eight core checks and seven complete-fit checks, including balanced cached timing events and guarded timing restarts.
- Independent review is clean. The interleaved-worklist execution-order issue and relative-path reopening issue have regression tests.
- WDL checks (28), image configuration (8), CI routing (10), and real Cromwell 85 syntax/type/localization checks pass. The four-mode real PLINK preparation-to-fit smoke and input/provenance checks pass.
- CPU runtime tests, the complete statistical regression suite, and the complete PGEN gene regression suite pass on the final PR source.
