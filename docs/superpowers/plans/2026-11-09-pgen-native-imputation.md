# Native PGEN Mean Imputation Implementation Plan

> **For agentic workers:** Use test-driven development and independent review.

**Goal:** Use pgenlibr mean imputation per feature without genotype caching.

**Architecture:** Use native mean imputation on the first PGEN read for ordinary
fits, after selecting the exact feature samples. Mark these genotypes as already
imputed. Active full-cohort ancestry MAF retains a raw read before its native
imputation function. Additional-genotype inputs keep their existing preparation.

**Tech Stack:** R, pgenlibr 0.7.0, existing SuSiE regression tests and image CI.

**Spec:** `docs/superpowers/specs/2026-11-09-pgen-native-imputation-design.md`

## Constraints and review focus

- No genotype cache, new runtime fallback, WDL change, or local Docker build.
- Native means use the exact fit samples and preserve sample and variant order.
- Full-cohort MAF precedes imputation when requested; normal fits decode directly with native imputation.
- Both missing policies and single-dimension matrices retain values and shape.
- Real hardcall matrices, complete fits, and RNG state agree with the old path.
- ReadList uses stored dosages; the intended input has integer hardcall dosages.
- Temporary native handles close after success and failure.
- Additional-genotype and other input paths retain their existing preparation.

## Task 1: Native imputation and fit integration

- [x] Add real PGEN tests; observe failure before implementation.
- [x] Save the absolute input path in `openPgenReader()` in `R/utils/PgenFunctions.R`.
- [x] Add `imputePgenHardcalls(reader, mat, missing_value = -1)` in `R/utils/PgenGeneFunctions.R`. Read missing variants on a temporary exact-sample reader in bounded blocks.
- [x] Add optional `genotype_imputer` and `genotypes_imputed` to `finemapPhenotype()`; skip repeated imputation for the first-read path.
- [x] Run native imputation, reader, gene parity, and existing statistical regressions.

## Task 2: CI, documentation, and review

- [x] Add the new regression to quick installed-image CI and verify local routing.
- [x] Document native per-feature imputation and its dosage input.
- [x] Run image configuration checks, preparation-to-fit smoke, and R lint.
- [x] Request independent review and resolve material findings.

## Execution record

- Reuse this feature checkout and preserve earlier upstream preparation edits and unrelated files. Do not reset or stage them.
- RED: the new real PGEN regression fails because `imputePgenHardcalls` does not yet exist.
- User clarification: verify equality in tests, then use the native result directly with no new fallback.
- User clarification: remove entirely missing variants before native imputation; keep existing empty-window handling and outcome statuses.
- User clarification: impute on the first read for ordinary PGEN-only fits instead of always rereading missing variants. Raw reads remain only for inputs that require full-cohort filtering or additional-genotype merging.
- [x] Add a first-read regression, observe unused-argument failure, then implement native `meanimpute` and exact `sample_ids` in `readPgenWindow()`.
- [x] Confirm the ordinary reader makes one native `ReadList(TRUE)` call and zero `ReadIntList()` calls. Remove only native-reported all-missing variant indices and retry; other errors propagate.
- PASS: reader suite (13 cases) and native imputation (16 cases), including reordered selected samples, integer stored dosages, first-read call counts, all-missing removal, duplicate physical records, multiple blocks, and handle cleanup.
- PASS: complete PGEN-versus-R fitted objects, decoded tables, RNG state, full-cohort MAF, and additional-genotype input tests.
- PASS: existing SuSiE regression suite; real pinned PLINK preparation, worker, and aggregation smoke; R lint; image configuration (8) and CI routing (9).
- Independent review found one relative-path reopening error. A new directory-change regression failed before the absolute-path fix, then passed. No material review finding remains.
- Independent review of the revised first-read path found no actionable issue. An independent shuffled-sample/window-offset fixture matched the R reference and removed the correct all-missing global variant indices.
- The real fit integration test confirms one native read, zero hardcall reads, and no repeated R imputation in an ordinary window. Entire fitted objects, decoded tables, and RNG state match the R reference.
- GitHub image CI is wired but has not been run for these local edits. No local Docker image was built, no Terra job was submitted, and no speed gain is claimed.
