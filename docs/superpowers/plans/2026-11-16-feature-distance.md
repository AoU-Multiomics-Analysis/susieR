# Feature Distance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the signed distance to each fine-mapped variant's feature.

**Architecture:** Pass authoritative feature chromosome and midpoint with region
metadata. Compute distance during compact formatting. Preserve it in batch and
combined outputs, with numeric missing values for older inputs.

**Tech Stack:** R, tidyverse, arrow, existing GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-16-feature-distance-design.md`

## Global Constraints

- Column: `distance_to_feature_bp`, numeric, variant position minus feature midpoint.
- Use the existing phenotype midpoint and chromosome; do not infer it from region.
- Preserve fitting, physical rows, other scientific values and full/LBF/CS schemas.
- No WDL change, local Docker build, or Terra submission.

## Review Focus

- Asymmetric regions: use phenotype position, not the region center.
- Repeated variant labels: preserve physical rows and cached output.
- Additional variants on other chromosomes: numeric missing distance.
- Old metadata and compact files: missing values, without invented coordinates.
- Empty output: numeric column and unchanged checkpointed defaults.

### Task 1: Format distance and pass feature coordinates

**Files:** `R/utils/SusieOutputFunctions.R`, `R/utils/InitFunctions.R`,
`R/utils/ImportFunctions.R`, `R/utils/PgenGeneFunctions.R`, output and PGEN gene tests.

**Interfaces:** Optional `feature_chromosome` and `feature_position` in `region_df`.
Compact output appends numeric `distance_to_feature_bp` before convergence.

- [x] Add signed/zero/cross-chromosome/missing/empty/Parquet regression cases; run RED.
- [x] Pass coordinates from phenotype metadata; compute the compact column; run GREEN.
- [x] Update exact schema tests and add native PGEN coordinate assertions.

### Task 2: Combine and document distance

**Files:** `R/utils/PgenBatchFunctions.R`, `tests/test_pgen_batches.R`, `docs/SusieR2.md`.

**Interfaces:** Recorded numeric distances pass through unchanged. Absent columns
receive numeric NA, including empty all-old inputs.

- [x] Add mixed/all-old/empty aggregation cases; run RED then implement and run GREEN.
- [x] Document sign, feature reference point, SV reference point and missing values.

### Task 3: Verify and update PR

- [x] Run local output, batch, lint and repository validation checks.
- [x] Get independent review and address material findings.
- [ ] Commit and push the feature to PR #83; update its title and description.
- [ ] Check GitHub Actions, including installed PGEN tests, and report results.

## Verification record

The output regression failed before the formatter change because the distance
column was absent. Output, batch/manifest/CLI, repository validation, and R lint
checks passed after implementation. Loader midpoint checks passed locally. Native
PGEN gene checks require the installed CI image; no local image was built.

Independent review found no material issues and separately passed output and
batch suites plus the diff check.
