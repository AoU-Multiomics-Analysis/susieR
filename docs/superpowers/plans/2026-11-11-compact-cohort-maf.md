# Compact cohort MAF implementation plan

> Use test-driven development and independent review for the implementation.

**Goal:** Add fit-sample `cohort_maf` only to compact SuSiE Parquets.

**Architecture:** Capture uncentered dosage MAF in predictor order and retain
it as fit/cache metadata. Append it to compact formatting and carry it through
the existing direct aggregation. The full output schema stays unchanged.

**Tech stack:** R, tidyverse, matrixStats, Arrow, pgenlibr; GitHub image CI.

**Spec:** `docs/superpowers/specs/2026-11-11-compact-cohort-maf-design.md`.

- [x] Add failing dosage, fit-path, cache, compact export, and aggregation cases.
- [x] Implement `cohortMafFromDosages(mat)` for the complete uncentered modeled
  dosage matrix; finite values within [0, 2] produce MAF in [0, 0.5]. Preserve
  all-missing legacy rows as NA using missingness recorded before imputation.
- [x] Store `cohort_maf` on ordinary/reusable fits and whole/rolling cache blocks;
  select by positions and pass through to prepared fits.
- [x] Format compact MAF from the fitted vector, with NA for older fits, empty
  outputs, and missing old compact inputs. Keep full/LBF/CS schemas unchanged.
- [ ] Verify native fit values and output row order with small fixtures. Cover
  fractional doses, ALT frequency above 0.5, missingness, sample subsets,
  duplicate labels, cache eviction, old data, and cached output formatting.
- [ ] Review and run focused local checks; push and verify both image builds
  and installed tests in GitHub. Merge and verify publication after checks pass.

Local dosage, cache-state, compact output, full aggregation, batch/manifest,
and lint checks pass. Independent review found no material issues. The complete
local R runner requires missing `eQTLUtils`, and native PGEN checks require
missing `pgenlibr`; the GitHub images provide these dependencies. One obsolete
QC metadata test failed against the base commit as well. Its expectations now
follow the direct-merge contract from PR #80 and verify unchanged table hashes.
