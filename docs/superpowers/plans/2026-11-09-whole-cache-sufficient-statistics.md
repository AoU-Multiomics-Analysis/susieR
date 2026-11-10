# Whole-cache sufficient statistics implementation plan

> **For agentic workers:** Use the existing isolated worktree and implement each part with a failing targeted test first. Root owns the fitter and integration. Separate agents own cache storage and fitter tests.

**Goal:** Automatically reuse one XtX for eligible whole-batch PGEN fits.

**Architecture:** The cache plans memory and builds a tiled cross-product once. Each prepared window receives an aligned subset. The fixed SuSiE model selects the SS fitter when supplied supported statistics and otherwise uses its current individual-data fitter.

**Tech stack:** R, pinned susieR 0.16.6, pgenlibr, existing quick GitHub image checks.

**Spec:** `docs/superpowers/specs/2026-11-09-whole-cache-sufficient-statistics-design.md`.

**Constraints:** Only whole caches; typed existing WDL inputs unchanged; no local Docker build; no full regression run; scientific outputs within 1e-8 with exact CS members; keep native purity sampling and guarded fallback.

**Review focus:** Physical variant indices with duplicate names; excluded and guarded columns; memory limits at the boundary; source version fallback; large-CS purity and RNG; saved-fit raw scales and fitted values.

1. Add failing cache planner/storage tests. Implement eligibility and lazy block-pair XtX construction, counters, memory summary, and cleanup in PgenCacheFunctions.R. Run only cache storage tests.
2. Add failing fitter tests for scientific agreement, invalid statistics, legacy/source fallback, and native CS use. Implement the guarded SS helper in SusiePerformance.R and optional statistics input in fitSusieModel. Restore individual-data fitted fields and compute marginal scores.
3. Pass aligned statistics from prepared PGEN windows after existing guards. Add eligibility logging at batch setup. Run the targeted cache and PIP output suites, with no full pipeline run.
4. Review the combined change independently. Fix findings, run the required focused checks, and report measured limits. Any GitHub image check must use the existing quick path with full_checks=false.
