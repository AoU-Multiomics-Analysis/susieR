# PGEN Input Manifest and Optional Samples Implementation Plan

> **For agentic workers:** Use subagent-driven-development or executing-plans, with failing regression tests before implementation.

**Goal:** Accept one PGEN TSV manifest and infer the genotype/phenotype sample intersection when SampleList is absent.

**Architecture:** The existing planner validates and normalizes the manifest; WDL promotes its path fields to typed Files. The batch opens the full PGEN cohort and gives its IDs to the loader only for automatic sample inference.

**Tech Stack:** WDL 1.0, R/dplyr/readr, miniwdl, WOMtool 85, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-06-pgen-input-manifest-design.md`.

## Constraints and review focus

- No workflow file writers, hidden unlocalized input paths, JSON CLI wrappers, extra VMs, or genotype copies in planning.
- Preserve manifest order, exact path text, shared PSAM paths, typed genotype File fields, and safe quoting.
- Keep a supplied sample list's behavior; infer PSAM-order ID intersection only when it is absent.
- Preserve full-cohort MAF, phenotype missing-value handling, covariate coverage, and separate gene outputs.
- Run only local tests and GitHub Actions; do not submit Terra jobs or build Docker locally.

## Tasks

1. [ ] Add manifest-parser/CLI failing tests; implement validation and normalized output in PgenBatchFunctions.R and plan_pgen_susie.R. Keep the old CLI route.
2. [ ] Add inferred-sample failing tests; implement guarded loader inference, optional raw covariate-header validation, earlier full reader opening, optional SampleList, and provenance in ImportFunctions.R and run_pgen_susie_batch.R.
3. [ ] Add WDL type/localization tests; replace the array input with PgenManifest, construct typed structs from normalized planner output, make SampleList optional, and update WOMtool metadata expectations.
4. [ ] Exercise both changes through the real installed CLI/PGEN pipeline. Verify exact supplied-list fits and inferred-versus-explicit intersection fits. Update examples and documentation.
5. [ ] Run regression, batch, WDL parser/localization, image config, lint, and repository checks; review and update the GitHub pull request.
