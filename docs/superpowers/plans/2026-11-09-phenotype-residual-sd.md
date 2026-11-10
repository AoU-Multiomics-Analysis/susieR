# Phenotype residual SD implementation plan

> **For agentic workers:** Use the shared existing worktree. Implement the independent export and aggregation tasks in parallel, then review the combined change.

**Goal:** Export the SD needed to express an existing per-copy effect in phenotype residual SD units.

**Architecture:** Record the sample SD at the common fit boundary, preserve it in native saved fits and variant tables, and carry it through the compact and streaming full aggregators. Missing legacy metadata remains missing.

**Tech Stack:** R, tidyverse, Arrow, existing PGEN tests and GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-09-phenotype-residual-sd-design.md`

## Global constraints

The model, output row scope, WDL interfaces, and existing scientific values must remain unchanged. Use the exact numeric name `phenotype_residual_SD`. Do not infer missing values from error variance. No local Docker build, full scientific run, or Terra submission.

## Review focus

- Verify the value uses the exact final phenotype vector, not raw phenotype or fitted noise SD.
- Check all genotype cache modes retain the same phenotype SD definition.
- Check empty files, old fits, old Parquets, and mixed old/new aggregation.
- Ensure one scale per outcome even when full output spans Arrow batches.
- Ensure new metadata does not change scientific comparisons or output PIPs.

## Task 1: Fit metadata and exports

- [ ] Add a regression check for the sample SD of a known fitted phenotype vector, and for full/summary serialization, cached extraction, and missing old-fit metadata. Run the check before implementation and confirm the missing-field failure.
- [ ] Save `stats::sd(drop(y))` in the common fit result and add the value to full and summary extraction/empty schemas.
- [ ] Run the focused output tests and relevant small fit tests. Check scientific fields remain unchanged.

## Task 2: Aggregation

- [ ] Add regression cases for retaining new values, filling old values with numeric NA, mixed files, and rejecting invalid or inconsistent scales. Confirm failure before implementation.
- [ ] Extend compact and full streaming aggregation with only this named optional schema addition. Keep all other schema and ownership checks.
- [ ] Run the full-aggregation fixtures and existing planning/aggregation checks.

## Task 3: Delivery

- [ ] Document the field definition, conversion, and legacy limits.
- [ ] Review the combined diff and run targeted checks.
- [ ] Commit and push the reviewed change. Run GitHub Actions image CI with `full_checks=false`; check the actual logs before reporting success.
