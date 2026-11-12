# Manifest-only Gene Outputs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Stop the main SuSiE PGEN workflow from copying or combining full gene Parquet files after fitting; return their locations and RDS locations through manifests.

**Architecture:** Retain source artifacts as fitting-task Files. Pass their URIs to aggregation as Strings. Expose compact Parquets and chromosome manifests at workflow level, without direct gene-artifact File arrays.

**Tech Stack:** WDL 1.0, Cromwell 85, miniwdl, existing R aggregation CLI.

**Spec:** `docs/superpowers/specs/2026-11-12-manifest-only-gene-outputs-design.md`

## Global Constraints

Preserve all existing scientific and genotype behavior. No new workflow-scope file writers, local Docker builds, cloud jobs, or changes to the standalone full-only workflow. Retain source files with `delete_intermediate_output_files=false`.

## Review Focus

- Full Parquet files must not enter aggregate-task localization.
- Full Parquet and RDS sources must still be uploaded by fitting tasks.
- Final output File collection must exclude the six removed large/direct output fields.
- Manifest URI strings must retain original worker locations without opening gene files.
- Existing standalone full-only functionality must remain valid.

## Task 1: Change and verify the main workflow output boundary

**Files:** `workflows/SusieR2.wdl`, `tests/test_pgen_wdl.py`, `tests/test_full_aggregate_wom_render.scala`, `docs/SusieR2.md`.

**Interface:** Preserve `SusieParquet`, `CredibleSetsParquet`, and `ChromosomeGeneManifests`. Remove `AggregateFullSuSiEParquet`, `GeneRds`, `GeneParquets`, `GeneTsvs`, `GeneLogs`, and `GeneProvenance` from main workflow outputs. Preserve fitting-task gene outputs except the redundant `FullSusieParquetFiles` glob. Remove `full_susie_files` from `AggregatePgenSusie`; leave the optional standalone CLI unchanged.

- [x] Add and run a failing workflow File-collection regression before changing production WDL.
- [x] Remove full-file aggregation wiring and direct per-gene workflow outputs.
- [x] Update existing tests for the new main-workflow boundary while retaining standalone full-only tests.
- [x] Document manifest access, removed Terra output mappings, source uploads, and intermediate-file retention.
- [x] Run the full PGEN WDL test file, miniwdl checks, manifest CLI regression, Cromwell parser/render checks when available, and repository validation.
- [x] Review the focused diff and report the checks and lack of a Terra execution.

## Verification results

- The new workflow File-collection regression failed against the previous WDL for all five per-gene output groups and the combined full-Parquet output.
- All 34 PGEN WDL tests passed after the change.
- The manifest CLI regression passed with opaque, nonexistent gene URIs and byte-identical compact scientific outputs.
- Cromwell 85 validation, input extraction, and command-localization tests passed with Java 11. Both main and standalone aggregation were checked.
- Repository validation passed: all WDLs, YAML, example JSON, R parsing, and diff whitespace checks.
- Independent review found no issues in the production or test changes.
- No local Docker build or cloud job was run. The complete workflow has not been tested on Terra.
