# Automatic PGEN Batch Memory Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans. Complete each task with a failing test, implementation, and verification.

**Goal:** Request enough memory for each Terra SuSiE batch from its metadata and report memory use at gene boundaries.

**Architecture:** Preparation estimates memory from the prepared PVAR, PSAM and exact feature windows. Typed WDL outputs select each worker's request. The worker collects unused R objects before reading the next gene and logs process and cgroup memory at stage boundaries.

**Tech Stack:** R, pgenlibr, tidyverse, WDL 1.0, GitHub Actions.

**Spec:** docs/superpowers/specs/2026-11-09-pgen-batch-memory-design.md

## Global Constraints

- Terra managed Cromwell; no workflow-scope file writers or file paths hidden in argument JSON.
- No local Docker builds and no Terra submissions.
- Positive BatchMemoryGB remains a fixed override; zero selects automatic sizing.
- Keep the existing cache contents and scientific calculations.
- Use decimal GB consistently for resource requests.

## Review Focus

- Unsorted worklists must keep numeric requests aligned with glob-sorted shards.
- Sample restrictions must not exceed the metadata estimate's cache allowance.
- Explicit cache caps and forced modes must fit the final task request.
- A killed process must have useful stage and memory information in prior log lines.
- Cleanup must preserve the prepared cache, tables, model results, and random state.

### Task 1: Metadata memory model

Files: R/utils/PgenMemoryFunctions.R, tests/test_pgen_memory.R.

- [x] Write and run failing tests for window size, cache modes, fixed overrides, explicit caps, duplicate physical rows, and invalid values.
- [x] Implement planPgenBatchMemory(reader, worklist, cache_mode, cache_memory_gb, memory_gb, cache_eligible). Use an 8 GB fitting floor, five dense matrix payloads for larger fits, two sample-square guard payloads, 2 GB overhead and 25% extra capacity. Account for the existing cache budget without circular sizing.
- [x] Run the memory tests.

### Task 2: Preparation and Terra resource wiring

Files: R/utils/PgenPreparationFunctions.R, R/scripts/prepare_pgen_susie.R, workflows/SusieR2.wdl, examples/inputs/SusieR2.inputs.json, preparation/WDL/WOM tests, container and CI file lists, docs/SusieR2.md.

- [x] Write and run failing tests for per-batch outputs, reversed worklists, input validation and WDL resource selection.
- [x] Calculate memory while the prepared shard reader is open. Add the breakdown to preparation provenance, write sorted integer JSON output data and a TSV report, and pass typed settings through named CLI arguments.
- [x] Set workflow memory from the matching preparation output. Pass that allocation to the worker's cache planner.
- [x] Install the new module and run its tests in the existing GitHub image checks. Update examples and documentation.
- [x] Run preparation, Python WDL/localization, real WOM, smoke, image-boundary and existing batch checks.

### Task 3: Between-gene cleanup and memory evidence

Files: R/utils/PgenRuntime.R, R/utils/PgenGeneFunctions.R, R/scripts/run_pgen_susie_batch.R, runtime and batch input tests.

- [x] Write and run failing tests for memory parsing and cleanup ordering.
- [x] Add safe Linux RSS, peak RSS and cgroup readings. Collect unused objects after each completed gene has returned its compact outputs and before the next genotype read.
- [x] Log memory at stage starts/completions and before/after cleanup. Preserve RNG and cached genotypes.
- [x] Run runtime, cache-fit parity and full applicable regressions.

### Final verification

- [x] Obtain an independent review and resolve material findings with regression tests.
- [x] Check whitespace and summarize the exact test results. State that the complete workflow has not run on Terra.

Validation record: memory (11), runtime (9), preparation (8), cache planning and fit agreement, WDL/localization (33), WOMtool 85, four-mode CLI smoke, input/provenance checks, R lint and repository checks passed. An old -1 sentinel-reference mismatch was corrected; the complete pipeline rerun then passed all layouts, both missing codes, and native/additional-feature Tabix comparisons. Installed-image CI exposed optimizer sensitivity to input rounding. The fixture now checks adjusted inputs at 1e-13, maximum absolute PIP error at 1e-8, and full fits/tables at 1e-7; variant order, credible sets, random state, convergence, iteration count, and safeguard fits retain exact checks. No Terra run or local Docker build was performed.
