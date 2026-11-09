# Upstream PGEN preparation implementation plan

> **For agentic workers:** Use test-driven implementation with independent file ownership and a final review. Preserve the prior local empty-selection fix.

**Goal:** Apply native PLINK QC once per input and localize only batch-sized triplets in SuSiE workers.

**Architecture:** A preparation task per unique input performs QC and extracts its planned batch windows. A small report binds each shard to shared QC, its worklist and sample cohort. Workers consume these files directly.

**Tech Stack:** R, PLINK 2, WDL 1.0, miniwdl and Cromwell 85.

**Spec:** `docs/superpowers/specs/2026-11-09-upstream-pgen-preparation-design.md`

## Global constraints

- WDL 1.0, managed Terra; no workflow-scope file writers.
- Keep File inputs typed through command rendering; use named CLI arguments and newline-delimited lists for file arrays.
- Preserve samples, phase, dosage and retained metadata; QC includes nonfounders.
- PLINK native MAC/geno filters are authorized. Record their semantics accurately.
- No local Docker build or Terra job. Only contained source changes trigger image builds.
- Preserve unrelated workspace files and the prior empty-selection fix.

## Review focus

- Multiple chromosome aliases share one input, while input manifest order differs from batch order.
- Shard records and samples match the worklist and the QC cohort.
- Empty selections fail clearly and do not create false row IDs.
- Prepared workers never require or execute PLINK, and wrong shard/provenance inputs fail before fitting.
- Disabled QC still allows upstream extraction and records that filters were disabled.

## Tasks

- [x] Native QC: add a failing real-PLINK test, move filter decisions into PLINK flags, preserve report and metadata validation, verify the QC tests.
- [x] Preparation and runner: add real multi-batch fixtures; implement `preparePgenBatches(opt)` and `prepare_pgen_susie.R` with named inputs `pgen`, `pvar`, `psam`, `worklists_file`, `cis_distance`, `mac`, `min_call_rate`, `plink2`, `plink_memory_mb`, `output_dir`. Output `prepared/batches/batch_*.pgen`, `.pvar`, `.psam`, `.preparation.json`; shared QC, when enabled, lives in `prepared/qc`. Add `--prepared-batch-manifest` to the worker; validate and open the shard directly. Test tampering and no-PLINK prepared execution.
- [x] Planner and WDL: emit `PgenBatchIndices` (global indices grouped by source) and `BatchWithinPgenIndices`. Select typed worklists in a declaration-only scatter. Replace workflow FilterPgen calls with `PreparePgenBatches`; select per-batch outputs using source and within-source indices. Preserve whole-cohort output meaning and add separate batch outputs. Test expression order and localization.
- [x] Integration: include preparation source in the pinned image and exact CI path filters; extend installed-code smoke checks, static WDL checks and real Cromwell command rendering. Run required regression tests and independent review. Record local validation and the untested Terra state.

## Progress

- Design approved by the user's explicit upstream-extraction and native-PLINK requests. Existing branch is used; no new worktree is needed for this continuation.
- Existing QC empty-selection fix has real-PLINK RED/GREEN evidence and remains in the working tree.
- Native QC: mixed-dosage filter-authority test RED then GREEN; all 11 QC checks pass. PLINK counts remain descriptive, `FILTER_MAC` records native dosage MAC, and metadata records the native filter semantics.
- Preparation: four batches across two chromosomes share one QC call. Disjoint and overlapping windows, nonfounders, disabled QC and report tampering checks pass. Reports store an ordered sample digest rather than a repeated sample list.
- Workflow: 25 Python checks plus real Cromwell 85 parsing, input metadata and command rendering pass. Array File lists contain mapped local paths. Global batch order is preserved even when manifest source order differs.
- Integration: real fit smoke passes with a nonexistent PLINK executable in the worker; legacy MAC extraction, batch planning and full SuSiE regression checks pass. Image source filters and quick test routing include preparation. No image build or Terra run has occurred.
- Final review is complete. Its findings are fixed: shared QC requires a companion report; extracted metadata is checked and restored before hashes are recorded; nonempty preparation directories are rejected. The source metadata scan permits unselected sentinel positions while shard coordinates remain strict.
- Pinned PLINK 2.0.0-a.6.9 passed native QC, four-batch preparation and a real converged fit/aggregate smoke test. WDL checks include 25 Python cases, WOMtool 85 syntax/types, input metadata and actual command localization. The preparation regression includes metadata precision, duplicate IDs, aliases, gzip PSAM, and source row order. Image configuration, quick/full test routing and R lint pass.
- Local implementation is ready. Publishing and image execution remain separate integration steps; no Terra run is authorized or claimed.
