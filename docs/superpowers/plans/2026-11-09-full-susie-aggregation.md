# Full SuSiE aggregation implementation plan

> **For agentic workers:** Use the current isolated worktree. Apply test-first development for each implementation part. Separate agents own R aggregation and WDL wiring; root owns quick CI integration, documentation, and final verification.

**Goal:** Produce one validated full SuSiE Parquet from completed per-gene outputs, with no refit.

**Architecture:** Stream native full-result files through Arrow into one temporary Parquet, enrich ownership from the completed outcome index, then finalize. Main and aggregation-only WDLs provide explicit localized File arrays to a named CLI.

**Spec:** `docs/superpowers/specs/2026-11-09-full-susie-aggregation-design.md`.

**Constraints:** WDL 1.0, Terra File localization, task-command list creation, native values preserved, bounded memory, no local Docker build, no full run, no Terra submission.

**Review focus:** Empty/native schema differences; outside-CS NA values; physical variant identities; missing, duplicate, or foreign files; large output memory; incomplete writes; file-list shell quoting; metadata URIs never opened; no claim that absent no-CS source rows exist.

1. Add failing Arrow fixtures and CLI tests. Implement aggregatePgenFullSusie(outcome_index_file, full_susie_files, output_file, batch_size=16384L). Add optional --full-susie-files to ordinary aggregation and an explicit --full-only mode using --outcome-index.
2. Add failing WDL localization tests. Pass typed full-file arrays in SusieR2.wdl, expose AggregateFullSuSiEParquet, and add AggregateFullSusie.wdl with one aggregation task only.
3. Link the small R suite to quick installed-image checks. Verify that a failure stops the runner. Document input/output sizes, output coverage, and the remaining scale and statistical limits.
4. Run focused R, CLI, WDL, localization, lint, and build-boundary checks. Obtain independent review. Push the branch and run quick GitHub image checks with full_checks=false; preserve the no-full-run constraint.
