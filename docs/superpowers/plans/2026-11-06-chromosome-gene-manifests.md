# Chromosome Gene Manifest Implementation Plan

> **For agentic workers:** Use subagent-driven-development or executing-plans to implement each task with its regression checks.

**Goal:** Return one gene-output TSV per processed chromosome, with all nine artifact locations.

**Architecture:** Extend the existing aggregate task. Read its validated outcome index; treat completed output URIs as metadata strings that are never opened. Preserve all typed input Files and existing scientific outputs.

**Tech Stack:** R with dplyr/readr, WDL 1.0, miniwdl, WOMtool 85, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-11-06-chromosome-gene-manifests-design.md`.

## Global constraints

- No workflow-scope file writers, input File serialization before localization, guessed cloud paths, or JSON argument wrappers.
- URI strings are output-manifest metadata only; files opened by scripts remain typed File inputs.
- No extra VM task, genotype-file copies, local Docker builds, or Terra submissions.

## Review focus

- Multiple outcomes, no credible sets, and nonconvergence must not drop a gene.
- Shuffled or duplicate artifact arrays must not change gene ownership.
- A variant chromosome must not replace the gene chromosome.
- Cloud URI metadata must survive localization unchanged and must never be opened.
- All existing scientific tables and fits must remain unchanged.

## Tasks

1. [ ] Add failing manifest-builder tests, then implement `writeChromosomeGeneManifests()` in `R/utils/PgenBatchFunctions.R`. Add tests to the batch runner; validate all identity and nine-role invariants.
2. [ ] Add WDL type/command tests, then extend `AggregatePgenSusie` with five metadata URI arrays, command-local list arguments, and a chromosome-manifest File array output. Preserve parser and localization guards.
3. [ ] Extend `aggregate_pgen_susie.R` with all-or-none metadata-list CLI options. Invoke the builder after successful aggregation. Extend installed CLI tests with exact old scientific-table comparisons and manifest ownership checks. Update documentation.
4. [ ] Run the full batch suite, end-to-end PGEN pipeline, WDL parser/localization checks, image configuration, lint, and repository validation. Review the final diff, then publish a PR with the GitHub image smoke test.
