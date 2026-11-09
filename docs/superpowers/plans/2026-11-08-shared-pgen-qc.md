# Shared PGEN genotype QC implementation plan

**Goal:** Calculate MAC and AN once per selected input PGEN triplet and reuse the filtered files in every SuSiE batch.

**Spec:** The user requested MAC and AN filtering, selected a minimum called fraction, and set its default to 95%. MAC remains 0 by default. For diploid autosomal hard calls, AN is twice the number of called samples, MAC is the smaller allele count, and call rate is AN / (2 × input sample count). Include nonfounders. Calculate counts from PGEN; do not trust PVAR INFO annotations.

**Architecture:** The planner groups selected PGEN/PVAR/PSAM triplets and writes a QC index for each batch. One conditional QC task per unique triplet calculates hard-call counts, filters native PGEN records, and returns explicit File outputs. Batches use those files with no repeated MAC calculation. The QC report is structured result data, not an argument configuration. Aggregation writes a reusable output manifest from output-location metadata.

**Constraints:** WDL 1.0 for managed Terra Cromwell; no workflow-scope file writers; File inputs stay typed until command rendering; flat task File inputs and named CLI arguments; no local Docker build or Terra submission. Only contained R source changes trigger automatic image builds.

**Defaults:** `MAC = 0`, `MinCallRate = 0.95`. Both must be valid nonnegative thresholds; call rate must be at most 1. Setting both to 0 skips shared QC and uses the original files. Counts use the complete input PGEN cohort before phenotype selection. Existing standalone batch `--mac` behavior remains available; workflow batches pass 0 and attach shared QC provenance.

## Review focus

- A call rate exactly equal to 0.95 passes; missing calls are not imputed before counting.
- Minor reference alleles, missing/stale INFO counts, nonfounders, and dosage-only missing hard calls use the correct hard-call counts.
- Duplicate and missing raw IDs cannot select the wrong record; native phase and dosage data and retained PVAR records are preserved.
- Repeated triplets across chromosomes run QC once; different companion files remain distinct. Batches receive the matching prepared files and report.
- Invalid thresholds, all supported variants removed, or a report that does not match the prepared files fail clearly. Unsupported records outside analysed windows retain existing handling.

## Tasks

1. Add `PgenQcFunctions.R`, a filter CLI, and real PLINK tests. First demonstrate failing tests. Use row-unique temporary IDs for count and extraction commands; verify retained order before restoring original metadata. Write compressed count results and a report with thresholds, counts, and source/output hashes.
2. Add planner metadata for unique selected triplets and batch indices. Validate thresholds before planning. Add an optional batch QC report, verify its output hashes before fitting, and record effective shared filters without repeating them.
3. Wire the conditional QC task into WDL, with readable localized inputs and explicit outputs. Export effective input files, QC results, and a reusable headered manifest. Test shared-file routing and disabled QC with the real Cromwell 85 renderer as well as static workflow writer checks.
4. Include the new R source in the pinned image and its exact CI source filters. Run QC and one-fit smoke checks in GitHub Actions. Check defaults, CLI validation, hashes, count boundaries, metadata, and current statistical regression tests. Obtain an independent review, open and attach a PR, and merge after checks pass under the user's existing instruction.

The complete workflow must not be described as Terra-validated unless it runs on Terra. No cloud analysis job is authorized by this implementation request.
