# Upstream PGEN preparation

## Intent and approved approach

The user requested upstream batch extraction so each SuSiE job localizes only its batch PGEN/PVAR/PSAM triplet. Apply variant QC once per complete input cohort. The user also selected PLINK's native `--mac`, `--geno` and `--nonfounders` filters; their input dosages are hard calls. Keep planning as a small metadata-only task.

## Design

- Run one preparation task per unique selected input triplet. Localize that original triplet once. Run shared QC there, then extract the union of cis windows for each assigned planned batch.
- Create one triplet per existing batch, with 50 genes per batch by default. Do not create a shard for each gene. Report batch and variant counts and output bytes so preparation cost can be assessed.
- Keep full input samples, including nonfounders, through QC and extraction. Phenotype sample selection remains in the fit task.
- Pass each fit task its own typed PGEN/PVAR/PSAM File inputs, worklist and preparation report. Prepared mode must not call PLINK or repeat filters.
- Preserve standalone runner extraction for existing callers. Add an explicit prepared-batch report option that selects the new mode.
- Bind the preparation report to the exact worklist, shard hashes, ordered samples and shared QC report. Preserve current whole-input report validation for standalone callers.
- Preserve reusable whole-cohort PreparedPgen outputs and their chromosome manifest. Expose batch shards as separate outputs.
- Planner emits integer metadata linking global batch order to each source and its within-source batch index. WDL preserves File types and uses those indices to select typed outputs.

## Approaches

The user selected preparation before fitting. Extraction inside each fit repeats full-file localization. A separate extraction task for each batch would also repeat full-file localization. Group QC and all extractions for one input in one preparation task instead. This serializes extraction within an input while retaining parallel preparation across different inputs.

## Constraints and validation

Use WDL 1.0 for managed Terra Cromwell. No workflow-scope file writers. Render safely quoted named arguments and task-local newline-delimited File lists. JSON preparation reports are structured provenance outputs, not argument wrappers. Keep the existing empty-selection QC fix.

Test real PLINK QC and extraction, multiple batches and shared input triplets, nonfounders, interval and sample order, malformed provenance, and prepared mode without PLINK. Test cloud-to-local File mapping with the real Cromwell renderer. Extend the existing GitHub Actions image smoke test and build source filters. Do not build Docker locally or submit Terra jobs. Report that the complete workflow has not been tested on Terra.
