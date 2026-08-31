# Checkpointed Window Univariate SuSiE Design

## Goal

Add a new univariate SuSiE workflow to the existing `susieR` repository. Run one prepared variant window per WDL task. Fine-map each phenotype already linked to that window. Save each complete SuSiE fit to Google Cloud Storage (GCS) before the task starts the next phenotype. Resume an interrupted preemptible attempt without reading every earlier fit.

The new workflow must not change the behavior of `R/scripts/susie.R` or any existing WDL workflow.

## Scope

The first version has the following scope:

- Use all usable variants in the supplied prepared window.
- Do not create a phenotype-centered cis interval.
- Trust the prepared phenotype manifest as the authoritative membership list.
- Process only the linked phenotypes present in that manifest.
- Run one univariate SuSiE model for each phenotype.
- Process phenotypes in ascending association P-value order.
- Save the complete fitted SuSiE object for every completed model.
- Use one WDL workflow invocation and one preemptible task for each window.
- Produce per-window result tables and a fit index.

The first version does not select cis phenotypes, change the prepare-window workflow, create a genome-wide cross-window aggregate, or change existing `susieR` entry points. A later prepare-window change will add the required `p_value` column before production use.

## Prepared Input Contract

The analysis workflow consumes outputs from the separate prepare-window workflow. The prepared phenotype manifest must contain these columns:

| Column | Meaning |
|---|---|
| `window_id` | Prepared variant-window identifier. |
| `phenotype_id` | Unique molecular-trait identifier within the window. |
| `modality` | Covariate modality for the phenotype. |
| `phenotype_file` | Prepared phenotype-data file name. |
| `p_value` | Association P value used for processing priority. |

The wrapper must fail before model fitting when:

- A required column is absent.
- A required value is missing or empty.
- A P value is not finite or is outside zero through one.
- A phenotype ID occurs more than once.
- The manifest contains more than one window ID.
- The manifest window ID does not equal the WDL `window_id` input.
- A manifest phenotype is absent from the prepared phenotype data.
- A manifest modality has no applicable covariate input.

Sort the authoritative manifest by `p_value`, `modality`, and `phenotype_id`. The resulting order is part of the analysis fingerprint. The wrapper must not read a global association table or add, remove, or reclassify phenotypes.

Other inputs are:

- The prepared wide dosage file for one window.
- The prepared phenotype data file.
- One or more covariate files and matching modality labels.
- An optional sample allowlist.
- The GCS checkpoint root.
- WDL runtime settings, including memory, CPU, disk, and preemptible attempts.
- The existing univariate SuSiE settings: `L = 10`, residual-variance estimation enabled, prior-variance estimation enabled, `scaled_prior_variance = 0.1`, univariate Z-score calculation enabled, and `min_abs_corr = 0.5`. Other SuSiE behavior comes from the package version in the pinned image.

## Repository Structure

Add an isolated workflow path:

```text
R/scripts/run_checkpointed_window_susie.R
R/utils/CheckpointedWindowSusieFunctions.R
workflows/CheckpointedWindowSusie.wdl
containers/CheckpointedWindowSusie/Dockerfile
.github/workflows/checkpointed-window-susie-image.yml
examples/inputs/CheckpointedWindowSusie.inputs.json
tests/test_checkpointed_window_susie.R
tests/test_checkpointed_window_susie_wdl.sh
```

The wrapper sources existing import, memory, and SuSiE helper files. The wrapper uses a preloaded prepared-window genotype matrix when it calls the existing phenotype-level fine-mapping functions. This prevents the existing cis-oriented dosage extractor from selecting variants around a linked trans phenotype.

The new helper file owns only the new input validation, deterministic ordering, checkpoint storage, fit validation, manifest handling, and window-output assembly. Existing entry points remain unchanged.

## Model Data Flow

The wrapper performs these steps:

1. Validate all inputs and the prepared manifest.
2. Sort the manifest by P value, modality, and phenotype ID.
3. Calculate content hashes for scientific inputs and covariates.
4. Build the analysis fingerprint from the content hashes, ordered manifest, model settings, wrapper version, and pinned container identity.
5. Load the prepared phenotype data and the complete prepared-window dosage matrix.
6. Apply the optional sample allowlist and align samples across phenotype, genotype, and applicable covariate inputs.
7. Select shared covariates and the covariates assigned to the phenotype modality.
8. Apply the repository's existing missing-dosage handling, variant filtering, rank-inverse-normal phenotype transformation, covariate residualization, and genotype standardization.
9. Fit one SuSiE model with the repository's existing scientific defaults.
10. Validate, serialize, and commit the phenotype result to GCS.
11. Update the window manifest in GCS.
12. Continue with the next ordered phenotype.
13. Assemble the compatible per-window Parquet outputs and fit index after all phenotypes have a terminal status.
14. Mark the window manifest `COMPLETE` and upload it after the final outputs are valid.

The full prepared-window variant set is the candidate set. The wrapper can remove variants only for existing data-quality rules, such as no variation after sample alignment or the configured MAF rule. The wrapper must not filter variants by distance from the linked phenotype.

## Covariates and Sample Alignment

The WDL accepts parallel arrays of covariate files and modality labels. The label `shared` applies to every phenotype. A phenotype also receives each file labeled with its modality. Multiple files can use the same label. Combine applicable files in WDL array order, reject duplicate covariate IDs across combined files, and record the file order. Reject array-length mismatches, missing sample identifiers, and singular covariate designs that cannot be reduced safely.

Use the same aligned sample order for the phenotype vector, genotype matrix, and covariate matrix. Remove zero-variance covariate columns. Use a pivoted QR decomposition to identify linearly dependent columns, retain a deterministic full-rank basis in input-column order, and record removed columns. Fail input validation if the resulting design is not finite. Record the input and retained sample counts in each phenotype manifest. A deterministic lack of enough samples for the retained design is a committed `SKIPPED` result, not an unexpected software error.

## GCS Layout

Use a content-derived analysis ID and the supplied window ID:

```text
gs://<bucket>/<checkpoint-root>/<analysis-id>/<window-id>/
  window_manifest.json
  window_fit_index.tsv
  results/
    credible_sets.parquet
    lbf_variable.parquet
    full_susie.parquet
  phenotypes/
    <phenotype-key>/
      fit_manifest.json
      susie_fit.<sha256>.rds
      credible_sets.parquet
      lbf_variable.parquet
      full_susie.parquet
```

Derive `phenotype-key` from a SHA-256 hash of the window ID, modality, and phenotype ID. Store the original identifiers in each manifest. This rule prevents unsafe or ambiguous GCS object names.

The complete fitted `susie` object is the RDS payload. Add the ordered `variant_id` field before serialization. Model metadata remains in the JSON manifest so `readRDS()` returns the fitted object directly.

## Phenotype Commit Protocol

Use this transaction order for each phenotype:

1. Write the fit and result tables to local temporary paths.
2. Read the RDS back and run structural validation.
3. Read the Parquet files back and validate their schemas.
4. Calculate SHA-256 checksums and byte sizes.
5. Upload the content-addressed RDS and phenotype result tables.
6. Upload `fit_manifest.json` after all payload uploads succeed.
7. Update the local `window_manifest.json`.
8. Upload `window_manifest.json` last.

An uploaded payload without `fit_manifest.json` is an uncommitted orphan. A phenotype is committed only when its fixed-path phenotype manifest exists and validates. An interrupted attempt ignores an orphan and safely recomputes the phenotype.

GCS object replacement is the durability boundary. Local files do not count as checkpoints. Every progress update must be copied to the GCS window directory before the wrapper starts the next phenotype.

## Phenotype Manifest

Each `fit_manifest.json` records:

- Schema version.
- Analysis ID and window ID.
- Phenotype ID, phenotype key, modality, and processing index.
- Association P value.
- Terminal status.
- Input fingerprints and model settings.
- Container and package versions.
- Input and retained sample and variant counts.
- Fit convergence state.
- Payload GCS paths, SHA-256 checksums, and byte sizes.
- Start and completion timestamps.
- Expected exclusion reason when status is `SKIPPED`.

Valid terminal statuses are `CONVERGED`, `NONCONVERGED`, and `SKIPPED`.

## Window Manifest and Fast Resume

`window_manifest.json` is the authoritative resume cursor and final window inventory. It records:

- Schema version.
- Analysis ID and window ID.
- Ordered phenotype IDs and phenotype keys.
- Scientific settings and input fingerprints.
- Last committed zero-based index.
- The status and phenotype-manifest path for each committed item.
- Current window status: `RUNNING`, `COMPLETE`, or `FAILED`.
- Failure context when the window stops unexpectedly.
- Final result paths and checksums after completion.

On task start, download the exact window-manifest path. Do not list the window directory. Validate the analysis ID, input fingerprint, phenotype order, and latest committed phenotype manifest. Download and read the latest RDS. Resume at the next index only when all checks pass.

After the latest cursor check, probe the known fixed manifest path for the next phenotype. This handles an interruption after a phenotype commit but before the window-manifest update.

If the window manifest is absent or invalid, use binary search over the known ordered fixed phenotype-manifest paths to find the gap-free committed boundary. Validate the boundary manifest and RDS before resume. This recovery uses approximately `log2(n)` GCS existence checks and does not read all earlier manifests.

The workflow permits one active writer for each analysis-ID and window-ID pair. Cromwell preemptible attempts are sequential. Concurrent independent launches with the same pair are outside the first-version contract. The caller must not launch duplicate active workflows for the same pair.

## Fit Validation and Corruption Handling

A reusable checkpoint must pass all checks:

- The phenotype manifest is valid JSON with the expected schema version.
- The manifest analysis ID, window ID, phenotype ID, modality, and index match the current run.
- The RDS object exists at the recorded path.
- The downloaded RDS SHA-256 checksum and byte size match the manifest.
- `readRDS()` succeeds.
- The object inherits the expected SuSiE class.
- `variant_id`, `pip`, `alpha`, `mu`, and `mu2` dimensions agree.
- Variant identifiers are nonempty and unique.
- PIPs and required model arrays contain valid numeric values.
- The recorded convergence state matches the fit object.

If the current boundary checkpoint fails validation, treat that boundary index as uncommitted. Record the validation problem and recompute the same phenotype before processing later indexes. Do not delete a corrupt object automatically. A new content hash produces a new RDS path, and an atomic replacement of the fixed-path phenotype manifest points to the valid payload. Rewrite the window cursor at the same index only after the replacement checkpoint validates.

## Failure Policy

Commit and continue after deterministic input exclusions:

- No usable variants.
- Zero phenotype variance.
- No alternative allele after filtering.
- Too few aligned samples.

Save the full fit and continue when SuSiE returns a valid but non-converged object. Mark the phenotype `NONCONVERGED` in both manifests and output QC.

For an unexpected R error, invalid numerical state, schema error, GCS error, or serialization error:

1. Add the phenotype index, phenotype ID, error class, and concise message to the window manifest.
2. Set the window status to `FAILED`.
3. Upload the failed window manifest when GCS is available.
4. Exit nonzero.

The WDL backend then starts a new preemptible attempt up to the configured limit. `FAILED` records the state of the last attempt; it does not prohibit resume. A later attempt validates the saved cursor, changes the window status back to `RUNNING`, and resumes at the same uncommitted phenotype. All WDL command blocks and major R stages must emit explicit logging messages.

## WDL Design

`CheckpointedWindowSusie.wdl` contains one workflow and one main task. The workflow does not scatter. Terra or another outer scheduler launches one workflow per prepared window.

The task runtime exposes configurable memory, CPU, local SSD size, and preemptible-attempt count. The checkpoint root is a `String` because the target objects can be absent at workflow submission time. Production checkpoint roots must use `gs://`.

The task outputs are local final copies of:

- `window_manifest.json`.
- `window_fit_index.tsv`.
- `credible_sets.parquet`.
- `lbf_variable.parquet`.
- `full_susie.parquet`.

Cromwell uploads these normal WDL outputs after task success. The direct GCS checkpoint paths remain the source for retry recovery and the saved per-phenotype RDS files.

## Container Design

Create a separate checkpoint-runner image. Inherit from an immutable digest of the repository's existing micromamba-based SuSiE image. Add only the new wrapper and helper. Do not reinstall the scientific stack or use a mutable `main` tag as the production base.

The image must contain `gsutil`, R, the existing SuSiE utility files, the new wrapper, and the new checkpoint helper. The wrapper records the base-image digest and installed R package versions in every phenotype manifest.

Build and smoke-test the image only in GitHub Actions. Do not build the image locally for the smoke test.

## Result Tables

Use the existing SuSiE result-extraction helpers. Produce per-phenotype Parquet tables compatible with the current credible-set, variable-log-Bayes-factor, and full-variant schemas. The final window tables are row-bind merges in deterministic phenotype order.

`window_fit_index.tsv` contains one row per ordered phenotype with:

- Window ID.
- Processing index.
- Phenotype ID and modality.
- Association P value.
- Terminal status.
- Convergence state.
- RDS GCS path and checksum when a fit exists.
- Phenotype-manifest GCS path.
- Exclusion reason when applicable.

## Tests and Continuous Integration

Add unit tests for:

- Required manifest columns and value validation.
- Deterministic ordering and tie breaks.
- Analysis-fingerprint stability and invalidation.
- Safe phenotype-key generation.
- Exact window-variant use without phenotype-distance filtering.
- Modality-specific covariate selection.
- Sample alignment and expected exclusions.
- SuSiE object structural validation.
- Content checksum mismatch and unreadable RDS handling.
- Phenotype commit ordering.
- Window-manifest cursor updates.
- Resume from the next phenotype.
- Recovery when a phenotype commit exists but the window cursor lags.
- Logarithmic boundary recovery when the window manifest is absent.
- `NONCONVERGED`, `SKIPPED`, and unexpected-error behavior.
- Deterministic fit-index and window-table assembly.

Use a filesystem-backed checkpoint store in unit tests. The production store uses GCS. Both stores implement the same small interface for exact-path existence checks, uploads, downloads, and atomic manifest replacement.

Extend GitHub Actions to:

1. Run repository validation and R lint.
2. Run `miniwdl check` on the new WDL.
3. Build the checkpoint-runner image.
4. Confirm that the wrapper and helper exist in the image.
5. Confirm that all required R packages and `gsutil` are available.
6. Run the wrapper help command.
7. Run a small synthetic univariate fit.
8. Simulate interruption, resume, a corrupt RDS, an expected exclusion, and non-convergence.

The GitHub Actions smoke test must not require real workspace credentials or write to a real GCS bucket.

## Documentation

Update the repository README and workflow documentation with:

- The new one-window workflow purpose.
- The prepared input contract and required future `p_value` column.
- The GCS checkpoint layout.
- Resume and corruption behavior.
- The one-writer restriction.
- A template input JSON.
- Instructions for launching one workflow per prepared window on Terra.

## Acceptance Criteria

The design is complete when all of these conditions hold:

- Existing WDLs and `R/scripts/susie.R` have unchanged behavior.
- One task processes one already-filtered window manifest.
- Every fitted phenotype has a readable, checksummed SuSiE RDS in GCS.
- The GCS window manifest advances only after a phenotype commit.
- A preempted retry reads the window manifest and latest fit instead of all earlier fits.
- Boundary recovery does not require a full GCS listing.
- Expected exclusions and valid non-converged fits do not block later phenotypes.
- Unexpected errors fail the task and preserve the last valid cursor.
- Final Parquet outputs and the fit index match the committed phenotype set.
- The new container is pinned and smoke-tested in GitHub Actions.
- Local development does not build the Docker image for smoke testing.
