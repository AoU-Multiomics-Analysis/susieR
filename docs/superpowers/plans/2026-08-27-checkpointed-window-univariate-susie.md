# Checkpointed Window Univariate SuSiE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an isolated one-window univariate SuSiE workflow that saves and validates each phenotype fit in GCS and resumes ordered work after preemption.

**Architecture:** A new R controller reads an already-filtered prepared-window bundle, orders phenotypes by P value, fits each phenotype against every usable variant in the supplied window, and commits a full SuSiE RDS plus compatible Parquet tables before advancing a window cursor. A small checkpoint-store interface supports GCS in production and the local filesystem in tests. A new one-task WDL and a thin image leave `R/scripts/susie.R` and all existing workflows unchanged.

**Tech Stack:** WDL 1.0, R, tidyverse, data.table, arrow, jsonlite, digest, susieR, gsutil, micromamba-based OCI images, GitHub Actions, miniwdl.

**Spec:** `docs/superpowers/specs/2026-08-27-checkpointed-window-univariate-susie-design.md`

## Global Constraints

- Do not change the behavior or command-line contract of `R/scripts/susie.R` or an existing WDL.
- Treat `window_phenotypes.tsv` as the authoritative, already-filtered phenotype list. Do not read a global association table.
- Require manifest columns `window_id`, `phenotype_id`, `modality`, `phenotype_file`, and `p_value`.
- Sort by ascending `p_value`, then `modality`, then `phenotype_id`.
- Test every usable variant in the prepared window. Do not filter variants by distance from a phenotype.
- Use the existing model settings: `L = 10`, residual-variance estimation enabled, prior-variance estimation enabled, `scaled_prior_variance = 0.1`, univariate Z-score calculation enabled, and `min_abs_corr = 0.5`.
- Save a complete fitted `susie` object for each `CONVERGED` or `NONCONVERGED` phenotype.
- Upload phenotype payloads first, `fit_manifest.json` second, and `window_manifest.json` last.
- Read the exact window-manifest path and the latest committed fit on normal resume. Do not list all checkpoint objects.
- Commit deterministic exclusions as `SKIPPED`. Save valid non-converged fits as `NONCONVERGED`. Stop on unexpected errors.
- Run one workflow and one preemptible task per prepared window. Do not add a WDL scatter.
- Use tidyverse syntax for R data transformations and emit explicit WDL and R progress messages.
- Base the new image on `ghcr.io/aou-multiomics-analysis/susier@sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764`.
- Build and smoke-test the image in GitHub Actions. Do not build the image locally for smoke testing.

## File Structure

- `R/utils/CheckpointedWindowSusieFunctions.R`: Pure manifest, fingerprint, sample/covariate, fit-validation, model, and output-schema helpers.
- `R/utils/CheckpointStore.R`: Exact-path filesystem and GCS storage adapters plus phenotype commit and resume-boundary logic.
- `R/scripts/run_checkpointed_window_susie.R`: CLI parsing and the ordered window controller.
- `tests/test_helpers.R`: Small base-R assertion and fixture helpers.
- `tests/test_checkpointed_window_manifest.R`: Manifest, key, settings, and fingerprint tests.
- `tests/test_checkpoint_store.R`: Commit ordering, fit validation, cursor, and boundary-recovery tests.
- `tests/fixtures/checkpointed_window/generate_fixture.R`: Deterministic dosage, phenotype, covariate, and manifest fixture generator.
- `tests/test_checkpointed_window_model.R`: Full-window model and output-schema tests.
- `tests/test_checkpointed_window_controller.R`: Interruption, resume, corruption, exclusion, and final-index integration tests.
- `workflows/CheckpointedWindowSusie.wdl`: One-window, one-task preemptible workflow.
- `tests/test_checkpointed_window_susie_wdl.sh`: Static WDL contract checks.
- `examples/inputs/CheckpointedWindowSusie.inputs.json`: Terra-style example inputs.
- `containers/CheckpointedWindowSusie/Dockerfile`: Thin pinned image containing the new runner.
- `.github/workflows/checkpointed-window-susie-image.yml`: Image build, local-store integration smoke tests, and image publish.
- `tests/test_checkpointed_window_container.sh`: Static image and CI contract checks.
- `.dockstore.yml`: New workflow registration.
- `.github/workflows/wdl-validation.yml`: Explicit check for the new descriptor.
- `README.md`, `docs/README.md`, `docs/inputs.md`, `docs/scripts.md`, `docs/workflows.md`, `docs/docker.md`: User and maintainer documentation.

---

### Task 1: Prepared manifest and analysis identity

**Files:**
- Create: `tests/test_helpers.R`
- Create: `tests/test_checkpointed_window_manifest.R`
- Create: `R/utils/CheckpointedWindowSusieFunctions.R`

**Interfaces:**
- Consumes: A tibble with `window_id`, `phenotype_id`, `modality`, `phenotype_file`, and `p_value`.
- Produces: `checkpointed_susie_settings() -> named list`, `checkpoint_phenotype_key(window_id, modality, phenotype_id) -> character(1)`, `validate_window_phenotype_manifest(manifest, window_id) -> ordered tibble`, `sha256_file(path) -> character(1)`, and `build_checkpoint_analysis_id(input_hashes, ordered_manifest, settings, container_digest, wrapper_hash) -> character(1)`.

- [ ] **Step 1: Write shared assertions and failing manifest tests**

Create `tests/test_helpers.R` with `expect_error_message()`, `expect_identical_value()`, and `run_named_test()` helpers. Create `tests/test_checkpointed_window_manifest.R` with these concrete cases:

```r
source("tests/test_helpers.R")
source("R/utils/CheckpointedWindowSusieFunctions.R")

manifest <- tibble::tribble(
  ~window_id, ~phenotype_id, ~modality, ~phenotype_file, ~p_value,
  "chr1_0_2000000", "trait_b", "splicing", "window.bed.gz", 1e-6,
  "chr1_0_2000000", "trait_a", "expression", "window.bed.gz", 1e-8,
  "chr1_0_2000000", "trait_c", "expression", "window.bed.gz", 1e-6
)

ordered <- validate_window_phenotype_manifest(manifest, "chr1_0_2000000")
stopifnot(identical(ordered$phenotype_id, c("trait_a", "trait_c", "trait_b")))
stopifnot(identical(ordered$processing_index, 0:2))
stopifnot(length(unique(ordered$phenotype_key)) == 3L)

expect_error_message(
  validate_window_phenotype_manifest(dplyr::select(manifest, -p_value), "chr1_0_2000000"),
  "missing required columns: p_value"
)
expect_error_message(
  validate_window_phenotype_manifest(dplyr::bind_rows(manifest, manifest[1, ]), "chr1_0_2000000"),
  "duplicate phenotype_id"
)
```

Also assert that P values outside `[0, 1]`, multiple window IDs, empty identifiers, and a requested-window mismatch fail. Assert that the same canonical inputs produce the same analysis ID and that a changed P value, file hash, setting, or container digest changes it.

- [ ] **Step 2: Run the manifest test and verify the red state**

Run: `Rscript tests/test_checkpointed_window_manifest.R`

Expected: FAIL because `R/utils/CheckpointedWindowSusieFunctions.R` does not exist.

- [ ] **Step 3: Implement canonical manifest and fingerprint helpers**

Use tidyverse transformations and explicit checks. The core implementation must follow this shape:

```r
checkpointed_susie_settings <- function() {
  list(
    L = 10L,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    scaled_prior_variance = 0.1,
    compute_univariate_zscore = TRUE,
    min_abs_corr = 0.5
  )
}

checkpoint_phenotype_key <- function(window_id, modality, phenotype_id) {
  digest::digest(
    paste(window_id, modality, phenotype_id, sep = "\n"),
    algo = "sha256",
    serialize = FALSE
  )
}

validate_window_phenotype_manifest <- function(manifest, window_id) {
  required <- c("window_id", "phenotype_id", "modality", "phenotype_file", "p_value")
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0L) {
    stop("Prepared manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  validated <- manifest %>%
    dplyr::transmute(
      window_id = as.character(.data$window_id),
      phenotype_id = as.character(.data$phenotype_id),
      modality = as.character(.data$modality),
      phenotype_file = as.character(.data$phenotype_file),
      p_value = suppressWarnings(as.numeric(.data$p_value))
    )

  identifier_columns <- c("window_id", "phenotype_id", "modality", "phenotype_file")
  has_empty_identifier <- purrr::map_lgl(
    validated[identifier_columns],
    ~ any(is.na(.x) | !nzchar(.x))
  ) %>%
    any()
  if (has_empty_identifier) {
    stop("Prepared manifest identifiers cannot be empty.", call. = FALSE)
  }
  if (anyNA(validated$p_value) || any(!is.finite(validated$p_value)) ||
      any(validated$p_value < 0 | validated$p_value > 1)) {
    stop("Prepared manifest p_value values must be finite and between zero and one.", call. = FALSE)
  }
  if (anyDuplicated(validated$phenotype_id)) {
    stop("Prepared manifest contains duplicate phenotype_id values.", call. = FALSE)
  }
  manifest_windows <- unique(validated$window_id)
  if (length(manifest_windows) != 1L || !identical(manifest_windows, window_id)) {
    stop("Prepared manifest window_id does not match the requested window.", call. = FALSE)
  }

  validated %>%
    dplyr::arrange(.data$p_value, .data$modality, .data$phenotype_id) %>%
    dplyr::mutate(
      processing_index = dplyr::row_number() - 1L,
      phenotype_key = purrr::pmap_chr(
        list(.data$window_id, .data$modality, .data$phenotype_id),
        checkpoint_phenotype_key
      )
    )
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
```

Implement `build_checkpoint_analysis_id()` with this canonical payload:

```r
build_checkpoint_analysis_id <- function(
    input_hashes,
    ordered_manifest,
    settings,
    container_digest,
    wrapper_hash
) {
  canonical_payload <- list(
    input_hashes = as.list(input_hashes[order(names(input_hashes))]),
    phenotypes = ordered_manifest %>%
      dplyr::select(
        .data$processing_index,
        .data$window_id,
        .data$phenotype_id,
        .data$modality,
        .data$p_value,
        .data$phenotype_key
      ),
    settings = settings[order(names(settings))],
    container_digest = container_digest,
    wrapper_hash = wrapper_hash
  )
  canonical_json <- jsonlite::toJSON(
    canonical_payload,
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    digits = NA
  )
  digest::digest(canonical_json, algo = "sha256", serialize = FALSE)
}
```

- [ ] **Step 4: Run the manifest test and R parser**

Run: `Rscript tests/test_checkpointed_window_manifest.R`

Expected: PASS with each named test reported once.

Run: `Rscript -e 'parse(file="R/utils/CheckpointedWindowSusieFunctions.R")'`

Expected: exit 0.

- [ ] **Step 5: Commit the manifest contract**

```bash
git add tests/test_helpers.R tests/test_checkpointed_window_manifest.R R/utils/CheckpointedWindowSusieFunctions.R
git commit -m "feat: add checkpointed window manifest contract"
```

---

### Task 2: Checkpoint store, fit validation, and fast resume

**Files:**
- Create: `R/utils/CheckpointStore.R`
- Create: `tests/test_checkpoint_store.R`
- Modify: `R/utils/CheckpointedWindowSusieFunctions.R`

**Interfaces:**
- Consumes: Analysis ID, ordered phenotype manifest, local artifact paths, and a GCS URI or local checkpoint directory.
- Produces: `new_checkpoint_store(root, gsutil = "gsutil") -> store list`, `window_checkpoint_paths(analysis_id, window_id) -> named list`, `phenotype_checkpoint_paths(analysis_id, window_id, phenotype_key, fit_sha256 = NULL) -> named list`, `validate_checkpointed_susie_fit(fit) -> TRUE`, `commit_phenotype_checkpoint(store, paths, local_artifacts, fit_manifest) -> manifest`, `new_window_run_manifest(analysis_id, window_id, ordered_manifest, input_hashes, settings) -> list`, `advance_window_run_manifest(window_manifest, fit_manifest) -> list`, `fail_window_run_manifest(window_manifest, condition) -> list`, and `resolve_resume_boundary(store, ordered_manifest, window_manifest = NULL) -> list(last_committed_index, next_index, recovered, window_manifest)`.

- [ ] **Step 1: Write failing local-store and resume tests**

Create `tests/test_checkpoint_store.R`. Use a temporary local store and a structurally valid fake fit:

```r
fake_fit <- structure(
  list(
    pip = c(0.1, 0.9),
    alpha = matrix(c(0.1, 0.9), nrow = 1L),
    mu = matrix(c(0.2, 0.4), nrow = 1L),
    mu2 = matrix(c(0.05, 0.2), nrow = 1L),
    variant_id = c("chr1_100_A_G", "chr1_200_C_T"),
    converged = TRUE
  ),
  class = "susie"
)
stopifnot(isTRUE(validate_checkpointed_susie_fit(fake_fit)))
```

The test must also:

- Reject a checksum mismatch, duplicate variant IDs, a PIP outside `[0, 1]`, and inconsistent matrix dimensions.
- Record upload calls in a mock store and assert that the RDS and Parquet payloads precede `fit_manifest.json`.
- Create 1,024 ordered phenotype keys with the first 731 fixed commit paths present.
- Call `resolve_resume_boundary()` with no window manifest and assert `last_committed_index == 730L`.
- Assert that the boundary search uses no more than 12 exact existence checks.
- Create a cursor that ends at index 9 and a valid index-10 commit, then assert resume advances to index 11 without fitting index 10 again.

- [ ] **Step 2: Run the checkpoint-store test and verify the red state**

Run: `Rscript tests/test_checkpoint_store.R`

Expected: FAIL because `R/utils/CheckpointStore.R` and the checkpoint functions do not exist.

- [ ] **Step 3: Implement the storage interface and exact path rules**

Use this explicit store contract:

```r
new_checkpoint_store <- function(root, gsutil = "gsutil") {
  is_gcs <- startsWith(root, "gs://")
  clean_root <- sub("/+$", "", root)
  object_uri <- function(relative_path) paste0(clean_root, "/", relative_path)
  run_gsutil <- function(arguments, operation, uri) {
    status <- system2(gsutil, arguments)
    if (!identical(status, 0L)) {
      stop("GCS ", operation, " failed for: ", uri, call. = FALSE)
    }
    invisible(TRUE)
  }

  list(
    root = clean_root,
    is_gcs = is_gcs,
    object_exists = function(relative_path) {
      if (is_gcs) {
        return(identical(system2(gsutil, c("-q", "stat", object_uri(relative_path))), 0L))
      }
      file.exists(file.path(clean_root, relative_path))
    },
    upload = function(local_path, relative_path) {
      if (is_gcs) {
        return(run_gsutil(c("-q", "cp", local_path, object_uri(relative_path)), "upload", object_uri(relative_path)))
      }
      destination <- file.path(clean_root, relative_path)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(local_path, destination, overwrite = TRUE)) {
        stop("Filesystem upload failed for: ", destination, call. = FALSE)
      }
      invisible(TRUE)
    },
    download = function(relative_path, local_path) {
      dir.create(dirname(local_path), recursive = TRUE, showWarnings = FALSE)
      if (is_gcs) {
        return(run_gsutil(c("-q", "cp", object_uri(relative_path), local_path), "download", object_uri(relative_path)))
      }
      source_path <- file.path(clean_root, relative_path)
      if (!file.copy(source_path, local_path, overwrite = TRUE)) {
        stop("Filesystem download failed for: ", source_path, call. = FALSE)
      }
      invisible(TRUE)
    },
    object_uri = object_uri
  )
}
```

For GCS, call `gsutil -q stat`, `gsutil -q cp <local> <uri>`, and `gsutil -q cp <uri> <local>` with `system2()`. Check every exit status and include the failing operation and exact object URI in the error. Do not use `gsutil ls`.

Return window-level paths from `window_checkpoint_paths()` and phenotype-level paths from `phenotype_checkpoint_paths()`:

```r
list(
  window_manifest = file.path(analysis_id, window_id, "window_manifest.json"),
  fit_index = file.path(analysis_id, window_id, "window_fit_index.tsv"),
  credible_sets = file.path(analysis_id, window_id, "results", "credible_sets.parquet"),
  lbf_variable = file.path(analysis_id, window_id, "results", "lbf_variable.parquet"),
  full_susie = file.path(analysis_id, window_id, "results", "full_susie.parquet")
)

list(
  fit_manifest = file.path(analysis_id, window_id, "phenotypes", phenotype_key, "fit_manifest.json"),
  fit_rds = file.path(analysis_id, window_id, "phenotypes", phenotype_key, paste0("susie_fit.", fit_sha256, ".rds")),
  credible_sets = file.path(analysis_id, window_id, "phenotypes", phenotype_key, "credible_sets.parquet"),
  lbf_variable = file.path(analysis_id, window_id, "phenotypes", phenotype_key, "lbf_variable.parquet"),
  full_susie = file.path(analysis_id, window_id, "phenotypes", phenotype_key, "full_susie.parquet")
)
```

Write JSON to a local temporary file with `jsonlite::write_json(auto_unbox = TRUE, pretty = TRUE, null = "null")`. Upload `fit_manifest.json` after all payloads. Upload `window_manifest.json` only through a separate cursor update.

- [ ] **Step 4: Implement fit validation and boundary recovery**

Add `validate_checkpointed_susie_fit()` to `CheckpointedWindowSusieFunctions.R`. Require class `susie`, nonempty unique `variant_id`, finite PIPs in `[0, 1]`, equal variant dimensions for `pip`, `alpha`, `mu`, and `mu2`, and a scalar logical `converged` value.

Implement binary search only over fixed `fit_manifest.json` paths. The committed path sequence is a gap-free prefix because the controller processes one phenotype at a time. After binary search, download and validate only the boundary manifest. For `CONVERGED` and `NONCONVERGED`, also download and validate the boundary RDS. For `SKIPPED`, require an allowed exclusion reason and no RDS payload. If the cursor points at an invalid boundary, return that index as `next_index` so the controller recomputes it.

- [ ] **Step 5: Run checkpoint tests and verify they pass**

Run: `Rscript tests/test_checkpoint_store.R`

Expected: PASS, including the `<= 12` exact-check assertion for 1,024 phenotypes.

Run: `Rscript tests/test_checkpointed_window_manifest.R`

Expected: PASS.

- [ ] **Step 6: Commit checkpoint storage and resume logic**

```bash
git add R/utils/CheckpointStore.R R/utils/CheckpointedWindowSusieFunctions.R tests/test_checkpoint_store.R
git commit -m "feat: add durable phenotype checkpoints"
```

---

### Task 3: Full-window univariate model adapter

**Files:**
- Create: `tests/fixtures/checkpointed_window/generate_fixture.R`
- Create: `tests/test_checkpointed_window_model.R`
- Modify: `R/utils/CheckpointedWindowSusieFunctions.R`

**Interfaces:**
- Consumes: A prepared wide dosage table, phenotype BED-style table, ordered manifest row, applicable covariate tables, and optional sample allowlist.
- Produces: `read_checkpointed_window_dosage(path) -> list(genotype, variant_info)`, `read_checkpointed_window_phenotypes(path, ordered_manifest) -> list(values, metadata)`, `read_checkpointed_covariates(paths, labels) -> named list`, `applicable_covariate_matrix(covariates, modality, sample_ids) -> matrix`, `load_checkpointed_window_inputs(config, ordered_manifest) -> list(genotype, variant_info, phenotypes, phenotype_metadata, covariates, sample_ids)`, `fit_checkpointed_window_phenotype(genotype, phenotype, covariates, variant_ids, settings) -> list(status, fit, reason, qc)`, `format_checkpointed_susie_tables(result, phenotype_record) -> named list of tibbles`, and `write_checkpointed_susie_tables(tables, output_dir) -> named file paths`.

- [ ] **Step 1: Generate a deterministic trans-window fixture**

Create a fixture generator with `set.seed(20260827)`, 40 samples, six chromosome-1 window variants, two linked phenotypes whose BED coordinates are on chromosomes 7 and 12, one shared covariate, and one modality-specific covariate. Write:

```text
tests/fixtures/checkpointed_window/window_dosage.tsv
tests/fixtures/checkpointed_window/window_phenotypes.bed.gz
tests/fixtures/checkpointed_window/window_phenotypes.tsv
tests/fixtures/checkpointed_window/shared_covariates.tsv
tests/fixtures/checkpointed_window/expression_covariates.tsv
tests/fixtures/checkpointed_window/keep_samples.tsv
```

Set the first phenotype to contain a strong contribution from window variants 2 and 5. Put `-1` in one dosage cell to exercise missing-dosage imputation. The manifest must already contain only the two linked phenotypes and must include finite P values.

- [ ] **Step 2: Write failing model and schema tests**

Create `tests/test_checkpointed_window_model.R` to assert:

```r
dosage <- read_checkpointed_window_dosage(fixture_path("window_dosage.tsv"))
stopifnot(nrow(dosage$genotype) == 6L)
stopifnot(all(is.finite(dosage$genotype)))

result <- fit_checkpointed_window_phenotype(
  genotype = dosage$genotype,
  phenotype = phenotype_values,
  covariates = covariate_matrix,
  variant_ids = dosage$variant_info$variant_id,
  settings = checkpointed_susie_settings()
)
stopifnot(result$status %in% c("CONVERGED", "NONCONVERGED"))
stopifnot(identical(result$fit$variant_id, dosage$variant_info$variant_id))
```

The phenotype chromosome must differ from the dosage chromosome. The equality above is the regression check that no phenotype-distance filter runs. Also test zero phenotype variance, no variable variants, too few aligned samples, modality-specific covariate selection, duplicate covariate IDs, and the exact required columns in all three output tables.

- [ ] **Step 3: Run the model test and verify the red state**

Run: `Rscript tests/fixtures/checkpointed_window/generate_fixture.R`

Expected: fixture files are created.

Run: `Rscript tests/test_checkpointed_window_model.R`

Expected: FAIL because the model adapter functions do not exist.

- [ ] **Step 4: Implement window input readers and covariate alignment**

Read the dosage with `data.table::fread()`. Require metadata columns `CHROM`, `POS`, `REF`, and `ALT`, accepting `#CHROM` as the chromosome-column spelling. Build `variant_id` as `CHROM_POS_REF_ALT`. Convert exact `-1` values to `NA_real_`, reject other finite values outside `[0, 2]`, mean-impute each variant, and retain the manifest variant order.

Define covariate files as covariate-by-sample tables. The first column is the covariate ID and remaining column names are sample IDs. Combine `shared` files and files matching the phenotype modality in WDL array order. Reject duplicate covariate IDs. Remove zero-variance columns, use pivoted QR to retain a deterministic full-rank basis, and return a sample-by-covariate matrix.

- [ ] **Step 5: Implement the full-window fit path**

Use tidyverse transformations for metadata and this numerical sequence for each phenotype:

```r
rank_inverse_normal <- function(values) {
  stats::qnorm((rank(values, ties.method = "average") - 0.5) / length(values))
}

design <- cbind(intercept = 1, covariates)
design_qr <- qr(design)
phenotype_residual <- qr.resid(design_qr, rank_inverse_normal(phenotype))
genotype_samples_by_variants <- t(genotype)
genotype_residual <- qr.resid(design_qr, genotype_samples_by_variants)

fit <- susieR::susie(
  genotype_residual,
  phenotype_residual,
  L = settings$L,
  estimate_residual_variance = settings$estimate_residual_variance,
  estimate_prior_variance = settings$estimate_prior_variance,
  scaled_prior_variance = settings$scaled_prior_variance,
  compute_univariate_zscore = settings$compute_univariate_zscore,
  min_abs_corr = settings$min_abs_corr,
  verbose = TRUE
)
fit$variant_id <- variant_ids
```

Filter variants only when they have no finite imputation mean or no variation after sample alignment. Return `SKIPPED` with a machine-readable reason for expected exclusions. Do not catch unexpected `susieR::susie()` errors.

Use existing `extractResults()` output and the empty constructors in `InitFunctions.R` where possible. Add a new formatter rather than modifying `R/scripts/susie.R`. Set `molecular_trait_id` to the manifest phenotype ID and `region` to the window ID. Preserve the current credible-set columns from `InitEmptyInCSVariantDf()`, the LBF columns from `InitEmptyLbfDf()`, and the full-variant columns from `InitEmptyVariantDf()`. Write those typed empty schemas for `SKIPPED` results and for valid fits with no credible set.

- [ ] **Step 6: Run model and earlier tests**

Run: `Rscript tests/test_checkpointed_window_model.R`

Expected: PASS with both linked phenotypes using all six prepared variants before quality filtering.

Run: `Rscript tests/test_checkpointed_window_manifest.R && Rscript tests/test_checkpoint_store.R`

Expected: both scripts pass.

- [ ] **Step 7: Commit the full-window model adapter**

```bash
git add R/utils/CheckpointedWindowSusieFunctions.R tests/fixtures/checkpointed_window tests/test_checkpointed_window_model.R
git commit -m "feat: fit linked phenotypes across prepared windows"
```

---

### Task 4: Ordered checkpointed controller and CLI

**Files:**
- Create: `R/scripts/run_checkpointed_window_susie.R`
- Create: `tests/test_checkpointed_window_controller.R`
- Modify: `R/utils/CheckpointedWindowSusieFunctions.R`
- Modify: `R/utils/CheckpointStore.R`

**Interfaces:**
- Consumes: `run_checkpointed_window(config, store = NULL, fit_function = fit_checkpointed_window_phenotype, interrupt_after_commits = Inf)` and the CLI flags defined below.
- Produces: `write_local_phenotype_artifacts(result, phenotype_record, directory) -> named file paths`, `build_fit_manifest(analysis_id, window_id, phenotype_record, fit_result, local_artifacts, input_hashes, settings) -> list`, `read_window_run_manifest_if_present(store, relative_path) -> list or NULL`, `upload_window_run_manifest(store, relative_path, window_manifest) -> TRUE`, `assemble_checkpointed_window_outputs(store, committed_records, output_dir) -> named file paths`, `commit_window_outputs(store, analysis_id, window_id, local_paths) -> named GCS records`, `complete_window_run_manifest(window_manifest, committed_window_outputs) -> list`, and a completed local output directory with `window_manifest.json`, `window_fit_index.tsv`, `credible_sets.parquet`, `lbf_variable.parquet`, and `full_susie.parquet`; durable per-phenotype GCS checkpoints after every commit.

- [ ] **Step 1: Write a failing interruption-and-resume integration test**

Create `tests/test_checkpointed_window_controller.R`. Use the filesystem checkpoint store and the fixture from Task 3. Inject a counting fit function. Run once with `interrupt_after_commits = 1L`, assert one committed phenotype manifest and a `RUNNING` window cursor, then run again without interruption and assert only the second phenotype is fitted.

Add these scenarios:

- Corrupt the latest committed RDS after the interrupted run. Assert the retry recomputes that boundary phenotype.
- Leave a valid second phenotype commit while the cursor still ends at the first. Assert the retry adopts the second commit without fitting it again.
- Return `SKIPPED` for the first phenotype and assert the cursor advances.
- Return a fake valid fit with `converged = FALSE` and assert the RDS is saved with `NONCONVERGED`.
- Throw `simpleError("synthetic unexpected failure")` and assert the window manifest is uploaded as `FAILED` and the controller rethrows.
- Assert the final fit index is in processing order and contains RDS URIs only for fitted phenotypes.

- [ ] **Step 2: Run the controller test and verify the red state**

Run: `Rscript tests/test_checkpointed_window_controller.R`

Expected: FAIL because `run_checkpointed_window()` and the CLI script do not exist.

- [ ] **Step 3: Implement the ordered controller**

Implement the controller in this order:

```r
run_checkpointed_window <- function(
    config,
    store = NULL,
    fit_function = fit_checkpointed_window_phenotype,
    interrupt_after_commits = Inf
) {
  window_state <- NULL
  window_paths <- NULL
  run_window_body <- function() {
    raw_manifest <- readr::read_tsv(config$window_phenotypes, show_col_types = FALSE)
    ordered_manifest <- validate_window_phenotype_manifest(raw_manifest, config$window_id)
    covariate_hash_names <- paste0(
      "covariate_",
      seq_along(config$covariate_files),
      "_",
      config$covariate_modalities
    )
    scientific_paths <- c(
      window_dosage = config$window_dosage,
      window_phenotypes = config$window_phenotypes,
      phenotype_data = config$phenotype_data,
      stats::setNames(config$covariate_files, covariate_hash_names)
    )
    if (!is.null(config$keep_samples)) {
      scientific_paths <- c(scientific_paths, keep_samples = config$keep_samples)
    }
    input_hashes <- purrr::map_chr(
      scientific_paths,
      sha256_file
    )
    settings <- checkpointed_susie_settings()
    analysis_id <- build_checkpoint_analysis_id(
    input_hashes,
    ordered_manifest,
    settings,
    Sys.getenv("CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST"),
    sha256_file(config$wrapper_path)
  )
    store <<- store %||% new_checkpoint_store(config$checkpoint_root)
    window_paths <<- window_checkpoint_paths(analysis_id, config$window_id)
    saved_window_manifest <- read_window_run_manifest_if_present(store, window_paths$window_manifest)
    resume <- resolve_resume_boundary(store, ordered_manifest, saved_window_manifest)
    window_state <<- resume$window_manifest %||% new_window_run_manifest(
    analysis_id,
    config$window_id,
    ordered_manifest,
    input_hashes,
    settings
  )
    window_state$status <<- "RUNNING"

    model_inputs <- load_checkpointed_window_inputs(config, ordered_manifest)
    commit_count <- 0L
    indexes <- ordered_manifest$processing_index[ordered_manifest$processing_index >= resume$next_index]
    for (index in indexes) {
    phenotype_record <- ordered_manifest %>% dplyr::filter(.data$processing_index == index)
    fit_result <- fit_function(
      genotype = model_inputs$genotype,
      phenotype = model_inputs$phenotypes[[phenotype_record$phenotype_id[[1L]]]],
      covariates = applicable_covariate_matrix(
        model_inputs$covariates,
        phenotype_record$modality[[1L]],
        model_inputs$sample_ids
      ),
      variant_ids = model_inputs$variant_info$variant_id,
      settings = settings
    )
    local_artifacts <- write_local_phenotype_artifacts(
      fit_result,
      phenotype_record,
      tempfile("checkpointed-phenotype-")
    )
    fit_manifest <- build_fit_manifest(
      analysis_id,
      config$window_id,
      phenotype_record,
      fit_result,
      local_artifacts,
      input_hashes,
      settings
    )
    fit_sha256 <- if (identical(fit_result$status, "SKIPPED")) {
      NULL
    } else {
      fit_manifest$payloads$susie_fit$sha256
    }
    committed <- commit_phenotype_checkpoint(
      store,
      phenotype_checkpoint_paths(
        analysis_id,
        config$window_id,
        phenotype_record$phenotype_key[[1L]],
        fit_sha256
      ),
      local_artifacts,
      fit_manifest
    )
      window_state <<- advance_window_run_manifest(window_state, committed)
      upload_window_run_manifest(store, window_paths$window_manifest, window_state)
      commit_count <- commit_count + 1L
      if (commit_count >= interrupt_after_commits) {
        stop(structure(list(message = "Synthetic interruption"), class = c("checkpoint_test_interrupt", "error", "condition")))
      }
    }

    final_paths <- assemble_checkpointed_window_outputs(
      store,
      window_state$committed,
      config$output_dir
    )
    committed_window_outputs <- commit_window_outputs(
      store,
      analysis_id,
      config$window_id,
      final_paths
    )
    window_state <<- complete_window_run_manifest(window_state, committed_window_outputs)
    upload_window_run_manifest(store, window_paths$window_manifest, window_state)
    final_paths
  }

  tryCatch(
    run_window_body(),
    checkpoint_test_interrupt = function(condition) stop(condition),
    error = function(condition) {
      if (!is.null(window_state) && !is.null(window_paths)) {
        failed_state <- fail_window_run_manifest(window_state, condition)
        try(upload_window_run_manifest(store, window_paths$window_manifest, failed_state), silent = TRUE)
      }
      stop(condition)
    }
  )
}
```

Define `%||%` in the helper as `function(x, y) if (is.null(x)) y else x`. `read_window_run_manifest_if_present()`, `upload_window_run_manifest()`, `build_fit_manifest()`, `commit_window_outputs()`, and `complete_window_run_manifest()` are implemented in this task and tested through the controller scenarios. `build_fit_manifest()` records schema version, timestamps, package versions, counts, convergence, status, checksums, byte sizes, and GCS object URIs.

The `checkpoint_test_interrupt` handler preserves a `RUNNING` cursor for the synthetic interruption test. The general error handler writes `FAILED` for real unexpected errors after state initialization.

After a phenotype commit and window-cursor upload, emit `message("[checkpointed-window] Committed ", phenotype_id, " at index ", index)`. Do not start the next phenotype before the GCS cursor upload returns success.

For a `FAILED` cursor from an earlier attempt, validate the saved boundary, set status to `RUNNING`, clear only the attempt failure field, and resume. Preserve prior failure history in an array.

- [ ] **Step 4: Implement the CLI without changing existing option files**

Use a separate `optparse` list in `run_checkpointed_window_susie.R`. Provide these flags:

```text
--window-id
--window-dosage
--window-phenotypes
--phenotype-data
--covariate-files
--covariate-modalities
--keep-samples
--checkpoint-root
--output-dir
```

Require equal nonzero covariate arrays. Accept comma-separated arrays with no empty members. Use `SUSIER_FUNCTIONS_PATH` for existing utility files and source the new installed helpers from the same directory. Put `main()` behind `if (sys.nframe() == 0L)` so tests can source the controller.

- [ ] **Step 5: Run controller, model, and CLI checks**

Run: `Rscript tests/test_checkpointed_window_controller.R`

Expected: PASS for interruption, corruption, lagging cursor, skip, non-convergence, and unexpected failure.

Run: `Rscript R/scripts/run_checkpointed_window_susie.R --help`

Expected: exit 0 and all nine flags are present.

Run: `Rscript tests/test_checkpointed_window_manifest.R && Rscript tests/test_checkpoint_store.R && Rscript tests/test_checkpointed_window_model.R`

Expected: all pass.

- [ ] **Step 6: Commit the controller**

```bash
git add R/scripts/run_checkpointed_window_susie.R R/utils/CheckpointedWindowSusieFunctions.R R/utils/CheckpointStore.R tests/test_checkpointed_window_controller.R
git commit -m "feat: resume ordered window fine-mapping"
```

---

### Task 5: One-window preemptible WDL

**Files:**
- Create: `workflows/CheckpointedWindowSusie.wdl`
- Create: `tests/test_checkpointed_window_susie_wdl.sh`
- Create: `examples/inputs/CheckpointedWindowSusie.inputs.json`
- Modify: `.dockstore.yml`
- Modify: `.github/workflows/wdl-validation.yml`

**Interfaces:**
- Consumes: `window_id`, `window_dosage`, `window_phenotypes`, `phenotype_data`, `covariate_files`, `covariate_modalities`, optional `keep_samples`, `checkpoint_root`, runtime settings, and a runner-image string.
- Produces: WDL outputs `WindowManifest`, `WindowFitIndex`, `SusieParquet`, `SusielbfParquet`, and `FullSusieParquet`.

- [ ] **Step 1: Write the failing WDL contract test**

Create `tests/test_checkpointed_window_susie_wdl.sh` with `set -euo pipefail`. Require:

```bash
miniwdl check workflows/CheckpointedWindowSusie.wdl
rg -q '^workflow CheckpointedWindowSusieWorkflow' workflows/CheckpointedWindowSusie.wdl
rg -q '^task RunCheckpointedWindowSusie' workflows/CheckpointedWindowSusie.wdl
rg -q 'preemptible: preemptible_attempts' workflows/CheckpointedWindowSusie.wdl
rg -q 'Rscript /opt/r/scripts/run_checkpointed_window_susie.R' workflows/CheckpointedWindowSusie.wdl
test "$(rg -c 'scatter[[:space:]]*\(' workflows/CheckpointedWindowSusie.wdl || true)" -eq 0
```

Also check the five WDL outputs, explicit start/completion log messages, all wrapper arguments, the Dockstore entry, and the example JSON keys.
Require a command-block guard that rejects a checkpoint root without the `gs://` prefix.

- [ ] **Step 2: Run the WDL test and verify the red state**

Run: `bash tests/test_checkpointed_window_susie_wdl.sh`

Expected: FAIL because the WDL does not exist.

- [ ] **Step 3: Implement the one-task WDL**

Use this workflow input contract:

```wdl
version 1.0

workflow CheckpointedWindowSusieWorkflow {
  input {
    String window_id
    File window_dosage
    File window_phenotypes
    File phenotype_data
    Array[File] covariate_files
    Array[String] covariate_modalities
    File? keep_samples
    String checkpoint_root
    String runner_image = "ghcr.io/aou-multiomics-analysis/susier/checkpointed-window:main"
    Int memory_gb = 16
    Int cpu = 1
    Int disk_gb = 500
    Int preemptible_attempts = 3
  }

  call RunCheckpointedWindowSusie {
    input:
      window_id = window_id,
      window_dosage = window_dosage,
      window_phenotypes = window_phenotypes,
      phenotype_data = phenotype_data,
      covariate_files = covariate_files,
      covariate_modalities = covariate_modalities,
      keep_samples = keep_samples,
      checkpoint_root = checkpoint_root,
      runner_image = runner_image,
      memory_gb = memory_gb,
      cpu = cpu,
      disk_gb = disk_gb,
      preemptible_attempts = preemptible_attempts
  }

  output {
    File WindowManifest = RunCheckpointedWindowSusie.WindowManifest
    File WindowFitIndex = RunCheckpointedWindowSusie.WindowFitIndex
    File SusieParquet = RunCheckpointedWindowSusie.SusieParquet
    File SusielbfParquet = RunCheckpointedWindowSusie.SusielbfParquet
    File FullSusieParquet = RunCheckpointedWindowSusie.FullSusieParquet
  }
}
```

The task command must start with `set -euo pipefail`, log the window ID and checkpoint root, and reject a production checkpoint root that does not start with `gs://`. Create `window_outputs`, call the wrapper once, test all five output files, and log completion. Use `~{sep="," covariate_files}` and `~{sep="," covariate_modalities}`. Add the optional keep-samples flag only when defined. The R controller keeps local-directory support for tests, but the WDL exposes only the GCS production path.

Use `runner_image`, `memory_gb`, `cpu`, `disk_gb`, and `preemptible_attempts` in the runtime block. Confirm the exact interpolation syntax with `miniwdl check`.

- [ ] **Step 4: Register and document machine-readable inputs**

Add this Dockstore entry:

```yaml
  - subclass: WDL
    primaryDescriptorPath: /workflows/CheckpointedWindowSusie.wdl
    name: "CheckpointedWindowSusie"
```

Add `miniwdl check CheckpointedWindowSusie.wdl` to `.github/workflows/wdl-validation.yml`.

Create an example JSON that uses `gs://example-bucket/prepared/chr1_0_2000000/` inputs and `gs://example-bucket/checkpoints/window-susie` as the checkpoint root. Set `preemptible_attempts` to `3`.

- [ ] **Step 5: Run WDL, YAML, and JSON validation**

Run: `bash tests/test_checkpointed_window_susie_wdl.sh`

Expected: PASS.

Run: `ruby -e 'require "yaml"; YAML.load_file(".dockstore.yml"); YAML.load_file(".github/workflows/wdl-validation.yml"); puts "YAML ok"'`

Expected: `YAML ok`.

Run: `ruby -e 'require "json"; JSON.parse(File.read("examples/inputs/CheckpointedWindowSusie.inputs.json")); puts "JSON ok"'`

Expected: `JSON ok`.

- [ ] **Step 6: Commit the WDL**

```bash
git add workflows/CheckpointedWindowSusie.wdl tests/test_checkpointed_window_susie_wdl.sh examples/inputs/CheckpointedWindowSusie.inputs.json .dockstore.yml .github/workflows/wdl-validation.yml
git commit -m "feat: add checkpointed window SuSiE workflow"
```

---

### Task 6: Thin pinned image and GitHub Actions smoke tests

**Files:**
- Create: `containers/CheckpointedWindowSusie/Dockerfile`
- Create: `.github/workflows/checkpointed-window-susie-image.yml`
- Create: `tests/test_checkpointed_window_container.sh`

**Interfaces:**
- Consumes: The wrapper and new helper files from Tasks 1-4 and the pinned existing SuSiE image.
- Produces: `ghcr.io/aou-multiomics-analysis/susier/checkpointed-window` with the installed runner at `/opt/r/scripts/run_checkpointed_window_susie.R`.

- [ ] **Step 1: Write a failing static container contract test**

Create `tests/test_checkpointed_window_container.sh`. Require the exact base digest, installed helper paths, `CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST`, the runner symlink, GitHub Actions path triggers, image build, wrapper `--help`, package checks, all four R tests, and no GCS credentials.

The digest check must be exact:

```bash
rg -Fq 'FROM ghcr.io/aou-multiomics-analysis/susier@sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764' containers/CheckpointedWindowSusie/Dockerfile
```

- [ ] **Step 2: Run the container contract test and verify the red state**

Run: `bash tests/test_checkpointed_window_container.sh`

Expected: FAIL because the Dockerfile and workflow do not exist.

- [ ] **Step 3: Create the thin micromamba-derived image**

Use this structure without package installation:

```dockerfile
FROM ghcr.io/aou-multiomics-analysis/susier@sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764

ENV CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST=sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764
ENV SUSIER_FUNCTIONS_PATH=/opt/r/lib

USER root
COPY R/utils/CheckpointedWindowSusieFunctions.R R/utils/CheckpointStore.R /opt/r/lib/
COPY R/scripts/run_checkpointed_window_susie.R /opt/r/scripts/
RUN ln -sf /opt/r/scripts/run_checkpointed_window_susie.R /tmp/run_checkpointed_window_susie.R
USER $MAMBA_USER

CMD ["bash"]
```

Do not add `micromamba install`, `R install.packages()`, or a source checkout. The pinned parent already contains the scientific environment and `gsutil`.

- [ ] **Step 4: Add the GitHub Actions image workflow**

Follow the existing Docker workflow pattern. Build `susier-checkpointed-window:smoke-test` with `load: true`. The smoke step must:

```bash
set -euo pipefail
docker run --rm susier-checkpointed-window:smoke-test test -f /opt/r/scripts/run_checkpointed_window_susie.R
docker run --rm susier-checkpointed-window:smoke-test test -f /opt/r/lib/CheckpointedWindowSusieFunctions.R
docker run --rm susier-checkpointed-window:smoke-test test -f /opt/r/lib/CheckpointStore.R
docker run --rm susier-checkpointed-window:smoke-test Rscript /opt/r/scripts/run_checkpointed_window_susie.R --help
docker run --rm susier-checkpointed-window:smoke-test Rscript -e 'library(tidyverse); library(data.table); library(arrow); library(jsonlite); library(digest); library(susieR); stopifnot(nzchar(Sys.which("gsutil")))'
docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpointed_window_manifest.R
docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpoint_store.R
docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpointed_window_model.R
docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpointed_window_controller.R
```

Push `ghcr.io/${{ github.repository }}/checkpointed-window` only on pushes to `main`. Give the job `packages: write`. Do not use GCS credentials or a real bucket.

- [ ] **Step 5: Run only static validation locally**

Run: `bash tests/test_checkpointed_window_container.sh`

Expected: PASS.

Run: `ruby -e 'require "yaml"; YAML.load_file(".github/workflows/checkpointed-window-susie-image.yml"); puts "YAML ok"'`

Expected: `YAML ok`.

Do not run `docker build` locally. The GitHub Actions job is the image build and smoke-test authority.

- [ ] **Step 6: Commit the image workflow**

```bash
git add containers/CheckpointedWindowSusie/Dockerfile .github/workflows/checkpointed-window-susie-image.yml tests/test_checkpointed_window_container.sh
git commit -m "ci: build checkpointed window SuSiE image"
```

---

### Task 7: Documentation and complete verification

**Files:**
- Modify: `README.md`
- Modify: `docs/README.md`
- Modify: `docs/inputs.md`
- Modify: `docs/scripts.md`
- Modify: `docs/workflows.md`
- Modify: `docs/docker.md`

**Interfaces:**
- Consumes: The completed CLI, WDL, checkpoint schema, image, and example input.
- Produces: User-facing instructions that define prepared-input ownership, GCS durability, resume behavior, outputs, and Terra launch rules.

- [ ] **Step 1: Write a failing documentation contract check**

Extend `tests/test_checkpointed_window_susie_wdl.sh` to require each of these exact terms in the documentation set:

```text
CheckpointedWindowSusieWorkflow
window_phenotypes.tsv
p_value
window_manifest.json
window_fit_index.tsv
NONCONVERGED
SKIPPED
one workflow per prepared window
```

Require the README workflow table and Docker image table to include the new descriptor and image.

- [ ] **Step 2: Run the documentation contract and verify the red state**

Run: `bash tests/test_checkpointed_window_susie_wdl.sh`

Expected: FAIL because the documentation does not yet describe the new workflow.

- [ ] **Step 3: Document the workflow and input ownership**

Update documentation with these explicit statements:

- The prepare workflow owns phenotype filtering.
- The analysis wrapper trusts the prepared manifest and requires `p_value`.
- The first version uses every usable variant in the supplied window and does not calculate a phenotype-centered interval.
- Terra launches one workflow per prepared window.
- The checkpoint root must be a writable `gs://` prefix.
- The RDS is the full fitted `susie` object, and `window_fit_index.tsv` locates each fit.
- GCS payloads precede the phenotype manifest, and the window manifest uploads last.
- Normal resume reads only the window cursor and latest fit.
- `SKIPPED` and `NONCONVERGED` do not block later phenotypes.
- The new runner image inherits the exact pinned scientific base digest and is built only in GitHub Actions.

Use ASD-STE100-style short sentences for technical instructions. Do not add plots.

- [ ] **Step 4: Run all non-container verification commands**

Run each command separately and inspect the complete output:

```bash
Rscript tests/fixtures/checkpointed_window/generate_fixture.R
Rscript tests/test_checkpointed_window_manifest.R
Rscript tests/test_checkpoint_store.R
Rscript tests/test_checkpointed_window_model.R
Rscript tests/test_checkpointed_window_controller.R
bash tests/test_checkpointed_window_susie_wdl.sh
bash tests/test_checkpointed_window_container.sh
Rscript tools/lint_r.R
tools/validate_repo.sh
git diff --check
```

Expected: every command exits 0. Do not claim the Docker image builds locally. Report the GitHub Actions image build as pending until the workflow runs on GitHub.

- [ ] **Step 5: Audit the implementation against the design acceptance criteria**

Read `docs/superpowers/specs/2026-08-27-checkpointed-window-univariate-susie-design.md` line by line. Confirm that each acceptance criterion has a code path and a test. Record any GitHub-only validation as such. Confirm `git diff -- R/scripts/susie.R workflows/susieR.wdl workflows/susieRonly.wdl` is empty.

- [ ] **Step 6: Commit documentation and final validation updates**

```bash
git add README.md docs/README.md docs/inputs.md docs/scripts.md docs/workflows.md docs/docker.md tests/test_checkpointed_window_susie_wdl.sh
git commit -m "docs: document checkpointed window fine-mapping"
```

- [ ] **Step 7: Request code review**

Use the `superpowers:requesting-code-review` skill. Give the reviewer the spec, this plan, the full commit range beginning after `8149aa8`, all verification output, and the explicit note that Docker build verification belongs to GitHub Actions.
