# Native PGEN mean imputation

The user approved native mean imputation for the current per-feature workflow.
Genotype caching and changes to upstream batch files are outside this change.

For ordinary PGEN-only fits, select the feature's exact samples before decoding
and use `pgenlibr::ReadList(meanimpute = TRUE)` on the first window read. Mark the
matrix as already imputed so the shared fitting helper does not impute again.
Reuse parsed PVAR metadata and close temporary subset handles after each read.

Active full-cohort ancestry MAF filtering still needs raw calls first. Keep
that raw read and apply native imputation after variant and sample selection.
Inputs with additional genotypes retain their existing merge and preparation.
The ordinary path must make no hardcall read and one native imputed read for
a window without entirely missing variants. No new genotype cache is added.

ReadList uses stored dosages. The user confirmed that the input dosages are hard
calls. Compare native values and complete fits with the current R calculation
on real hardcall PGEN fixtures before adopting the native call. The user asked
for no new runtime numerical comparison or imputation fallback.

The user subsequently requested removal of entirely missing variants. The
first native read can identify them with its specific all-missing error. Remove
only the reported global variant index, then retry. Propagate every other
error. Keep existing empty-window handling and outcome statuses. PGEN-only raw
filtering uses NA rather than -1, which prevents missing sentinels from canceling
real ALT counts. Compare both missing-code configurations against this intended
NA-based filter rule and compare complete fits on the same retained variants.

No WDL input, localization, container dependency, statistical setting, or cloud
execution change is required. Existing GitHub image CI must exercise the new
native imputation regression. Do not build a Docker image locally or submit a
Terra job. Do not claim a speed gain without measurement.
