# Cohort AN and MAC in compact SuSiE outputs

Add numeric `cohort_an` and `cohort_mac` beside `cohort_maf` in compact per-gene,
batch, and combined variant Parquets. Use the final samples in each fit, including
the additional-genotype sample intersection. For uncentered diploid ALT dosages,
exclude NA and the -1 missing code before imputation. AN is twice the observed
sample count. MAC is the smaller of the observed ALT dosage sum and AN minus that
sum. Fractional dosages can give a fractional MAC. All-missing retained rows have
AN=0 and MAC=0; nonmissing values outside [0,2] or nonfinite values give unknown
counts. Existing `cohort_maf` values and all scientific model inputs stay unchanged.

Use the observed genotype representation actually supplied to the fit, as for
`cohort_maf`: native-first/cached rows use stored dosages; unchanged fallback and
additional-feature paths use their hardcalls. A fallback row replaced by native
dosage uses raw native dosage counts. Counts from an unused stored representation
would describe different predictors from the existing MAF and fitted effects.

Record counts by physical predictor position. Carry them through ordinary,
reusable, native-first PGEN, and whole/rolling-cache paths. Do not join metadata
by variant label. Cache eviction and reordering must apply the same row masks.
Keep the current native mean-imputation call. Its API loses missingness, so an
opt-in bounded raw-dosage read supplies counts for newly loaded blocks. Cache hits
reuse those vectors. Avoid replacing dosage counts with hardcall counts.

Old fits and old compact files get numeric missing values for each absent count
independently. Empty compact files have numeric empty columns. Full variant, LBF,
credible-set, diagnostic and public extraction schemas keep their current fields.
Checkpointed constructor defaults stay unchanged. WDL inputs and cloud execution
stay unchanged. Use local R checks and existing GitHub image smoke tests; no local
Docker build or Terra job is needed.

Approach: extend the existing fit-sample MAF metadata path. Reconstructing AN or
MAC from an imputed matrix would count imputed alleles; hardcall-only counts would
lose observed fractional dosages. Neither supplies the requested observed counts.
