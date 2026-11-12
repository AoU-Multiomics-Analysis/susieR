# Cohort MAF in compact SuSiE outputs

The user selected the samples in each fit and requested `cohort_maf` only in
the compact SuSiE Parquets. Compute a numeric value for each physical predictor
before genotype centering and covariate adjustment. For diploid ALT dosage,
`p = mean(dosage) / 2` and `cohort_maf = min(p, 1 - p)`.

Use the final fit sample subset, including any additional-genotype intersection.
The ordinary path records the mean-imputed dosage matrix. Mean imputation
preserves the observed mean on the same samples and dosage representation.
Retain unknown MAF for entirely missing legacy rows and for values outside the
diploid dosage range. Do not infer MAF from adjusted genotypes, fitted effects,
genotype SD, ancestry summaries, or QC counts.

Store MAF vectors on the fit and prepared genotype caches by physical predictor
position. Preserve this alignment through window selection and cache eviction.
Repeated variant labels must not cause a join or row expansion. Add the vector
to compact output only. Per-gene, batch, and combined compact files carry the
column. Old fits or compact inputs with no recorded MAF receive numeric NA.
Full, LBF, and credible-set tables keep their current schemas.

This is metadata only. Keep genotypes, sample selection, model arguments, PIPs,
CSs, moments, RNG state, aggregation behavior, and WDL inputs unchanged. Use
small tests and the GitHub image checks. Do not run a full analysis, build a
local Docker image, or submit a Terra job.
