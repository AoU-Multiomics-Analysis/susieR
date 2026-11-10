# Whole-cache sufficient statistics

Historical design: superseded on 2026-11-09 at the user's request. The active
pipeline uses individual-data SuSiE in every cache mode and does not reuse XtX.

Use the shared sufficient-statistics approach selected after the three-trait synthetic benchmark. Enable it automatically for eligible whole-batch PGEN caches. The cache must keep one sample order and covariate design. Do not change scientific model settings or output scales.

Considered approaches: calculate a new cross-product for every trait; share one cross-product for a whole batch; or maintain cross-products as a rolling window changes. Use the shared whole-batch approach. It amortizes construction and has a clear cache lifetime. Keep rolling and uncached fitting on the current individual-data path.

Require at least two nonempty windows, sample count at least twice the largest window, union squared no larger than the sum of window sizes squared, and enough cache memory for the shared matrix and temporary matrices. Preserve the existing genotype-cache budget and record the separate sufficient-statistics estimate. Build the matrix from cache block pairs without a second sample-by-union matrix. Map columns by physical PGEN index.

Use unscaled centered statistics. Pass the current model settings explicitly, including prior variance 0.1 and residual variance lower bound var(y)/10000. Use native credible-set construction with X after fitting so purity sampling and RNG behavior stay consistent. Preserve marginal scores, raw-scale posterior summaries, sample fitted values, and saved-fit fields. Dense numerical guards, legacy fits, unsupported source versions, and unavailable statistics use the current fitter.

Log eligibility, construction, reuse, and the fit method. Add small regression tests to the existing quick image checks. Compare native PIPs, posterior means, posterior SDs, marginal scores, CS members and purity, sample and variant order, convergence, RNG state, memory fallback, and cache cleanup. Do not run the full regression workflow or build Docker locally. No WDL or Terra changes are required.
