# Final source check

These runs include the compatibility fallback and the covariate-explained
variant diagnostics. `source_sha256.tsv` identifies the tested production files.
Full run files remain in the ignored `results/10000x20000/final_pr/` directory.

Each run used one synthetic gene, 10,000 samples, 20,000 variants, 20 covariates
plus an intercept, three causal effects, and 10 SuSiE effects. R 4.3.1,
susieR 0.16.6, OpenBLAS 0.3.23, and one thread match the historical benchmark.
Input and model seeds are 20261105 and 20261106. Runs were sequential on the
local machine, with filesystem caching enabled. These are single pipeline
measurements, not repeated runtime estimates.

| Measurement | Original dense baseline | Final Tabix pipeline | Final PGEN prototype |
|---|---:|---:|---:|
| R entry point | 184.615 s | 75.334 s | 48.884 s |
| Whole-worker peak RSS | 8.623 GiB | 6.882 GiB | 6.670 GiB |

The original time excludes 25.570 seconds of benchmark-only matrix checks.
The final runs contain no such matrix-check step. The final Tabix pipeline is
2.45 times faster than the original. The PGEN prototype is 3.78 times faster.
PGEN remains benchmark code; production inputs remain Tabix dosage files.

The Tabix run used the original fixture with mixed `NA` and `-1` missing codes.
The PGEN run used the matching hardcalls with one missing-value convention.
The preparation check found zero changes to the variant filter and checked
all 200 million PGEN hardcalls against the source values. The phenotype,
covariate, sample, and model settings are the same.

All 22 checks in `pgen_checks.tsv` passed. The complete final fits and decoded
outputs match between the two backends, including diagnostic metadata and RNG
states. All pre-existing fitted fields match the previous production run
exactly. A separate check confirmed that all previous columns in all three
Parquet outputs match that run exactly; the primary table adds only the logical
credible-set flag.

`fit_summary.tsv` records zero PIP difference from the earlier production fit
and a maximum PIP difference of `1.40e-10` from the original dense baseline.
Both final runs retained the original credible-set membership and converged in
four iterations with three credible sets. This large fixture has no flagged
variants. The regression suite separately tests flagged variants, all fit
paths, per-phenotype lists, and credible-set flag propagation. Those fallback
cases preserve every previous native fitted field exactly.

The complete R regression suite, checkpointed-runner compatibility tests, R
lint, repository checks, and WDL static checks passed locally. The complete
workflow has not been tested on Terra for this change. No local Docker build
or cloud job was run. Image smoke tests run in GitHub Actions.

The copied check table uses `NA` for empty detail fields and a relative label
for the previous production run. Published files contain no local absolute
paths. See the [benchmark commands](../../README.md) for reproduction.
