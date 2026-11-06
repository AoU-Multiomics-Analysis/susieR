# Real IKZF1 cis benchmark

The optimized pipeline produced the same scientific results as the native SuSiE baseline in four local runs. This test used one IKZF1 expression outcome, 8,895 samples, and 155 expression covariates plus an intercept. The normal genotype filter retained 23,376 of the 23,393 supplied cis variants. Both versions used the same inputs, sample alignment, covariate adjustment, imputation, model settings, and random seed.

The baseline is commit `f0f7bcf3b9a35e758ccef96934d376ee31c5d9e4`. It uses native SuSiE fitting with the existing fast marginal-score helper. The optimized version adds the fit caches, shared marginal statistics, and shared output extraction in this change. The optimizer and convergence settings are unchanged. This comparison measures these changes; it does not measure the total gain from the earliest dosage-based pipeline.

## Results and timing

The run order was baseline, optimized, optimized, baseline. Each run used a fresh R process, R 4.3.1, susieR 0.16.6, and one OpenBLAS thread on macOS. The table gives the mean of two runs per version.

| Measurement | Baseline | Optimized | Runtime or memory reduction |
| --- | ---: | ---: | ---: |
| Fit and marginal z-scores | 98.032 s | 53.351 s | 45.6% |
| Marginal z-scores only | 2.671 s | 0.330 s | 87.7% |
| Gene step, including input verification | 121.903 s | 78.233 s | 35.8% |
| Complete worker | 130.904 s | 87.334 s | 33.3% |
| Worker peak RSS | 7.147 GiB | 6.449 GiB | 9.8% |

The gene step was 1.56 times faster. The fit and marginal scores were 1.84 times faster. SuSiE fitting accounted for about 68% of the optimized gene time. Each fit converged in 15 iterations and returned five credible sets. Both versions used the fast covariate adjustment; neither used its dense fallback. No variants were flagged as covariate-explained.

The fit timer starts after the benchmark hashes the adjusted inputs. The gene timer includes PGEN window reading, preparation, fitting, marginal scores, output writing, artifact hashing, and the same input-verification trace. It excludes experiment construction, reader opening, batch PLINK extraction, and R startup. The complete-worker timer includes setup and startup. Peak RSS covers the complete worker. These measurements do not measure cloud runtime.

## Exact checks

All 150 final checks passed. All seven scientific artifacts were byte-for-byte identical across all four runs: the full RDS, three Parquet files, covariate-variant diagnostic TSV, variant-position TSV, and outcome-index TSV. Strict comparisons also checked the complete decoded RDS, three decoded Parquet tables, ordered fit inputs, block hashes of the full adjusted genotype matrix and response, and random-number states before and after fitting. No model fields were removed before comparison.

A separate full input check compared the PGEN with all 208,080,735 original genotype cells. Sample order, variant positions and alleles, hardcalls, and all 132,310 missing values matched. No allele flip was applied.

The first pair had a benchmark observer error: it checked for the cached-statistics argument inside a tibble data mask. Those observations remain unmeasured. The second pair used a corrected observer and confirmed that cached input statistics were supplied to the optimized marginal calculation. This correction changed only the temporary benchmark worker. All scientific equality checks remained strict. Output extraction decreased from two calls to one per outcome.

Individual-level inputs and fitted objects are not included. The aggregate measurements are in [stage_timings.tsv](stage_timings.tsv) and [performance_summary.tsv](performance_summary.tsv). [source_sha256.tsv](source_sha256.tsv) identifies the four tested production R files. The benchmark did not build a local Docker image or submit a Terra job. The complete workflow has not been tested on Terra.
