# Synthetic scale benchmark

This benchmark uses one gene, 10,000 samples, 20,000 variants, 20 covariates
plus an intercept, and three simulated causal effects. Input generation uses
seed 20261105. Model runs use seed 20261106 and 10 SuSiE effects.

The [historical reports](reports/2026-11-06/README.md) contain local measurements
from 2026-11-06 UTC. The production Tabix entry point took 75.227 seconds,
compared with 184.615 seconds for the original dense calculation after its
input checks were subtracted: a 2.45-fold speed increase. Peak process RSS was
6.753 GiB, compared with 8.623 GiB. These were individual runs on one machine.
The reports precede the later guard for variants with almost zero residual
variance. They do not validate the final source revision.

The [final source check](reports/final-pr/README.md) includes that guard and the
diagnostic output columns. It took 75.334 seconds with Tabix and 48.884 seconds
with the PGEN prototype. All previous output columns match the earlier
production run exactly. See that report for source hashes and memory use.

The separate [PGEN prototype](PGEN.md) took 49.694 seconds, a 3.72-fold speed
increase relative to that original baseline. PGEN is benchmark code. Production inputs
remain Tabix dosage files. These measurements do not establish cloud runtime.
The complete workflow has not been tested on Terra for this change.

## Dependencies

Use Python 3 and an R library with the dependencies of
`containers/SusieR/Dockerfile`. The benchmark also uses `digest` and `jsonlite`.
The R code needs `data.table`, `devtools`, `dplyr`, `tidyr`, `stringr`, `purrr`,
`readr`, `tibble`, `arrow`, `optparse`, `Rfast`, `matrixStats`, `caret`,
`Rsamtools`, `SummarizedExperiment`, `susieR`, and `eQTLUtils`.

The measured runs used R 4.3.1, susieR 0.16.6, OpenBLAS 0.3.23, and one
thread. Install `eQTLUtils` from revision
`26242562a4e244334fd9691d03bc1ef4d2d6c1d9`. The PGEN test also needs
`pgenlibr` 0.7.0. Selected measured versions are in
[`software_versions.tsv`](reports/2026-11-06/software_versions.tsv).
This version record is not a complete environment lock file. Package and BLAS
differences can change numerical results and runtime.

For an existing library, the two pinned R packages can be installed with:

```sh
export R_LIBS_USER="${BENCH_LIBRARY:?Set BENCH_LIBRARY to an existing R package library}"
Rscript -e 'remotes::install_version("susieR", version = "0.16.6", upgrade = "never")'
Rscript -e 'remotes::install_github("kauralasoo/eQTLUtils@26242562a4e244334fd9691d03bc1ef4d2d6c1d9", upgrade = "never")'
```

All benchmark runs are local. No image build or cloud submission is required.

## Generate inputs and run current production code

Start in the repository root. Set `BENCH_LIBRARY` to the R library that contains
the dependencies. The commands create all inputs and outputs in a new temporary
directory. Keep that directory if you need the full fitted objects.

```sh
set -eu
bench_repo="$(git rev-parse --show-toplevel)"
bench_library="${BENCH_LIBRARY:?Set BENCH_LIBRARY to an existing R package library}"
export BENCH_LIBRARY="$bench_library"
export R_LIBS_USER="$bench_library"
export VROOM_THREADS=1
bench_root="$(mktemp -d "${TMPDIR:-/tmp}/susier-scale.XXXXXX")"
bench_inputs="$bench_root/inputs"
bench_output="$bench_root/results"
mkdir -p "$bench_output"

# The reference option uses the BLAS backend of the installed R.
bench_blas="${BENCH_BLAS:-reference}"

Rscript "$bench_repo/benchmarks/synthetic_scale/generate_inputs.R" \
  "$bench_inputs" 10000 20000

python3 "$bench_repo/benchmarks/synthetic_scale/run_process.py" \
  --worker "$bench_repo/benchmarks/synthetic_scale/production_worker.R" \
  --inputs "$bench_inputs" --output "$bench_output/production" \
  --mode production --library "$bench_library" --repo "$bench_repo" \
  --samples 10000 --variants 20000 --blas "$bench_blas" --threads 1
```

For the measured macOS OpenBLAS override, set `BENCH_BLAS=openblas` before
these commands. This option requires
`/opt/homebrew/opt/openblas/lib/libopenblas.dylib`. It changes only the worker
process. The default `reference` option works with the installed R backend;
its name does not assert that R uses the reference BLAS implementation.

`production_worker.R` copies current production helpers and the entry point.
It adds CLI arguments and timing hooks. It does not replace calculations.
The runner saves source copies, logs, stage times, and process resource use.
Input generation and output verification are outside the measured run.

## Recreate the historical comparisons

The historical `worker.R` modes require production source from commit
`dd9b8efc2c0f0a01b3456a9a25cc5ec4a6797267`. They cannot patch current production
source. Create a separate source directory with `git archive`. This does not
change the current checkout or create a Git worktree.

Continue with the variables and generated inputs from the previous example:

```sh
bench_baseline="$bench_root/baseline-source"
mkdir -p "$bench_baseline"
git -C "$bench_repo" archive dd9b8efc2c0f0a01b3456a9a25cc5ec4a6797267 | \
  tar -x -C "$bench_baseline"

for mode in baseline candidate factored; do
  python3 "$bench_repo/benchmarks/synthetic_scale/run_process.py" \
    --worker "$bench_repo/benchmarks/synthetic_scale/worker.R" \
    --inputs "$bench_inputs" --output "$bench_output/$mode" \
    --mode "$mode" --library "$bench_library" --repo "$bench_baseline" \
    --samples 10000 --variants 20000 --blas "$bench_blas" --threads 1
done

Rscript "$bench_repo/benchmarks/synthetic_scale/compare_outputs.R" \
  --baseline "$bench_output/baseline" \
  --candidate "$bench_output/candidate" \
  --output "$bench_output/candidate_comparison"

Rscript "$bench_repo/benchmarks/synthetic_scale/verify_production.R" \
  "$bench_output"
```

These commands create the fitted objects needed by the verifiers. They do not
require a fitted object from this repository. `verify_production.R` compares
current production with the historical factored calculation and checks its
agreement with the original dense calculation. A failure needs investigation;
new production changes can make an exact comparison with an old fit invalid.

The `candidate` mode preserves dense residualization and tests temporary
memory/output changes. Its recorded outputs match the original exactly.
The `factored` mode changes floating-point operation order. Its recorded
maximum PIP difference from the original was `1.40e-10`, with the same
credible-set membership. Exact comparison of these two complete fits is
expected to fail. For a diagnostic report, run `compare_outputs.R` with the
factored directory as the candidate; it returns exit status 1 for differences.

`summarize_results.R` requires both comparison directories. The historical
`factored_fast` mode and the isolated `fit_worker.R` runs are retained as
prototype code. They require their earlier source/cache setup and are outside
the reproduction commands above.

## Stored files

`reports/2026-11-06/` contains small text copies of recorded metrics and checks.
Its README states their source and limits. Generated `results/`, `inputs/`, and
Python caches are ignored. Full local run files remain on disk. New runs do not
replace the published historical reports.
