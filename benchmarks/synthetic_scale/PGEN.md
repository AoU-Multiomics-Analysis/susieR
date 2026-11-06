# Synthetic PGEN loading test

This test uses the existing one-gene fixture: 10,000 samples, 20,000 diploid
biallelic variants, 20 covariates plus an intercept, and three causal effects.
It uses hardcalls, not fractional dosages. The R reader is `pgenlibr` 0.7.0.

`prepare_pgen_inputs.R` reads the original Tabix fixture. It writes a matching
Tabix TSV and PGEN/PVAR/PSAM files. The PGEN has two chr2 variants before the
20,000 chr1 variants and two chr2 variants after them. These extra records test
that the reader uses indices from the complete PGEN. The chr1 query is the same
for both formats. File creation is outside all timed runs.

The original fixture uses both `-1` and `NA` for missing values. PGEN has one
missing hardcall state. Both new files use `NA`. The preparation report records
whether this changes the existing nonzero-variant filter, which runs before
imputation. This matters for rare or all-missing variants. The builder checks
every PGEN chr1 hardcall against the source values after this conversion.

The recorded reader test used three fresh R processes per format. Each pair
ran Tabix first and PGEN second. All calls used one thread. Files remained in
the filesystem cache; these were local
warm-cache tests. Both workers load the same packages before the reader clock
starts. The timed call includes metadata, region selection, file opening,
reading, sample ordering, transpose, and conversion to the common double matrix
with variants in rows. Peak RSS covers the whole worker, including package
loading and checks after reading. It is not a measurement of the reader alone.

An extra integer run measures the PGEN adapter without conversion to double.
It still includes transpose and metadata. The pipeline comparison uses the
double result to keep its input matrix identical to the Tabix result.

`pgen_pipeline_worker.R` runs copies of the current entry point and
helpers. It changes only the reader calls for the PGEN case, and adds the same
CLI and timing hooks to both cases. All preprocessing and model settings are
the same. The PGEN wrapper clips a query start below 1 to 1, since PVAR positions
are positive. This selects the same records as the Tabix query. This test uses
seed 20261106, OpenBLAS, and one thread. The verifier
compares the complete saved fits, decoded Parquet tables, variant summaries,
and random-number states. It can also compare the earlier production result.

The adapter is benchmark code. These tests do not add PGEN inputs to the
production CLI, image, or WDL. They do not measure chromosome-scale metadata,
cloud file localization, fractional dosages, phased calls, or other ploidy.

## Results

Small text reports from the 2026-11-06 UTC runs are published in
[`reports/2026-11-06/pgen/`](reports/2026-11-06/pgen/).
The [report notes](reports/2026-11-06/README.md) state their provenance.
These measurements precede the later guard for variants with almost zero
residual variance. New commands below run the current source.
The [final source check](reports/final-pr/README.md) includes the guard and
diagnostic outputs, with the original Tabix fixture and matching PGEN hardcalls.

| Measurement | Tabix TSV | PGEN |
|---|---:|---:|
| Region reader, median of 3 runs | 23.677 s | 1.661 s |
| Reader worker peak RSS, median of 3 runs | 5.315 GiB | 4.167 GiB |
| Complete R entry point, 1 run each | 72.114 s | 49.694 s |
| Complete pipeline peak RSS | 7.382 GiB | 7.592 GiB |

PGEN region loading was 14.25 times faster. It reduced complete entry point
time by 31.09%. Full pipeline peak memory did not improve in this pair of runs.
The later matrix operations still allocate double matrices.

All 200 million chr1 hardcalls passed exact roundtrip checks. All 7 reader runs
had identical block checks, including sample and variant names. The 21 pipeline
checks passed: complete saved fits, PIPs, credible sets, variant order,
convergence, decoded output tables, position summaries, environment fields,
and random-number states. Maximum PIP difference was zero. Both runs converged
in 4 iterations with 3 credible sets. The PGEN fit also exactly matched the
earlier production fit from the original mixed-missing-value fixture.

The source had 100,000 `NA` cells and 100,000 `-1` cells. Conversion to one
missing-value convention changed zero filter decisions in this fixture. This
does not establish that conversion preserves filter decisions for other data;
the small tests include rare/all-missing rows where the existing filter differs.

The extra integer adapter run took 1.318 s. Its matrix used 0.747 GiB, compared
with 1.492 GiB for double output. Its whole-worker peak RSS was 2.590 GiB.
The pipeline comparison above uses double output in both cases.

The complete local run files and generated inputs are excluded from version
control. Use the generator to recreate inputs. The commands below do not need
saved fitted objects or files from an earlier run.

## Run

Start in the repository root. Set `BENCH_LIBRARY` to an existing R package
library with the [benchmark dependencies](README.md#dependencies).
Install `pgenlibr` 0.7.0 in that library. Use `eQTLUtils` revision
`26242562a4e244334fd9691d03bc1ef4d2d6c1d9`.

The commands create a fresh fixture. New reader repetitions alternate the
backend order. This differs from the recorded measurements, which always
used Tabix first. All runs keep the files in the filesystem cache.

```sh
set -eu
bench_repo="$(git rev-parse --show-toplevel)"
bench_library="${BENCH_LIBRARY:?Set BENCH_LIBRARY to an existing R package library}"
export BENCH_LIBRARY="$bench_library"
export R_LIBS_USER="$bench_library"
export VROOM_THREADS=1
bench_root="$(mktemp -d "${TMPDIR:-/tmp}/susier-pgen.XXXXXX")"
bench_inputs="$bench_root/original-inputs"
bench_output="$bench_root/results"
mkdir -p "$bench_output"
bench_blas="${BENCH_BLAS:-reference}"

Rscript "$bench_repo/benchmarks/synthetic_scale/generate_inputs.R" \
  "$bench_inputs" 10000 20000

Rscript \
  "$bench_repo/benchmarks/synthetic_scale/prepare_pgen_inputs.R" \
  "$bench_inputs" "$bench_output/inputs"

for repetition in 1 2 3; do
  if [ "$((repetition % 2))" -eq 1 ]; then
    set -- tabix_read pgen_read
  else
    set -- pgen_read tabix_read
  fi
  for backend in "$@"; do
    python3 "$bench_repo/benchmarks/synthetic_scale/run_process.py" \
      --worker "$bench_repo/benchmarks/synthetic_scale/pgen_read_worker.R" \
      --inputs "$bench_output/inputs" \
      --output "$bench_output/${backend}_${repetition}" --mode "$backend" \
      --library "$bench_library" --repo "$bench_repo" \
      --samples 10000 --variants 20000 --blas "$bench_blas" --threads 1
  done
done

python3 "$bench_repo/benchmarks/synthetic_scale/run_process.py" \
  --worker "$bench_repo/benchmarks/synthetic_scale/pgen_read_worker.R" \
  --inputs "$bench_output/inputs" --output "$bench_output/pgen_integer" \
  --mode pgen_integer --library "$bench_library" --repo "$bench_repo" \
  --samples 10000 --variants 20000 --blas "$bench_blas" --threads 1

Rscript "$bench_repo/benchmarks/synthetic_scale/verify_pgen_reads.R" "$bench_output"

for backend in tabix_pipeline pgen_pipeline; do
  python3 "$bench_repo/benchmarks/synthetic_scale/run_process.py" \
    --worker "$bench_repo/benchmarks/synthetic_scale/pgen_pipeline_worker.R" \
    --inputs "$bench_output/inputs" --output "$bench_output/$backend" \
    --mode "$backend" --library "$bench_library" --repo "$bench_repo" \
    --samples 10000 --variants 20000 --blas "$bench_blas" --threads 1
done

Rscript "$bench_repo/benchmarks/synthetic_scale/verify_pgen_pipeline.R" \
  "$bench_output"
```

For an optional comparison with a production run you created separately, pass
its directory as the verifier's second argument:

```sh
Rscript "$bench_repo/benchmarks/synthetic_scale/verify_pgen_pipeline.R" \
  "$bench_output" "$legacy_production_dir"
```

That directory must contain `synthetic_susie.rds`. This additional comparison
checks exact equality and can fail if the production source differs.

The default BLAS option uses the installed R backend. To use the measured
macOS override, set `BENCH_BLAS=openblas` before the setup commands. It requires
`/opt/homebrew/opt/openblas/lib/libopenblas.dylib` and applies only to the worker
process. The runner saves source copies, logs, stage times, and process resource
use for each run. The complete workflow has not been tested on Terra for this
change.
