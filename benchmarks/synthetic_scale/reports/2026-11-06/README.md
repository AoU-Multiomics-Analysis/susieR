# Recorded local benchmark

These text reports were copied from local runs collected on 2026-11-06 UTC
(2026-11-05 in America/Los_Angeles). The complete local records are in the
ignored `results/10000x20000/` directory. That directory contains generated
fits, decoded output tables, source copies, logs, and process records.

Each pipeline run used one synthetic gene, 10,000 samples, 20,000 variants,
20 covariates plus an intercept, three causal effects, and 10 SuSiE effects.
The machine used R 4.3.1, OpenBLAS 0.3.23, and one thread. See
`software_versions.tsv` for selected versions and seeds. `input_sha256.tsv`
records hashes of the original generated inputs; it omits local path names.

The original source was commit
`dd9b8efc2c0f0a01b3456a9a25cc5ec4a6797267`. The historical worker created
the `candidate` and `factored` source copies from that commit. The `production`
and PGEN runs used locally modified source after the performance changes.
They ran before the later production guard for variants with almost zero
residual variance. An additional covariate-header fix was made after the
production source copy. These reports do not establish performance or complete
workflow validation for the final source revision.

SHA-256 hashes identify the saved production helper copies used in this run:

```text
000686e90fccd6cade1f33b6b2a9e2dd8a099dde5dcd1a1e475e5de622e8fb54  SusieFunctions.R
a2936271ba6b47292ec14bd40b193353aa2a360df75fb3edd1b5c1d57ed732f8  ImportFunctions.R
64d56163f6dbb2cdf1c4e279352842c9e05e3e4cb767cb85593d27f4b67d3207  SusiePerformance.R
```

| Published files | Original local records | Purpose |
|---|---|---|
| `benchmark_results.tsv`, `numerical_checks.tsv` | Same names at the run root | Dense baseline and historical candidates |
| `production_results.tsv`, `production_checks.tsv` | Same names at the run root | Earlier production run and output checks |
| `baseline_stages.tsv`, `baseline_resource_usage.tsv` | `baseline/stages.tsv`, `baseline/resource_usage.tsv` | Original timing and whole-process resource use |
| `production_stages.tsv`, `production_resource_usage.tsv` | `production/stages.tsv`, `production/resource_usage.tsv` | Production timing and whole-process resource use |
| `covariate_parser_checks.tsv`, `fast_z_small_checks.tsv`, `fast_z_model_checks.tsv` | Same names at the run root | Historical prototype checks, separate from final regression tests |
| `pgen/*` | Same names under `pgen/` | Matched reader/pipeline tests and preparation checks |

The dense baseline entry point took 210.185 seconds, including 25.570 seconds
of matrix fingerprint work. Its adjusted time was 184.615 seconds. Production
took 75.227 seconds with no fingerprint work: a 2.45-fold speed increase.
Whole-process peak RSS fell from 8.623 to 6.753 GiB in these individual runs.
The production fit had the same credible-set membership as the original
baseline. Its maximum PIP difference from the baseline was `1.40e-10`.

The matched Tabix and PGEN entry points took 72.114 and 49.694 seconds. This
pair gives a 31.09% time reduction. PGEN gave a 3.72-fold speed increase relative
to the original adjusted dense baseline. It is a separate benchmark prototype.
Its whole-pipeline peak RSS did not improve in the matched pair. Three reader
repetitions ran Tabix first, then PGEN. They used files in the filesystem cache.
Pipeline times are single runs; they are not repeated runtime estimates.

Values in the copied TSV files are unchanged. The PGEN preparation report's
source directory and the optional legacy comparison's directory were replaced
with descriptive placeholders. The software table records versions from the
saved environment/session records. These edits omit absolute user and temporary
paths. Empty check-detail fields are stored as `NA` to keep explicit TSV cells.
No input matrices, fitted objects, binary files, logs, or symlinks are
published here.

Use the [benchmark commands](../../README.md) to generate new inputs and outputs.
Use the [PGEN commands](../../PGEN.md) for the reader prototype. Neither command
sequence requires a saved fitted object from this repository. The reports cover
local execution. The complete workflow has not been tested on Terra for this
change.
