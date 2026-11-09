#!/usr/bin/env bash
# Run in the GitHub-built image with only tests mounted, using installed code.
set -euo pipefail
full_checks=false
if (( $# == 1 )) && [[ "$1" == "--full" ]]; then
  full_checks=true
elif (( $# != 0 )); then
  echo "Usage: bash tests/test_pgen_container.sh [--full]" >&2
  exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 VROOM_THREADS=1
echo 'CHECK audited SuSiE source enables fit optimizations'
Rscript - <<'RS'
helpers <- new.env(parent = globalenv())
sys.source(file.path(Sys.getenv("SUSIER_FUNCTIONS_PATH"), "SusiePerformance.R"),
           envir = helpers)
stopifnot("Installed PGEN image must enable audited SuSiE fit optimizations" =
            isTRUE(helpers$.susieFitSourceSupported(asNamespace("susieR"))))
RS
echo 'CHECK installed PGEN tools'
plink2 --version
bgzip --version
tabix --version
echo 'TEST task CPU detection and thread settings'
Rscript tests/test_pgen_runtime.R .
echo 'TEST optional MAC with real PLINK extraction'
Rscript tests/test_pgen_mac.R .
echo 'TEST shared PGEN MAC and AN call-rate QC'
Rscript tests/test_pgen_qc.R .
echo 'SMOKE one PLINK extraction, one SuSiE fit, and aggregation'
Rscript tests/test_pgen_smoke.R .
if [[ "$full_checks" != "true" ]]; then
  echo 'COMPLETE quick installed PGEN smoke test'
  exit 0
fi
echo 'FULL installed-code regression and parity checks'
for script in plan_pgen_susie.R filter_pgen_susie.R run_pgen_susie_batch.R aggregate_pgen_susie.R; do
  Rscript "${SUSIER_SCRIPTS_PATH}/${script}" --help
done
echo 'TEST real PGEN reader and per-gene parity'
Rscript tests/run_pgen_reader_tests.R .
Rscript tests/run_pgen_gene_tests.R .
echo 'TEST planning and strict aggregation'
Rscript tests/run_pgen_batch_tests.R .
echo 'TEST existing statistical and output behavior'
Rscript tests/run_tests.R .
echo 'TEST PLINK extraction through installed batch and aggregate CLIs'
Rscript tests/test_pgen_pipeline.R .
echo 'COMPLETE full installed PGEN regression tests'
