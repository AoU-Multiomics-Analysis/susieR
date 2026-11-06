#!/usr/bin/env bash
# Run in the GitHub-built image with only tests mounted, using installed code.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 VROOM_THREADS=1
echo 'CHECK installed PGEN tools and scripts'
plink2 --version
bgzip --version
tabix --version
for script in plan_pgen_susie.R run_pgen_susie_batch.R aggregate_pgen_susie.R; do
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
echo 'COMPLETE installed PGEN pipeline smoke tests'
