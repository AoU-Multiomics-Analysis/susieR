#!/usr/bin/env bash
set -euo pipefail

wdl="workflows/CheckpointedWindowSusie.wdl"
inputs="examples/inputs/CheckpointedWindowSusie.inputs.json"

miniwdl check "$wdl"
rg -q '^workflow CheckpointedWindowSusieWorkflow' "$wdl"
rg -q '^task RunCheckpointedWindowSusie' "$wdl"
rg -q 'preemptible: preemptible_attempts' "$wdl"
rg -q 'Rscript /opt/r/scripts/run_checkpointed_window_susie.R' "$wdl"
scatter_count="$(rg -c 'scatter[[:space:]]*\(' "$wdl" || true)"
test "${scatter_count:-0}" -eq 0

for output_name in \
  WindowManifest \
  WindowFitIndex \
  SusieParquet \
  SusielbfParquet \
  FullSusieParquet; do
  rg -q "File $output_name" "$wdl"
done

for message in \
  'START window' \
  'VALIDATE checkpoint root' \
  'RUN runner' \
  'CHECK outputs' \
  'COMPLETE window'; do
  rg -Fq "$message" "$wdl"
done

for argument in \
  --window-id \
  --window-dosage \
  --window-phenotypes \
  --phenotype-data \
  --covariate-files \
  --covariate-modalities \
  --keep-samples \
  --checkpoint-root \
  --output-dir; do
  rg -Fq -- "$argument" "$wdl"
done

rg -Fq '~{sep="," covariate_files}' "$wdl"
rg -Fq '~{sep="," covariate_modalities}' "$wdl"
rg -q 'defined\(keep_samples\)' "$wdl"
rg -Fq 'checkpoint_root="~{checkpoint_root}"' "$wdl"
rg -Fq 'case "$checkpoint_root" in' "$wdl"
rg -Fq 'gs://*)' "$wdl"

rg -Fq 'primaryDescriptorPath: /workflows/CheckpointedWindowSusie.wdl' .dockstore.yml
rg -Fq 'name: "CheckpointedWindowSusie"' .dockstore.yml
rg -Fq 'miniwdl check CheckpointedWindowSusie.wdl' .github/workflows/wdl-validation.yml

for key in \
  CheckpointedWindowSusieWorkflow.window_id \
  CheckpointedWindowSusieWorkflow.window_dosage \
  CheckpointedWindowSusieWorkflow.window_phenotypes \
  CheckpointedWindowSusieWorkflow.phenotype_data \
  CheckpointedWindowSusieWorkflow.covariate_files \
  CheckpointedWindowSusieWorkflow.covariate_modalities \
  CheckpointedWindowSusieWorkflow.keep_samples \
  CheckpointedWindowSusieWorkflow.checkpoint_root \
  CheckpointedWindowSusieWorkflow.runner_image \
  CheckpointedWindowSusieWorkflow.memory_gb \
  CheckpointedWindowSusieWorkflow.cpu \
  CheckpointedWindowSusieWorkflow.disk_gb \
  CheckpointedWindowSusieWorkflow.preemptible_attempts; do
  rg -Fq "\"$key\"" "$inputs"
done

rg -Fq '"CheckpointedWindowSusieWorkflow.preemptible_attempts": 3' "$inputs"
