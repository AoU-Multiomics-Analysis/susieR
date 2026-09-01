#!/usr/bin/env bash
set -euo pipefail

wdl="workflows/CheckpointedWindowSusie.wdl"
inputs="examples/inputs/CheckpointedWindowSusie.inputs.json"
documentation=(
  README.md
  docs/README.md
  docs/inputs.md
  docs/scripts.md
  docs/workflows.md
  docs/docker.md
)

require_count() {
  local pattern="$1"
  local expected="$2"
  local actual
  actual="$(rg -c -- "$pattern" "$wdl" || true)"
  actual="${actual:-0}"
  if [ "$actual" -ne "$expected" ]; then
    printf 'Expected %s matches for %s; found %s.\n' \
      "$expected" "$pattern" "$actual" >&2
    exit 1
  fi
}

miniwdl check "$wdl"
require_count '^workflow ' 1
require_count '^task ' 1
require_count '^[[:space:]]*call RunCheckpointedWindowSusie[[:space:]]*\{' 1
require_count 'Rscript /opt/r/scripts/run_checkpointed_window_susie.R' 1
require_count 'scatter[[:space:]]*\(' 0
require_count '^[[:space:]]*test -s window_outputs/' 5

checked_output_paths="$(
  awk '/^[[:space:]]*test -s window_outputs\// { print $3 }' "$wdl" | sort
)"
expected_checked_output_paths=$'window_outputs/credible_sets.parquet\nwindow_outputs/full_susie.parquet\nwindow_outputs/lbf_variable.parquet\nwindow_outputs/window_fit_index.tsv\nwindow_outputs/window_manifest.json'
if [ "$checked_output_paths" != "$expected_checked_output_paths" ]; then
  printf 'Checked output paths do not match the required set.\n' >&2
  exit 1
fi

rg -q '^workflow CheckpointedWindowSusieWorkflow' "$wdl"
rg -q '^task RunCheckpointedWindowSusie' "$wdl"
rg -q 'Rscript /opt/r/scripts/run_checkpointed_window_susie.R' "$wdl"
default_runner_image='ghcr.io/aou-multiomics-analysis/susier/checkpointed-window@sha256:6d390ae3e186f2500abf8585e9ab16261ad7f7651ca7dc82e57497eed71aa228'
require_count "^[[:space:]]+String runner_image = \"${default_runner_image}\"$" 1
require_count '^[[:space:]]+String runner_image$' 1

workflow_outputs="$(
  awk '
    /^workflow CheckpointedWindowSusieWorkflow[[:space:]]*\{/ { in_workflow = 1; next }
    /^task RunCheckpointedWindowSusie[[:space:]]*\{/ { exit }
    in_workflow && /^  output \{/ { in_output = 1; next }
    in_output && /^  \}/ { exit }
    in_output && /^[[:space:]]+File[[:space:]]/ { print $2 }
  ' "$wdl"
)"
expected_workflow_outputs=$'WindowManifest\nWindowFitIndex\nSusieParquet\nSusielbfParquet\nFullSusieParquet'
if [ "$workflow_outputs" != "$expected_workflow_outputs" ]; then
  printf 'Workflow outputs do not match the required ordered set.\n' >&2
  exit 1
fi

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
  require_count "$argument" 1
done

for runtime_input in \
  runner_image \
  memory_gb \
  cpu \
  disk_gb \
  preemptible_attempts; do
  rg -q "^[[:space:]]+${runtime_input}[[:space:]]=[[:space:]]${runtime_input},?\$" "$wdl"
  rg -q "^[[:space:]]+(String|Int)[[:space:]]+${runtime_input}\$" "$wdl"
done

rg -q '^[[:space:]]*docker: runner_image$' "$wdl"
rg -Fq 'memory: "${memory_gb} GB"' "$wdl"
rg -q '^[[:space:]]*cpu: cpu$' "$wdl"
rg -Fq 'disks: "local-disk ${disk_gb} SSD"' "$wdl"
rg -q '^[[:space:]]*preemptible: preemptible_attempts$' "$wdl"

rg -Fq 'File window_id_values = write_lines([window_id])' "$wdl"
rg -Fq 'File checkpoint_root_values = write_lines([checkpoint_root])' "$wdl"
rg -Fq 'File covariate_files_values = write_lines(covariate_files)' "$wdl"
rg -Fq 'File covariate_modalities_values = write_lines(covariate_modalities)' "$wdl"
rg -Fq 'File keep_samples_values = write_lines(select_all([keep_samples]))' "$wdl"
rg -Fq 'if IFS= read -r keep_samples < "~{keep_samples_values}" && [[ -n "$keep_samples" ]]; then' "$wdl"
rg -Fq 'File runner_image_values = write_lines([runner_image])' "$wdl"
rg -Fq 'CHECKPOINTED_SUSIE_RUNTIME_IMAGE="$runner_image"' "$wdl"
rg -Fq 'runner image must end with @sha256 and 64 lowercase hex characters' "$wdl"

command_body="$(awk '/^[[:space:]]*command <<</, /^[[:space:]]*>>>/' "$wdl")"
if printf '%s\n' "$command_body" | rg -q '~\{(window_id|checkpoint_root|covariate_modalities)\}'; then
  printf 'Raw user String interpolation is not allowed in the Bash command.\n' >&2
  exit 1
fi
if rg -n \
  '^[[:space:]]*[[:alpha:]_][[:alnum:]_]*=.*~\{(window_id|checkpoint_root|covariate_modalities)\}' \
  "$wdl"; then
  printf 'Raw user String interpolation is not allowed in Bash assignments.\n' >&2
  exit 1
fi
if rg -n \
  '^[[:space:]]*(echo|printf)[^\n]*~\{(window_id|checkpoint_root|covariate_modalities)\}' \
  "$wdl"; then
  printf 'Raw user String interpolation is not allowed in Bash log messages.\n' >&2
  exit 1
fi

guard_line="$(rg -n -F 'case "$checkpoint_root" in' "$wdl" | cut -d: -f1)"
runner_line="$(rg -n -F 'Rscript /opt/r/scripts/run_checkpointed_window_susie.R' "$wdl" | cut -d: -f1)"
test "$guard_line" -lt "$runner_line"

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
  CheckpointedWindowSusieWorkflow.memory_gb \
  CheckpointedWindowSusieWorkflow.cpu \
  CheckpointedWindowSusieWorkflow.disk_gb \
  CheckpointedWindowSusieWorkflow.preemptible_attempts; do
  rg -Fq "\"$key\"" "$inputs"
done

rg -Fq '"CheckpointedWindowSusieWorkflow.preemptible_attempts": 3' "$inputs"
if rg -Fq '"CheckpointedWindowSusieWorkflow.runner_image"' "$inputs"; then
  printf 'The example must demonstrate the default runner image.\n' >&2
  exit 1
fi

documentation_text="$(cat "${documentation[@]}")"
for term in \
  CheckpointedWindowSusieWorkflow \
  window_phenotypes.tsv \
  p_value \
  window_manifest.json \
  window_fit_index.tsv \
  NONCONVERGED \
  SKIPPED \
  'one workflow per prepared window'; do
  if ! printf '%s\n' "$documentation_text" | rg -Fq "$term"; then
    printf 'Documentation is missing required term: %s\n' "$term" >&2
    exit 1
  fi
done

rg -Fq '| [`workflows/CheckpointedWindowSusie.wdl`](workflows/CheckpointedWindowSusie.wdl) | `CheckpointedWindowSusieWorkflow` |' README.md
rg -Fq '| `containers/CheckpointedWindowSusie/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/checkpointed-window` |' docs/docker.md
