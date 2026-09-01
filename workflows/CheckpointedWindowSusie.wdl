version 1.0

workflow CheckpointedWindowSusieWorkflow {
  input {
    String window_id
    File window_dosage
    File window_phenotypes
    File phenotype_data
    Array[File] covariate_files
    Array[String] covariate_modalities
    File? keep_samples
    String checkpoint_root
    String runner_image = "ghcr.io/aou-multiomics-analysis/susier/checkpointed-window@sha256:554a3af851a96bc22d5579e1b4c06f5f2d7ae90f583c53449c8289d228b588f0"
    Int memory_gb = 16
    Int cpu = 1
    Int disk_gb = 500
    Int preemptible_attempts = 3
  }

  call RunCheckpointedWindowSusie {
    input:
      window_id = window_id,
      window_dosage = window_dosage,
      window_phenotypes = window_phenotypes,
      phenotype_data = phenotype_data,
      covariate_files = covariate_files,
      covariate_modalities = covariate_modalities,
      keep_samples = keep_samples,
      checkpoint_root = checkpoint_root,
      runner_image = runner_image,
      memory_gb = memory_gb,
      cpu = cpu,
      disk_gb = disk_gb,
      preemptible_attempts = preemptible_attempts
  }

  output {
    File WindowManifest = RunCheckpointedWindowSusie.WindowManifest
    File WindowFitIndex = RunCheckpointedWindowSusie.WindowFitIndex
    File SusieParquet = RunCheckpointedWindowSusie.SusieParquet
    File SusielbfParquet = RunCheckpointedWindowSusie.SusielbfParquet
    File FullSusieParquet = RunCheckpointedWindowSusie.FullSusieParquet
  }
}

task RunCheckpointedWindowSusie {
  input {
    String window_id
    File window_dosage
    File window_phenotypes
    File phenotype_data
    Array[File] covariate_files
    Array[String] covariate_modalities
    File? keep_samples
    String checkpoint_root
    String runner_image
    Int memory_gb
    Int cpu
    Int disk_gb
    Int preemptible_attempts
  }

  File window_id_values = write_lines([window_id])
  File checkpoint_root_values = write_lines([checkpoint_root])
  File covariate_modalities_values = write_lines(covariate_modalities)
  File runner_image_values = write_lines([runner_image])

  command <<<
    set -euo pipefail

    IFS= read -r window_id < "~{window_id_values}"
    IFS= read -r checkpoint_root < "~{checkpoint_root_values}"
    IFS= read -r runner_image < "~{runner_image_values}"

    covariate_modalities=()
    while IFS= read -r covariate_modality; do
      covariate_modalities+=("$covariate_modality")
    done < "~{covariate_modalities_values}"

    comma_join() {
      local IFS=,
      printf '%s' "$*"
    }
    covariate_files_csv="~{sep=',' covariate_files}"
    covariate_modalities_csv="$(comma_join "${covariate_modalities[@]}")"

    keep_samples_args=()
    keep_samples_value="~{default="" keep_samples}"
    if [[ -n "$keep_samples_value" ]]; then
      keep_samples_args=(--keep-samples "$keep_samples_value")
    fi

    echo "START window ${window_id}"
    echo "VALIDATE checkpoint root ${checkpoint_root}"
    case "$checkpoint_root" in
      gs://*) ;;
      *)
        echo "ERROR checkpoint root must start with gs://" >&2
        exit 1
        ;;
    esac
    if [[ ! "$runner_image" =~ @sha256:[0-9a-f]{64}$ ]]; then
      echo "ERROR runner image must end with @sha256 and 64 lowercase hex characters" >&2
      exit 1
    fi
    export CHECKPOINTED_SUSIE_RUNTIME_IMAGE="$runner_image"

    mkdir -p window_outputs

    echo "RUN runner ${window_id}"
    Rscript /opt/r/scripts/run_checkpointed_window_susie.R \
      --window-id "$window_id" \
      --window-dosage "~{window_dosage}" \
      --window-phenotypes "~{window_phenotypes}" \
      --phenotype-data "~{phenotype_data}" \
      --covariate-files "$covariate_files_csv" \
      --covariate-modalities "$covariate_modalities_csv" \
      "${keep_samples_args[@]}" \
      --checkpoint-root "$checkpoint_root" \
      --output-dir window_outputs

    echo "CHECK outputs ${window_id}"
    test -s window_outputs/window_manifest.json
    test -s window_outputs/window_fit_index.tsv
    test -s window_outputs/credible_sets.parquet
    test -s window_outputs/lbf_variable.parquet
    test -s window_outputs/full_susie.parquet
    echo "COMPLETE window ${window_id}"
  >>>

  runtime {
    docker: runner_image
    memory: "${memory_gb} GB"
    cpu: cpu
    disks: "local-disk ${disk_gb} SSD"
    preemptible: preemptible_attempts
  }

  output {
    File WindowManifest = "window_outputs/window_manifest.json"
    File WindowFitIndex = "window_outputs/window_fit_index.tsv"
    File SusieParquet = "window_outputs/credible_sets.parquet"
    File SusielbfParquet = "window_outputs/lbf_variable.parquet"
    File FullSusieParquet = "window_outputs/full_susie.parquet"
  }
}
