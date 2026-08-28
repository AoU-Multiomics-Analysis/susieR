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
    String runner_image = "ghcr.io/aou-multiomics-analysis/susier/checkpointed-window:main"
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

  command <<<
    set -euo pipefail

    checkpoint_root="~{checkpoint_root}"
    echo "START window ~{window_id}"
    echo "VALIDATE checkpoint root ${checkpoint_root}"
    case "$checkpoint_root" in
      gs://*) ;;
      *)
        echo "ERROR checkpoint root must start with gs://" >&2
        exit 1
        ;;
    esac

    mkdir -p window_outputs

    keep_samples_args=()
    if ~{if defined(keep_samples) then "true" else "false"}; then
      keep_samples_args=(--keep-samples "~{keep_samples}")
    fi

    echo "RUN runner ~{window_id}"
    Rscript /opt/r/scripts/run_checkpointed_window_susie.R \
      --window-id "~{window_id}" \
      --window-dosage "~{window_dosage}" \
      --window-phenotypes "~{window_phenotypes}" \
      --phenotype-data "~{phenotype_data}" \
      --covariate-files "~{sep="," covariate_files}" \
      --covariate-modalities "~{sep="," covariate_modalities}" \
      "${keep_samples_args[@]}" \
      --checkpoint-root "$checkpoint_root" \
      --output-dir window_outputs

    echo "CHECK outputs ~{window_id}"
    test -s window_outputs/window_manifest.json
    test -s window_outputs/window_fit_index.tsv
    test -s window_outputs/credible_sets.parquet
    test -s window_outputs/lbf_variable.parquet
    test -s window_outputs/full_susie.parquet
    echo "COMPLETE window ~{window_id}"
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
