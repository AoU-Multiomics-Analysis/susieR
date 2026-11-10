version 1.0

# Aggregate saved full-result files without genotype preparation or model fits.
workflow AggregateFullSusie {
  input {
    File outcome_index
    Array[File] full_susie_files
    String docker_image = "ghcr.io/aou-multiomics-analysis/susier/pgen:main"
    Int memory_gb = 8
    Int disk_gb = 50
  }

  call AggregateFullSusieFiles {
    input:
      outcome_index = outcome_index,
      full_susie_files = full_susie_files,
      docker_image = docker_image,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    File FullSuSiEParquet = AggregateFullSusieFiles.FullSuSiEParquet
    File TaskStdout = AggregateFullSusieFiles.TaskStdout
    File TaskStderr = AggregateFullSusieFiles.TaskStderr
  }
}

task AggregateFullSusieFiles {
  input {
    File outcome_index
    Array[File] full_susie_files
    String docker_image
    Int memory_gb
    Int disk_gb
  }

  command <<<
    set -euo pipefail
    echo 'START saved full-result aggregation'
    require_readable() {
      local path
      for path in "$@"; do
        case "$path" in
          *://*) printf 'ERROR localization: unresolved URI: %s\n' "$path" >&2; exit 1 ;;
        esac
        if [[ -z "$path" || "$path" == *$'\n'* || "$path" == *$'\r'* || ! -f "$path" || ! -r "$path" ]]; then
          printf 'ERROR input must be a readable local file: %s\n' "$path" >&2
          exit 1
        fi
      done
    }
    require_readable '~{sub(outcome_index, "'", "'\"'\"'")}'
    # Preserve File types until command-time localization and file-list creation.
    cp '~{write_lines(full_susie_files)}' full_susie_files.txt
    while IFS= read -r path; do require_readable "$path"; done < full_susie_files.txt
    Rscript /opt/r/scripts/aggregate_pgen_susie.R \
      --full-only \
      --outcome-index '~{sub(outcome_index, "'", "'\"'\"'")}' \
      --full-susie-files full_susie_files.txt \
      --output-dir aggregate
    echo 'COMPLETE saved full-result aggregation'
  >>>

  output {
    File FullSuSiEParquet = "aggregate/full_susie.parquet"
    File TaskStdout = stdout()
    File TaskStderr = stderr()
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}
