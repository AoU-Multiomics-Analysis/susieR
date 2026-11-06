version 1.0

# WDL 1.0 has no keys(Map). An ordered array keeps chromosome IDs and Files
# together. Planner indices refer to positions in this input array.
struct ChromosomePgen {
  String chromosome
  File pgen
  File pvar
  File psam
}

workflow SusieR2 {
  input {
    Array[ChromosomePgen]+ ChromosomePgens
    File PhenotypeBed
    File QtlResults
    File Covariates
    File SampleList
    String LegacyMissingCode
    File? GeneMap
    File? Ancestry
    File? VariantList
    File? AdditionalGenotypes
    Int BatchSize = 50
    Int CisDistance = 1000000
    Int BaseSeed = 1
    Boolean MemoryDebug = false
    Float MAF = 0.0
    Boolean SelectTopPhenotypePerCluster = false
    String TopPhenotypePvalueColumn = "qval"
    String DockerImage = "ghcr.io/aou-multiomics-analysis/susier/pgen:main"
    Int Cpu = 1
    Int BatchMemoryGB = 16
    Int BatchDiskGB = 100
    Int PlinkMemoryMB = 1024
    Int PlanMemoryGB = 8
    Int PlanDiskGB = 50
    Int AggregateMemoryGB = 8
    Int AggregateDiskGB = 50
    Int PreemptibleAttempts = 3
    Int MaxRetries = 1
  }

  scatter (chromosome_pgen in ChromosomePgens) {
    String chromosome = chromosome_pgen.chromosome
  }
  call PlanPgenSusie {
    input:
      chromosomes = chromosome,
      phenotype_bed = PhenotypeBed,
      qtl_results = QtlResults,
      gene_map = GeneMap,
      batch_size = BatchSize,
      cis_distance = CisDistance,
      select_top_phenotype_per_cluster = SelectTopPhenotypePerCluster,
      top_phenotype_pvalue_column = TopPhenotypePvalueColumn,
      docker_image = DockerImage,
      memory_gb = PlanMemoryGB,
      disk_gb = PlanDiskGB,
      preemptible_attempts = PreemptibleAttempts,
      max_retries = MaxRetries
  }

  scatter (i in range(length(PlanPgenSusie.Worklists))) {
    call RunPgenSusieBatch {
      input:
        pgen_files = ChromosomePgens[PlanPgenSusie.ChromosomeIndices[i]],
        worklist = PlanPgenSusie.Worklists[i],
        phenotype_bed = PlanPgenSusie.PhenotypeBeds[i],
        qtl_results = QtlResults,
        covariates = Covariates,
        sample_list = SampleList,
        missing_code = LegacyMissingCode,
        ancestry = Ancestry,
        variant_list = VariantList,
        additional_genotypes = AdditionalGenotypes,
        cis_distance = CisDistance,
        maf = MAF,
        base_seed = BaseSeed,
        memory_debug = MemoryDebug,
        select_top_phenotype_per_cluster = SelectTopPhenotypePerCluster,
        top_phenotype_pvalue_column = TopPhenotypePvalueColumn,
        plink_memory_mb = PlinkMemoryMB,
        docker_image = DockerImage,
        cpu = Cpu,
        memory_gb = BatchMemoryGB,
        disk_gb = BatchDiskGB,
        preemptible_attempts = PreemptibleAttempts,
        max_retries = MaxRetries
    }
  }

  call AggregatePgenSusie {
    input:
      expected_outcomes = PlanPgenSusie.ExpectedOutcomes,
      summary_files = RunPgenSusieBatch.Summary,
      credible_set_files = RunPgenSusieBatch.CredibleSets,
      index_files = RunPgenSusieBatch.OutcomeIndex,
      diagnostic_files = RunPgenSusieBatch.Diagnostics,
      docker_image = DockerImage,
      memory_gb = AggregateMemoryGB,
      disk_gb = AggregateDiskGB,
      preemptible_attempts = PreemptibleAttempts,
      max_retries = MaxRetries
  }

  output {
    File SusieParquet = AggregatePgenSusie.Summary
    File CredibleSetsParquet = AggregatePgenSusie.CredibleSets
    File OutcomeIndex = AggregatePgenSusie.OutcomeIndex
    File CovariateExplainedVariants = AggregatePgenSusie.Diagnostics
    File ExpectedOutcomes = PlanPgenSusie.ExpectedOutcomes
    File BatchIndex = PlanPgenSusie.BatchIndex
    Array[File] GeneRds = flatten(RunPgenSusieBatch.GeneRds)
    Array[File] GeneParquets = flatten(RunPgenSusieBatch.GeneParquets)
    Array[File] GeneTsvs = flatten(RunPgenSusieBatch.GeneTsvs)
    Array[File] GeneLogs = flatten(RunPgenSusieBatch.GeneLogs)
    Array[File] GeneProvenance = flatten(RunPgenSusieBatch.GeneProvenance)
    Array[File] BatchLogs = RunPgenSusieBatch.TaskStdout
    Array[File] BatchErrors = RunPgenSusieBatch.TaskStderr
    Array[File] PlinkLogs = flatten(RunPgenSusieBatch.PlinkLogs)
    File PlanLog = PlanPgenSusie.TaskStderr
    File AggregateLog = AggregatePgenSusie.TaskStderr
  }
}

task PlanPgenSusie {
  input {
    Array[String]+ chromosomes
    File phenotype_bed
    File qtl_results
    File? gene_map
    Int batch_size
    Int cis_distance
    Boolean select_top_phenotype_per_cluster
    String top_phenotype_pvalue_column
    String docker_image
    Int memory_gb
    Int disk_gb
    Int preemptible_attempts
    Int max_retries
  }
  command <<<
    set -euo pipefail
    echo 'START PGEN batch planning'
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
    require_readable '~{sub(phenotype_bed, "'", "'\"'\"'")}' '~{sub(qtl_results, "'", "'\"'\"'")}'
    args=(
      --phenotype-bed '~{sub(phenotype_bed, "'", "'\"'\"'")}'
      --qtl-results '~{sub(qtl_results, "'", "'\"'\"'")}'
      --chromosomes-file chromosomes.txt
      --batch-size '~{batch_size}'
      --cis-distance '~{cis_distance}'
      --top-phenotype-pvalue-column '~{sub(top_phenotype_pvalue_column, "'", "'\"'\"'")}'
      --output-dir plan
    )
    if [[ -n '~{sub(select_first([gene_map, ""]), "'", "'\"'\"'")}' ]]; then
      args+=(--gene-map '~{sub(select_first([gene_map, ""]), "'", "'\"'\"'")}')
    fi
    if [[ ~{true="true" false="false" select_top_phenotype_per_cluster} == true ]]; then
      args+=(--select-top-phenotype-per-cluster)
    fi
    for (( j=0; j<${#args[@]}; j++ )); do
      if [[ "${args[j]}" == --gene-map ]]; then require_readable "${args[j+1]}"; fi
    done
    printf '%s\n' '~{sub(sep("\n", chromosomes), "'", "'\"'\"'")}' > chromosomes.txt
    Rscript /opt/r/scripts/plan_pgen_susie.R "${args[@]}"
    echo 'COMPLETE PGEN batch planning'
  >>>
  output {
    Array[File] Worklists = read_lines("plan/batch_worklists.txt")
    Array[File] PhenotypeBeds = read_lines("plan/batch_phenotype_beds.txt")
    Array[String] Chromosomes = read_lines("plan/batch_chromosomes.txt")
    Array[Int] ChromosomeIndices = read_json("plan/batch_chromosome_indices.json")
    File ExpectedOutcomes = "plan/expected_outcomes.tsv"
    File BatchIndex = "plan/batch_index.tsv"
    File TaskStdout = stdout()
    File TaskStderr = stderr()
  }
  runtime {
    docker: docker_image
    cpu: 1
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

task RunPgenSusieBatch {
  input {
    ChromosomePgen pgen_files
    File worklist
    File phenotype_bed
    File qtl_results
    File covariates
    File sample_list
    String missing_code
    File? ancestry
    File? variant_list
    File? additional_genotypes
    Int cis_distance
    Float maf
    Int base_seed
    Boolean memory_debug
    Boolean select_top_phenotype_per_cluster
    String top_phenotype_pvalue_column
    Int plink_memory_mb
    String docker_image
    Int cpu
    Int memory_gb
    Int disk_gb
    Int preemptible_attempts
    Int max_retries
  }
  command <<<
    set -euo pipefail
    echo 'START PGEN SuSiE batch'
    export SUSIER_RUNTIME_IMAGE='~{sub(docker_image, "'", "'\"'\"'")}'
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 VROOM_THREADS=1
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
    require_readable '~{sub(pgen_files.pgen, "'", "'\"'\"'")}' '~{sub(pgen_files.pvar, "'", "'\"'\"'")}' '~{sub(pgen_files.psam, "'", "'\"'\"'")}' '~{sub(worklist, "'", "'\"'\"'")}' '~{sub(phenotype_bed, "'", "'\"'\"'")}' '~{sub(qtl_results, "'", "'\"'\"'")}' '~{sub(covariates, "'", "'\"'\"'")}' '~{sub(sample_list, "'", "'\"'\"'")}'
    args=(
      --pgen '~{sub(pgen_files.pgen, "'", "'\"'\"'")}'
      --pvar '~{sub(pgen_files.pvar, "'", "'\"'\"'")}'
      --psam '~{sub(pgen_files.psam, "'", "'\"'\"'")}'
      --worklist '~{sub(worklist, "'", "'\"'\"'")}'
      --phenotype-bed '~{sub(phenotype_bed, "'", "'\"'\"'")}'
      --qtl-results '~{sub(qtl_results, "'", "'\"'\"'")}'
      --covariates '~{sub(covariates, "'", "'\"'\"'")}'
      --sample-list '~{sub(sample_list, "'", "'\"'\"'")}'
      '--missing-code=~{sub(missing_code, "'", "'\"'\"'")}'
      --cis-distance '~{cis_distance}'
      --maf '~{maf}'
      --base-seed '~{base_seed}'
      --plink-memory-mb '~{plink_memory_mb}'
      --top-phenotype-pvalue-column '~{sub(top_phenotype_pvalue_column, "'", "'\"'\"'")}'
      --output-dir batch
    )
    if [[ -n '~{sub(select_first([ancestry, ""]), "'", "'\"'\"'")}' ]]; then
      args+=(--ancestry '~{sub(select_first([ancestry, ""]), "'", "'\"'\"'")}')
    fi
    if [[ -n '~{sub(select_first([variant_list, ""]), "'", "'\"'\"'")}' ]]; then
      args+=(--variant-list '~{sub(select_first([variant_list, ""]), "'", "'\"'\"'")}')
    fi
    if [[ -n '~{sub(select_first([additional_genotypes, ""]), "'", "'\"'\"'")}' ]]; then
      args+=(--additional-genotypes '~{sub(select_first([additional_genotypes, ""]), "'", "'\"'\"'")}')
    fi
    if [[ ~{true="true" false="false" select_top_phenotype_per_cluster} == true ]]; then
      args+=(--select-top-phenotype-per-cluster)
    fi
    if [[ ~{true="true" false="false" memory_debug} == true ]]; then
      args+=(--memory-debug)
    fi
    for (( j=0; j<${#args[@]}; j++ )); do
      case "${args[j]}" in
        --ancestry|--variant-list|--additional-genotypes) require_readable "${args[j+1]}" ;;
      esac
    done
    Rscript /opt/r/scripts/run_pgen_susie_batch.R "${args[@]}"
    echo 'COMPLETE PGEN SuSiE batch'
  >>>
  output {
    File Summary = "batch/batch.susie.parquet"
    File CredibleSets = "batch/batch.credible_sets.parquet"
    File OutcomeIndex = "batch/batch.outcome_index.tsv"
    File Diagnostics = "batch/batch.covariate_explained_variants.tsv"
    Array[File] GeneRds = glob("batch/genes/*/*_susie.rds")
    Array[File] GeneParquets = glob("batch/genes/*/*.parquet")
    Array[File] GeneTsvs = glob("batch/genes/*/*.tsv")
    Array[File] GeneLogs = glob("batch/genes/*/*.log")
    Array[File] GeneProvenance = glob("batch/genes/*/*.provenance.json")
    Array[File] PlinkLogs = glob("batch/*.log")
    File TaskStdout = stdout()
    File TaskStderr = stderr()
  }
  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

task AggregatePgenSusie {
  input {
    File expected_outcomes
    Array[File] summary_files
    Array[File] credible_set_files
    Array[File] index_files
    Array[File] diagnostic_files
    String docker_image
    Int memory_gb
    Int disk_gb
    Int preemptible_attempts
    Int max_retries
  }
  command <<<
    set -euo pipefail
    echo 'START PGEN batch aggregation'
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
    require_readable '~{sub(expected_outcomes, "'", "'\"'\"'")}'
    printf '%s\n' '~{sub(sep("\n", summary_files), "'", "'\"'\"'")}' > summary_files.txt
    while IFS= read -r path; do require_readable "$path"; done < summary_files.txt
    printf '%s\n' '~{sub(sep("\n", credible_set_files), "'", "'\"'\"'")}' > credible_set_files.txt
    while IFS= read -r path; do require_readable "$path"; done < credible_set_files.txt
    printf '%s\n' '~{sub(sep("\n", index_files), "'", "'\"'\"'")}' > index_files.txt
    while IFS= read -r path; do require_readable "$path"; done < index_files.txt
    printf '%s\n' '~{sub(sep("\n", diagnostic_files), "'", "'\"'\"'")}' > diagnostic_files.txt
    while IFS= read -r path; do require_readable "$path"; done < diagnostic_files.txt
    Rscript /opt/r/scripts/aggregate_pgen_susie.R \
      --expected-outcomes '~{sub(expected_outcomes, "'", "'\"'\"'")}' \
      --summary-files summary_files.txt \
      --credible-set-files credible_set_files.txt \
      --index-files index_files.txt \
      --diagnostic-files diagnostic_files.txt \
      --output-dir aggregate
    echo 'COMPLETE PGEN batch aggregation'
  >>>
  output {
    File Summary = "aggregate/susie.parquet"
    File CredibleSets = "aggregate/credible_sets.parquet"
    File OutcomeIndex = "aggregate/outcome_index.tsv"
    File Diagnostics = "aggregate/covariate_explained_variants.tsv"
    File TaskStdout = stdout()
    File TaskStderr = stderr()
  }
  runtime {
    docker: docker_image
    cpu: 1
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}
