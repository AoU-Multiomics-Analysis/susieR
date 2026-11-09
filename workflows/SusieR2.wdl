version 1.0

# QC indices refer to unique selected PGEN/PVAR/PSAM triplets. The workflow
# promotes their paths to typed Files before QC or batch localization.
struct ChromosomePgen {
  String chromosome
  File pgen
  File pvar
  File psam
}

workflow SusieR2 {
  input {
    File PgenManifest
    File PhenotypeBed
    File QtlResults
    File Covariates
    File? SampleList
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
    Int MAC = 0
    Float MinCallRate = 0.95
    Boolean SelectTopPhenotypePerCluster = false
    String TopPhenotypePvalueColumn = "qval"
    String DockerImage = "ghcr.io/aou-multiomics-analysis/susier/pgen:main"
    Int Cpu = 1
    Int BatchMemoryGB = 16
    Int BatchDiskGB = 100
    Int PlinkMemoryMB = 1024
    Int PlanMemoryGB = 8
    Int PlanDiskGB = 50
    Int QcMemoryGB = 8
    Int QcDiskGB = 100
    Int AggregateMemoryGB = 8
    Int AggregateDiskGB = 50
    Int PreemptibleAttempts = 3
    Int MaxRetries = 1
  }

  call PlanPgenSusie {
    input:
      pgen_manifest = PgenManifest,
      phenotype_bed = PhenotypeBed,
      qtl_results = QtlResults,
      gene_map = GeneMap,
      batch_size = BatchSize,
      cis_distance = CisDistance,
      mac = MAC,
      min_call_rate = MinCallRate,
      select_top_phenotype_per_cluster = SelectTopPhenotypePerCluster,
      top_phenotype_pvalue_column = TopPhenotypePvalueColumn,
      docker_image = DockerImage,
      memory_gb = PlanMemoryGB,
      disk_gb = PlanDiskGB,
      preemptible_attempts = PreemptibleAttempts,
      max_retries = MaxRetries
  }

  # The planner opens only the small manifest. These declarations make its
  # genotype paths explicit File inputs. Identical selected triplets share QC.
  Array[Array[String]] pgen_rows = read_tsv(PlanPgenSusie.UniquePgenManifest)
  scatter (pgen_row in pgen_rows) {
    ChromosomePgen chromosome_pgen = object {
      chromosome: pgen_row[0],
      pgen: pgen_row[1],
      pvar: pgen_row[2],
      psam: pgen_row[3]
    }
    if (MAC > 0 || MinCallRate > 0) {
      call FilterPgen {
        input:
          pgen = chromosome_pgen.pgen,
          pvar = chromosome_pgen.pvar,
          psam = chromosome_pgen.psam,
          mac = MAC,
          min_call_rate = MinCallRate,
          plink_memory_mb = PlinkMemoryMB,
          docker_image = DockerImage,
          memory_gb = QcMemoryGB,
          disk_gb = QcDiskGB,
          preemptible_attempts = PreemptibleAttempts,
          max_retries = MaxRetries
      }
    }
    File prepared_pgen = select_first([FilterPgen.Pgen, chromosome_pgen.pgen])
    File prepared_pvar = select_first([FilterPgen.Pvar, chromosome_pgen.pvar])
    File prepared_psam = select_first([FilterPgen.Psam, chromosome_pgen.psam])
    ChromosomePgen prepared_files = object {
      chromosome: chromosome_pgen.chromosome,
      pgen: prepared_pgen,
      pvar: prepared_pvar,
      psam: prepared_psam
    }
  }
  Array[ChromosomePgen] ChromosomePgens = prepared_files

  # zip rejects unequal output counts before any batch starts localization.
  scatter (batch in zip(zip(PlanPgenSusie.Worklists, PlanPgenSusie.PhenotypeBeds), zip(PlanPgenSusie.BatchPgenIndices, PlanPgenSusie.Chromosomes))) {
    # Keep a named Int for the chromosome lookup in Dockstore's Cromwell parser.
    Int chromosome_index = batch.right.left
    # These locations are output-manifest metadata, never local input paths.
    String pgen_chromosome = batch.right.right
    String pgen_uri = ChromosomePgens[chromosome_index].pgen
    String pvar_uri = ChromosomePgens[chromosome_index].pvar
    String psam_uri = ChromosomePgens[chromosome_index].psam
    call RunPgenSusieBatch {
      input:
        pgen = ChromosomePgens[chromosome_index].pgen,
        pvar = ChromosomePgens[chromosome_index].pvar,
        psam = ChromosomePgens[chromosome_index].psam,
        worklist = batch.left.left,
        phenotype_bed = batch.left.right,
        qtl_results = QtlResults,
        covariates = Covariates,
        sample_list = SampleList,
        missing_code = LegacyMissingCode,
        ancestry = Ancestry,
        variant_list = VariantList,
        additional_genotypes = AdditionalGenotypes,
        genotype_qc_report = FilterPgen.Report[chromosome_index],
        cis_distance = CisDistance,
        maf = MAF,
        mac = 0,
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
      # Completed output locations are string metadata for manifests only.
      # The aggregate task never opens or localizes these URI values.
      gene_rds_uris = flatten(RunPgenSusieBatch.GeneRds),
      gene_parquet_uris = flatten(RunPgenSusieBatch.GeneParquets),
      gene_tsv_uris = flatten(RunPgenSusieBatch.GeneTsvs),
      gene_log_uris = flatten(RunPgenSusieBatch.GeneLogs),
      gene_provenance_uris = flatten(RunPgenSusieBatch.GeneProvenance),
      pgen_chromosomes = pgen_chromosome,
      pgen_uris = pgen_uri,
      pvar_uris = pvar_uri,
      psam_uris = psam_uri,
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
    File PreparedPgenManifest = AggregatePgenSusie.PreparedPgenManifest
    Array[File] PreparedPgens = prepared_pgen
    Array[File] PreparedPvars = prepared_pvar
    Array[File] PreparedPsams = prepared_psam
    Array[File] GenotypeQcCounts = select_all(FilterPgen.Counts)
    Array[File] GenotypeQcReports = select_all(FilterPgen.Report)
    Array[File] GenotypeQcLogs = select_all(FilterPgen.TaskStdout)
    Array[File] GenotypeQcErrors = select_all(FilterPgen.TaskStderr)
    Array[File] GeneRds = flatten(RunPgenSusieBatch.GeneRds)
    Array[File] GeneParquets = flatten(RunPgenSusieBatch.GeneParquets)
    Array[File] GeneTsvs = flatten(RunPgenSusieBatch.GeneTsvs)
    Array[File] GeneLogs = flatten(RunPgenSusieBatch.GeneLogs)
    Array[File] GeneProvenance = flatten(RunPgenSusieBatch.GeneProvenance)
    Array[File] ChromosomeGeneManifests = AggregatePgenSusie.ChromosomeGeneManifests
    Array[File] BatchLogs = RunPgenSusieBatch.TaskStdout
    Array[File] BatchErrors = RunPgenSusieBatch.TaskStderr
    Array[File] PlinkLogs = flatten(RunPgenSusieBatch.PlinkLogs)
    File PlanLog = PlanPgenSusie.TaskStderr
    File AggregateLog = AggregatePgenSusie.TaskStderr
  }
}

task PlanPgenSusie {
  input {
    File pgen_manifest
    File phenotype_bed
    File qtl_results
    File? gene_map
    Int batch_size
    Int cis_distance
    Int mac = 0
    Float min_call_rate = 0.95
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
    require_readable '~{sub(pgen_manifest, "'", "'\"'\"'")}' '~{sub(phenotype_bed, "'", "'\"'\"'")}' '~{sub(qtl_results, "'", "'\"'\"'")}'
    args=(
      --phenotype-bed '~{sub(phenotype_bed, "'", "'\"'\"'")}'
      --qtl-results '~{sub(qtl_results, "'", "'\"'\"'")}'
      --pgen-manifest '~{sub(pgen_manifest, "'", "'\"'\"'")}'
      --batch-size '~{batch_size}'
      --cis-distance '~{cis_distance}'
      '--mac=~{mac}'
      '--min-call-rate=~{min_call_rate}'
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
    Rscript /opt/r/scripts/plan_pgen_susie.R "${args[@]}"
    echo 'CHECK planned batch files before cloud output collection'
    require_readable plan/batch_worklists.txt plan/batch_phenotype_beds.txt
    require_readable plan/unique_pgens.tsv plan/batch_pgen_indices.json
    if [[ ! -s plan/batch_worklists.txt || ! -s plan/batch_phenotype_beds.txt ]]; then
      echo 'ERROR planner must list at least one worklist and phenotype BED' >&2
      exit 1
    fi
    shopt -s nullglob
    worklists=(plan/batches/batch_*.tsv)
    phenotype_beds=(plan/batches/batch_*.bed)
    if [[ ! -f "${worklists[0]-}" || ! -f "${phenotype_beds[0]-}" ]]; then
      echo 'ERROR planner produced no worklists or phenotype BEDs' >&2
      exit 1
    fi
    listed_worklists=()
    listed_beds=()
    while IFS= read -r path || [[ -n "$path" ]]; do listed_worklists+=("$path"); done < plan/batch_worklists.txt
    while IFS= read -r path || [[ -n "$path" ]]; do listed_beds+=("$path"); done < plan/batch_phenotype_beds.txt
    if ! [[ ${#worklists[@]} -eq ${#phenotype_beds[@]} && ${#worklists[@]} -eq ${#listed_worklists[@]} && ${#phenotype_beds[@]} -eq ${#listed_beds[@]} ]]; then
      echo 'ERROR planner batch file counts do not match the ordered output lists' >&2
      exit 1
    fi
    for (( i=0; i<${#worklists[@]}; i++ )); do
      if [[ "${worklists[i]%.tsv}" != "${phenotype_beds[i]%.bed}" || "${worklists[i]}" != "${listed_worklists[i]}" || "${phenotype_beds[i]}" != "${listed_beds[i]}" ]]; then
        printf 'ERROR planner batch files or output order differ at batch index %s\n' "$i" >&2
        exit 1
      fi
    done
    require_readable "${worklists[@]}" "${phenotype_beds[@]}"
    echo 'COMPLETE PGEN batch planning'
  >>>
  output {
    # Explicit globs let the Google backend collect files before the VM exits.
    # Both patterns follow the same batch-ID order checked in the command.
    Array[File] Worklists = glob("plan/batches/batch_*.tsv")
    Array[File] PhenotypeBeds = glob("plan/batches/batch_*.bed")
    Array[String] Chromosomes = read_lines("plan/batch_chromosomes.txt")
    Array[Int] ChromosomeIndices = read_json("plan/batch_chromosome_indices.json")
    File ExpectedOutcomes = "plan/expected_outcomes.tsv"
    File BatchIndex = "plan/batch_index.tsv"
    File NormalizedManifest = "plan/chromosome_pgens.tsv"
    File UniquePgenManifest = "plan/unique_pgens.tsv"
    Array[Int] BatchPgenIndices = read_json("plan/batch_pgen_indices.json")
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

task FilterPgen {
  input {
    File pgen
    File pvar
    File psam
    Int mac = 0
    Float min_call_rate = 0.95
    Int plink_memory_mb
    String docker_image
    Int memory_gb
    Int disk_gb
    Int preemptible_attempts
    Int max_retries
  }
  command <<<
    set -euo pipefail
    echo 'START shared PGEN MAC and call-rate QC'
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
    require_readable '~{sub(pgen, "'", "'\"'\"'")}' '~{sub(pvar, "'", "'\"'\"'")}' '~{sub(psam, "'", "'\"'\"'")}'
    Rscript /opt/r/scripts/filter_pgen_susie.R \
      --pgen '~{sub(pgen, "'", "'\"'\"'")}' \
      --pvar '~{sub(pvar, "'", "'\"'\"'")}' \
      --psam '~{sub(psam, "'", "'\"'\"'")}' \
      '--mac=~{mac}' \
      '--min-call-rate=~{min_call_rate}' \
      --plink-memory-mb '~{plink_memory_mb}' \
      --output-dir qc
    require_readable qc/filtered.pgen qc/filtered.pvar qc/filtered.psam qc/variant_qc.tsv.gz qc/qc.json
    echo 'COMPLETE shared PGEN MAC and call-rate QC'
  >>>
  output {
    File Pgen = "qc/filtered.pgen"
    File Pvar = "qc/filtered.pvar"
    File Psam = "qc/filtered.psam"
    File Counts = "qc/variant_qc.tsv.gz"
    File Report = "qc/qc.json"
    Array[File] PlinkLogs = glob("qc/*.log")
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
    # Separate File inputs let Cromwell map paths before sub() escapes them.
    File pgen
    File pvar
    File psam
    File worklist
    File phenotype_bed
    File qtl_results
    File covariates
    File? sample_list
    String missing_code
    File? ancestry
    File? variant_list
    File? additional_genotypes
    File? genotype_qc_report
    Int cis_distance
    Float maf
    Int mac = 0
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
    require_readable '~{sub(pgen, "'", "'\"'\"'")}' '~{sub(pvar, "'", "'\"'\"'")}' '~{sub(psam, "'", "'\"'\"'")}' '~{sub(worklist, "'", "'\"'\"'")}' '~{sub(phenotype_bed, "'", "'\"'\"'")}' '~{sub(qtl_results, "'", "'\"'\"'")}' '~{sub(covariates, "'", "'\"'\"'")}'
    args=(
      --pgen '~{sub(pgen, "'", "'\"'\"'")}'
      --pvar '~{sub(pvar, "'", "'\"'\"'")}'
      --psam '~{sub(psam, "'", "'\"'\"'")}'
      --worklist '~{sub(worklist, "'", "'\"'\"'")}'
      --phenotype-bed '~{sub(phenotype_bed, "'", "'\"'\"'")}'
      --qtl-results '~{sub(qtl_results, "'", "'\"'\"'")}'
      --covariates '~{sub(covariates, "'", "'\"'\"'")}'
      '--missing-code=~{sub(missing_code, "'", "'\"'\"'")}'
      --cis-distance '~{cis_distance}'
      --maf '~{maf}'
      '--mac=~{mac}'
      --base-seed '~{base_seed}'
      --plink-memory-mb '~{plink_memory_mb}'
      --top-phenotype-pvalue-column '~{sub(top_phenotype_pvalue_column, "'", "'\"'\"'")}'
      --output-dir batch
    )
    if [[ -n '~{sub(select_first([sample_list, ""]), "'", "'\"'\"'")}' ]]; then
      args+=(--sample-list '~{sub(select_first([sample_list, ""]), "'", "'\"'\"'")}')
    fi
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
    if [[ -n '~{sub(select_first([genotype_qc_report, ""]), "'", "'\"'\"'")}' ]]; then
      args+=(--genotype-qc-report '~{sub(select_first([genotype_qc_report, ""]), "'", "'\"'\"'")}')
    fi
    for (( j=0; j<${#args[@]}; j++ )); do
      case "${args[j]}" in
        --sample-list|--ancestry|--variant-list|--additional-genotypes|--genotype-qc-report) require_readable "${args[j+1]}" ;;
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
    # These completed output paths are metadata, not local input filenames.
    Array[String] gene_rds_uris
    Array[String] gene_parquet_uris
    Array[String] gene_tsv_uris
    Array[String] gene_log_uris
    Array[String] gene_provenance_uris
    Array[String] pgen_chromosomes
    Array[String] pgen_uris
    Array[String] pvar_uris
    Array[String] psam_uris
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
    # Command-time writers receive localized Files. Direct File placeholders
    # also let Cromwell map each engine-generated list into the container.
    cp '~{write_lines(summary_files)}' summary_files.txt
    while IFS= read -r path; do require_readable "$path"; done < summary_files.txt
    cp '~{write_lines(credible_set_files)}' credible_set_files.txt
    while IFS= read -r path; do require_readable "$path"; done < credible_set_files.txt
    cp '~{write_lines(index_files)}' index_files.txt
    while IFS= read -r path; do require_readable "$path"; done < index_files.txt
    cp '~{write_lines(diagnostic_files)}' diagnostic_files.txt
    while IFS= read -r path; do require_readable "$path"; done < diagnostic_files.txt
    # Preserve output URI metadata unchanged. Do not open or localize its paths.
    cp '~{write_lines(gene_rds_uris)}' gene_rds_uris.txt
    cp '~{write_lines(gene_parquet_uris)}' gene_parquet_uris.txt
    cp '~{write_lines(gene_tsv_uris)}' gene_tsv_uris.txt
    cp '~{write_lines(gene_log_uris)}' gene_log_uris.txt
    cp '~{write_lines(gene_provenance_uris)}' gene_provenance_uris.txt
    cp '~{write_lines(pgen_chromosomes)}' pgen_chromosomes.txt
    cp '~{write_lines(pgen_uris)}' pgen_uris.txt
    cp '~{write_lines(pvar_uris)}' pvar_uris.txt
    cp '~{write_lines(psam_uris)}' psam_uris.txt
    Rscript /opt/r/scripts/aggregate_pgen_susie.R \
      --expected-outcomes '~{sub(expected_outcomes, "'", "'\"'\"'")}' \
      --summary-files summary_files.txt \
      --credible-set-files credible_set_files.txt \
      --index-files index_files.txt \
      --diagnostic-files diagnostic_files.txt \
      --gene-rds-uris-file gene_rds_uris.txt \
      --gene-parquet-uris-file gene_parquet_uris.txt \
      --gene-tsv-uris-file gene_tsv_uris.txt \
      --gene-log-uris-file gene_log_uris.txt \
      --gene-provenance-uris-file gene_provenance_uris.txt \
      --pgen-chromosomes-file pgen_chromosomes.txt \
      --pgen-uris-file pgen_uris.txt \
      --pvar-uris-file pvar_uris.txt \
      --psam-uris-file psam_uris.txt \
      --output-dir aggregate
    echo 'COMPLETE PGEN batch aggregation'
  >>>
  output {
    File Summary = "aggregate/susie.parquet"
    File CredibleSets = "aggregate/credible_sets.parquet"
    File OutcomeIndex = "aggregate/outcome_index.tsv"
    File Diagnostics = "aggregate/covariate_explained_variants.tsv"
    File PreparedPgenManifest = "aggregate/pgen_manifest.tsv"
    Array[File] ChromosomeGeneManifests = glob("aggregate/gene_manifests/chr*.gene_outputs.tsv")
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
