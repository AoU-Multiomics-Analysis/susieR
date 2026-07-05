version 1.0


task susieR {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        File PhenotypeBed
        Int CisDistance
        String OutputPrefix
        Int memory
        Int NumPrempt
        Float? MAF
        Boolean MatchPhenotypeIDSubstring = false
        Boolean ReuseGenotypeMatrix = false
        Boolean SelectTopPhenotypePerCluster = false
        String TopPhenotypePerClusterPvalueColumn = "qval"
        File? VariantList
        File? AncestryFile
        File? AdditionalGenotypesBed
    }
    String phenotype_match_mode = if MatchPhenotypeIDSubstring then "contains" else "exact"

    command <<<

        echo "Phenotype match mode: ~{phenotype_match_mode}"
        echo "Output prefix selected for this run: ~{OutputPrefix}"
        zcat ~{PhenotypeBed} | head -n 1 > header.txt

        report_phenotype_match_failure() {
            echo "Requested phenotype ID(s): $1" >&2
            zcat ~{PhenotypeBed} \
                | awk -F'\t' '
                    FNR == 1 && $4 == "phenotype_id" {next}
                    $4 != "" && !seen[$4]++ {
                        count++
                        if (count <= 20) preview[count] = $4
                    }
                    END {
                        printf "Available phenotype IDs in PhenotypeBed: %d unique ID(s)\n", count > "/dev/stderr"
                        if (count > 0) {
                            print "First available phenotype IDs:" > "/dev/stderr"
                            for (i = 1; i <= count && i <= 20; i++) print preview[i] > "/dev/stderr"
                            if (count > 20) printf "... (%d additional phenotype IDs not shown)\n", count - 20 > "/dev/stderr"
                        }
                    }'
            zcat ~{PhenotypeBed} \
                | awk -v needle="$1" 'FNR == 1 && $4 == "phenotype_id" {next} index($4, needle) > 0' \
                > phenotype_substring_matches.bed
            substring_match_count=$(wc -l < phenotype_substring_matches.bed | tr -d ' ')
            if [ "$substring_match_count" -gt 0 ]; then
                echo "Substring fallback found ${substring_match_count} PhenotypeBed row(s) containing the requested ID." >&2
                echo "First matching phenotype IDs from substring fallback:" >&2
                awk -F'\t' '!seen[$4]++ {print $4; if (++n == 20) exit}' phenotype_substring_matches.bed >&2
                echo "If these are intended, rerun with MatchPhenotypeIDSubstring=true." >&2
            else
                echo "Substring fallback found 0 PhenotypeBed rows containing the requested ID." >&2
            fi
        }

        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat ~{PhenotypeBed} \
                | awk -v needle="~{OutputPrefix}" 'FNR == 1 && $4 == "phenotype_id" {next} index($4, needle) > 0' \
                > input_gene.txt
        else
            zcat ~{PhenotypeBed} \
                | awk -v phenotype_id="~{OutputPrefix}" 'FNR == 1 && $4 == "phenotype_id" {next} $4 == phenotype_id' \
                > input_gene.txt
        fi
        if [ ! -s input_gene.txt ]; then
            echo "No rows in PhenotypeBed matched the requested phenotype IDs" >&2
            report_phenotype_match_failure "~{OutputPrefix}"
            exit 1
        fi
        head -n 1 input_gene.txt | awk -F'\t' 'BEGIN{OFS="\t"} {$4="skip"; print}' > skip.txt        
        cat header.txt input_gene.txt skip.txt > input_gene.bed  

        Rscript /tmp/susie.R ~{if defined(MAF) then "--MAF ~{MAF}  " else ""} ~{if ReuseGenotypeMatrix then "--reuse_genotype_matrix true  " else ""} ~{if SelectTopPhenotypePerCluster then "--select_top_phenotype_per_cluster true --top_phenotype_pvalue_column " else ""}~{if SelectTopPhenotypePerCluster then TopPhenotypePerClusterPvalueColumn else ""} ~{if defined(AncestryFile) then "--AncestryMetadata ~{AncestryFile}  "  else ""} ~{if defined(VariantList) then "--VariantList ~{VariantList}  "  else ""}  ~{if defined(AdditionalGenotypesBed) then "--AdditionalGenotypesBed ~{AdditionalGenotypesBed}  "  else ""} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix input_gene.bed \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance} 

    >>>

    runtime {
        docker: 'ghcr.io/aou-multiomics-analysis/susier:main'
        memory: "${memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        preemptible: "${NumPrempt}"
        cpu: "1"
    }

    output {
        File SusieParquet = "${OutputPrefix}.parquet"
        File lbfParquet = "${OutputPrefix}.lbf_variable.parquet"
        File FullSusieParquet = "${OutputPrefix}.full_susie.parquet"
        File SusieObject = "${OutputPrefix}_susie.rds"
    }
}


workflow SusieROnlyWorkflow {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        File PhenotypeBed
        Int CisDistance
        Int memory
        Int NumPrempt
        String OutputPrefix
        Float? MAF
        Boolean MatchPhenotypeIDSubstring = false
        Boolean ReuseGenotypeMatrix = false
        Boolean SelectTopPhenotypePerCluster = false
        String TopPhenotypePerClusterPvalueColumn = "qval"
        File? VariantList
        File? AncestryFile
        File? AdditionalGenotypesBed
    }
    call susieR {
        input:
            GenotypeDosages = GenotypeDosages,
            GenotypeDosageIndex = GenotypeDosageIndex,
            QTLCovariates = QTLCovariates,
            TensorQTLPermutations = TensorQTLPermutations,
            SampleList = SampleList,
            PhenotypeBed = PhenotypeBed ,
            CisDistance = CisDistance,
            OutputPrefix = OutputPrefix,
            memory = memory,
            NumPrempt = NumPrempt,
            MAF = MAF,
            MatchPhenotypeIDSubstring = MatchPhenotypeIDSubstring,
            ReuseGenotypeMatrix = ReuseGenotypeMatrix,
            SelectTopPhenotypePerCluster = SelectTopPhenotypePerCluster,
            TopPhenotypePerClusterPvalueColumn = TopPhenotypePerClusterPvalueColumn,
            VariantList = VariantList,
            AncestryFile = AncestryFile,
            AdditionalGenotypesBed = AdditionalGenotypesBed
        }
        output {
            File SusieParquet = susieR.SusieParquet
            File SusielbfParquet = susieR.lbfParquet
            File FullSusieParquet = susieR.FullSusieParquet
            File SusieObject = susieR.SusieObject
        }
}
