version 1.0

task ExtractMultiPhenotypeInputs {
    input {
        File PhenotypeBed
        File TensorQTLPermutations
        String PhenotypeID
        Array[String]? PhenotypeIDs
        Boolean MatchPhenotypeIDSubstring = false
        Boolean AddSkipRow = true
        Int NumPrempt = 2
    }

    Array[String] phenotype_ids = select_first([PhenotypeIDs, [PhenotypeID]])
    File PhenotypeIDFile = write_lines(phenotype_ids)
    String phenotype_match_mode = if defined(PhenotypeIDs) then "exact-list" else if MatchPhenotypeIDSubstring then "contains" else "exact"

    command <<<
        echo "Extracting phenotype inputs"
        echo "Phenotype match mode: ~{phenotype_match_mode}"
        echo "Requested phenotype IDs:"
        cat "~{PhenotypeIDFile}"

        zcat "~{PhenotypeBed}" | head -n 1 > phenotype_header.txt

        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat "~{PhenotypeBed}" \
                | awk -v needle="~{PhenotypeID}" 'FNR == 1 && $4 == "phenotype_id" {next} index($4, needle) > 0' \
                > phenotype_rows.bed
        else
            zcat "~{PhenotypeBed}" \
                | awk 'NR==FNR {ids[$1]=1; next} FNR == 1 && $4 == "phenotype_id" {next} ($4 in ids)' "~{PhenotypeIDFile}" - \
                > phenotype_rows.bed
        fi

        if [ ! -s phenotype_rows.bed ]; then
            echo "No rows in PhenotypeBed matched the requested phenotype IDs" >&2
            exit 1
        fi

        awk -F'\t' '{print $4}' phenotype_rows.bed | sort -u > matched_phenotype_ids.txt
        echo "Matched phenotype IDs:"
        cat matched_phenotype_ids.txt

        if ~{AddSkipRow}; then
            head -n 1 phenotype_rows.bed | awk -F'\t' 'BEGIN{OFS="\t"} {$4="skip"; print}' > skip.txt
            cat phenotype_header.txt phenotype_rows.bed skip.txt > "~{PhenotypeID}.phenotype_matrix.bed"
        else
            cat phenotype_header.txt phenotype_rows.bed > "~{PhenotypeID}.phenotype_matrix.bed"
        fi

        cat phenotype_header.txt phenotype_rows.bed | bgzip -c - > "~{PhenotypeID}.phenotypes.bed.gz"

        zcat "~{TensorQTLPermutations}" | head -n 1 > tensorqtl_header.txt
        zcat "~{TensorQTLPermutations}" \
            | awk 'NR==FNR {ids[$1]=1; next} FNR == 1 && $1 == "phenotype_id" {next} ($1 in ids)' matched_phenotype_ids.txt - \
            > tensorqtl_rows.txt

        if [ ! -s tensorqtl_rows.txt ]; then
            echo "No rows in TensorQTLPermutations matched the extracted phenotype IDs" >&2
            exit 1
        fi

        cat tensorqtl_header.txt tensorqtl_rows.txt > "~{PhenotypeID}.tensorqtl.txt"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/susier:main"
        disks: "local-disk 10 SSD"
        preemptible: "${NumPrempt}"
        memory: "1GB"
        cpu: "1"
    }

    output {
        File PhenotypeMatrix = "~{PhenotypeID}.phenotype_matrix.bed"
        File SubsetPhenotypeBed = "~{PhenotypeID}.phenotypes.bed.gz"
        File SubsetPermutationPvals = "~{PhenotypeID}.tensorqtl.txt"
        File MatchedPhenotypeIDs = "matched_phenotype_ids.txt"
    }
}

workflow ExtractMultiPhenotypeInputsWorkflow {
    input {
        File PhenotypeBed
        File TensorQTLPermutations
        String PhenotypeID
        Array[String]? PhenotypeIDs
        Boolean MatchPhenotypeIDSubstring = false
        Boolean AddSkipRow = true
        Int NumPrempt = 2
    }

    call ExtractMultiPhenotypeInputs {
        input:
            PhenotypeBed = PhenotypeBed,
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
            PhenotypeIDs = PhenotypeIDs,
            MatchPhenotypeIDSubstring = MatchPhenotypeIDSubstring,
            AddSkipRow = AddSkipRow,
            NumPrempt = NumPrempt
    }

    output {
        File PhenotypeMatrix = ExtractMultiPhenotypeInputs.PhenotypeMatrix
        File SubsetPhenotypeBed = ExtractMultiPhenotypeInputs.SubsetPhenotypeBed
        File SubsetPermutationPvals = ExtractMultiPhenotypeInputs.SubsetPermutationPvals
        File MatchedPhenotypeIDs = ExtractMultiPhenotypeInputs.MatchedPhenotypeIDs
    }
}
