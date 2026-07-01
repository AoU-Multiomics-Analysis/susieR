version 1.0

task ExtractMultiPhenotypeInputs {
    input {
        File PhenotypeBed
        File TensorQTLPermutations
        String PhenotypeID
        Boolean MatchPhenotypeIDSubstring = false
        Boolean SelectTopPhenotypePerCluster = false
        String TopPhenotypePerClusterPvalueColumn = "qval"
        Boolean AddSkipRow = true
        Int NumPrempt = 2
    }

    String phenotype_match_mode = if MatchPhenotypeIDSubstring then "contains" else "exact"

    command <<<
        echo "Extracting phenotype inputs"
        echo "Phenotype match mode: ~{phenotype_match_mode}"
        echo "Requested phenotype ID: ~{PhenotypeID}"

        zcat "~{PhenotypeBed}" | head -n 1 > phenotype_header.txt

        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat "~{PhenotypeBed}" \
                | awk -v needle="~{PhenotypeID}" 'FNR == 1 && $4 == "phenotype_id" {next} index($4, needle) > 0' \
                > phenotype_rows.bed
        else
            zcat "~{PhenotypeBed}" \
                | awk -v phenotype_id="~{PhenotypeID}" 'FNR == 1 && $4 == "phenotype_id" {next} $4 == phenotype_id' \
                > phenotype_rows.bed
        fi

        if [ ! -s phenotype_rows.bed ]; then
            echo "No rows in PhenotypeBed matched the requested phenotype IDs" >&2
            exit 1
        fi

        awk -F'\t' '{print $4}' phenotype_rows.bed | sort -u > candidate_phenotype_ids.txt

        zcat "~{TensorQTLPermutations}" | head -n 1 > tensorqtl_header.txt
        zcat "~{TensorQTLPermutations}" \
            | awk 'NR==FNR {ids[$1]=1; next} FNR == 1 && $1 == "phenotype_id" {next} ($1 in ids)' candidate_phenotype_ids.txt - \
            > tensorqtl_rows.txt

        if [ ! -s tensorqtl_rows.txt ]; then
            echo "No rows in TensorQTLPermutations matched the extracted phenotype IDs" >&2
            exit 1
        fi

        if [ "~{SelectTopPhenotypePerCluster}" = "true" ]; then
            echo "Selecting strongest significant phenotype per LeafCutter cluster"
            awk -v requested_pvalue="~{TopPhenotypePerClusterPvalueColumn}" '
                BEGIN {FS=OFS="\t"}
                FNR == NR {
                    for (i = 1; i <= NF; i++) {
                        if ($i == requested_pvalue) {pcol = i; pcol_name = $i}
                        if ($i == "pval_beta") {pval_beta = i}
                        if ($i == "pval_perm") {pval_perm = i}
                        if ($i == "pval_nominal") {pval_nominal = i}
                        if ($i == "pval_true_df") {pval_true_df = i}
                        if ($i == "qval") {qval = i; qval_col = i}
                    }
                    if (!pcol && pval_beta) {pcol = pval_beta; pcol_name = "pval_beta"}
                    if (!pcol && pval_perm) {pcol = pval_perm; pcol_name = "pval_perm"}
                    if (!pcol && pval_nominal) {pcol = pval_nominal; pcol_name = "pval_nominal"}
                    if (!pcol && pval_true_df) {pcol = pval_true_df; pcol_name = "pval_true_df"}
                    if (!pcol && qval) {pcol = qval; pcol_name = "qval"}
                    if (!pcol) {
                        print "No usable p-value column found for top phenotype selection" > "/dev/stderr"
                        exit 2
                    }
                    print "Top phenotype per cluster p-value column: " pcol_name > "/dev/stderr"
                    next
                }
                function cluster_id(pid, fields, n, i) {
                    n = split(pid, fields, ":")
                    for (i = 1; i <= n; i++) {
                        if (fields[i] ~ /^clu_/) {return fields[i]}
                    }
                    return pid
                }
                function is_missing(value) {
                    return value == "" || value == "NA" || value == "NaN" || value == "nan"
                }
                {
                    if (qval_col && (is_missing($qval_col) || $qval_col + 0 >= 0.05)) {
                        next
                    }
                    cid = cluster_id($1)
                    missing = is_missing($pcol)
                    pvalue = missing ? 0 : $pcol + 0
                    if (!(cid in seen) ||
                        (best_missing[cid] && !missing) ||
                        (best_missing[cid] == missing && !missing && pvalue < best_pvalue[cid]) ||
                        (best_missing[cid] == missing && (missing || pvalue == best_pvalue[cid]) && $1 < best_id[cid])) {
                        seen[cid] = 1
                        best_id[cid] = $1
                        best_pvalue[cid] = pvalue
                        best_missing[cid] = missing
                    }
                }
                END {
                    if (pcol) {
                        for (cid in best_id) {print best_id[cid]}
                    }
                }
            ' tensorqtl_header.txt tensorqtl_rows.txt > selected_phenotype_ids.txt

            if [ ! -s selected_phenotype_ids.txt ]; then
                echo "No representative phenotypes were selected" >&2
                exit 1
            fi
            awk 'NR==FNR {ids[$1]=1; next} ($1 in ids)' selected_phenotype_ids.txt tensorqtl_rows.txt > tensorqtl_rows.selected.txt
            awk 'NR==FNR {ids[$1]=1; next} ($4 in ids)' selected_phenotype_ids.txt phenotype_rows.bed > phenotype_rows.selected.bed
            mv tensorqtl_rows.selected.txt tensorqtl_rows.txt
            mv phenotype_rows.selected.bed phenotype_rows.bed
            if [ ! -s tensorqtl_rows.txt ] || [ ! -s phenotype_rows.bed ]; then
                echo "Representative phenotype selection produced no overlapping TensorQTL/BED rows" >&2
                exit 1
            fi
        fi

        awk -F'\t' '{print $4}' phenotype_rows.bed | sort -u > matched_phenotype_ids.txt
        echo "Matched phenotype IDs after optional cluster selection:"
        cat matched_phenotype_ids.txt

        if ~{AddSkipRow}; then
            head -n 1 phenotype_rows.bed | awk -F'\t' 'BEGIN{OFS="\t"} {$4="skip"; print}' > skip.txt
            cat phenotype_header.txt phenotype_rows.bed skip.txt > "~{PhenotypeID}.phenotype_matrix.bed"
        else
            cat phenotype_header.txt phenotype_rows.bed > "~{PhenotypeID}.phenotype_matrix.bed"
        fi

        cat phenotype_header.txt phenotype_rows.bed | bgzip -c - > "~{PhenotypeID}.phenotypes.bed.gz"

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
        Boolean MatchPhenotypeIDSubstring = false
        Boolean SelectTopPhenotypePerCluster = false
        String TopPhenotypePerClusterPvalueColumn = "qval"
        Boolean AddSkipRow = true
        Int NumPrempt = 2
    }

    call ExtractMultiPhenotypeInputs {
        input:
            PhenotypeBed = PhenotypeBed,
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
            MatchPhenotypeIDSubstring = MatchPhenotypeIDSubstring,
            SelectTopPhenotypePerCluster = SelectTopPhenotypePerCluster,
            TopPhenotypePerClusterPvalueColumn = TopPhenotypePerClusterPvalueColumn,
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
