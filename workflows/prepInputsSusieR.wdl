version 1.0

task PrepInputs {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        String PhenotypeID
        Boolean MatchPhenotypeIDSubstring = false
        Boolean SelectTopPhenotypePerCluster = false
        String TopPhenotypePerClusterPvalueColumn = "pval_beta"
        File PhenotypeBed
        File TensorQTLPermutations
        Int NumPrempt
        Int WindowSize = 1000000
    }
    String phenotype_match_mode = if MatchPhenotypeIDSubstring then "contains" else "exact"

    command <<<
        echo "Extracting headers from files"
        headerPermutations=$(zcat "~{TensorQTLPermutations}" | head -n 1)
        headerBed=$(zcat "~{PhenotypeBed}" | head -n 1)
        echo "Phenotype match mode: ~{phenotype_match_mode}"
        echo "Phenotype ID selected for this run: ~{PhenotypeID}"
        
 
        echo "Bed file header:"
        echo "$headerBed"

        echo "TensorQTL file header"
        echo "$headerPermutations"
        
        #zcat ~{GenotypeDosages} |  awk 'NR==1 { if ($0 ~ /^#/) print; else print "#" $0; exit }'  > dosage_header.txt
        zcat "~{GenotypeDosages}" |  awk 'NR==1 {print $0; exit }'  > dosage_header.txt

        echo "Subsetting bed file"
        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat "~{PhenotypeBed}" \
                | awk -v needle="~{PhenotypeID}" 'BEGIN{OFS="\t"} FNR == 1 && $4 == "phenotype_id" {next} index($4, needle) > 0 {$2=$2-~{WindowSize}; $3=$3+~{WindowSize}; if($2<1) $2=1; print}' \
                > feature.bed
        else
            zcat "~{PhenotypeBed}" \
                | awk -v phenotype_id="~{PhenotypeID}" 'BEGIN{OFS="\t"} FNR == 1 && $4 == "phenotype_id" {next} $4 == phenotype_id {$2=$2-~{WindowSize}; $3=$3+~{WindowSize}; if($2<1) $2=1; print}' \
                > feature.bed
        fi
        if [ ! -s feature.bed ]; then
            echo "No rows in PhenotypeBed matched the requested phenotype IDs" >&2
            exit 1
        fi
        
        echo "Subsetting TensorQTL file"
        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat "~{TensorQTLPermutations}" \
                | awk -v needle="~{PhenotypeID}" 'FNR == 1 && $1 == "phenotype_id" {next} index($1, needle) > 0' \
                > feature.txt
        else
            zcat "~{TensorQTLPermutations}" \
                | awk -v phenotype_id="~{PhenotypeID}" 'FNR == 1 && $1 == "phenotype_id" {next} $1 == phenotype_id' \
                > feature.txt
        fi
        if [ ! -s feature.txt ]; then
            echo "No rows in TensorQTLPermutations matched the requested phenotype IDs" >&2
            exit 1
        fi
        echo "$headerPermutations" > temp_header_perm.txt

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
            ' temp_header_perm.txt feature.txt > selected_phenotype_ids.txt

            if [ ! -s selected_phenotype_ids.txt ]; then
                echo "No representative phenotypes were selected" >&2
                exit 1
            fi
            awk 'NR==FNR {ids[$1]=1; next} ($1 in ids)' selected_phenotype_ids.txt feature.txt > feature.selected.txt
            awk 'NR==FNR {ids[$1]=1; next} ($4 in ids)' selected_phenotype_ids.txt feature.bed > feature.selected.bed
            mv feature.selected.txt feature.txt
            mv feature.selected.bed feature.bed
            if [ ! -s feature.txt ] || [ ! -s feature.bed ]; then
                echo "Representative phenotype selection produced no overlapping TensorQTL/BED rows" >&2
                exit 1
            fi
        fi

        #echo $headerBed > temp_header.txt
        zcat "~{PhenotypeBed}" | head -n 1 > temp_header.txt
        cat temp_header.txt feature.bed | bgzip -c - > ~{PhenotypeID}.bed.bgz
        #tabix ~{PhenotypeID}.bed.bgz

        echo "Merging overlapping phenotype windows for dosage extraction"
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' feature.bed \
            | sort -k1,1V -k2,2n -k3,3n \
            | awk 'BEGIN{OFS="\t"} NR == 1 {chr=$1; start=$2; end=$3; next} $1 == chr && $2 <= end {if ($3 > end) end=$3; next} {print chr,start,end; chr=$1; start=$2; end=$3} END {if (NR > 0) print chr,start,end}' \
            > dosage_regions.bed
        echo "Matched phenotype rows after optional cluster selection: $(wc -l < feature.bed)"
        echo "Merged dosage extraction regions: $(wc -l < dosage_regions.bed)"

        cat temp_header_perm.txt feature.txt > ~{PhenotypeID}.tensorqtl.txt
        echo "Matched TensorQTL rows: $(wc -l < feature.txt)"

        echo "Subsetting dose file"
        #(cat dosage_header.txt; tabix ~{GenotypeDosages} -R ~{PhenotypeID}.bed.bgz) | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
        tabix "~{GenotypeDosages}" -R dosage_regions.bed > dose.tmp.tsv
        echo "Extracted dosage rows before de-duplication: $(wc -l < dose.tmp.tsv)"

        #(head -n 1 dose.tmp.tsv && tail -n +2 dose.tmp.tsv | sort -k1,1V -k2,2n) \
        #| bgzip -c > ~{PhenotypeID}.dose.tsv.gz    
        (cat dosage_header.txt && sort -k1,1V -k2,2n dose.tmp.tsv | awk '!seen[$0]++') \
            | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
        echo "Extracted dosage rows after de-duplication: $(zcat ~{PhenotypeID}.dose.tsv.gz | tail -n +2 | wc -l)"
        #tabix  ~{GenotypeDosages} -R ~{PhenotypeID}.bed.bgz | bgzip -c - > ~{PhenotypeID}.dose.tsv.gz
        tabix -s1 -b2 -e2 -S1 "~{PhenotypeID}.dose.tsv.gz"
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/susier:main"
        disks: "local-disk 500 SSD"
        preemptible: "${NumPrempt}"
        memory: "2GB"
        cpu: "1"
    }
    
    output {
        File DosageHeader = "dosage_header.txt"
        File SubsetBed = "~{PhenotypeID}.bed.bgz"
        #File SubsetBedIndex = "~{PhenotypeID}.bed.bgz.tbi" 
        File SubsetPermutationPvals = "~{PhenotypeID}.tensorqtl.txt"
        File SubsetDosages = "~{PhenotypeID}.dose.tsv.gz"
        File SubsetDosagesIndex = "~{PhenotypeID}.dose.tsv.gz.tbi"
    }

}


workflow PrepSusieRWorkflow {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File TensorQTLPermutations
        File PhenotypeBed
        Int NumPrempt
        String PhenotypeID
        Boolean MatchPhenotypeIDSubstring = false
        Boolean SelectTopPhenotypePerCluster = false
        String TopPhenotypePerClusterPvalueColumn = "pval_beta"
        Int WindowSize = 1000000

    }

    call PrepInputs {
        input:
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
            MatchPhenotypeIDSubstring = MatchPhenotypeIDSubstring,
            SelectTopPhenotypePerCluster = SelectTopPhenotypePerCluster,
            TopPhenotypePerClusterPvalueColumn = TopPhenotypePerClusterPvalueColumn,
            GenotypeDosages = GenotypeDosages,
            GenotypeDosageIndex = GenotypeDosageIndex,
            PhenotypeBed = PhenotypeBed,
            NumPrempt = NumPrempt,
            WindowSize = WindowSize
    }
   
    #call MergeSusie {
    #    input:
    #        SusieOutput = susieR.SusieParquet,
    #        OutputPrefix = OutputPrefix 
    #
    #} 
    output {
        File SubsetBed = PrepInputs.SubsetBed
        #File SubsetBedIndex = PrepInputs.SubsetBedIndex
        File SubsetDosages = PrepInputs.SubsetDosages
        File SubsetPermutationPvals = PrepInputs.SubsetPermutationPvals
        File SubsetDosagesIndex = PrepInputs.SubsetDosagesIndex
    }
}
