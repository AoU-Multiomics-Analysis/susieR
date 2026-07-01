version 1.0

task PrepInputs {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        String PhenotypeID
        Array[String]? PhenotypeIDs
        Boolean MatchPhenotypeIDSubstring = false
        File PhenotypeBed
        File TensorQTLPermutations
        Int NumPrempt
        Int WindowSize = 1000000
    }
    Array[String] phenotype_ids = select_first([PhenotypeIDs, [PhenotypeID]])
    File PhenotypeIDFile = write_lines(phenotype_ids)
    String phenotype_match_mode = if defined(PhenotypeIDs) then "exact-list" else if MatchPhenotypeIDSubstring then "contains" else "exact"

    command <<<
        echo "Extracting headers from files"
        headerPermutations=$(zcat "~{TensorQTLPermutations}" | head -n 1)
        headerBed=$(zcat "~{PhenotypeBed}" | head -n 1)
        echo "Phenotype match mode: ~{phenotype_match_mode}"
        echo "Phenotype IDs selected for this run:"
        cat "~{PhenotypeIDFile}"
        
 
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
                | awk 'BEGIN{OFS="\t"} NR==FNR {ids[$1]=1; next} FNR == 1 && $4 == "phenotype_id" {next} ($4 in ids) {$2=$2-~{WindowSize}; $3=$3+~{WindowSize}; if($2<1) $2=1; print}' "~{PhenotypeIDFile}" - \
                > feature.bed
        fi
        if [ ! -s feature.bed ]; then
            echo "No rows in PhenotypeBed matched the requested phenotype IDs" >&2
            exit 1
        fi
        
        #echo $headerBed > temp_header.txt
        zcat "~{PhenotypeBed}" | head -n 1 > temp_header.txt
        cat temp_header.txt feature.bed | bgzip -c - > ~{PhenotypeID}.bed.bgz
        #tabix ~{PhenotypeID}.bed.bgz

        echo "Subsetting TensorQTL file"
        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat "~{TensorQTLPermutations}" \
                | awk -v needle="~{PhenotypeID}" 'FNR == 1 && $1 == "phenotype_id" {next} index($1, needle) > 0' \
                > feature.txt
        else
            zcat "~{TensorQTLPermutations}" \
                | awk 'NR==FNR {ids[$1]=1; next} FNR == 1 && $1 == "phenotype_id" {next} ($1 in ids)' "~{PhenotypeIDFile}" - \
                > feature.txt
        fi
        if [ ! -s feature.txt ]; then
            echo "No rows in TensorQTLPermutations matched the requested phenotype IDs" >&2
            exit 1
        fi
        echo "$headerPermutations" > temp_header_perm.txt
        cat temp_header_perm.txt feature.txt > ~{PhenotypeID}.tensorqtl.txt

        echo "Subsetting dose file"
        #(cat dosage_header.txt; tabix ~{GenotypeDosages} -R ~{PhenotypeID}.bed.bgz) | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
        tabix "~{GenotypeDosages}" -R "~{PhenotypeID}.bed.bgz" > dose.tmp.tsv

        #(head -n 1 dose.tmp.tsv && tail -n +2 dose.tmp.tsv | sort -k1,1V -k2,2n) \
        #| bgzip -c > ~{PhenotypeID}.dose.tsv.gz    
        (cat dosage_header.txt && sort -k1,1V -k2,2n dose.tmp.tsv | awk '!seen[$0]++') \
            | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
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
        Array[String]? PhenotypeIDs
        Boolean MatchPhenotypeIDSubstring = false
        Int WindowSize = 1000000

    }

    call PrepInputs {
        input:
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
            PhenotypeIDs = PhenotypeIDs,
            MatchPhenotypeIDSubstring = MatchPhenotypeIDSubstring,
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
