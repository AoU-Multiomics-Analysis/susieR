version 1.0

task PrepInputs {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        String PhenotypeID
        Boolean MatchPhenotypeIDSubstring = false
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
        : > feature_sites.tsv
        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat "~{PhenotypeBed}" \
                | awk -v needle="~{PhenotypeID}" -v window_size=~{WindowSize} 'BEGIN{OFS="\t"} FNR == 1 && $4 == "phenotype_id" {next} index($4, needle) > 0 {feature_pos=$2; window_start=$2-window_size; window_end=$3+window_size; if(window_start<1) window_start=1; print $4,$1,feature_pos,window_start,window_end,window_end-window_start+1 >> "feature_sites.tsv"; $2=window_start; $3=window_end; print}' \
                > feature.bed
        else
            zcat "~{PhenotypeBed}" \
                | awk -v phenotype_id="~{PhenotypeID}" -v window_size=~{WindowSize} 'BEGIN{OFS="\t"} FNR == 1 && $4 == "phenotype_id" {next} $4 == phenotype_id {feature_pos=$2; window_start=$2-window_size; window_end=$3+window_size; if(window_start<1) window_start=1; print $4,$1,feature_pos,window_start,window_end,window_end-window_start+1 >> "feature_sites.tsv"; $2=window_start; $3=window_end; print}' \
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

        #echo $headerBed > temp_header.txt
        zcat "~{PhenotypeBed}" | head -n 1 > temp_header.txt
        cat temp_header.txt feature.bed | bgzip -c - > ~{PhenotypeID}.bed.bgz
        #tabix ~{PhenotypeID}.bed.bgz

        echo "Merging overlapping phenotype windows for dosage extraction"
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' feature.bed \
            | sort -k1,1V -k2,2n -k3,3n \
            | awk 'BEGIN{OFS="\t"} NR == 1 {chr=$1; start=$2; end=$3; next} $1 == chr && $2 <= end {if ($3 > end) end=$3; next} {print chr,start,end; chr=$1; start=$2; end=$3} END {if (NR > 0) print chr,start,end}' \
            > dosage_regions.bed
        echo "Matched phenotype rows: $(wc -l < feature.bed)"
        echo "Merged dosage extraction regions: $(wc -l < dosage_regions.bed)"

        cat temp_header_perm.txt feature.txt > ~{PhenotypeID}.tensorqtl.txt
        echo "Matched TensorQTL rows: $(wc -l < feature.txt)"

        echo "Subsetting dose file"
        #(cat dosage_header.txt; tabix ~{GenotypeDosages} -R ~{PhenotypeID}.bed.bgz) | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
        tabix "~{GenotypeDosages}" -R dosage_regions.bed > dose.tmp.tsv
        echo "Extracted dosage rows before de-duplication: $(wc -l < dose.tmp.tsv)"

        #(head -n 1 dose.tmp.tsv && tail -n +2 dose.tmp.tsv | sort -k1,1V -k2,2n) \
        #| bgzip -c > ~{PhenotypeID}.dose.tsv.gz    
        sort -k1,1V -k2,2n dose.tmp.tsv | awk '!seen[$0]++' > dose.sorted.tsv
        echo "Extracted dosage rows after de-duplication: $(wc -l < dose.sorted.tsv)"
        printf "phenotype_id\tfeature_chromosome\tfeature_position\trequested_window_start\trequested_window_end\trequested_window_size_bp\ttotal_variants_extracted\tupstream_variants\tdownstream_variants\tvariants_at_feature_coordinate\tmin_abs_distance_bp\tmax_abs_distance_bp\tmin_signed_distance_bp\tmax_signed_distance_bp\n" > "~{PhenotypeID}.prep_variant_position_summary.tsv"
        awk 'BEGIN{FS=OFS="\t"} NR==FNR {variant_count++; variant_chr[variant_count]=$1; variant_pos[variant_count]=$2+0; next} {phenotype_id=$1; feature_chr=$2; feature_pos=$3+0; window_start=$4+0; window_end=$5+0; window_size=$6+0; total=0; upstream=0; downstream=0; at_feature=0; min_abs="NA"; max_abs="NA"; min_signed="NA"; max_signed="NA"; for (i=1; i<=variant_count; i++) {if (variant_chr[i] != feature_chr || variant_pos[i] < window_start || variant_pos[i] > window_end) next; distance=variant_pos[i]-feature_pos; abs_distance=distance < 0 ? -distance : distance; total++; if (distance < 0) upstream++; else if (distance > 0) downstream++; else at_feature++; if (min_abs == "NA" || abs_distance < min_abs) min_abs=abs_distance; if (max_abs == "NA" || abs_distance > max_abs) max_abs=abs_distance; if (min_signed == "NA" || distance < min_signed) min_signed=distance; if (max_signed == "NA" || distance > max_signed) max_signed=distance} print phenotype_id,feature_chr,feature_pos,window_start,window_end,window_size,total,upstream,downstream,at_feature,min_abs,max_abs,min_signed,max_signed}' dose.sorted.tsv feature_sites.tsv >> "~{PhenotypeID}.prep_variant_position_summary.tsv"
        (cat dosage_header.txt && cat dose.sorted.tsv) | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
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
        File PrepVariantPositionSummary = "~{PhenotypeID}.prep_variant_position_summary.tsv"
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
        Int WindowSize = 1000000

    }

    call PrepInputs {
        input:
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
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
        File PrepVariantPositionSummary = PrepInputs.PrepVariantPositionSummary
    }
}
