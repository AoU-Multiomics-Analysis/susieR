version 1.0

task PrepInputs {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        String PhenotypeID
        File PhenotypeBed
        File TensorQTLPermutations
        Int NumPrempt
    }
    command <<<
        echo "Extracting headers from files"
        headerPermutations=$(zcat ~{TensorQTLPermutations} | head -n 1)
        headerBed=$(zcat "~{PhenotypeBed}" | head -n 1)
        
 
        echo "Bed file header:"
        echo $headerBed

        echo "TensorQTL file header"
        echo $headerPermutations
        
        #zcat ~{GenotypeDosages} |  awk 'NR==1 { if ($0 ~ /^#/) print; else print "#" $0; exit }'  > dosage_header.txt
        zcat ~{GenotypeDosages} |  awk 'NR==1 {print $0; exit }'  > dosage_header.txt

        echo "Subsetting bed file"
        zcat ~{PhenotypeBed} | grep "~{PhenotypeID}" \
            | awk 'BEGIN{OFS="\t"} {$2=$2-1000000; $3=$3+1000000; if($2<1) $2=1; print}' \
            > feature.bed
        
        #echo $headerBed > temp_header.txt
        zcat ~{PhenotypeBed} | head -n 1 > temp_header.txt
        cat temp_header.txt feature.bed | bgzip -c - > ~{PhenotypeID}.bed.bgz
        #tabix ~{PhenotypeID}.bed.bgz

        echo "Subsetting TensorQTL file"
        zcat ~{TensorQTLPermutations} | grep "~{PhenotypeID}" > feature.txt
        echo $headerPermutations > temp_header_perm.txt
        cat temp_header_perm.txt feature.txt > ~{PhenotypeID}.tensorqtl.txt

        echo "Subsetting dose file"
        #(cat dosage_header.txt; tabix ~{GenotypeDosages} -R ~{PhenotypeID}.bed.bgz) | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
        tabix ~{GenotypeDosages}  -R ~{PhenotypeID}.bed.bgz -h | bgzip -c > ~{PhenotypeID}.dose.tsv.gz


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


workflow susieR_workflow {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        File PhenotypeBed
        Int CisDistance
        Int NumPrempt
        String PhenotypeID
    }

    call PrepInputs {
        input:
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
            GenotypeDosages = GenotypeDosages,
            GenotypeDosageIndex = GenotypeDosageIndex,
            PhenotypeBed = PhenotypeBed,
            NumPrempt = NumPrempt
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
        File SubsetDosagesIndex = PrepInputs.SubsetDosagesIndex
    }
}
