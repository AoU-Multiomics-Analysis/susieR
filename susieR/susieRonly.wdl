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
        File? VariantList
        File? AncestryFile
    }

    command <<<

        zcat ~{PhenotypeBed} | head -n 1 > header.txt
        zcat ~{PhenotypeBed} | grep ~{OutputPrefix} > input_gene.txt
        awk -F'\t' 'BEGIN{OFS="\t"} { $4="skip"; print }' input_gene.txt > skip.txt
        cat header.txt input_gene.txt skip.txt > input_gene.bed  

        Rscript /tmp/susie.R ~{if defined(MAF) then "--MAF ~{MAF}  " else ""} ~{if defined(AncestryFile) then "--AncestryMetadata ~{AncestryFile}  "  else ""} ~{if defined(VariantList) then "--VariantList ~{VariantList}  "  else ""} \
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
        Int memory
        Int NumPrempt
        String OutputPrefix
        Float? MAF
        File? VariantList
        File? AncestryFile 

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
            VariantList = VariantList,
            AncestryFile = AncestryFile
        }
        output {
            File SusieParquet = susieR.SusieParquet
            File SusielbfParquet = susieR.lbfParquet
            File FullSusieParquet = susieR.FullSusieParquet
        }
}
