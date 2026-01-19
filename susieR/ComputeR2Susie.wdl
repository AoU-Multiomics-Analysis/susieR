version 1.0


task ComputeR2 {
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
        File VariantList
        File CVMetadata
    }

    command <<<
        Rscript /tmp/ComputeR2Susie.R  \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix ~{PhenotypeBed} \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance} \
            --cv_meta ~{CVMetadata} \
            --variant_list ~{VariantList}

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
        File SusieR2 = "${OutputPrefix}_R2_CV.tsv"
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
        String PhenotypeID
        File VariantList 
        File CVMetadata
    }

    call ComputeR2 {
        input:
            GenotypeDosages = GenotypeDosages,
            CVMetadata = CVMetadata,
            VariantList = VariantList,
            GenotypeDosageIndex = GenotypeDosageIndex,
            QTLCovariates = QTLCovariates,
            TensorQTLPermutations = TensorQTLPermutations,
            SampleList = SampleList,
            PhenotypeBed = PhenotypeBed ,
            CisDistance = CisDistance,
            OutputPrefix = PhenotypeID,
            memory = memory,
            NumPrempt = NumPrempt

        }
    
    output {
            File R2Susie =  ComputeR2.SusieR2
            }
}
