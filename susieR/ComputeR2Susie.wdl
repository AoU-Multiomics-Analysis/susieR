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
        Int n_folds 
        Float train_test_split
        File AncestryMetadata
    }

    command <<<
        Rscript /tmp/ComputeR2Susie.R  \
            --AncestryMetadata ~{AncestryMetadata} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix ~{PhenotypeBed} \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance} \
            --train_test_split ~{train_test_split} \
            --n_folds ~{n_folds}

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
        Int n_folds
        Float train_test_split
        Int NumPrempt
        String PhenotypeID
        File AncestryMetadata 
    }

    call ComputeR2 {
        input:
            AncestryMetadata = AncestryMetadata,
            GenotypeDosages = GenotypeDosages,
            train_test_split = train_test_split,
            n_folds = n_folds,
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
