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
        File susie_rscript
        Int memory
        Int NumPrempt
        Float MAF
    }

    command <<<
        Rscript ~{susie_rscript} \
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
        docker: 'quay.io/kfkf33/susier:v24.01.2'
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
        File susie_rscript
        Int memory
        Int NumPrempt
        String OutputPrefix
        String PhenotypeID
        Float MAF
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

    call susieR {
        input:
            GenotypeDosages = PrepInputs.SubsetDosages,
            GenotypeDosageIndex = PrepInputs.SubsetDosagesIndex,
            QTLCovariates = QTLCovariates,
            TensorQTLPermutations = PrepInputs.SubsetPermutationPvals,
            SampleList = SampleList,
            PhenotypeBed = PhenotypeBed ,
            CisDistance = CisDistance,
            OutputPrefix = PhenotypeID,
            susie_rscript = susie_rscript,
            memory = memory,
            NumPrempt = NumPrempt,
            MAF = MAF

        }
    
    #call MergeSusie {
    #    input:
    #        SusieOutput = susieR.SusieParquet,
    #        OutputPrefix = OutputPrefix 
    #
    #} 
    output {
        File SusieParquet = susieR.SusieParquet
        File SusielbfParquet = susieR.lbfParquet
        File FullSusieParquet = susieR.FullSusieParquet
        File SubsetBed = PrepInputs.SubsetBed
        #File SubsetBedIndex = PrepInputs.SubsetBedIndex
        File SubsetDosages = PrepInputs.SubsetDosages
        File SubsetDosagesIndex = PrepInputs.SubsetDosagesIndex
    }
}
