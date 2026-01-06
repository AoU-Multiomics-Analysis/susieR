version 1.0

task RunSusieR {
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
        File AncestryData
    }

    command <<<
        Rscript ~{susie_rscript} \
            --AncestryMetadata ~{AncestryData} \
            --MAF ~{MAF} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix ~{PhenotypeBed} \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance} \

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
        File MAFTsv = "MAF.tsv"
    }
}


workflow susieR {
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
        File AncestryData
    }
    call RunSusieR {
        input:
            GenotypeDosages = GenotypeDosages,
            GenotypeDosageIndex = GenotypeDosageIndex,
            QTLCovariates = QTLCovariates,
            TensorQTLPermutations = TensorQTLPermutations,
            SampleList = SampleList,
            PhenotypeBed = PhenotypeBed ,
            CisDistance = CisDistance,
            OutputPrefix = OutputPrefix,
            susie_rscript = susie_rscript,
            memory = memory,
            NumPrempt = NumPrempt,
            MAF = MAF,
            AncestryData = AncestryData
        }
 
    output {
        File SusieParquet = RunSusieR.SusieParquet         
        File lbfParquet = RunSusieR.lbfParquet         
        File FullSusieParquet = RunSusieR.FullSusieParquet
        File MAFTsv = RunSusieR.MAFTsv
    }

}
