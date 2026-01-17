version 1.0

#task splitPhenotypeBed {
#    input {
#        File TensorQTLPermutations
#    }

    #String baseName = basename(PhenotypeBed, ".gz")

#    command <<<
#        zcat ~{TensorQTLPermutations} | awk '$18 < 0.05' | head -n 100  > significant_qtls.txt 
#        awk 'NR==1 {header=$0; next} {out=$1".txt"; print header > out; print >> out}' significant_qtls.txt 
#    >>>
#
#    output {
#        Array[File] splitFiles = glob("*.txt")
#    }
#    runtime {
#        docker: "quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
#        disks: "local-disk 500 SSD"
#        memory: "2GB"
#        cpu: "1"
#    }
#}

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
        File? VariantList
    }

    command <<<
        Rscript ~{susie_rscript} \
            --MAF ~{MAF} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix ~{PhenotypeBed} \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            ~{if defined(VariantList) then "--VariantList {VariantList}  "     else ""} --cisdistance ~{CisDistance} 

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
        Float MAF
        File? VariantList

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
            susie_rscript = susie_rscript,
            memory = memory,
            NumPrempt = NumPrempt,
            MAF = MAF,
            VariantList = VariantList
        }
        output {
        File SusieParquet = susieR.SusieParquet
        File SusielbfParquet = susieR.lbfParquet
        File FullSusieParquet = susieR.FullSusieParquet
        }
}
