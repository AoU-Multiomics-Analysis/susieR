task susieR {

    File GenotypeDosages 
    File QTLCovariates 
    File TensorQTLPermutations
    File SampleList
    File PhenotypeBed
    Int CisDistance
    String OutputPrefix
    File susie_rscript
    Int memory

command{

    Rscript ${susie_rscript} \
        --genotype_matrix ${GenotypeDosages} \
        --sample_meta ${SampleList} \
        --phenotype_list ${TensorQTLPermutations} \
        --expresion_matrix ${PhenotypeBed} \
        --covariates ${QTLCovariates} \
        --out_prefix ${OutputPrefix} \
        --cisdistance ${CisDistance}

    } 

runtime {
        docker: 'quay.io/kfkf33/susier:v24.01.1'        
        memory: "${memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        cpu: "4"
        zones: ["us-central1-c"]
    }

output {
    File SusieParquet = "${OutputPrefix}.parquet" 
    File lbfParquet = "${OutputPrefix}.lbf_variable.parquet"
    File FullSusieParquet = "${OutputPrefix}.full_susie.parquet"

    }
}

workflow susieR_workflow {
    call susieR 
}
