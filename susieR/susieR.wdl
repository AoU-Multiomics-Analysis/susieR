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

task PrepInputs {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        String PhenotypeID
        File PhenotypeBed
        File TensorQTLPermutations
    }
    command <<<
    # grab headers from all files
    headerPermutations=$(zcat header=$(head -n 1 ~{TensorQTLPermutations}) | head -n 1)
    headerBed=$(zcat header=$(head -n 1 ~{PhenotypeBed}) | head -n 1)
    
    # create subset featured bed  
    zcat ~{PhenotypeBed} | grep ~{PhenotypeID} > feature.bed
    cat $headerBed feature.bed > ~{PhenotypeID}.bed
  
    # create subset tensorQTL file
    zcat ~{TensorQTLPermutations} | grep ~{PhenotypeID} > feature.txt
    cat $headerPermutations feature.txt > ~{PhenotypeID}.tensorQTL.txt
   
    # create subset dose file
    tabix ~{GenotypeDosages} ~{PhenotypeID}.bed > ~{PhenotypeID}.dose.tsv.gz
    tabix ~{PhenotypeID}.dose.tsv.gz 
    >>>
    
    runtime {
        docker: "quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
        disks: "local-disk 500 SSD"
        memory: "2GB"
        cpu: "1"
    }
    
    output {
        File SubsetBed = "~{PhenotypeID}.bed" 
        File SubsetPermutationPvals = "~{PhenotypeID}.tensorqtl.txt"
        File SubsetDosages = "~{PhenotypeID}.dose.tsv.gz"
        File SubsetDosagesIndex = "~{PhenotypeID}.dose.tsv.gz.tbi"
    }

}


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
        String PhenotypeID
    }

    command <<<
        Rscript ~{susie_rscript} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix ~{PhenotypeBed} \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance} \

    >>>

    runtime {
        docker: 'quay.io/kfkf33/susier:v24.01.1'
        memory: "${memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        cpu: "1"
    }

    output {
        File SusieParquet = "${OutputPrefix}.parquet"
        File lbfParquet = "${OutputPrefix}.lbf_variable.parquet"
        File FullSusieParquet = "${OutputPrefix}.full_susie.parquet"
    }
}


#task MergeSusie {
#    input {
#    Array[File] SusieOutput
#    Int memory
#    String OutputPrefix
    #}
#    
#    command <<<
#    for file in ~{sep='\n' SusieOutput}; do
#    echo $file >> filelist.txt
#    done
#
#    Rscript merge_susie.R \ 
#       --FilePaths filelist.txt \
#       --OutputPrefix ~{OutputPrefix}
#   >>>
#
#runtime {
#        docker: 'quay.io/kfkf33/susier:v24.01.1'
#        memory: "${memory}GB"
#        disks: "local-disk 500 SSD"
#        bootDiskSizeGb: 25
#        cpu: "1"
#    }
#
#
#    output {
#    File MergedSusieParquet = "${OutputPrefix}_SusieMerged.parquet" 
#    File MergedSusieTsv = "${OutputPrefix}_SusieMerged.tsv.gz" 
#
#    }
#
#}

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
        String OutputPrefix
    }

    call PrepInputs {
        input:
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
            GenotypeDosages = GenotypeDosages,
            GenotypeDosageIndex = GenotypeDosageIndex,
            PhenotypeBed = PhenotypeBed
    }

    call susieR {
        input:
            GenotypeDosages = PrepInputs.SubsetDosages,
            GenotypeDosageIndex = PrepInputs.SubsetDosagesIndex,
            QTLCovariates = QTLCovariates,
            TensorQTLPermutations = PrepInputs.SubsetPermutationPvals,
            SampleList = SampleList,
            PhenotypeBed = PrepInputs.SubsetBed ,
            CisDistance = CisDistance,
            OutputPrefix = PhenotypeID,
            susie_rscript = susie_rscript,
            memory = memory
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
    }
}
