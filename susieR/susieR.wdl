version 1.0

task splitPhenotypeBed {
    input {
        File TensorQTLPermutations
        #File PhenotypeBed
        #Int numSplits
    }

    #String baseName = basename(PhenotypeBed, ".gz")

    command <<<
        # Decompress the BED file
        gunzip -c ~{PhenotypeBed} > ~{baseName}
        
        # filter permutations to significant QTLs
        #zcat ~{TensorQTLPermutations} | awk '$18 < 0.05' | awk '{print $1}' > feature_list.txt
        zcat ~{TensorQTLPermutations} | awk '$18 < 0.05' | head -n 100  > significant_qtls.txt 
        awk 'NR==1 {header=$0; next} {out=$1".txt"; print header > out; print >> out}' significant_qtls.txt 
        # Extract the header line
        #header=$(head -n 1 ~{baseName})
        #header=$(head -n 1 ~{TensorQTLPermutations})


        # Get the total number of lines excluding the header
        #total_lines=$(wc -l < ~{baseName})
        #lines_per_file=2

        # Split the file into parts, excluding the header
        #tail -n +2 ~{baseName} | grep -Ff feature_list.txt    | split -l ${lines_per_file} -a 8 - ~{baseName}.part_
        #tail -n +2 significant_qtls.txt | grep -Ff feature_list.txt    | split -l ${lines_per_file} -a 8 - ~{baseName}.part_


        # Add the header to each split file and compress with bgzip
        #for file in ~{baseName}.part_*; do
        #    (echo "${header}" && cat "${file}") > "${file}.with_header"
        #    mv "${file}.with_header" "${file}"
        #    #bgzip "${file}"
        #done
    >>>

    #output {
    #    Array[File] splitFiles = glob("${baseName}.part_*")
   # }
    output {
        Array[File] split_files = glob("*.txt")
    }
    runtime {
        docker: "quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
        disks: "local-disk 500 SSD"
        memory: "2GB"
        cpu: "1"
    }
}

task susieR {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        File PhenotypeBedPart
        Int CisDistance
        String OutputPrefix
        File susie_rscript
        Int memory
    }

    command <<<
        Rscript ~{susie_rscript} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix ~{PhenotypeBedPart} \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance}
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


task MergeSusie {
    input {
    Array[File] SusieOutput
    Int memory
    String OutputPrefix
    }
    
    command <<<
    for file in ~{sep='\n' SusieOutput}; do
    echo $file >> filelist.txt
    done

    Rscript merge_susie.R \ 
       --FilePaths filelist.txt \
       --OutputPrefix ~{OutputPrefix}
   >>>

runtime {
        docker: 'quay.io/kfkf33/susier:v24.01.1'
        memory: "${memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        cpu: "1"
    }


    output {
    File MergedSusieParquet = "${OutputPrefix}_SusieMerged.parquet" 
    File MergedSusieTsv = "${OutputPrefix}_SusieMerged.tsv.gz" 

    }

}

workflow susieR_workflow {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        #File PhenotypeBed
        Int CisDistance
        File susie_rscript
        Int memory
        #Int numSplits
        String OutputPrefix
    }

    call splitPhenotypeBed {
        input:
            TensorQTLPermutations = TensorQTLPermutations,
            #numSplits = numSplits
    }

    scatter (phenotype_id in splitPhenotypeBed.splitFiles) {
        String phenotype_id = basename(phenotype_id, ".txt")
        call susieR {
            input:
                GenotypeDosages = GenotypeDosages,
                GenotypeDosageIndex = GenotypeDosageIndex,
                QTLCovariates = QTLCovariates,
                TensorQTLPermutations = part,
                SampleList = SampleList,
                PhenotypeBed = PhenotypeBed ,
                CisDistance = CisDistance,
                OutputPrefix = "~{phenotype_id}",
                susie_rscript = susie_rscript,
                memory = memory
        }
    }
    
    call MergeSusie {
        input:
            SusieOutput = susieR.SusieParquet,
            OutputPrefix = OutputPrefix 

    } 
    output {
        File SusieParquet = MergeSusie.MergedSusieParquet
        File SusieTsv = MergeSusie.MergedSusieTsv
        #Array[File] SusieParquets = susieR.SusieParquet
        #Array[File] lbfParquets = susieR.lbfParquet
        #Array[File] FullSusieParquets = susieR.FullSusieParquet
    }
}
