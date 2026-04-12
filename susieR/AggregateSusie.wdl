version 1.0



task AggregateSusie{
    input{
        File SusieParquetsFOFN
        Int Memory
        String OutputPrefix
        Int NumThreads 
        #String AggregateMode # should be either lbf or pip
    }

    command <<<
    export GSUTIL_PARALLEL_PROCESS_COUNT=32
    export GSUTIL_PARALLEL_THREAD_COUNT=8

    awk '{print $1}' ~{SusieParquetsFOFN} | grep -v '^$' > file_paths.txt 

    mkdir -p localized
    gsutil -m cp -I localized/ < file_paths.txt 

    # Write the new local file paths into filelist.txt
    ls -1 "$(pwd)/localized/"* > filelist.txt
    Rscript /MergeSusie.R --FilePaths filelist.txt  --OutputPrefix ~{OutputPrefix} 
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/susier/postanalysis:main"
        disks: "local-disk 500 SSD"
        memory: "~{Memory}GB"
        cpu: "~{NumThreads}"
    }
     output {
        File MergedSusieParquet = "${OutputPrefix}_SusieMerged.parquet" 
        File MergedSusieTsv = "${OutputPrefix}_SusieMerged.tsv.gz" 
    } 

}

task AnnotateSusie {
    input {
        File SusieTSV 
        File GencodeGTF
        String OutputPrefix
        Int Memory
        File AnnotationENCODE 
        File AnnotationFANTOM5 
        File AnnotationGnomad 
        File AnnotationPhyloP 
        File VATData
    }
    command <<<
    Rscript /AnnotateSusie.R \
        --OutputPrefix ~{OutputPrefix} \
        --GencodeGTF ~{GencodeGTF} \
        --SusieTSV ~{SusieTSV} \
        --phyloPBigWig ~{AnnotationPhyloP} \
        --FANTOM5 ~{AnnotationFANTOM5} \
        --gnomadConstraint ~{AnnotationGnomad} \
        --ENCODEcCRES ~{AnnotationENCODE} \
        --VAT ~{VATData}
    >>>
   runtime {
        docker: "ghcr.io/aou-multiomics-analysis/susier/postanalysis:main"
        disks: "local-disk 500 SSD"
        memory: "~{Memory}GB"
        cpu: "1"
    }


    output {
        File AnnotatedSusieParquetOut = "~{OutputPrefix}_SusieMerged.annotated.tsv" 
    }

}



workflow AggregateSusieWorkflow {
    input {
        File SusieParquetsFOFN
        Int Memory 
        String OutputPrefix
        Int NumThreads
        File GencodeGTF 
        File AnnotationPhyloP 
        File AnnotationENCODE 
        File AnnotationFANTOM5 
        File AnnotationGnomad
        File VATData  
    }
    
    call AggregateSusie {
        input:
            SusieParquetsFOFN = SusieParquetsFOFN,
            OutputPrefix = OutputPrefix,
            Memory = Memory,
            NumThreads = NumThreads
    }

    call AnnotateSusie {
        input:
            SusieTSV = AggregateSusie.MergedSusieTsv,
            GencodeGTF = GencodeGTF,
            OutputPrefix = OutputPrefix,
            Memory = Memory,
            AnnotationPhyloP = AnnotationPhyloP,
            AnnotationENCODE = AnnotationENCODE,
            AnnotationFANTOM5 = AnnotationFANTOM5,
            AnnotationGnomad = AnnotationGnomad,
            VATData = VATData
    } 

    output {
        File AnnotatedMergedSusieParquet = AnnotateSusie.AnnotatedSusieParquetOut
        #File MergedSusieParquet = AggregateSusie.MergedSusieParquet
        #File MergedSusieTsv = AggregateSusie.MergedSusieTsv
    }

}
