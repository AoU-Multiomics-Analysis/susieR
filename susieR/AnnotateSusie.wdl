version 1.0

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
        File AnnotatedSusieTsv = "~{OutputPrefix}_SusieMerged.annotated.tsv"
    }
}

workflow AnnotateSusieWorkflow {
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

    call AnnotateSusie {
        input:
            SusieTSV = SusieTSV,
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
        File AnnotatedSusieTsv = AnnotateSusie.AnnotatedSusieTsv
    }
}
