version 1.0

import "AggregateSusieTask.wdl" as AggregateSusieTask
import "AnnotateSusie.wdl" as AnnotateSusieTask
import "ComputeAncestrySkew.wdl" as AncestrySkew

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
        Int AncestrySkewVariantsPerShard = 10000
        Float AncestrySkewPipThreshold = 0.9
        String AncestrySkewAdmixedSubpops = "oth"
    }
    
    call AggregateSusieTask.AggregateSusie as AggregateSusie {
        input:
            SusieParquetsFOFN = SusieParquetsFOFN,
            OutputPrefix = OutputPrefix,
            Memory = Memory,
            NumThreads = NumThreads
    }

    call AnnotateSusieTask.AnnotateSusie as AnnotateSusie {
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

    call AncestrySkew.ComputeAncestrySkew as ComputeAncestrySkew {
        input:
            AnnotationData = AnnotateSusie.AnnotatedSusieTsv,
            OutputFile = OutputPrefix + "_SusieMerged.annotated.ancestry_skew.tsv.gz",
            VariantsPerShard = AncestrySkewVariantsPerShard,
            PipThreshold = AncestrySkewPipThreshold,
            AdmixedSubpops = AncestrySkewAdmixedSubpops,
            KeepInputColumns = true
    }

    output {
        File AnnotatedMergedSusieTsv = ComputeAncestrySkew.Output
        File AnnotatedMergedSusieWithoutAncestrySkewTsv = AnnotateSusie.AnnotatedSusieTsv
        #File MergedSusieParquet = AggregateSusie.MergedSusieParquet
        #File MergedSusieTsv = AggregateSusie.MergedSusieTsv
    }

}
