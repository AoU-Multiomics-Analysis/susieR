version 1.0

task AggregateSusie {
    input {
        File SusieParquetsFOFN
        Int Memory
        String OutputPrefix
        Int NumThreads
        #String AggregateMode # should be either lbf or pip
    }

    command <<<
    mkdir -p input_files
    mkdir -p localized

    export GSUTIL_PARALLEL_PROCESS_COUNT=32
    export GSUTIL_PARALLEL_THREAD_COUNT=8

    awk '{print $1}' ~{SusieParquetsFOFN} | grep -v '^$' > file_paths.txt

    echo "=== file_paths.txt count ==="
    wc -l file_paths.txt

    xargs -a file_paths.txt -n 100 sh -c 'gsutil -m cp "$@" localized/' sh

    echo "=== localized file count ==="
    find localized -maxdepth 1 -type f | wc -l

    find "$(pwd)/localized" -maxdepth 1 -type f | sort > filelist.txt
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

workflow AggregateSusieTaskWorkflow {
    input {
        File SusieParquetsFOFN
        Int Memory
        String OutputPrefix
        Int NumThreads
    }

    call AggregateSusie {
        input:
            SusieParquetsFOFN = SusieParquetsFOFN,
            OutputPrefix = OutputPrefix,
            Memory = Memory,
            NumThreads = NumThreads
    }

    output {
        File MergedSusieParquet = AggregateSusie.MergedSusieParquet
        File MergedSusieTsv = AggregateSusie.MergedSusieTsv
    }
}
