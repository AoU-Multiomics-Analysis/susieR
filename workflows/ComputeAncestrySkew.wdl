version 1.0

workflow ComputeAncestrySkew {
    input {
        File AnnotationData
        String OutputFile
        Int VariantsPerShard
        Float PipThreshold = 0.9
        String AdmixedSubpops = "oth"
        Boolean KeepInputColumns = false
    }
    
    call ShardVariants {
        input:
            AnnotationData = AnnotationData,
            rows_per_shard = VariantsPerShard
    }

    scatter (shard in ShardVariants.shard_files) {
        call GetAncestrySkew {
            input:
                AnnotationData = shard,
                PipThreshold = PipThreshold,
                AdmixedSubpops = AdmixedSubpops,
                KeepInputColumns = KeepInputColumns
        }
    }
    
    call AggregateAncestrySkew {
        input:
            shard_outputs =  GetAncestrySkew.AncestrySkewOutput,
            OutputFile = OutputFile
    }
    output {
        File Output = AggregateAncestrySkew.AggregatedOutput
    }
}


task ShardVariants {
    input {
        File AnnotationData
        Int rows_per_shard
    }
    command <<<
      python3 <<PY
    import csv
    import json
    import os
    import gzip

    infile = "~{AnnotationData}"
    rows_per_shard = ~{rows_per_shard}

    os.makedirs("shards", exist_ok=True)

    def is_gzip(path):
        with open(path, "rb") as probe:
            return probe.read(2) == b"\x1f\x8b"

    # Open plain TSV or gzip-compressed TSV, even if the suffix is not .gz.
    if is_gzip(infile):
        fh = gzip.open(infile, "rt")
    else:
        fh = open(infile, "r")

    reader = csv.DictReader(fh, delimiter="\t")
    header = reader.fieldnames
    if header is None:
        raise ValueError("Input file has no header")

    shard_idx = 0
    rows_in_current_shard = 0
    out = None
    writer = None
    shard_paths = []

    def open_shard(idx):
        path = f"shards/shard_{idx:04d}.tsv"
        handle = open(path, "w", newline="")
        w = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        w.writeheader()
        return path, handle, w

    for row in reader:
        if writer is None or rows_in_current_shard >= rows_per_shard:
            if out is not None:
                out.close()
            shard_path, out, writer = open_shard(shard_idx)
            shard_paths.append(shard_path)
            shard_idx += 1
            rows_in_current_shard = 0

        writer.writerow(row)
        rows_in_current_shard += 1

    fh.close()

    if out is not None:
        out.close()

    with open("shard_manifest.json", "w") as f:
        json.dump(shard_paths, f)
    PY
    >>>    
    output {
      Array[File] shard_files = glob("shards/*.tsv")
      File manifest = "shard_manifest.json"
    }
    runtime {
        docker: "python:3.10"
        memory: "32G"
        cpu: 1
        disks: "local-disk 2500 SSD"
  }
}

task GetAncestrySkew {
    input {
        File AnnotationData
        Float PipThreshold
        String AdmixedSubpops
        Boolean KeepInputColumns
    }
    String shard_base = basename(AnnotationData, ".tsv")

    command <<< 

    Rscript /ComputeAncestrySkew.R \
       --AnnotationData ~{AnnotationData} \
       --OutputPrefix ~{shard_base} \
       --PipThreshold ~{PipThreshold} \
       --AdmixedSubpops "~{AdmixedSubpops}" \
       ~{if KeepInputColumns then "--KeepInputColumns" else ""}
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/susier/postanalysis:main"
        memory: "96G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
    
    output {
        File AncestrySkewOutput = shard_base + ".AncestrySkew.tsv.gz" 
    }
}

task AggregateAncestrySkew {

    input {
        Array[File] shard_outputs
        String OutputFile
    }

    command <<<
    set -euo pipefail

    first=1

    for f in ~{sep=' ' shard_outputs}; do
        if [ $first -eq 1 ]; then
            zcat "$f"
            first=0
        else
            zcat "$f" | tail -n +2
        fi
    done | gzip > "~{OutputFile}"
    >>>
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/susier/postanalysis:main"
        memory: "32G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
    
    output {
        File AggregatedOutput = OutputFile
    }
}
