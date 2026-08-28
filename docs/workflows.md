# WDL Workflows

The WDL descriptors live in `workflows/`. Each descriptor has a unique workflow name so Dockstore, Terra, and Cromwell output namespaces are easier to distinguish.

| Descriptor | Workflow | Use |
|---|---|---|
| `workflows/susieR.wdl` | `SusieRWorkflow` | Prepare phenotype-specific inputs and run fine-mapping. |
| `workflows/susieRonly.wdl` | `SusieROnlyWorkflow` | Run fine-mapping when inputs are already prepared. |
| `workflows/ExtractMultiPhenotypeInputs.wdl` | `ExtractMultiPhenotypeInputsWorkflow` | Extract multi-phenotype BED and TensorQTL inputs without subsetting dosages. |
| `workflows/prepInputsSusieR.wdl` | `PrepSusieRWorkflow` | Prepare phenotype-specific input files only. |
| `workflows/CheckpointedWindowSusie.wdl` | `CheckpointedWindowSusieWorkflow` | Fine-map one prepared window with GCS checkpoints. |
| `workflows/ComputeR2Susie.wdl` | `ComputeR2SusieWorkflow` | Run cross-validation R2 evaluation. |
| `workflows/AggregateSusieTask.wdl` | `AggregateSusieTaskWorkflow` | Merge sharded Susie Parquet outputs. |
| `workflows/AnnotateSusie.wdl` | `AnnotateSusieWorkflow` | Annotate a merged Susie TSV. |
| `workflows/ComputeAncestrySkew.wdl` | `ComputeAncestrySkew` | Compute ancestry skew for annotated Susie variants. |
| `workflows/AggregateSusie.wdl` | `AggregateSusieWorkflow` | Run aggregate, annotate, and ancestry skew together. |

Template input JSONs live in [`../examples/inputs/`](../examples/inputs/). They use placeholder paths and values; replace those with workspace-specific files before submitting.

## Fine-Mapping Workflows

### `workflows/susieR.wdl` - Full Pipeline

Runs both input preparation and fine-mapping in a single workflow. First calls `PrepInputs` to subset all input files to the `CisDistance` region around the target phenotype, then calls `susieR` to perform fine-mapping on those subsetted files. Representative intron-per-cluster selection is applied inside `susie.R` during fine-mapping, not by the WDL prep shell code.

| Input | Type | Description |
|---|---|---|
| `GenotypeDosages` | File | Tabix-indexed genotype dosage file, which can be generated from VCF using `bcftools dose`. |
| `GenotypeDosageIndex` | File | Tabix `.tbi` index for the dosage file. |
| `TensorQTLPermutations` | File | Permutation p-values output from tensorQTL. |
| `PhenotypeBed` | File | BED file for the gene or phenotype to be fine-mapped. Prep and fine-mapping center windows on the midpoint of columns 2-3. |
| `CisDistance` | Int | Window size in bp added to each side of the phenotype BED interval midpoint for both prep extraction and fine-mapping. |
| `PhenotypeID` | String | Legacy single phenotype ID. When substring matching is used, this becomes the output prefix and gene ID to match within splice-junction phenotype IDs. |
| `MatchPhenotypeIDSubstring` | Boolean | If `true`, select all phenotype IDs containing `PhenotypeID`. This supports splice-junction IDs that embed the gene ID. |
| `ReuseGenotypeMatrix` | Boolean | If `true`, reuse one residualized genotype matrix when selected phenotype windows merge into a single region. |
| `SelectTopPhenotypePerCluster` | Boolean | If `true`, fine-map one representative FDR-passing intron per parsed LeafCutter `clu_*` cluster instead of every matched intron. |
| `TopPhenotypePerClusterPvalueColumn` | String | Preferred TensorQTL p-value column used to choose the representative intron. Defaults to `qval`. |
| `QTLCovariates` | File | Covariate table used in QTL calling. |
| `SampleList` | File | Sample IDs used in fine-mapping. Requires a header. |
| `memory` | Int | Memory in GB for the fine-mapping task. |
| `NumPrempt` | Int | Number of preemptible retries. |
| `MAF` | Float | Optional minor allele frequency cutoff. Defaults to `0`, which includes all variants. |

| Output | Description |
|---|---|
| `SusieParquet` | Variants in credible sets. |
| `SusielbfParquet` | Log-Bayes factors per credible set. |
| `FullSusieParquet` | Fine-mapping statistics for all tested variants. |
| `VariantPositionSummary` | TSV summary of tested variant counts, upstream/downstream counts, min/max variant distance from the feature, and requested cis-window size. |
| `SubsetBed` | Subsetted BED file for the phenotype region. |
| `SubsetDosages` | Subsetted genotype dosage file. |
| `SubsetDosagesIndex` | Index for the subsetted dosage file. |
| `PrepVariantPositionSummary` | TSV summary of de-duplicated dosage variants extracted during prep, with per-feature upstream/downstream counts and requested prep-window size. |

### `workflows/susieRonly.wdl` - Fine-Mapping Only

Runs only the fine-mapping step. Use this when inputs are already prepared, such as from `prepInputsSusieR.wdl`, or when re-running fine-mapping with different parameters.

Optional inputs in addition to the shared fine-mapping inputs:

| Input | Type | Description |
|---|---|---|
| `MAF` | Float? | Minor allele frequency cutoff. |
| `MatchPhenotypeIDSubstring` | Boolean | If `true`, select all phenotype IDs containing `OutputPrefix`. |
| `ReuseGenotypeMatrix` | Boolean | If `true`, reuse one residualized genotype matrix when selected phenotype windows merge into a single region. |
| `SelectTopPhenotypePerCluster` | Boolean | If `true`, reduce selected phenotypes to the strongest FDR-passing intron per parsed LeafCutter `clu_*` cluster. |
| `TopPhenotypePerClusterPvalueColumn` | String | Preferred TensorQTL p-value column used to choose the representative intron. |
| `PreparedWindowSize` | Int | Optional guard for pre-subset dosage inputs. Set to the prep `WindowSize`; values below `CisDistance` stop the task before fine-mapping. Defaults to `-1`, which disables the shell guard. |
| `VariantList` | File? | Single-column file of variants formatted as `chr_pos_ref_alt` to restrict analysis. |
| `AncestryFile` | File? | Ancestry metadata for per-population MAF filtering. |
| `AdditionalGenotypesBed` | File? | Additional genotype BED file. |

Outputs are `SusieParquet`, `SusielbfParquet`, `FullSusieParquet`, and `VariantPositionSummary`.

### `workflows/CheckpointedWindowSusie.wdl` - Checkpointed Prepared Window

This workflow runs one task for one prepared variant window. Terra launches one workflow per prepared window. Do not use a WDL scatter for this workflow.

The prepare workflow owns phenotype filtering. The analysis task trusts
`window_phenotypes.tsv` as the prepared manifest. The manifest requires
`window_id`, `phenotype_id`, `modality`, `phenotype_file`, and `p_value`. The
task orders phenotypes by `p_value`, modality, and phenotype ID.

The current prepare workflows do not yet create this production manifest. This
is a deployment gap. Before production use, add a preparation adapter that
writes the final linked and filtered phenotype set. The adapter must write the
five required columns. It must use one `phenotype_file` value that matches the
`phenotype_data` file base name.

The first version tests every usable variant in the supplied window. It does
not calculate a phenotype-centered interval. Set `checkpoint_root` to a
writable `gs://` prefix. One active writer may use an analysis ID and window ID
pair.

Set `runner_image` to an immutable published `@sha256` digest. The WDL has no
runner-image default. It rejects a mutable tag before it starts R. The runtime
hashes the installed runner wrapper and both checkpoint helper files. The
analysis ID changes when one of these source files changes. Window and
phenotype manifests record these hashes, the runner-image identity, and the
pinned base-image identity.

The controller first aligns dosage, phenotype headers, and the optional sample
allowlist. For each phenotype, it then aligns shared covariates and covariates
for that phenotype's modality. It applies the finite-value phenotype mask
after this alignment. A gap in one modality-specific covariate file does not
remove a sample from other modalities. The controller uses an exact cache for
each retained-sample mask and covariate design. The cache identity covers
ordered overlap and retained sample IDs, retained variants, retained
covariate columns, and design values. The controller imputes, filters,
transposes, and residualizes genotype once for each cache. It then transforms
and residualizes only the phenotype before each SuSiE fit.

The controller releases a cache after its last ordered phenotype. It releases
the raw genotype after it prepares the last distinct cache. Phenotype order
does not change to reduce memory. Interleaved cache groups can therefore be in
memory at the same time. As a first estimate, allow at least eight bytes for
each sample-by-retained-variant value in each active cache, plus R and SuSiE
overhead. Increase `memory_gb` when a window has many distinct retained-sample
masks or covariate designs.

A prior IKZF1 prepared-window run completed in 1307.44 seconds and used 8.01
GB. This result is evidence from one input. It is not a resource guarantee.

Each committed phenotype has an RDS with the full fitted `susie` object. GCS
payloads upload before the phenotype manifest. The task uploads
`window_manifest.json` last. Normal resume reads only the window cursor and
the latest fit RDS. Light hydration validates every phenotype manifest and
checks that each saved fit exists, but it does not download every prior RDS.
If the cursor is missing, the task probes fixed manifest paths. It does not
list the full GCS prefix.

If boundary validation fails, the task does not delete the corrupt object. It
records the index, phenotype ID, phenotype key, fixed manifest path, concise
reason, timestamp, and attempt in `recovery_history`. It then recomputes that
phenotype. Later attempts preserve this history.

Final assembly validates the checksum, byte count, and table schema of all
three phenotype Parquet payloads. If a deterministic interior payload error is
found, the controller copies the exact invalid object to a deterministic
`recovery_evidence` path before it changes the cursor or recomputes a fit. The
recovery record contains the source path, evidence path, and evidence URI. If
the source object is missing, the record marks it as missing and does not
create evidence. The controller then recomputes that phenotype and its suffix
in the same invocation. It does not delete the source or evidence object. GCS
authentication, DNS, quota, copy, and service errors stop the task; they do
not cause recomputation.

The WDL returns `window_manifest.json`, `window_fit_index.tsv`, and three
Parquet tables. `window_fit_index.tsv` identifies each fitted RDS.
`SKIPPED` and `NONCONVERGED` are terminal phenotype states. They do not block
later phenotypes. Unexpected failures preserve the last valid cursor and fail
the task.

The checkpoint root uses this layout:

```text
gs://<bucket>/<checkpoint-root>/<analysis-id>/<window-id>/
  window_manifest.json
  window_fit_index.tsv
  results/
  phenotypes/<phenotype-key>/
    fit_manifest.json
    susie_fit.<sha256>.rds
```

Start from
[`examples/inputs/CheckpointedWindowSusie.inputs.json`](../examples/inputs/CheckpointedWindowSusie.inputs.json).
Replace every example path with a Terra workspace path. Submit one input JSON
for each prepared window. Replace the runner-image placeholder with the digest
of a published checkpointed-window image.

### `workflows/ExtractMultiPhenotypeInputs.wdl` - Phenotype Extraction Only

Extracts phenotype rows and matching TensorQTL permutation rows for multi-phenotype fine-mapping cases, without touching genotype dosages. Use this to materialize all splice-junction phenotypes for a gene before running fine-mapping against an already-prepared dosage file. Representative intron-per-cluster selection is handled by the fine-mapping workflows, not this extraction helper.

| Input | Type | Description |
|---|---|---|
| `PhenotypeBed` | File | BED file containing phenotype rows. |
| `TensorQTLPermutations` | File | TensorQTL permutation output. |
| `PhenotypeID` | String | Legacy single phenotype ID, or gene/output prefix when selecting multiple embedded phenotype IDs. |
| `MatchPhenotypeIDSubstring` | Boolean | If `true`, select all phenotype IDs containing `PhenotypeID`. |
| `AddSkipRow` | Boolean | Add the legacy `skip` row used by `susie.R` phenotype matrices. Defaults to `true`. |
| `NumPrempt` | Int | Number of preemptible retries. |

| Output | Description |
|---|---|
| `PhenotypeMatrix` | Uncompressed BED-style phenotype matrix containing the header, matched phenotype rows, and optional `skip` row. |
| `SubsetPhenotypeBed` | Gzipped BED containing the header and matched phenotype rows only. |
| `SubsetPermutationPvals` | TensorQTL permutation rows for the matched phenotype IDs. |
| `MatchedPhenotypeIDs` | Exact phenotype IDs discovered from the phenotype BED. |

### `workflows/prepInputsSusieR.wdl` - Input Preparation Only

Subsets genotype dosages, the phenotype BED file, and TensorQTL permutation results to the region surrounding `PhenotypeID`, or all phenotype IDs containing `PhenotypeID` when `MatchPhenotypeIDSubstring` is true. Prep windows are centered on the midpoint of each matched phenotype BED interval before adding `WindowSize` on each side. This is the first step of `workflows/susieR.wdl` exposed as a standalone workflow, useful for preparing inputs once before running fine-mapping multiple times. It always prepares all matched rows; representative intron-per-cluster selection happens later in fine-mapping.

| Output | Description |
|---|---|
| `SubsetBed` | Subsetted BED file for the phenotype region. |
| `SubsetDosages` | Subsetted genotype dosage file. |
| `SubsetDosagesIndex` | Tabix index for the subsetted dosage file. |
| `PrepVariantPositionSummary` | TSV summary of de-duplicated dosage variants extracted during prep, with per-feature upstream/downstream counts and requested prep-window size. |

### `workflows/ComputeR2Susie.wdl` - Cross-Validation R2

Runs the `ComputeR2Susie.R` script to evaluate fine-mapping predictive R2 via cross-validation. Requires a `CVMetadata` file generated by `PrepSusieCVPCs.R`.

| Input | Type | Description |
|---|---|---|
| `CVMetadata` | File | Cross-validation fold and PC metadata generated by `PrepSusieCVPCs.R`. |
| `VariantList` | File | Variants to include in the 1% AF threshold analysis. |

Output: `<PhenotypeID>_SusiePredictions.tsv`, with per-fold predicted and observed expression values.

## Aggregate and Annotation Workflows

The post-fine-mapping workflows are split into standalone WDLs so each step can be launched independently or as one combined workflow through Dockstore or Terra.

| Descriptor | Workflow | Purpose |
|---|---|---|
| `workflows/AggregateSusieTask.wdl` | `AggregateSusieTaskWorkflow` | Localizes sharded Susie Parquet outputs from a FOFN and merges them into one Parquet plus one gzipped TSV. |
| `workflows/AnnotateSusie.wdl` | `AnnotateSusieWorkflow` | Annotates an existing merged Susie TSV with GENCODE, ENCODE, FANTOM5, gnomAD constraint, phyloP, and VAT data. |
| `workflows/ComputeAncestrySkew.wdl` | `ComputeAncestrySkew` | Computes ancestry skew on an annotated Susie TSV. |
| `workflows/AggregateSusie.wdl` | `AggregateSusieWorkflow` | Runs aggregate, annotate, and ancestry-skew analysis together by importing the three standalone post-analysis WDLs. |

### Aggregate-Only Inputs

| Input | Type | Description |
|---|---|---|
| `SusieParquetsFOFN` | File | Text file listing Susie Parquet shard paths to merge. |
| `OutputPrefix` | String | Prefix for merged output files. |
| `Memory` | Int | Memory in GB for the aggregate task. |
| `NumThreads` | Int | CPU threads for the aggregate task. |

### Aggregate-Only Outputs

| Output | Description |
|---|---|
| `MergedSusieParquet` | `<OutputPrefix>_SusieMerged.parquet`. |
| `MergedSusieTsv` | `<OutputPrefix>_SusieMerged.tsv.gz`. |

### Annotate-Only Inputs

| Input | Type | Description |
|---|---|---|
| `SusieTSV` | File | Merged Susie TSV, usually from `AggregateSusieTaskWorkflow.MergedSusieTsv`. |
| `GencodeGTF` | File | GENCODE GTF annotation file. |
| `AnnotationENCODE` | File | ENCODE cCRE annotation file. |
| `AnnotationFANTOM5` | File | FANTOM5 annotation file. |
| `AnnotationGnomad` | File | gnomAD constraint annotation file. |
| `AnnotationPhyloP` | File | phyloP bigWig annotation file. |
| `VATData` | File | VAT annotation data. |
| `AllelicFoldChangeData` | File? | Optional allelic fold-change table with `pid`, `sid`, `log2_aFC`, `log2_aFC_lower`, and `log2_aFC_upper`; `sid` may be formatted as `chr:pos_ref_alt`. |
| `OutputPrefix` | String | Prefix for annotated output. |
| `Memory` | Int | Memory in GB for the annotation task. |

Annotate-only output: `AnnotatedSusieTsv`, written as `<OutputPrefix>_SusieMerged.annotated.tsv`.

The combined `AggregateSusieWorkflow` keeps the same aggregate and annotation inputs, adds optional ancestry-skew controls (`AncestrySkewVariantsPerShard`, `AncestrySkewPipThreshold`, and `AncestrySkewAdmixedSubpops`), and returns both the ancestry-skew annotated TSV and the annotation-only TSV. When `AllelicFoldChangeData` is supplied, annotation appends `log2_aFC`, `log2_aFC_lower`, and `log2_aFC_upper` by matching normalized `gene_id`/`pid` and `variant`/`sid` keys.

### Ancestry Skew Inputs

| Input | Type | Default | Description |
|---|---|---|---|
| `AnnotationData` | File | required | Annotated variant table. Plain TSV and gzip-compressed TSV inputs are supported. |
| `OutputFile` | String | required | Name for the aggregated output file, usually ending in `.tsv.gz`. |
| `VariantsPerShard` | Int | required | Number of input rows per shard. |
| `PipThreshold` | Float | `0.9` | Variants with `pip >= PipThreshold` are included. |
| `AdmixedSubpops` | String | `"oth"` | Comma-separated GVS subpopulation labels to remove for the no-admixed skew calculation. |
| `KeepInputColumns` | Boolean | `false` | Keep all input annotation columns and append ancestry skew columns. |

### Ancestry Skew Outputs

The aggregated output is a gzip-compressed TSV with one row per variant passing the PIP threshold. In standalone mode, the output keeps identifiers plus `gvs_max_subpop`, recomputed MAF-based `gvs_max_af`, `gvs_odds_ratio`, `gvs_p_value`, `gvs_no_admixed_odds_ratio`, and `gvs_no_admixed_p_value`. When `KeepInputColumns` is true, those compact ancestry-skew columns are appended to the full input annotation table. Fisher tests add a pseudocount of one to each valid count cell before estimating odds ratios and p-values.

## Validation

Run local descriptor validation with:

```bash
tools/validate_repo.sh
```
