# Common Inputs

These inputs are shared across the main fine-mapping workflows.

| Input | Description |
|---|---|
| `GenotypeDosages` | Genotype dosage file. Can be generated from a VCF using `bcftools dose`. |
| `GenotypeDosageIndex` | Tabix `.tbi` index for the dosage file. |
| `TensorQTLPermutations` | Permutation p-values output from tensorQTL. |
| `PhenotypeBed` | BED file for the gene or phenotype to be fine-mapped. The fine-mapping scripts use the midpoint of columns 2-3 as the feature coordinate. |
| `CisDistance` | Window in bp added to each side of the phenotype BED interval midpoint for fine-mapping. In the full workflow, this same value is used for prep extraction. |
| `WindowSize` | Standalone prep-only window in bp added to each side of the phenotype BED interval midpoint. |
| `PreparedWindowSize` | Optional fine-mapping-only guard for pre-subset dosage inputs. Set this to the prep `WindowSize` to fail early when `CisDistance` asks for variants outside the prepared dosage window. |
| `PhenotypeID` | Legacy single phenotype ID. When substring matching is used, this is used as the output prefix and gene ID to match within splice-junction phenotype IDs. |
| `MatchPhenotypeIDSubstring` | If `true`, select all phenotype IDs containing `PhenotypeID`. Use for splice-junction IDs that embed the gene ID. |
| `ReuseGenotypeMatrix` | If `true`, reuse one residualized genotype matrix when selected phenotype windows merge into a single region. |
| `SelectTopPhenotypePerCluster` | If `true`, reduce matched, FDR-passing LeafCutter phenotypes to the strongest intron per parsed `clu_*` cluster before fine-mapping. Defaults to `false` to preserve all-intron behavior. |
| `TopPhenotypePerClusterPvalueColumn` | Preferred TensorQTL p-value column for selecting the representative intron per cluster. Defaults to `qval`, with fallbacks in the code. |
| `QTLCovariates` | Covariate table used in QTL calling. Use the same file given to tensorQTL. |
| `SampleList` | List of sample IDs used in fine-mapping. Requires a header. |

## Optional Inputs

| Input | Description |
|---|---|
| `MAF` | Optional float MAF cutoff for variants. Defaults to `0`. Requires `AncestryMetadata` when greater than `0`. MAF is calculated per population and a variant must pass the cutoff in at least one population. Individuals not assigned to any population are excluded from the MAF calculation. |
| `AncestryMetadata` | Ancestry metadata file. Requires a column named `ancestry_pred_oth` for population assignment. |
| `VariantList` | Single-column file of variants formatted as `chr_pos_ref_alt`. Restricts analysis to listed variants and takes precedence over `MAF` filtering. |

## Checkpointed Window Inputs

`CheckpointedWindowSusieWorkflow` uses one prepared variant window. The
prepare workflow owns phenotype filtering. The analysis wrapper trusts the
prepared manifest. It does not read a global association table.

Set `window_phenotypes` to `window_phenotypes.tsv`. The file must contain one
row for each phenotype in the prepared window.

The current prepare workflows do not yet write this production manifest. This
is a deployment gap. A preparation step must create the file before production
use. That step must use the final linked and filtered phenotype set. It must
write one row for each phenotype. It must not add phenotypes that the analysis
wrapper must filter again.

| Column | Description |
|---|---|
| `window_id` | Prepared window identifier. It must equal the WDL `window_id`. |
| `phenotype_id` | Unique phenotype identifier in this window. |
| `modality` | Covariate modality for the phenotype. |
| `phenotype_file` | Prepared phenotype-data file name. All rows must use one value. Its base name must equal the `phenotype_data` base name. |
| `p_value` | Finite association P value from zero through one. The wrapper requires this column. |

The wrapper orders rows by `p_value`, `modality`, and `phenotype_id`. The first
version uses every usable variant in the supplied window. It does not calculate
a phenotype-centered interval.

Set `checkpoint_root` to a writable `gs://` prefix. Do not use a local path in
the WDL. Give parallel `covariate_files` and `covariate_modalities` arrays. Use
the label `shared` for covariates that apply to every phenotype.

Set `runner_image` to a published image digest. Use this form:

```text
ghcr.io/aou-multiomics-analysis/susier/checkpointed-window@sha256:<64-lowercase-hex-characters>
```

The WDL rejects tags such as `:main`. Replace
`REPLACE_WITH_PUBLISHED_DIGEST` in the example JSON. Do not use the pinned base
image digest as the runner-image digest. They identify different images.

### Prepared dosage schema

The dosage TSV must contain `CHROM` or `#CHROM`, `POS`, `REF`, and `ALT`. All
other columns are sample IDs. Sample IDs must be unique. Dosage values must be
numeric from zero through two. The exact value `-1` means missing. The runner
imputes missing dosage values inside each exact retained-sample cache.

### Prepared phenotype schema

The phenotype table must contain one chromosome column, `start`, `end`, and
`phenotype_id`. Accepted chromosome names are `chromosome`, `chrom`, `CHROM`,
`#CHROM`, `chr`, and `#chr`. All other columns are sample IDs. A phenotype can
contain a nonfinite value for one sample. The runner removes that sample only
for that phenotype. It does not remove the sample from other phenotypes.

### Covariate and allowlist schemas

Each covariate TSV is a covariate-by-sample table. The first column contains
unique covariate IDs. The other column names are unique sample IDs. All values
must be finite. The runner combines `shared` files and files for the phenotype
modality in WDL array order. It removes constant and dependent covariates for
each retained-sample mask.

The optional sample allowlist must have exactly one named column. Its sample
IDs must be unique and nonempty.

### GCS access

The task identity must be able to read, create, and replace objects under
`checkpoint_root`. Exact-path `stat` and download operations require object
read access. Cursor and manifest replacement also require the permissions that
the bucket uses for object replacement. The runner does not require a GCS
prefix listing. Use the smallest bucket or prefix scope that your platform
supports.
