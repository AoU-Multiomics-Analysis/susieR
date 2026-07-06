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
