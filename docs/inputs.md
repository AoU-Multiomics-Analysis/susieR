# Common Inputs

These inputs are shared across the main fine-mapping workflows.

| Input | Description |
|---|---|
| `GenotypeDosages` | Genotype dosage file. Can be generated from a VCF using `bcftools dose`. |
| `GenotypeDosageIndex` | Tabix `.tbi` index for the dosage file. |
| `TensorQTLPermutations` | Permutation p-values output from tensorQTL. |
| `PhenotypeBed` | BED file for the gene or phenotype to be fine-mapped. |
| `CisDistance` | Window in bp added to each side of the TSS for fine-mapping. |
| `PhenotypeID` | Legacy single phenotype ID. When substring matching is used, this is used as the output prefix and gene ID to match within splice-junction phenotype IDs. |
| `MatchPhenotypeIDSubstring` | If `true`, select all phenotype IDs containing `PhenotypeID`. Use for splice-junction IDs that embed the gene ID. |
| `ReuseGenotypeMatrix` | If `true`, reuse one residualized genotype matrix when selected phenotype windows merge into a single region. |
| `QTLCovariates` | Covariate table used in QTL calling. Use the same file given to tensorQTL. |
| `SampleList` | List of sample IDs used in fine-mapping. Requires a header. |

## Optional Inputs

| Input | Description |
|---|---|
| `MAF` | Float MAF cutoff for variants. Requires `AncestryMetadata`. MAF is calculated per population and a variant must pass the cutoff in at least one population. Individuals not assigned to any population are excluded from the MAF calculation. |
| `AncestryMetadata` | Ancestry metadata file. Requires a column named `ancestry_pred_oth` for population assignment. |
| `VariantList` | Single-column file of variants formatted as `chr_pos_ref_alt`. Restricts analysis to listed variants and takes precedence over `MAF` filtering. |
