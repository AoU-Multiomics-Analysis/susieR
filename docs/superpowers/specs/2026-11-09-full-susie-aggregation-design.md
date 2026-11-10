# Full SuSiE aggregation

The user requested one aggregated version of the existing per-gene full SuSiE Parquet outputs, added to the PGEN workflow. This must support a later SV/SNV effect comparison without refitting models. Keep every source column, including native PIPs, posterior moments, X_column_scale_factors, alpha, mu, and mu2. Preserve rows outside credible sets and their missing CS fields.

Use a separate streaming helper after the existing compact aggregation. Read only explicitly localized files, matched by exact gene-key basename against the validated outcome index. Add ownership columns. Validate schemas, gene hashes, outcome ownership, variant identities, positive scale factors, and duplicates. Allow the known zero-row schema differences. Some source files are empty when no credible set was retained; aggregation must not reconstruct absent rows.

Write a single Parquet file through bounded Arrow batches. Do not bind all genes into one R table. Write a temporary output, close it, and rename it only after validation succeeds. The full matrix can be much larger than the compact summary, so aggregation local disk must cover localized inputs and output.

The main WDL passes Array[File] inputs and creates file lists during task command rendering. URI metadata for output manifests remains separate. Expose AggregateFullSuSiEParquet while keeping the compact outputs. Add an aggregation-only WDL and explicit CLI mode to reuse completed full outputs and their outcome index without PLINK or SuSiE calls.

No phenotype residual SD field is added in this request. Full outputs support genotype-scaled effects, but a comparison standardized on both genotype and phenotype still needs the phenotype SD after covariate adjustment. Neither the compact nor current full schema records it. The fitted model error SD is not a substitute.

Use small test fixtures, real Arrow streaming, static WDL checks, and command localization checks. Test the image in quick GitHub Actions mode only. Do not build Docker locally, run a full scientific workflow, or submit a Terra job. State that the complete new workflows have not run on Terra.
