# Export phenotype residual SD

The user wants `phenotype_residual_SD` in the SuSiE Parquet so effects can be converted to phenotype SD units without another fit.

Compute the sample SD of the exact residual phenotype vector passed to `fitSusieModel()`. Save the value as fitted-object metadata. Repeat it on the corresponding rows in the full and credible-set summary Parquets, including their aggregated versions. Keep the fitted means, posterior SDs, PIPs, credible sets, and model settings unchanged.

## Approaches

The conversation selected exporting the scale rather than refitting or changing effect units. Compute it at the common fit boundary in every genotype cache mode. Recovering it from model error variance is invalid. Older saved fits and Parquets use a numeric missing value; aggregation must not invent a scale.

## Constraints

- Use the exact name `phenotype_residual_SD` and a numeric column.
- Keep typed empty schemas and legacy saved-output support.
- Preserve record-batch streaming and strict ownership/schema validation.
- Reject supplied negative, infinite, nonnumeric, or inconsistent within-outcome scales. Zero is a valid recorded SD but cannot be used as a divisor.
- Keep existing WDL localization, task interfaces, and WDL 1.0 unchanged.
- Use small local tests and GitHub Actions quick image checks. Do not build Docker locally or submit Terra jobs.
