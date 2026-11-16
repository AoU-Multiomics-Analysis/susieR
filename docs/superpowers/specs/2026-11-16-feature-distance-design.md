# Distance from a fine-mapped variant to its feature

Add numeric `distance_to_feature_bp` to the compact per-gene, batch, and combined
variant Parquets. Use variant position minus the existing feature midpoint
(`phenotype_pos`), in base pairs, when both are on the same chromosome. Negative
values indicate a lower genomic coordinate. Do not use strand direction or SV
length. The position in the variant ID is the reference point for an SV.

Use feature coordinates from phenotype metadata. Do not infer the midpoint from
a cis region, which can be clipped or asymmetric. Accept chromosome labels with
or without the `chr` prefix. Cross-chromosome variants and absent feature
coordinates receive numeric missing values. Old compact files also receive
numeric missing values during aggregation; preserve recorded distances.

Keep fitting, scientific values, physical row order, and full/LBF/CS schemas.
Keep checkpointed output defaults. Do not change WDL or submit Terra jobs.
Use local R checks and GitHub Actions; do not build a local Docker image.
