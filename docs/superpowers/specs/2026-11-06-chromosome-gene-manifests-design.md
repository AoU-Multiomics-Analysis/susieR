# Chromosome gene-output manifests

The user approved one TSV per chromosome, one row per gene, with all nine gene files. The workflow still accepts one or several chromosome PGEN bundles. Gene fitting and its outputs stay unchanged.

The selected approach extends the existing aggregation task. Its validated outcome index provides the original gene ID, gene key, chromosome, and batch. Completed gene-output locations are passed separately as string metadata, by artifact group. These paths are never opened or localized by the aggregate script. Files that the task reads remain typed WDL Files. This applies the existing instruction that cloud paths used only in output manifests are string metadata.

Alternatives were a separate manifest task, which adds a VM start, and a post-run resolver, which cannot return the TSVs from this workflow. Neither is needed. The implementation choice is to reuse aggregation.

Each row has chromosome, gene_id, gene_key, batch_id, and URI columns for the RDS, summary Parquet, LBF Parquet, full Parquet, covariate diagnostic TSV, variant-position TSV, outcome-index TSV, log, and provenance JSON. A gene with several outcomes still has one row. Zero-credible-set and nonconverged outcomes remain represented. Missing, duplicate, foreign, or incorrectly associated file metadata is an error.

The aggregate command creates newline-delimited metadata lists during command rendering. No workflow-scope file writers or guessed Cromwell paths are used. The workflow returns `Array[File] ChromosomeGeneManifests`, with files named `chr<N>.gene_outputs.tsv` for processed chromosomes. Existing inputs and scientific tables stay unchanged.

Tests must cover multiple chromosomes, repeated outcomes for a gene, shuffled URI arrays, original IDs, all nine artifact roles, invalid ownership, zero credible sets, nonconvergence, exact old table outputs, command-time File localization, unchanged URI metadata, and both WDL parsers. No local Docker build or Terra job is required. The complete workflow has not been tested on Terra.
