# Manifest-only per-gene SuSiE workflow outputs

The user requested that `SusieR2` copy and combine only the compact batch summary and credible-set Parquet files. Per-gene full Parquet files and SuSiE RDS fits must be available through manifests rather than copied into aggregation or a final workflow-output directory.

## Approach

Use the existing chromosome gene-output manifests. Keep the fitting task's per-gene File outputs so Cromwell uploads and retains their source objects. Remove the main workflow's full-table aggregation input, command argument, and combined full-table output. Remove the redundant dedicated full-Parquet task glob and the five direct workflow-level per-gene File arrays. Continue passing all per-gene output locations to aggregation as String metadata.

This preserves the original fitting artifacts and avoids another full-file transfer. Returning the large files as workflow Files would still select them for final-output copying. Removing their task File declarations would lose the source objects at VM shutdown. An explicit standalone full-only aggregation workflow remains available and is outside this change.

## Scope and constraints

- Use the existing pipeline worktree, which contains the current full aggregation implementation.
- Preserve WDL 1.0, command logs, typed compact File inputs, and task-command-only list writers.
- Preserve compact scientific outputs, fitting, cohort MAF, prepared genotypes, and existing diagnostic outputs.
- Return `ChromosomeGeneManifests` with all nine URI roles, including `full_parquet_uri` and `rds_uri`.
- Keep `delete_intermediate_output_files=false`; manifest strings do not prevent deletion of intermediate source files.
- Existing Terra output mappings for the six removed outputs must use `ChromosomeGeneManifests` instead.
- Do not build Docker locally or submit Terra jobs. Report local checks separately from Terra validation.

## Verification

Evaluate the real WDL workflow outputs with cloud-path fixtures and inspect the Files selected for final copying. Check that the worker still exports full files and RDS fits, while aggregation receives only compact data Files plus URI metadata. Run the existing manifest CLI test with nonexistent gene artifact URIs, the complete PGEN WDL tests, parser/type checks, and repository validation.
