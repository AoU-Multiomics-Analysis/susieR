# Combine supplied SuSiE tables directly

The user requested direct file merging and removed the requirement for
aggregation safeguards. This replaces the stricter approach in the preceding
repeated-label fix. Do not inspect the Terra submission or run new fits.

- Read the explicit local inputs and concatenate each table in supplied file
  order. Keep every row, including repeated paths and repeated result rows.
- Remove scientific consistency checks from compact and full aggregation:
  outcome ownership, status, CS membership/counts, variant labels, PIPs, scales,
  and metadata agreement do not decide whether rows can be combined.
- Keep supplied fields and values. Add missing legacy metadata without
  overwriting known values. Where index metadata is ambiguous, use missing
  values rather than selecting one conflicting value.
- Let table binding combine column sets and add missing cells. Keep full
  aggregation in bounded Arrow batches, with a temporary output until writes
  complete. File readability, localization, and write failures still apply.
- Keep typed WDL inputs, CLI names, output paths, native export behavior, and
  the predictor-index join unchanged. No new WDL or container recipe is needed.
- Replace old guard-error tests with small direct-merge fixtures. Verify exact
  row order, values, repeated paths, metadata, mixed/empty sources, and CLI I/O.
  Run the normal image checks in GitHub Actions, then merge after they pass.

## Verification

The compact/full merge cases fail under the previous implementation and pass
with direct concatenation. The compact, full, manifest, and CLI suites pass
locally. These include the normal CLI with gene and PGEN metadata, duplicate
failed-status index rows, and a readable placeholder expected-outcomes input.
Mixed missing/text scale columns, source ownership, empty native schemas,
stream-boundary order, and cleanup after a corrupt input are covered. R lint
and diff checks pass. Independent review found no unresolved material issues.

The gene manifest writer and the CLI's expected-outcome/PGEN agreement checks
were also simplified so they do not reintroduce the removed scientific checks.
Artifact URIs remain output-only text. Missing or conflicting locations remain
missing. No Terra submission, full scientific run, or local Docker build is used.
