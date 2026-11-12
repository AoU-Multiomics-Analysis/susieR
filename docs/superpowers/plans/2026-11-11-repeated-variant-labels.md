# Preserve repeated variant labels during aggregation

The user reports `summary contains duplicate rows` and wants aggregation to
continue. The current guard treats a repeated outcome/variant/CS label as a
repeated row, even when values differ. The reader already supports distinct
physical predictors with equal labels. The exporter also joins CS membership
by label, which can expand rows through a many-to-many join.

## Decisions

- Preserve every summary, diagnostic, and full-output variant row and its input
  order. Report repeated labels with source, outcome, and count. Do not select,
  combine, or drop effect/PIP values based on labels or numerical equality.
- Keep ownership, schema, failed-status, credible-set identity, metadata, purity,
  and definite membership-error checks. Keep unique-label CS counts exact.
- For repeated-label CSs in old exports, reject a distinct-label count above
  cs_size or a total row count below cs_size. Otherwise retain rows and warn
  that exact physical membership cannot be verified from legacy labels.
- Keep streamed full aggregation bounded and atomic. Repeated outcome/variant
  labels across Arrow batches are warnings, not errors.
- For future exports, use predictor column index for the posterior/CS join,
  then remove the temporary index so the public output schema stays unchanged.
  This avoids row expansion without refitting or changing native CSs or PIPs.
- Actual failed batch files are not yet available. Terra access returned 401;
  test the supported repeated-label cases and clearly state this limit.

## Work

1. Add red compact/full fixtures for identical and differing repeated rows,
   physical-member counts, cross-batch repeats, and retained hard failures.
2. Add a red native-fit export case with repeated predictor labels. Verify one
   row per physical predictor and no extra CS membership rows.
3. Implement the aggregation and index-based export fixes independently.
4. Review, run small local checks, and verify the installed code in GitHub
   Actions. Existing batch/full outputs can be reaggregated without refitting.
   Do not launch a full scientific run or a Terra job.

## Local verification

- Both aggregation fixtures failed on the original duplicate guards before
  the fix. They pass after the fix, including identical and differing rows,
  stream-boundary repeats, and retained definite membership errors.
- The native-fit exporter fixture failed on the original label join. It passes
  with the index join, including equal labels inside and outside CSs. A
  unique-label baseline has identical public tables before and after the fix.
- The output regression runner, R lint, 13 container-routing checks, and
  `git diff --check` pass. The compact aggregation suite now runs in the default
  GitHub image checks.
- The broad local R runner cannot start because `eQTLUtils` is absent. No local
  dependency install or Docker build was used. Image checks run in GitHub.
