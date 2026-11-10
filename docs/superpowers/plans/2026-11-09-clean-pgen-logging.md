# Clean PGEN task logging

The user wants readable pipeline logs. Use concise batch setup and feature
progress in the main Terra log. Keep the full fit output in the existing gene
log. Use concise progress as the default.

## Design

- Keep SuSiE `verbose = TRUE`, model arguments, samples, RNG, and output values.
- For each fit, route stdout and ordinary messages into the existing gene log.
  Use one append-mode connection for output and conditions. Preserve warnings
  as original conditions, including `options(warn = 2)` behavior.
- Restore caller output and message handling on success or failure. Do not
  install a message sink. Close and flush the log before hashing artifacts.
- Main output shows batch/setup/resource/cache information, the active gene and
  fit, and one fit result with status, credible-set count, elapsed/stage times,
  and warning count. Keep detailed seeds, stages, variants, hints, and iterations
  in the gene log. Show nonconvergence and failure clearly.
- Remove duplicate completion messages and use consistent stage labels. Keep
  shared setup messages useful and suppress only package startup chatter.
- Scope is PGEN batch/gene logging; WDL task interfaces and outputs stay intact.

## Implementation

1. Add failing scoped-capture tests for stdout, messages, warnings, `warn = 2`,
   nested caller sinks, and errors. Add fit-progress/detail-log tests.
2. Add scoped capture and concise gene/fit messages. Keep timing and log hashes.
3. Clean batch/resource headings and remove duplicate completion output.
4. Run focused logging, gene/cache, output, and small CLI smoke checks. Review
   the combined diff, push, and verify GitHub image checks. No local Docker build,
   full scientific dataset, or Terra job is required.
