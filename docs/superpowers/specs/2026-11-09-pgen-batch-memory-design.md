# Automatic PGEN batch memory

Use the prepared batch's sample count, physical variant records, and feature windows to request Terra task memory before a SuSiE worker starts. Keep all file inputs typed and preserve sequential feature execution.

The preparation task already opens each output triplet. Calculate the largest feature window and batch union there without decoding genotypes. Reserve at least 8 GB for fitting, increase that reserve with matrix size, add cache storage, allow for the existing dense numerical safeguard, and add fixed overhead and a percentage buffer. Record all terms and the final integer GB request. These are planning estimates, not measured peak RSS.

Use five dense matrix payloads plus retained model arrays as the fitting allowance. Reserve two `N` by `N` double matrices for the dense numerical safeguard. Use at least 2 GB for other data, with a count-based increase for retained batch tables. Add 25% capacity. If the implicit cache limit can be used, reserve 20% of the final allocation and solve the two percentage terms together. Additional genotypes require a positive fixed override because PVAR does not describe them.

Collect unused R objects after a gene has saved its full results and returned compact tables. Retain the whole-batch cache. For rolling mode, evict columns outside the next window and collect unused blocks before decoding incoming columns. Log stage RSS, peak RSS and available cgroup counters. A shell exit trap records failure status and available cgroup memory/OOM counters after a killed R process. This supplies evidence; it does not claim that every `Killed` message proves OOM.

`BatchMemoryGB = 0` selects automatic sizing. A positive value preserves a fixed request for all batches. The preparation task emits one ordered integer array as structured output data; the workflow uses the matching element for each batch's runtime request and `--job-memory-gb`. No additional PGEN localization is needed.

Whole and automatic cache modes reserve the batch union; rolling reserves the largest window; none reserves no genotype cache. An explicit cache budget must be covered by the job request. Cache mode selection continues to use its existing memory and reuse rules. Metadata-only sizing may overestimate memory when later sample or variant restrictions reduce the data.

Alternatives considered: a single fixed request does not adapt to batch size; a separate sizing task adds localization and scheduling. Preparation has the required metadata already, so calculation there is the smallest change.

Verification covers matrix-size scaling, cache modes and explicit budgets, invalid settings, physical duplicate variants, per-batch ordering, fixed overrides, task command localization, and Terra's workflow-scope file-writing restrictions. Add the estimator to installed-code GitHub Actions checks. Do not build Docker locally or submit Terra jobs.
