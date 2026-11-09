# Adaptive prepared-genotype caching

The user approved automatic per-batch selection of whole-batch, rolling-window,
and no caching. Cache mean-imputed, centered, covariate-adjusted genotype columns.
Do not cache fitted models or share phenotype residuals between features.

## Policy

Use physical PGEN row indices from the PVAR and exact worklist windows in their
actual processing order. Let S be all requested variant counts, U the distinct
union count, and D the first-window count plus each subsequent window's new
variants relative to the previous window. Whole-batch saved fraction is 1-U/S;
rolling saved fraction is 1-D/S.

The automatic cache budget is 20% of the declared task memory. An explicit cache
memory limit overrides it. Use decimal GB consistently with the CLI inputs.
Estimate stored columns as 8*N*P bytes, plus four preparation blocks of at most
256 variants each. This reserves cache preparation buffers, not the full fit's
working memory. Peak-memory validation remains necessary.

Choose whole when its estimate fits and saved fraction is at least 0.5. Otherwise
choose rolling when its largest-window estimate fits and saved fraction is at
least 0.5. Otherwise select none. Explicit modes whole and rolling bypass only
the reuse threshold, and must fit their budget. Mode none preserves the current
per-feature path. Invalid modes and memory limits fail before fitting.

The first cache implementation supports the native-first PGEN-only path with a
fixed sample cohort and covariate design. Active ancestry MAF, a variant list,
or additional genotypes select none in automatic mode. Forced caching of these
inputs fails clearly. Preserve those inputs' existing statistical preparation.

## Storage and preparation

Read each newly required physical variant once with native mean imputation.
Remember all-missing and zero-alt exclusions. Store prepared values in bounded
blocks keyed by physical indices, preserving duplicate biological variant keys.
Whole mode prepares the union once. Rolling mode releases columns outside the
next window, then prepares only new columns. Assemble a fit matrix in original
PGEN order without changing outcome order or seeds. Release all blocks at batch
completion and on failure.

Covariate adjustment is columnwise with the same sample order and covariate
design. Compute the fast factorized projection for cache blocks and record the
existing near-covariate-explained diagnostics. If a selected genotype or phenotype
requires the existing dense numerical calculation, fit that feature through its
existing per-feature path. This preserves the established full-window guard;
it is not a new mean-imputation fallback. Other outcomes can keep using the cache.

## Integration and evidence

Expose GenotypeCacheMode and GenotypeCacheMemoryGB in WDL 1.0 and named CLI
arguments. WDL passes the declared task memory explicitly. Do not change File
localization or introduce workflow-scope file-writing. Return a per-batch JSON
decision/metrics report as a typed File output; this is structured output data,
not an argument wrapper. Record the selected mode in gene provenance and logs.

Compare none, whole, and rolling against real PGEN fixtures and complete SuSiE
fits, including dense and response-triggered guards. Compare scientific values
with tight tolerance for ordinary blockwise floating-point differences and
exact values for guarded recomputation. Verify original sample and variant order,
RNG state, output schemas, new/reused/evicted counts, memory decisions, and cleanup.
Run local R and WDL checks and wire smoke tests to GitHub Actions. Do not build
Docker locally or submit Terra jobs. Do not claim Terra validation or a speed
gain without its corresponding execution evidence.
