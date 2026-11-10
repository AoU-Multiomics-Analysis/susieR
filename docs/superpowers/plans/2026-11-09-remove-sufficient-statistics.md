# Remove sufficient-statistics fitting

The user selected individual-data SuSiE for all fits. Remove the sufficient-statistics
fitter and shared XtX cache. Keep whole and rolling genotype-column caches,
ordinary fit optimizations, scientific settings, and phenotype residual SD export.

The triggering legacy-image failure was a test compatibility error: the new
phenotype SD regression required the pinned sufficient-statistics implementation
in an image that deliberately uses native individual-data fallback. It was not
evidence of an inaccurate sufficient-statistics fit. The feature removal follows
the user's stated preference.

1. Add a regression that eligible overlapping whole-cache windows invoke the
   ordinary fitter. Run it against the old code and confirm failure.
2. Remove the shared-matrix planner, storage, fit routing, counters, and obsolete
   feature-specific tests. Keep strict genotype/cache/PIP/mean/SD/CS checks.
3. Keep phenotype SD tests active for ordinary and legacy fits, output files,
   cached extraction, old fits, and empty tables.
4. Update active documentation and mark the prior implementation plans historical.
5. Run small local regressions, review the diff, push, and check both legacy and
   PGEN image workflows in GitHub Actions. Do not build Docker locally, run a full
   scientific dataset, or submit a Terra job.
