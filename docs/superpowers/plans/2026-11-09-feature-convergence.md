# Mark credible sets from nonconverged features

The user wants to identify credible sets produced by a molecular-trait fit that
did not converge. Append logical `feature_converged` to the credible-set table
and the compact/full variant tables. Keep all scientific rows and values.

Formatting uses the native fit's convergence value for each molecular trait.
Cached extraction must not require the new field. Missing convergence in an old
fit remains unknown. Empty tables retain a logical column, while shared
checkpointed constructors keep their existing defaults.

Aggregation uses the validated outcome index as the authoritative source. This
lets old Parquets receive TRUE/FALSE without a refit. A supplied nonmissing flag
must be logical and must match the index. Normalize only this optional schema
addition; preserve other schema, ownership, streaming, and atomic-output checks.

1. Add failing mixed-feature export and aggregation tests.
2. Annotate exports and aggregation, including old/new/empty inputs and conflicting flags.
3. Update empty-schema checks, the small end-to-end smoke, and documentation.
4. Review, open a separate convergence PR, and verify both image workflows.
   Do not run a full scientific dataset or submit a Terra job.
