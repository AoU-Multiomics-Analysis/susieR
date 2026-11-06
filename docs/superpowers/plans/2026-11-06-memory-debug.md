# Optional memory diagnostics

User approved disabling detailed memory debugging by default. This is a small
logging/performance change to PR #54, not a statistical change.

- Use logical R option `susier.memory_debug`, default FALSE. Check it directly
  in the core fitting function so default use does not require debug helpers.
- Gate the two detailed memory-report blocks in finemapPhenotype, including
  object-size scans, allocation estimates and their diagnostic collections.
  Retain the three explicit collections after large-object removal.
- Expose opt-in PGEN CLI `--memory-debug`, legacy CLI `--memory_debug true`, and
  Terra WDL Boolean `MemoryDebug = false`. Preserve ordinary progress/error logs.
- Test default-off/opt-in logging and that expensive reporters are not called
  when off. Compare complete fitted objects and RNG state for normal and legacy
  fallback fixtures. Test installed CLI/WDL propagation and restore R options
  when a batch function returns or fails.
- Measure elapsed time and peak memory using the existing 10k x 20k fixture if
  available. Numerical fit identity is required; do not claim an unmeasured
  speedup. Run existing regression checks, then GitHub installed-image smoke.
- Keep unrelated user files, model settings, output formatting, and statistical
  preparation unchanged. No local Docker build or Terra job.
