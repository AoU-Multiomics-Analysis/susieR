# Batch timing and CPU use implementation plan

> **For agentic workers:** Use superpowers:executing-plans to implement each task and verify its results.

**Goal:** Report elapsed time for each gene and outcome, and use the CPU capacity available to each SuSiE batch.

**Architecture:** Keep gene fits sequential. Measure loading, preparation, and fitting through optional callbacks outside the scientific fit objects. Store outcome times in the existing outcome index and gene totals in the existing log and provenance. Detect CPU capacity with cgroup-aware tooling and apply a common bounded thread count to PLINK, BLAS, and supported data readers.

**Tech stack:** R, PLINK 2, WDL 1.0, micromamba, GitHub Actions.

## Constraints and decisions

- Keep the existing preemptible settings. The user withdrew the request to change them.
- Preserve scientific model arguments, RNG behavior, output ownership, and Terra File localization.
- Measure elapsed seconds. Include marginal z-score calculation in fine-mapping time.
- Report batch setup separately; do not charge shared input loading to each gene.
- Keep timing out of fit objects. Avoid a new gene artifact role.
- Use task CPU limits, not only the host CPU count. Keep one gene fit active at a time.
- Build and test images in GitHub Actions. Do not build Docker locally or submit Terra jobs.

## Task 1: Timing

- [ ] Add failing real-PGEN checks for stage times, failure paths, multiple outcomes, fit/RNG parity, and aggregation.
- [ ] Add optional stage callbacks in `SusieFunctions.R` and `PgenGeneFunctions.R`.
- [ ] Add timing columns to the outcome index, gene totals to logs/provenance, and separate batch setup logs.
- [ ] Preserve numeric timing columns during aggregation. Verify them in the installed image smoke test.

## Task 2: CPU use

- [ ] Add failing checks for bounded CPU selection, thread application, PLINK arguments, and provenance.
- [ ] Detect the task CPU allowance, configure BLAS and supported readers, and use the selected count for PLINK.
- [ ] Pin any added R packages and update image source boundaries only as needed.
- [ ] Test quota-aware detection in GitHub Actions without a local image build.

## Final checks

- [ ] Run the relevant R suites, pipeline tests, WDL checks, and image configuration checks.
- [ ] Obtain an independent code review and resolve material findings.
- [ ] Open a pull request, wait for checks, merge, and verify publication of the default image.
- [ ] State that the complete workflow has not been tested on Terra.
