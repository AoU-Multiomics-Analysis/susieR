# susieR WDL

[![WDL validation](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/wdl-validation.yml)
[![Repo validation](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/repo-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/repo-validation.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/r-lint.yml)
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/docker-image.yml)
[![Post-analysis Docker Image CI](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/PostAnalysisImage.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/susieR/actions/workflows/PostAnalysisImage.yml)

This repository provides WDL workflows and R scripts for running [susieR](https://github.com/stephenslab/susieR) fine-mapping on Terra. The workflows are designed for Terra data tables and support Bayesian fine-mapping of QTL signals identified by tensorQTL.

## Documentation

| Need | Start here |
|---|---|
| Workflow entry points and outputs | [`docs/workflows.md`](docs/workflows.md) |
| R script behavior and command-line arguments | [`docs/scripts.md`](docs/scripts.md) |
| Docker images and rebuild triggers | [`docs/docker.md`](docs/docker.md) |
| Shared WDL inputs | [`docs/inputs.md`](docs/inputs.md) |
| Documentation index | [`docs/README.md`](docs/README.md) |
| Template input JSONs | [`examples/inputs/`](examples/inputs/) |

## Repository Layout

```text
.
|-- R/
|   |-- scripts/        # R entrypoint scripts copied into Docker images
|   `-- utils/          # Shared R helper functions
|-- workflows/          # WDL workflow descriptors
|-- containers/         # Docker image definitions
|-- examples/inputs/    # Template WDL input JSONs
|-- docs/               # Focused user and maintainer documentation
|-- tools/              # Local validation and lint helpers
`-- .github/workflows/  # CI workflows
```

## Main Workflows

| Descriptor | Workflow | Use |
|---|---|---|
| [`workflows/susieR.wdl`](workflows/susieR.wdl) | `SusieRWorkflow` | Prepare phenotype-specific inputs and run fine-mapping. |
| [`workflows/susieRonly.wdl`](workflows/susieRonly.wdl) | `SusieROnlyWorkflow` | Run fine-mapping when inputs are already prepared. |
| [`workflows/ExtractMultiPhenotypeInputs.wdl`](workflows/ExtractMultiPhenotypeInputs.wdl) | `ExtractMultiPhenotypeInputsWorkflow` | Extract multi-phenotype BED and TensorQTL inputs without subsetting dosages. |
| [`workflows/prepInputsSusieR.wdl`](workflows/prepInputsSusieR.wdl) | `PrepSusieRWorkflow` | Prepare phenotype-specific input files only. |
| [`workflows/ComputeR2Susie.wdl`](workflows/ComputeR2Susie.wdl) | `ComputeR2SusieWorkflow` | Run cross-validation R2 evaluation. |
| [`workflows/AggregateSusieTask.wdl`](workflows/AggregateSusieTask.wdl) | `AggregateSusieTaskWorkflow` | Merge sharded Susie Parquet outputs. |
| [`workflows/AnnotateSusie.wdl`](workflows/AnnotateSusie.wdl) | `AnnotateSusieWorkflow` | Annotate a merged Susie TSV. |
| [`workflows/ComputeAncestrySkew.wdl`](workflows/ComputeAncestrySkew.wdl) | `ComputeAncestrySkew` | Compute ancestry skew for annotated Susie variants. |
| [`workflows/AggregateSusie.wdl`](workflows/AggregateSusie.wdl) | `AggregateSusieWorkflow` | Run aggregate, annotate, and ancestry skew together. |

## Post-Analysis

The aggregate, annotate, and ancestry skew steps are available as standalone workflows and as one combined `AggregateSusieWorkflow`. Ancestry skew now lives in this repository and is validated with the rest of the susieR WDLs.

## Local Checks

Run the repository validation and R lint checks before opening workflow or script changes:

```bash
tools/validate_repo.sh
Rscript tools/lint_r.R
```

The WDL validation checks all workflow descriptors in `workflows/`.
