# Docker Images

This repository builds three Docker images from the `containers/` directory.

| Dockerfile | Image | Contains |
|---|---|---|
| `containers/SusieR/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier` | Fine-mapping and CV scripts `susie.R`, `ComputeR2Susie.R`, and `PrepSusieCVPCs.R` plus shared helpers from `R/utils/`. |
| `containers/PostAnalysis/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/postanalysis` | Post-analysis scripts `AnnotateSusie.R`, `MergeSusie.R`, and `ComputeAncestrySkew.R`. |
| `containers/CheckpointedWindowSusie/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/checkpointed-window` | Checkpointed-window runner and its checkpoint helpers. |

The fine-mapping image sets `SUSIER_FUNCTIONS_PATH=/opt/r/lib`, and the R scripts default to that location when sourcing shared helpers.

All images install their owned scripts under `/opt/r/scripts`. Compatibility symlinks are kept at `/` and `/tmp` for WDLs and existing Terra method configurations that reference those paths.

## GitHub Actions

The Docker workflows build smoke-test images on pushes and pull requests to `main`. They push images to GitHub Container Registry (`ghcr.io`) only when their workflow permits it.

| Workflow | Dockerfile | Image | Rebuilds when |
|---|---|---|---|
| `.github/workflows/docker-image.yml` | `containers/SusieR/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier` | Fine-mapping image dependencies, owned fine-mapping scripts, or `R/utils/**` change. |
| `.github/workflows/PostAnalysisImage.yml` | `containers/PostAnalysis/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/postanalysis` | Post-analysis image dependencies, `R/scripts/AnnotateSusie.R`, `R/scripts/MergeSusie.R`, or `R/scripts/ComputeAncestrySkew.R` change. |
| `.github/workflows/checkpointed-window-susie-image.yml` | `containers/CheckpointedWindowSusie/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/checkpointed-window` | Checkpointed runner, checkpoint helpers, tests, or image workflow changes. |

The post-analysis image copies the scripts called by the aggregate, annotate, and ancestry skew WDL tasks.

Each Docker workflow builds a local smoke-test image before pushing. The smoke tests verify that expected script paths exist and that the image can load the R packages needed by its scripts.

The checkpointed-window runner inherits the exact pinned scientific base digest
`sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764`.
GitHub Actions builds and smoke-tests this image. Do not build it locally for a
smoke test. The workflow uses no GCS credentials and writes no real GCS data.
It pushes the checkpointed image only for a push to `main`.
