# Docker Images

This repository builds two Docker images from the `containers/` directory.

| Dockerfile | Image | Contains |
|---|---|---|
| `containers/SusieR/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier` | Fine-mapping scripts from `R/scripts/` plus shared helpers from `R/utils/`. |
| `containers/PostAnalysis/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/postanalysis` | Post-analysis scripts `AnnotateSusie.R`, `MergeSusie.R`, and `ComputeAncestrySkew.R`. |

The fine-mapping image sets `SUSIER_FUNCTIONS_PATH=/opt/r/lib`, and the R scripts default to that location when sourcing shared helpers.

## GitHub Actions

The Docker workflows build and push images to GitHub Container Registry (`ghcr.io`) on pushes and pull requests to `main`.

| Workflow | Dockerfile | Image | Rebuilds when |
|---|---|---|---|
| `.github/workflows/docker-image.yml` | `containers/SusieR/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier` | Fine-mapping image dependencies, `R/scripts/**`, or `R/utils/**` change. |
| `.github/workflows/PostAnalysisImage.yml` | `containers/PostAnalysis/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/postanalysis` | Post-analysis image dependencies, `R/scripts/AnnotateSusie.R`, `R/scripts/MergeSusie.R`, or `R/scripts/ComputeAncestrySkew.R` change. |

The post-analysis image copies the scripts called by the aggregate, annotate, and ancestry skew WDL tasks.
