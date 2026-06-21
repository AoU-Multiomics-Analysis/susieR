# Docker Images

This repository builds two Docker images from the `containers/` directory.

| Dockerfile | Image | Contains |
|---|---|---|
| `containers/SusieR/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier` | Fine-mapping scripts from `R/scripts/` plus shared helpers from `R/utils/`. |
| `containers/PostAnalysis/Dockerfile` | `ghcr.io/aou-multiomics-analysis/susier/postanalysis` | Post-analysis scripts `AnnotateSusie.R` and `MergeSusie.R`. |

The fine-mapping image sets `SUSIER_FUNCTIONS_PATH=/opt/r/lib`, and the R scripts default to that location when sourcing shared helpers.

GitHub Actions rebuild the images when the relevant Dockerfile or copied R scripts change.
