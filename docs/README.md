# Documentation

This directory holds focused notes for using and maintaining the susieR workflows.

| Document | Use |
|---|---|
| [`workflows.md`](workflows.md) | WDL entry points, workflow names, inputs, outputs, ancestry skew, and validation guidance. |
| [`scripts.md`](scripts.md) | R entrypoint scripts, expected inputs, and generated outputs. |
| [`docker.md`](docker.md) | Docker image layout, image names, and CI rebuild triggers. |
| [`inputs.md`](inputs.md) | Common WDL inputs and optional fine-mapping controls. |

Template WDL input JSONs live in [`../examples/inputs`](../examples/inputs). They use placeholder paths and values that should be replaced before submitting workflows.

For local repository checks, run:

```bash
tools/validate_repo.sh
Rscript tools/lint_r.R
```
