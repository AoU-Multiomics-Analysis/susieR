# Documentation

This directory holds focused notes for using and maintaining the susieR workflows.

| Document | Use |
|---|---|
| [`workflows.md`](workflows.md) | WDL entry points, workflow names, and validation guidance. |
| [`docker.md`](docker.md) | Docker image layout, image names, and build triggers. |

Template WDL input JSONs live in [`../examples/inputs`](../examples/inputs). They use placeholder paths and values that should be replaced before submitting workflows.

For a local repository health check, run:

```bash
tools/validate_repo.sh
```
