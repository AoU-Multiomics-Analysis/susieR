# WDL Workflows

The WDL descriptors live in `workflows/`. Each descriptor has a unique workflow name so Dockstore, Terra, and Cromwell output namespaces are easier to distinguish.

| Descriptor | Workflow | Use |
|---|---|---|
| `workflows/susieR.wdl` | `SusieRWorkflow` | Prepare phenotype-specific inputs and run fine-mapping. |
| `workflows/susieRonly.wdl` | `SusieROnlyWorkflow` | Run fine-mapping when inputs are already prepared. |
| `workflows/prepInputsSusieR.wdl` | `PrepSusieRWorkflow` | Prepare phenotype-specific input files only. |
| `workflows/ComputeR2Susie.wdl` | `ComputeR2SusieWorkflow` | Run cross-validation R2 evaluation. |
| `workflows/AggregateSusieTask.wdl` | `AggregateSusieTaskWorkflow` | Merge sharded Susie Parquet outputs. |
| `workflows/AnnotateSusie.wdl` | `AnnotateSusieWorkflow` | Annotate a merged Susie TSV. |
| `workflows/AggregateSusie.wdl` | `AggregateSusieWorkflow` | Run aggregate, annotate, and ancestry skew together. |

Template input JSONs live in `examples/inputs/`. They use placeholder paths and values; replace those with workspace-specific files before submitting.

Run local descriptor validation with:

```bash
tools/validate_repo.sh
```
