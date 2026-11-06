#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

require_output() {
  local wdl="$1"
  local scope="$2"
  local declaration="$3"
  local outputs
  outputs="$(
    awk -v expected_scope="$scope" '
      /^[[:space:]]*(task|workflow)[[:space:]]/ {
        scope = $1 " " $2
        in_output = 0
      }
      scope == expected_scope && /^[[:space:]]*output[[:space:]]*\{/ {
        in_output = 1
        next
      }
      in_output && /^[[:space:]]*\}/ { in_output = 0 }
      in_output {
        sub(/^[[:space:]]+/, "")
        sub(/[[:space:]]+$/, "")
        print
      }
    ' "$wdl"
  )"
  if ! printf '%s\n' "$outputs" | awk -v expected="$declaration" '
    $0 == expected { found = 1 }
    END { exit !found }
  '; then
    printf '%s is missing the required %s output: %s\n' \
      "$wdl" "$scope" "$declaration" >&2
    exit 1
  fi
}

for wdl in workflows/susieR.wdl workflows/susieRonly.wdl; do
  if [ "$wdl" = workflows/susieR.wdl ]; then
    workflow_name="SusieRWorkflow"
  else
    workflow_name="SusieROnlyWorkflow"
  fi

  require_output "$wdl" 'task susieR' \
    'File CovariateExplainedVariants = "${OutputPrefix}.covariate_explained_variants.tsv"'
  require_output "$wdl" "workflow $workflow_name" \
    'File CovariateExplainedVariants = susieR.CovariateExplainedVariants'

  workflow_text="$(
    awk '
      /^[[:space:]]*workflow[[:space:]]/ { in_workflow = 1 }
      /^[[:space:]]*task[[:space:]]/ { in_workflow = 0 }
      in_workflow { sub(/#.*/, ""); print }
    ' "$wdl"
  )"
  if printf '%s\n' "$workflow_text" | awk '
    /write_(lines|tsv|map|json)[[:space:]]*[(]/ { print NR ":" $0; found = 1 }
    END { exit !found }
  '; then
    printf '%s contains a workflow-scope file-writing function.\n' "$wdl" >&2
    exit 1
  fi
done

printf 'Covariate-explained WDL output checks passed.\n'
