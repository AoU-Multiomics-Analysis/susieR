#!/usr/bin/env bash
# Validate the complete descriptor with Cromwell's parser; do not run tasks.
set -euo pipefail
repo_dir="${1:-.}"
: "${WOMTOOL_JAR:?Set WOMTOOL_JAR to the pinned WOMtool 85 JAR}"
java_bin="${JAVA_BIN:-java}"
descriptor="$repo_dir/workflows/SusieR2.wdl"
echo 'CHECK PGEN WDL syntax and types with WOMtool 85'
"$java_bin" -jar "$WOMTOOL_JAR" validate "$descriptor"
echo 'CHECK PGEN input metadata extraction'
metadata_dir=$(mktemp -d)
trap 'rm -rf "$metadata_dir"' EXIT
"$java_bin" -jar "$WOMTOOL_JAR" inputs "$descriptor" > "$metadata_dir/inputs.json"
python3 - "$metadata_dir/inputs.json" <<'PY'
import json
import sys

with open(sys.argv[1]) as stream:
    inputs = json.load(stream)
required = {"ChromosomePgens", "PhenotypeBed", "QtlResults", "Covariates",
            "SampleList", "LegacyMissingCode"}
missing = {"SusieR2." + name for name in required} - inputs.keys()
assert not missing, "WOMtool did not extract required inputs: " + str(sorted(missing))
print("PGEN WDL validation and input metadata checks passed")
PY
