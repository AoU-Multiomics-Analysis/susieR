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
required = {"PgenManifest", "PhenotypeBed", "QtlResults", "Covariates", "LegacyMissingCode"}
missing = {"SusieR2." + name for name in required} - inputs.keys()
assert not missing, "WOMtool did not extract required inputs: " + str(sorted(missing))
actual_required = {name for name, description in inputs.items() if "optional" not in description}
assert actual_required == {"SusieR2." + name for name in required}, "Unexpected required workflow inputs"
assert inputs["SusieR2.SampleList"].startswith("File?"), "SampleList must be optional"
assert "SusieR2.ChromosomePgens" not in inputs, "The old array input must be replaced by the manifest"
print("PGEN WDL validation and input metadata checks passed")
PY
echo 'CHECK PGEN command localization with Cromwell 85'
"$java_bin" -Dscala.usejavacp=true -cp "$WOMTOOL_JAR" \
  scala.tools.nsc.MainGenericRunner "$repo_dir/tests/test_pgen_wom_render.scala" "$repo_dir"
