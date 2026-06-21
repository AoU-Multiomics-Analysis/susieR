#!/usr/bin/env bash
set -euo pipefail

for wdl in workflows/*.wdl; do
  miniwdl check "$wdl"
done

ruby -e 'require "yaml"; Dir[".github/workflows/*.yml", ".dockstore.yml"].each { |f| YAML.load_file(f) }; puts "YAML ok"'

ruby -e 'require "json"; Dir["examples/inputs/*.json"].each { |f| JSON.parse(File.read(f)) }; puts "Example JSON ok"'

Rscript -e 'for (f in c(list.files("R/scripts", pattern="[.]R$", full.names=TRUE), list.files("R/utils", pattern="[.]R$", full.names=TRUE))) parse(file=f); cat("R parse ok\n")'

git diff --check
