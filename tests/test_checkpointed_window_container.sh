#!/usr/bin/env bash
set -euo pipefail

dockerfile="containers/CheckpointedWindowSusie/Dockerfile"
workflow=".github/workflows/checkpointed-window-susie-image.yml"

require_file() {
  local path="$1"
  if [ ! -f "$path" ]; then
    printf 'Required file is missing: %s\n' "$path" >&2
    exit 1
  fi
}

require_literal() {
  local literal="$1"
  local path="$2"
  if ! rg -Fq -- "$literal" "$path"; then
    printf 'Required literal is missing from %s: %s\n' "$path" "$literal" >&2
    exit 1
  fi
}

require_absent() {
  local pattern="$1"
  local path="$2"
  if rg -n -i -- "$pattern" "$path"; then
    printf 'Forbidden content is present in %s: %s\n' "$path" "$pattern" >&2
    exit 1
  fi
}

prohibited_install_pattern='(^[[:space:]]*RUN[[:space:]]+|[;&|][[:space:]]*)(apt-get|apt|conda|mamba|micromamba|pip|pip3)([[:space:]]+-{1,2}[^[:space:]]+)*[[:space:]]+install([[:space:]]|$)|(^[[:space:]]*RUN[[:space:]]+|[;&|][[:space:]]*)(python|python3)([[:space:]]+-{1,2}[^[:space:]]+)*[[:space:]]+-m[[:space:]]+pip([[:space:]]+-{1,2}[^[:space:]]+)*[[:space:]]+install([[:space:]]|$)|install\.packages\(|remotes::install|git[[:space:]]+clone'

require_prohibited_install_fixture() {
  local fixture="$1"
  if ! printf '%s\n' "$fixture" | rg -iq -- "$prohibited_install_pattern"; then
    printf 'Prohibited-install fixture was not rejected: %s\n' "$fixture" >&2
    exit 1
  fi
}

require_allowed_install_text_fixture() {
  local fixture="$1"
  if printf '%s\n' "$fixture" | rg -iq -- "$prohibited_install_pattern"; then
    printf 'Allowed text fixture was rejected as an install command: %s\n' "$fixture" >&2
    exit 1
  fi
}

step_block() {
  local step_name="$1"
  awk -v step_name="$step_name" '
    $0 == "      - name: " step_name { in_step = 1; next }
    in_step && /^      - name: / { exit }
    in_step { print }
  ' "$workflow"
}

require_publish_condition() {
  local step_name="$1"
  local block
  block="$(step_block "$step_name")"
  if [ -z "$block" ]; then
    printf 'Required workflow step is missing: %s\n' "$step_name" >&2
    exit 1
  fi
  if ! printf '%s\n' "$block" | rg -Fqx "        if: github.event_name == 'push' && github.ref == 'refs/heads/main'"; then
    printf 'Publish-sensitive workflow step lacks the push-to-main condition: %s\n' "$step_name" >&2
    exit 1
  fi
}

require_unconditional_step() {
  local step_name="$1"
  local block
  block="$(step_block "$step_name")"
  if [ -z "$block" ]; then
    printf 'Required workflow step is missing: %s\n' "$step_name" >&2
    exit 1
  fi
  if printf '%s\n' "$block" | rg -q '^        if:'; then
    printf 'Workflow step must run for pull requests, manual runs, and main pushes: %s\n' "$step_name" >&2
    exit 1
  fi
}

require_file "$dockerfile"
require_file "$workflow"

rg -Fq 'FROM ghcr.io/aou-multiomics-analysis/susier@sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764' containers/CheckpointedWindowSusie/Dockerfile
require_literal 'ENV CHECKPOINTED_SUSIE_BASE_IMAGE_DIGEST=sha256:07f9ddcb00391cceb6d5432144e38b16358b7a6ca7766ae3bc1b8b4aa3bac764' "$dockerfile"
require_literal 'ENV SUSIER_FUNCTIONS_PATH=/opt/r/lib' "$dockerfile"
require_literal 'COPY R/utils/CheckpointedWindowSusieFunctions.R R/utils/CheckpointStore.R /opt/r/lib/' "$dockerfile"
require_literal 'COPY R/scripts/run_checkpointed_window_susie.R /opt/r/scripts/' "$dockerfile"
require_literal 'RUN ln -sf /opt/r/scripts/run_checkpointed_window_susie.R /tmp/run_checkpointed_window_susie.R' "$dockerfile"
require_literal 'USER root' "$dockerfile"
require_literal 'USER $MAMBA_USER' "$dockerfile"
require_literal 'CMD ["bash"]' "$dockerfile"
require_absent "$prohibited_install_pattern" "$dockerfile"

for prohibited_install_fixture in \
  'RUN apt-get -y install curl' \
  'RUN pip --no-cache-dir install requests' \
  'RUN conda -y install r-base' \
  'RUN mamba --yes install r-base' \
  'RUN micromamba -y install r-base' \
  'RUN python3 -m pip --no-cache-dir install requests'; do
  require_prohibited_install_fixture "$prohibited_install_fixture"
done

for allowed_install_text_fixture in \
  "RUN echo 'apt-get -y install curl'" \
  "RUN printf '%s\\n' 'pip --no-cache-dir install requests'"; do
  require_allowed_install_text_fixture "$allowed_install_text_fixture"
done

require_literal '.github/workflows/checkpointed-window-susie-image.yml' "$workflow"
require_literal 'containers/CheckpointedWindowSusie/Dockerfile' "$workflow"
require_literal 'R/scripts/run_checkpointed_window_susie.R' "$workflow"
require_literal 'R/utils/CheckpointedWindowSusieFunctions.R' "$workflow"
require_literal 'R/utils/CheckpointStore.R' "$workflow"
require_literal 'tests/test_checkpointed_window_manifest.R' "$workflow"
require_literal 'tests/test_checkpoint_store.R' "$workflow"
require_literal 'tests/test_checkpointed_window_model.R' "$workflow"
require_literal 'tests/test_checkpointed_window_controller.R' "$workflow"
require_literal 'packages: write' "$workflow"
require_literal 'images: ghcr.io/${{ github.repository }}/checkpointed-window' "$workflow"
require_literal 'uses: actions/checkout@v7' "$workflow"
require_literal 'uses: docker/build-push-action@v7' "$workflow"
require_literal 'file: containers/CheckpointedWindowSusie/Dockerfile' "$workflow"
require_literal 'load: true' "$workflow"
require_literal 'tags: susier-checkpointed-window:smoke-test' "$workflow"
require_literal "if: github.event_name == 'push' && github.ref == 'refs/heads/main'" "$workflow"
require_literal 'push: true' "$workflow"
require_literal 'set -euo pipefail' "$workflow"
require_literal 'uses: docker/login-action@v4' "$workflow"
require_literal 'password: ${{ secrets.GITHUB_TOKEN }}' "$workflow"

login_line="$(rg -n -F '      - name: Log in to ghcr.io with docker' "$workflow" | cut -d: -f1)"
smoke_build_line="$(rg -n -F '      - name: Build docker image for smoke tests' "$workflow" | cut -d: -f1)"
if [ -z "$login_line" ] || [ -z "$smoke_build_line" ] || [ "$login_line" -ge "$smoke_build_line" ]; then
  printf 'GHCR login must occur before the smoke-image build.\n' >&2
  exit 1
fi
require_unconditional_step 'Log in to ghcr.io with docker'
require_publish_condition 'Extract Docker metadata'
require_publish_condition 'Build and push docker image'

for smoke_command in \
  'docker run --rm susier-checkpointed-window:smoke-test test -f /opt/r/scripts/run_checkpointed_window_susie.R' \
  'docker run --rm susier-checkpointed-window:smoke-test test -f /opt/r/lib/CheckpointedWindowSusieFunctions.R' \
  'docker run --rm susier-checkpointed-window:smoke-test test -f /opt/r/lib/CheckpointStore.R' \
  'docker run --rm susier-checkpointed-window:smoke-test Rscript /opt/r/scripts/run_checkpointed_window_susie.R --help' \
  "docker run --rm susier-checkpointed-window:smoke-test Rscript -e 'library(tidyverse); library(data.table); library(arrow); library(jsonlite); library(digest); library(susieR); stopifnot(nzchar(Sys.which(\"gsutil\")))'" \
  'docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpointed_window_manifest.R' \
  'docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpoint_store.R' \
  'docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpointed_window_model.R' \
  'docker run --rm -v "$PWD:/work" -w /work susier-checkpointed-window:smoke-test Rscript tests/test_checkpointed_window_controller.R'; do
  require_literal "$smoke_command" "$workflow"
done

require_absent 'GOOGLE_APPLICATION_CREDENTIALS|GCS_[A-Z_]*CREDENTIAL|gcloud[[:space:]]+auth|google-github-actions/auth|service_account|workload_identity|gsutil[[:space:]]+(cp|ls|rsync|rm|mb)|gs://' "$workflow"

printf 'Checkpointed window container contract passed.\n'
