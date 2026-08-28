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

# Permitted Dockerfile RUN syntax is shell form or unescaped exec form with optional BuildKit RUN flags.
# Unknown or escaped exec syntax fails closed.
docker_run_invokes_install() {
  awk '
    function tokenize(command, tokens,    key, i, character, token, quote, escaped, count) {
      for (key in tokens) {
        delete tokens[key]
      }
      token = ""
      quote = ""
      escaped = 0
      count = 0
      for (i = 1; i <= length(command); i++) {
        character = substr(command, i, 1)
        if (escaped) {
          token = token character
          escaped = 0
          continue
        }
        if (quote != "") {
          if (quote == "\"" && character == "\\") {
            escaped = 1
          } else if (character == quote) {
            quote = ""
          } else {
            token = token character
          }
          continue
        }
        if (character == "\"" || character == "\047") {
          quote = character
        } else if (character == "\\") {
          escaped = 1
        } else if (character ~ /[[:space:]]/) {
          if (token != "") {
            tokens[++count] = token
            token = ""
          }
        } else if (character == "#" && token == "") {
          break
        } else {
          token = token character
        }
      }
      if (token != "") {
        tokens[++count] = token
      }
      return count
    }

    function tokens_invoke_install(tokens, count,    token_start, tool, token_index) {
      token_start = 1
      while (token_start <= count && (tokens[token_start] ~ /^--/ || tokens[token_start] ~ /^[[:alpha:]_][[:alnum:]_]*=/)) {
        token_start++
      }
      tool = tokens[token_start]
      if (tool == "apt" || tool == "apt-get" || tool == "conda" || tool == "mamba" || tool == "micromamba" || tool == "pip" || tool == "pip3") {
        for (token_index = token_start + 1; token_index <= count; token_index++) {
          if (tokens[token_index] == "install") {
            return 1
          }
        }
      }
      if ((tool == "python" || tool == "python3") && tokens[token_start + 1] == "-m" && tokens[token_start + 2] == "pip") {
        for (token_index = token_start + 3; token_index <= count; token_index++) {
          if (tokens[token_index] == "install") {
            return 1
          }
        }
      }
      return 0
    }

    function exec_form_installs(command, tokens,    key, i, character, token, in_string, count, closed, trailing) {
      for (key in tokens) {
        delete tokens[key]
      }
      sub(/^[[:space:]]*/, "", command)
      if (substr(command, 1, 1) != "[") {
        return 0
      }
      token = ""
      in_string = 0
      count = 0
      closed = 0
      for (i = 2; i <= length(command); i++) {
        character = substr(command, i, 1)
        if (in_string) {
          if (character == "\\") {
            return 1
          } else if (character == "\"") {
            tokens[++count] = token
            token = ""
            in_string = 0
          } else {
            token = token character
          }
          continue
        }
        if (character == "\"") {
          in_string = 1
        } else if (character == "]") {
          trailing = substr(command, i + 1)
          if (trailing !~ /^[[:space:]]*$/) {
            return 1
          }
          closed = 1
          break
        } else if (character != "," && character !~ /[[:space:]]/) {
          return 1
        }
      }
      if (in_string || !closed || count == 0) {
        return 1
      }
      return tokens_invoke_install(tokens, count)
    }

    function simple_command_installs(command, tokens,    count, trimmed) {
      trimmed = command
      sub(/^[[:space:]]*/, "", trimmed)
      while (trimmed ~ /^--/) {
        if (trimmed !~ /^--[[:alnum:]][[:alnum:]-]*=[^[:space:]]+[[:space:]]+/) {
          return 1
        }
        sub(/^--[[:alnum:]][[:alnum:]-]*=[^[:space:]]+[[:space:]]+/, "", trimmed)
      }
      if (trimmed ~ /^\[[[:space:]]*\"/) {
        return exec_form_installs(trimmed, tokens)
      }
      count = tokenize(trimmed, tokens)
      return tokens_invoke_install(tokens, count)
    }

    function run_body_installs(body,    i, character, quote, escaped, segment) {
      quote = ""
      escaped = 0
      segment = ""
      for (i = 1; i <= length(body); i++) {
        character = substr(body, i, 1)
        if (escaped) {
          segment = segment character
          escaped = 0
          continue
        }
        if (quote != "") {
          segment = segment character
          if (quote == "\"" && character == "\\") {
            escaped = 1
          } else if (character == quote) {
            quote = ""
          }
          continue
        }
        if (character == "\"" || character == "\047") {
          quote = character
          segment = segment character
        } else if (character == "\\") {
          escaped = 1
          segment = segment character
        } else if (character == ";" || character == "&" || character == "|") {
          if (simple_command_installs(segment)) {
            return 1
          }
          segment = ""
          if ((character == "&" || character == "|") && substr(body, i + 1, 1) == character) {
            i++
          }
        } else {
          segment = segment character
        }
      }
      return simple_command_installs(segment)
    }

    function inspect_run_line(line) {
      if (line !~ /^[[:space:]]*[Rr][Uu][Nn][[:space:]]+/) {
        return 0
      }
      sub(/^[[:space:]]*[Rr][Uu][Nn][[:space:]]+/, "", line)
      return run_body_installs(line)
    }

    {
      line = $0
      if (continued_line != "") {
        line = continued_line line
        continued_line = ""
      }
      if (line ~ /\\[[:space:]]*$/) {
        sub(/\\[[:space:]]*$/, " ", line)
        continued_line = line
        next
      }
      if (inspect_run_line(line)) {
        found = 1
        exit
      }
    }

    END {
      if (!found && continued_line != "" && inspect_run_line(continued_line)) {
        found = 1
      }
      exit(found ? 0 : 1)
    }
  '
}

require_no_prohibited_install() {
  local path="$1"
  if docker_run_invokes_install < "$path"; then
    printf 'Forbidden Docker RUN install command is present in %s.\n' "$path" >&2
    exit 1
  fi
}

require_prohibited_install_fixture() {
  local fixture="$1"
  if ! printf '%s\n' "$fixture" | docker_run_invokes_install; then
    printf 'Prohibited-install fixture was not rejected: %s\n' "$fixture" >&2
    exit 1
  fi
}

require_allowed_install_text_fixture() {
  local fixture="$1"
  if printf '%s\n' "$fixture" | docker_run_invokes_install; then
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
require_no_prohibited_install "$dockerfile"
require_absent 'install\.packages\(|remotes::install|git[[:space:]]+clone' "$dockerfile"

for prohibited_install_fixture in \
  'RUN apt-get -y install curl' \
  'RUN apt -o Acquire::Retries=3 install curl' \
  'RUN apt-get -o Acquire::Retries=3 install curl' \
  'RUN pip --no-cache-dir install requests' \
  'RUN pip --cache-dir /tmp install requests' \
  'RUN python -m pip --cache-dir /tmp install requests' \
  'RUN python3 -m pip --cache-dir /tmp install requests' \
  'RUN conda -y install r-base' \
  'RUN conda --rc-file /tmp/condarc install r-base' \
  'RUN mamba --yes install r-base' \
  'RUN mamba --rc-file /tmp/condarc install r-base' \
  'RUN micromamba -y install r-base' \
  'RUN micromamba --rc-file /tmp/condarc install r-base'; do
  require_prohibited_install_fixture "$prohibited_install_fixture"
done

for prohibited_exec_install_fixture in \
  'RUN ["apt", "-o", "Acquire::Retries=3", "install", "curl"]' \
  'RUN ["apt-get", "-y", "install", "curl"]' \
  'RUN ["pip", "--cache-dir", "/tmp", "install", "requests"]' \
  'RUN ["pip3", "--cache-dir", "/tmp", "install", "requests"]' \
  'RUN ["python", "-m", "pip", "--cache-dir", "/tmp", "install", "requests"]' \
  'RUN ["python3", "-m", "pip", "--cache-dir", "/tmp", "install", "requests"]' \
  'RUN ["conda", "--rc-file", "/tmp/condarc", "install", "r-base"]' \
  'RUN ["mamba", "--rc-file", "/tmp/condarc", "install", "r-base"]' \
  'RUN ["micromamba", "--rc-file", "/tmp/condarc", "install", "r-base"]'; do
  require_prohibited_install_fixture "$prohibited_exec_install_fixture"
done

for prohibited_buildkit_exec_install_fixture in \
  'RUN --mount=type=cache ["pip3", "install", "requests"]' \
  'RUN --network=none ["conda", "--rc-file", "/tmp/condarc", "install", "r-base"]' \
  'RUN ["pip3", install, "requests"]' \
  'RUN ["pip3", "\u0069nstall", "requests"]'; do
  require_prohibited_install_fixture "$prohibited_buildkit_exec_install_fixture"
done

for allowed_install_text_fixture in \
  '# RUN apt-get -o Acquire::Retries=3 install curl' \
  '# RUN ["apt-get", "-y", "install", "curl"]' \
  "ENV INSTALL_COMMAND='[\"apt-get\", \"-y\", \"install\", \"curl\"]'" \
  "RUN echo 'apt-get -y install curl'" \
  "RUN echo '[\"apt-get\", \"-y\", \"install\", \"curl\"]'" \
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

require_literal 'Run installed end-to-end smoke test' "$workflow"
require_literal 'Rscript /opt/r/scripts/run_checkpointed_window_susie.R' "$workflow"
require_literal 'diff /tmp/first/window_fit_index.tsv /tmp/resumed/window_fit_index.tsv' "$workflow"

require_absent 'GOOGLE_APPLICATION_CREDENTIALS|GCS_[A-Z_]*CREDENTIAL|gcloud[[:space:]]+auth|google-github-actions/auth|service_account|workload_identity|gsutil[[:space:]]+(cp|ls|rsync|rm|mb)|gs://' "$workflow"

printf 'Checkpointed window container contract passed.\n'
