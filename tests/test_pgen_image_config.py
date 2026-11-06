"""Check image inputs and CI build/publication boundaries without Docker."""

import fnmatch
from pathlib import Path
import re
import shlex
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
IMAGE_DIR = ROOT / "containers/SusieR2"
WORKFLOW = ROOT / ".github/workflows/susier-v2-image.yml"
UTILS = {
    "ImportFunctions.R", "InitFunctions.R", "MemoryHelpers.R", "OptParser.R",
    "SusiePerformance.R", "SusieFunctions.R", "PgenFunctions.R",
    "SusieOutputFunctions.R", "PgenGeneFunctions.R", "PgenBatchFunctions.R",
}
SCRIPTS = {
    "plan_pgen_susie.R", "run_pgen_susie_batch.R", "aggregate_pgen_susie.R",
    "susie.R",
}
CONTAINED_SOURCE = {f"R/utils/{name}" for name in UTILS} | {
    f"R/scripts/{name}" for name in SCRIPTS
}


def load_yaml(path):
    # BaseLoader keeps GitHub's `on` key as a string, rather than YAML 1.1 bool.
    return yaml.load(path.read_text(), Loader=yaml.BaseLoader)


class ImageContractTests(unittest.TestCase):
    def test_package_recording_uses_supported_micromamba_json_output(self):
        """Catch the unsupported list --explicit flag that stopped image CI."""
        dockerfile = (IMAGE_DIR / "Dockerfile").read_text().replace("\\\n", " ")
        records = re.findall(r"micromamba list\s+([^>]+)>\s*([^\s&]+)", dockerfile)
        self.assertEqual(len(records), 1, "Record the installed package metadata once")
        arguments, output = records[0]
        self.assertIn("--json", shlex.split(arguments),
                      "micromamba 1.5.3 supports JSON package metadata")
        self.assertNotIn("--explicit", shlex.split(arguments),
                         "micromamba 1.5.3 list does not support --explicit")
        self.assertEqual(output, "/opt/conda/susier-v2-conda-packages.json")

    def test_config_job_installs_wdl_test_dependency_before_tests(self):
        """Catch a missing WDL Python module in the lightweight config job."""
        workflow = load_yaml(WORKFLOW)
        commands = "\n".join(step.get("run", "") for step in
                             workflow["jobs"]["config"]["steps"])
        install = commands.find("miniwdl==1.14.2")
        check = commands.find("python3 tests/test_pgen_wdl.py")
        self.assertGreaterEqual(install, 0, "The WDL test requires pinned miniwdl")
        self.assertGreater(check, install, "Install miniwdl before WDL tests")

    def test_environment_pins_runtime_and_source_build_dependencies(self):
        """Catch a missing native library or an unconstrained direct dependency."""
        path = IMAGE_DIR / "environment.yml"
        self.assertTrue(path.is_file(), "Pinned v2 environment is missing")
        env = load_yaml(path)
        self.assertEqual(env["channels"], ["conda-forge", "bioconda", "nodefaults"])
        dependencies = env["dependencies"]
        self.assertTrue(all(isinstance(item, str) and re.fullmatch(
            r"[A-Za-z0-9_.:-]+=[A-Za-z0-9_.-]+", item
        ) for item in dependencies), "Every direct Conda dependency needs an exact version")
        pins = dict(item.rsplit("=", 1) for item in dependencies)
        for name in {
            "r-base", "r-tidyverse", "r-data.table", "r-optparse", "r-arrow",
            "r-devtools", "r-rfast", "r-matrixstats", "r-digest", "r-jsonlite",
            "r-assertthat", "r-caret", "r-igraph", "r-lobstr", "r-rcpp",
            "r-cpp11", "r-cpp11armadillo", "r-mixsqp", "r-reshape", "r-crayon",
            "cxx-compiler", "make", "pkg-config", "plink2", "htslib",
            "zlib", "zstd", "libdeflate",
            "bioconductor-genomicranges", "bioconductor-summarizedexperiment",
            "bioconductor-rsamtools",
        }:
            self.assertIn(name, pins, f"Required image dependency is missing: {name}")

    def test_native_headers_are_checked_before_source_installation(self):
        """Catch missing compression development headers before R compilation."""
        dockerfile = (IMAGE_DIR / "Dockerfile").read_text()
        source_install = dockerfile.index("RUN Rscript -e")
        for header in ("zlib.h", "zstd.h", "libdeflate.h"):
            check = f"test -r /opt/conda/include/{header}"
            self.assertIn(check, dockerfile, f"Check required native header: {header}")
            self.assertLess(dockerfile.index(check), source_install)

    def test_image_copies_complete_explicit_runtime_source(self):
        """Catch missing CLI helpers and broad copies of unrelated source files."""
        path = IMAGE_DIR / "Dockerfile"
        self.assertTrue(path.is_file(), "V2 Dockerfile is missing")
        dockerfile = path.read_text()
        instructions = dockerfile.replace("\\\n", " ").splitlines()
        copies = [shlex.split(line)[1:] for line in instructions if line.startswith("COPY ")]
        copied = {source for fields in copies for source in fields[:-1] if source.startswith("R/")}
        self.assertEqual(copied, CONTAINED_SOURCE)
        self.assertIn("ENV SUSIER_FUNCTIONS_PATH=/opt/r/lib", dockerfile)
        self.assertIn("ENV SUSIER_SCRIPTS_PATH=/opt/r/scripts", dockerfile)
        self.assertRegex(dockerfile, r"(?m)^FROM mambaorg/micromamba:\d+\.\d+\.\d+(?:@sha256:[0-9a-f]{64})?$")
        self.assertIn("--strict-channel-priority", dockerfile)
        for source in ("susieR_0.16.6.tar.gz", "pgenlibr_0.7.0.tar.gz",
                       "26242562a4e244334fd9691d03bc1ef4d2d6c1d9"):
            self.assertIn(source, dockerfile)
        self.assertIn("dependencies = FALSE", dockerfile)
        self.assertIn('algo = "sha256"', dockerfile)

    def test_ci_rebuilds_only_for_contained_source(self):
        """Catch image builds caused by tests, WDL, docs or unrelated R files."""
        self.assertTrue(WORKFLOW.is_file(), "V2 image workflow is missing")
        workflow = load_yaml(WORKFLOW)
        events = workflow["on"]
        self.assertIn("workflow_dispatch", events)
        self.assertEqual(events["push"]["branches"], ["main"])
        self.assertEqual(events["pull_request"]["branches"], ["main", "codex/susier-performance"])
        for event in ("push", "pull_request"):
            filters = events[event]["paths"]
            self.assertEqual(set(filters), CONTAINED_SOURCE)
            for source in CONTAINED_SOURCE:
                self.assertTrue(any(fnmatch.fnmatchcase(source, pattern) for pattern in filters))
            for unrelated in ("tests/test_pgen_pipeline.R", "workflows/SusieR2.wdl",
                              "docs/SusieR2.md", "R/utils/CheckpointStore.R",
                              "containers/SusieR2/Dockerfile"):
                self.assertFalse(any(fnmatch.fnmatchcase(unrelated, pattern) for pattern in filters))

    def test_ci_tests_installed_code_before_main_only_publication(self):
        """Catch publication from a PR or publication before the installed smoke test."""
        self.assertTrue(WORKFLOW.is_file(), "V2 image workflow is missing")
        workflow = load_yaml(WORKFLOW)
        steps = workflow["jobs"]["image"]["steps"]
        smoke = next(i for i, step in enumerate(steps)
                     if "tests/test_pgen_container.sh" in step.get("run", ""))
        command = steps[smoke]["run"]
        self.assertIn('"$PWD/tests:/repo/tests:ro"', command)
        self.assertIn("--workdir /repo", command)
        self.assertNotIn('"$PWD:/repo', command, "Repository source must not shadow installed code")
        build = next(i for i, step in enumerate(steps)
                     if step.get("with", {}).get("load") == "true")
        self.assertLess(build, smoke)
        publications = [(i, step) for i, step in enumerate(steps)
                        if step.get("with", {}).get("push") == "true"]
        self.assertEqual(len(publications), 1)
        for i, step in publications:
            self.assertGreater(i, smoke)
            self.assertEqual(step["if"], "github.event_name == 'push' && github.ref == 'refs/heads/main'")
            self.assertIn("ghcr.io/aou-multiomics-analysis/susier/pgen", step["with"]["tags"])
            self.assertIn("ghcr.io/aou-multiomics-analysis/susier/pgen:main",
                          step["with"]["tags"].splitlines())
        for step in steps:
            if "docker/login-action@" in step.get("uses", ""):
                self.assertEqual(step["if"], "github.event_name == 'push' && github.ref == 'refs/heads/main'")


if __name__ == "__main__":
    unittest.main()
