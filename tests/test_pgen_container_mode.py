"""Check smoke/full routing without running models or building an image."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "tests/test_pgen_container.sh"
FULL_SUITES = (
    "tests/run_pgen_gene_tests.R",
    "tests/run_pgen_batch_tests.R",
    "tests/run_tests.R",
    "tests/test_pgen_pipeline.R",
)


class ContainerModeTests(unittest.TestCase):
    @staticmethod
    def recording_commands(directory, names):
        log = directory / "calls.jsonl"
        stub = """#!/usr/bin/env python3
import json, os, sys
from pathlib import Path
with open(os.environ["PGEN_COMMAND_LOG"], "a") as output:
    output.write(json.dumps([Path(sys.argv[0]).name, *sys.argv[1:]]) + "\\n")
if os.environ.get("FAIL_SMOKE") == "true" and sys.argv[1:2] == ["tests/test_pgen_smoke.R"]:
    raise SystemExit(42)
if os.environ.get("FAIL_MAC") == "true" and sys.argv[1:2] == ["tests/test_pgen_mac.R"]:
    raise SystemExit(43)
if os.environ.get("FAIL_QC") == "true" and sys.argv[1:2] == ["tests/test_pgen_qc.R"]:
    raise SystemExit(44)
if os.environ.get("FAIL_PREPARATION") == "true" and sys.argv[1:2] == ["tests/test_pgen_preparation.R"]:
    raise SystemExit(45)
if os.environ.get("FAIL_IMPUTATION") == "true" and sys.argv[1:2] == ["tests/run_pgen_reader_tests.R"]:
    raise SystemExit(46)
if os.environ.get("FAIL_CACHE") == "true" and sys.argv[1:2] == ["tests/run_pgen_cache_tests.R"]:
    raise SystemExit(47)
"""
        for name in names:
            executable = directory / name
            executable.write_text(stub)
            executable.chmod(0o755)
        return log

    def run_runner(self, arguments=(), fail_smoke=False, fail_mac=False, fail_qc=False,
                   fail_preparation=False, fail_imputation=False, fail_cache=False):
        with tempfile.TemporaryDirectory(prefix="pgen-runner-") as directory:
            directory = Path(directory)
            log = self.recording_commands(directory, ("Rscript", "plink2", "bgzip", "tabix"))
            environment = os.environ.copy()
            environment.update(
                PATH=str(directory) + os.pathsep + environment["PATH"],
                PGEN_COMMAND_LOG=str(log),
                FAIL_SMOKE="true" if fail_smoke else "false",
                FAIL_MAC="true" if fail_mac else "false",
                FAIL_QC="true" if fail_qc else "false",
                FAIL_PREPARATION="true" if fail_preparation else "false",
                FAIL_IMPUTATION="true" if fail_imputation else "false",
                FAIL_CACHE="true" if fail_cache else "false",
                SUSIER_FUNCTIONS_PATH="/installed functions",
                SUSIER_SCRIPTS_PATH="/installed scripts",
            )
            result = subprocess.run(
                ["bash", str(RUNNER), *arguments], cwd=ROOT, env=environment,
                capture_output=True, text=True, check=False,
            )
            calls = [json.loads(line) for line in log.read_text().splitlines()] if log.exists() else []
            return result, calls

    def test_default_runs_mac_and_smoke_without_full_suites(self):
        result, calls = self.run_runner()
        self.assertEqual(result.returncode, 0, result.stderr)
        scripts = [call[1] for call in calls if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
        self.assertEqual(scripts, ["tests/test_pgen_runtime.R", "tests/test_pgen_mac.R", "tests/test_pgen_qc.R",
                                  "tests/test_pgen_preparation.R",
                                  "tests/run_pgen_reader_tests.R",
                                  "tests/run_pgen_cache_tests.R", "tests/test_pgen_smoke.R"])
        self.assertFalse(any(suite in scripts for suite in FULL_SUITES))
        self.assertTrue(any(call[:2] == ["plink2", "--version"] for call in calls))

    def test_full_flag_retains_all_regression_suites(self):
        result, calls = self.run_runner(["--full"])
        self.assertEqual(result.returncode, 0, result.stderr)
        scripts = [call[1] for call in calls if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
        self.assertEqual(scripts.count("tests/test_pgen_mac.R"), 1)
        self.assertEqual(scripts.count("tests/test_pgen_runtime.R"), 1)
        self.assertEqual(scripts.count("tests/test_pgen_qc.R"), 1)
        self.assertEqual(scripts.count("tests/test_pgen_preparation.R"), 1)
        self.assertEqual(scripts.count("tests/test_pgen_smoke.R"), 1)
        self.assertEqual(scripts.count("tests/run_pgen_reader_tests.R"), 1)
        self.assertEqual(scripts.count("tests/run_pgen_cache_tests.R"), 1)
        self.assertLess(scripts.index("tests/test_pgen_mac.R"), scripts.index("tests/test_pgen_smoke.R"))
        self.assertLess(scripts.index("tests/test_pgen_preparation.R"), scripts.index("tests/test_pgen_smoke.R"))
        for suite in FULL_SUITES:
            self.assertEqual(scripts.count(suite), 1)

    def test_smoke_failure_stops_full_checks(self):
        result, calls = self.run_runner(["--full"], fail_smoke=True)
        self.assertEqual(result.returncode, 42, result.stderr)
        scripts = [call[1] for call in calls if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
        self.assertFalse(any(suite in scripts for suite in FULL_SUITES))

    def test_mac_failure_stops_smoke_and_full_checks_in_both_modes(self):
        for arguments in ([], ["--full"]):
            with self.subTest(arguments=arguments):
                result, calls = self.run_runner(arguments, fail_mac=True)
                self.assertEqual(result.returncode, 43, result.stderr)
                scripts = [call[1] for call in calls
                           if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
                self.assertEqual(scripts, ["tests/test_pgen_runtime.R", "tests/test_pgen_mac.R"])

    def test_unknown_or_extra_arguments_fail_before_running_tools(self):
        for arguments in (["--ful"], ["--full", "extra"], [""]):
            with self.subTest(arguments=arguments):
                result, calls = self.run_runner(arguments)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(calls, [])
                self.assertIn("usage", result.stderr.lower())

    def test_shared_qc_failure_stops_model_smoke_in_both_modes(self):
        for arguments in ([], ["--full"]):
            with self.subTest(arguments=arguments):
                result, calls = self.run_runner(arguments, fail_qc=True)
                self.assertEqual(result.returncode, 44, result.stderr)
                scripts = [call[1] for call in calls
                           if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
                self.assertEqual(scripts, ["tests/test_pgen_runtime.R", "tests/test_pgen_mac.R", "tests/test_pgen_qc.R"])

    def test_workflow_passes_full_flag_only_when_selected(self):
        workflow = yaml.load(
            (ROOT / ".github/workflows/susier-v2-image.yml").read_text(),
            Loader=yaml.BaseLoader,
        )
        step = next(
            step for step in workflow["jobs"]["image"]["steps"]
            if "tests/test_pgen_container.sh" in step.get("run", "")
        )
        self.assertIn("inputs.full_checks", step["env"]["FULL_CHECKS"])
        for full_checks in ("false", "true"):
            with self.subTest(full_checks=full_checks), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                log = self.recording_commands(directory, ("docker",))
                environment = os.environ.copy()
                environment.update(
                    PATH=str(directory) + os.pathsep + environment["PATH"],
                    PGEN_COMMAND_LOG=str(log), FULL_CHECKS=full_checks,
                )
                result = subprocess.run(
                    ["bash", "-c", step["run"]], cwd=ROOT, env=environment,
                    capture_output=True, text=True, check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                calls = [json.loads(line) for line in log.read_text().splitlines()]
                self.assertEqual(len(calls), 1)
                arguments = calls[0]
                runner = arguments.index("tests/test_pgen_container.sh")
                self.assertEqual(arguments[runner + 1:], ["--full"] if full_checks == "true" else [])

    def test_preparation_failure_stops_smoke_and_full_suites(self):
        for arguments in ([], ["--full"]):
            with self.subTest(arguments=arguments):
                result, calls = self.run_runner(arguments, fail_preparation=True)
                self.assertEqual(result.returncode, 45, result.stderr)
                scripts = [call[1] for call in calls
                           if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
                self.assertEqual(scripts, ["tests/test_pgen_runtime.R", "tests/test_pgen_mac.R", "tests/test_pgen_qc.R",
                                          "tests/test_pgen_preparation.R"])

    def test_native_imputation_failure_stops_smoke_in_both_modes(self):
        for arguments in ([], ["--full"]):
            with self.subTest(arguments=arguments):
                result, calls = self.run_runner(arguments, fail_imputation=True)
                self.assertEqual(result.returncode, 46, result.stderr)
                scripts = [call[1] for call in calls
                           if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
                self.assertEqual(scripts, ["tests/test_pgen_runtime.R", "tests/test_pgen_mac.R", "tests/test_pgen_qc.R",
                                          "tests/test_pgen_preparation.R",
                                          "tests/run_pgen_reader_tests.R"])

    def test_cache_failure_stops_smoke_in_both_modes(self):
        for arguments in ([], ["--full"]):
            with self.subTest(arguments=arguments):
                result, calls = self.run_runner(arguments, fail_cache=True)
                self.assertEqual(result.returncode, 47, result.stderr)
                scripts = [call[1] for call in calls
                           if call[0] == "Rscript" and len(call) > 1 and call[1] != "-"]
                self.assertEqual(scripts, ["tests/test_pgen_runtime.R", "tests/test_pgen_mac.R", "tests/test_pgen_qc.R",
                                          "tests/test_pgen_preparation.R",
                                          "tests/run_pgen_reader_tests.R",
                                          "tests/run_pgen_cache_tests.R"])


if __name__ == "__main__":
    unittest.main()
