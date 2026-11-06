#!/usr/bin/env python3
"""Run one synthetic benchmark worker and record process resource use."""

import argparse
import csv
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time


def positive_integer(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("The value must be a positive integer.")
    return number


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", required=True, type=Path)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=("baseline", "candidate", "factored", "factored_fast", "fit_baseline", "fit_quiet", "fit_fast_z", "production", "tabix_read", "pgen_read", "pgen_integer", "tabix_pipeline", "pgen_pipeline"))
    parser.add_argument("--library", required=True, type=Path)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=positive_integer)
    parser.add_argument("--variants", required=True, type=positive_integer)
    parser.add_argument("--blas", choices=("reference", "openblas", "accelerate"), default="reference")
    parser.add_argument("--threads", type=positive_integer, default=1)
    parser.add_argument("--save-fit-inputs", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--profile-r", type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["BENCH_LIBRARY"] = str(args.library.resolve())
    environment["R_MAX_VSIZE"] = "16G"
    if args.save_fit_inputs:
        environment["BENCH_FIT_CACHE"] = str(args.save_fit_inputs.resolve())
    else:
        environment.pop("BENCH_FIT_CACHE", None)
    environment["BENCH_PREPARE_ONLY"] = "true" if args.prepare_only else "false"
    if args.profile_r:
        environment["SUSIER_R_PROFILE"] = str(args.profile_r.resolve())
    else:
        environment.pop("SUSIER_R_PROFILE", None)
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "ARROW_NUM_THREADS",
        "RCPP_PARALLEL_NUM_THREADS",
        "VROOM_THREADS",
    ):
        environment[name] = str(args.threads)

    # Freeze the R source before launch. R can read a --file input incrementally.
    worker_snapshot = output / "worker_snapshot.R"
    shutil.copyfile(args.worker.resolve(), worker_snapshot)
    r_arguments = [
        str(args.inputs.resolve()),
        str(output),
        args.mode,
        str(args.repo.resolve()),
        str(args.samples),
        str(args.variants),
    ]
    if args.blas == "reference":
        command = ["Rscript", str(worker_snapshot), *r_arguments]
    else:
        if sys.platform != "darwin":
            raise RuntimeError("The optional BLAS overrides require macOS.")
        r_home = subprocess.check_output(["R", "RHOME"], text=True).strip()
        target = (
            Path("/opt/homebrew/opt/openblas/lib/libopenblas.dylib")
            if args.blas == "openblas"
            else Path(r_home) / "lib" / "libRblas.vecLib.dylib"
        )
        if not target.is_file():
            raise FileNotFoundError(target)
        shadow = output / "blas"
        shadow.mkdir(exist_ok=True)
        link = shadow / "libRblas.dylib"
        if link.is_symlink():
            link.unlink()
        link.symlink_to(target.resolve())
        environment["R_HOME"] = r_home
        environment["DYLD_LIBRARY_PATH"] = str(shadow)
        command = [str(Path(r_home) / "bin" / "exec" / "R"),
                   "--no-echo", "--no-restore", f"--file={worker_snapshot}",
                   "--args", *r_arguments]
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    process = None
    exit_code = 127
    with (output / "run.log").open("w", encoding="utf-8") as log:
        try:
            process = subprocess.Popen(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=environment,
                start_new_session=True,
            )
            print(f"R worker PID: {process.pid}", flush=True)
            try:
                while process.poll() is None:
                    time.sleep(0.25)
                exit_code = process.returncode
            except KeyboardInterrupt:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                exit_code = 130
        except OSError as error:
            log.write(f"Cannot start the R worker: {error}\n")
            print(f"Cannot start the R worker: {error}", file=sys.stderr)

    elapsed = time.monotonic() - started
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    # Darwin reports bytes. Linux and other supported Unix systems report KiB.
    peak_rss_bytes = int(usage_after.ru_maxrss)
    if sys.platform != "darwin":
        peak_rss_bytes *= 1024
    record = {
        "mode": args.mode,
        "blas_requested": args.blas,
        "threads_requested": args.threads,
        "samples": args.samples,
        "variants": args.variants,
        "pid": process.pid if process is not None else "",
        "exit_code": exit_code,
        "wall_seconds": f"{elapsed:.6f}",
        "user_cpu_seconds": f"{usage_after.ru_utime - usage_before.ru_utime:.6f}",
        "system_cpu_seconds": f"{usage_after.ru_stime - usage_before.ru_stime:.6f}",
        "peak_rss_bytes": peak_rss_bytes,
        "peak_rss_gib": f"{peak_rss_bytes / 1024 ** 3:.6f}",
        "platform": sys.platform,
    }
    with (output / "resource_usage.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=record.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(record)
    if exit_code < 0:
        return 128 + abs(exit_code)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
