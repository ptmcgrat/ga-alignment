#!/usr/bin/env python3
"""Pairwise alignment of panel assemblies to the reference with minimap2 (spec, section 4).

For each assembly named, one minimap2 run:
  forward (default):  minimap2 -c -x <preset> --cs -t <threads> <reference> <assembly>
                      (--preset asm5 asm10 gives one run per preset and assembly)
                      the assembly is the query, the reference the target
  --reverse:          minimap2 -c -x <preset> --cs -t <threads> <assembly> <reference>
                      the reciprocity check: the reference is the query
Genomes are located through ga_common (never a hand-built path) and read as they are (BGZF).

Outputs, in <scratch>/pairwise/<version>/:
  <assembly>.<preset>.paf            forward runs
  reverse/<assembly>.<preset>.paf    reverse runs
  logs/<run>.json, logs/<run>.stderr  per run: command, minimap2 version, inputs with their recorded sha256,
                                     start and end (UTC), wall and CPU seconds, peak memory, exit status;
                                     minimap2's own messages
  manifest.tsv                       one row per finished run, appended
A PAF is written under <name>.partial and renamed only when minimap2 exits 0; a failed run's output is kept
as <name>.failed-<timestamp> for inspection.

A dry run unless --go: it then prints each command and where its output would go. It never overwrites: if a
PAF it would write, or its .partial, already exists, it stops before running anything and says what to
remove. --threads is required (minimap2's -t); --max-jobs runs that many assemblies at once (default 1);
--exclude leaves assemblies out, for example those already aligned.
Exit status: 0 if every run succeeded, 1 if any failed, 2 if it refused to start.
"""

import argparse
import concurrent.futures
import datetime
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-alignment"
SUPPLEMENTARY = {"YH_LG10Y": "h2tg000014l"}  # held beside the panel; aligned after the 24 (spec, section 6)
PRESETS = ("asm5", "asm10", "asm20")
MANIFEST_COLUMNS = ["version", "run", "direction", "preset", "query_id", "target_id", "paf", "paf_bytes",
                    "paf_records", "minimap2_version", "query_sha256_recorded", "target_sha256_recorded",
                    "start_utc", "end_utc", "wall_s", "cpu_s", "peak_rss_gb", "exit_status", "command"]


# ---------------------------------------------------------------- the coordinator's library

class Deps:
    """What the runner needs from ga_common, in a form tests can stand in for."""

    def __init__(self, locate, accessions, reference_id):
        self.locate, self.accessions, self.reference_id = locate, accessions, reference_id


def load_deps():
    from ga_common.paths import Paths
    from ga_common import panel
    paths = Paths.load()
    pan = panel.read_panel(str(paths.config_dir / "genome_panel_P1.tsv"))
    accessions = {r["id"]: r["accession_genbank"] for r in pan.rows}
    accessions.update(SUPPLEMENTARY)

    def locate(assembly_id):
        return Path(paths.local("genome_fasta", assembly_id=assembly_id,
                                accession_version=accessions[assembly_id]))
    return Deps(locate, accessions, pan.reference()["id"])


def add_throttle_arguments(parser):
    try:
        from ga_common import throttle
        throttle.add_arguments(parser)
    except ImportError:  # the laptop's tests run without ga_common
        parser.add_argument("--threads", type=int, default=None)
        parser.add_argument("--max-jobs", dest="max_jobs", type=int, default=None)


# ---------------------------------------------------------------- helpers

def utc_now():
    return datetime.datetime.now(datetime.timezone.utc)


def recorded_sha256(fasta):
    """The genome file's sha256 as provenance.json recorded it (checked against the file by check_inputs.py)."""
    try:
        provenance = json.loads((Path(fasta).parent / "provenance.json").read_text())
        return next(f["sha256"] for f in provenance.get("files", []) if f.get("role") == "genome")
    except (OSError, ValueError, StopIteration, KeyError):
        return ""


def count_lines(path):
    with open(path, "rb") as handle:
        return sum(block.count(b"\n") for block in iter(lambda: handle.read(1 << 20), b""))


def exit_code(status):
    if hasattr(os, "waitstatus_to_exitcode"):
        return os.waitstatus_to_exitcode(status)
    return os.WEXITSTATUS(status) if os.WIFEXITED(status) else -os.WTERMSIG(status)


def minimap2_version(binary):
    try:
        done = subprocess.run([binary, "--version"], capture_output=True, text=True, timeout=60)
        return done.stdout.strip() if done.returncode == 0 else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def plan_runs(args, deps):
    """[(run name, direction, query id, target id, query FASTA, target FASTA, output PAF)] for the request."""
    reference = deps.reference_id
    names = list(deps.accessions) if args.all else args.assemblies
    unknown = [a for a in names if a not in deps.accessions]
    if unknown:
        raise SystemExit("unknown assembly id(s): {} (see the panel)".format(", ".join(unknown)))
    if args.all:
        names = [a for a in names if a != reference and a not in SUPPLEMENTARY]
    unknown = [a for a in args.exclude if a not in deps.accessions]
    if unknown:
        raise SystemExit("unknown assembly id(s) in --exclude: {}".format(", ".join(unknown)))
    names = [a for a in names if a not in args.exclude]
    if reference in names:
        raise SystemExit("{} is the reference; it is never aligned to itself".format(reference))
    if not names:
        raise SystemExit("no assemblies left to align")
    base = Path(args.scratch).expanduser() / "pairwise" / args.version
    runs = []
    for preset in args.preset:
        for assembly in names:
            if args.reverse:
                query, target, folder, direction = reference, assembly, base / "reverse", "reverse"
            else:
                query, target, folder, direction = assembly, reference, base, "forward"
            name = "{}.{}".format(assembly, preset) + (".reverse" if args.reverse else "")
            runs.append({"run": name, "direction": direction, "preset": preset, "query_id": query,
                         "target_id": target, "query": deps.locate(query), "target": deps.locate(target),
                         "paf": folder / "{}.{}.paf".format(assembly, preset)})
    return base, runs


def command_for(binary, threads, run):
    return [binary, "-c", "-x", run["preset"], "--cs", "-t", str(threads), str(run["target"]), str(run["query"])]


# ---------------------------------------------------------------- one run

def run_one(binary, threads, run, base, version, mm2_version, manifest_lock):
    paf, partial = run["paf"], Path(str(run["paf"]) + ".partial")
    logs = base / "logs"
    command = command_for(binary, threads, run)
    record = {"version": version, "run": run["run"], "direction": run["direction"], "preset": run["preset"],
              "query_id": run["query_id"], "target_id": run["target_id"], "paf": str(paf.relative_to(base)),
              "minimap2_version": mm2_version, "command": " ".join(command),
              "query": str(run["query"]), "target": str(run["target"]),
              "query_sha256_recorded": recorded_sha256(run["query"]),
              "target_sha256_recorded": recorded_sha256(run["target"])}
    start, clock = utc_now(), time.monotonic()
    record["start_utc"] = start.strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(partial, "wb") as out, open(logs / (run["run"] + ".stderr"), "wb") as err:
        proc = subprocess.Popen(command, stdout=out, stderr=err)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = exit_code(status)
    record["end_utc"] = utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")
    record["wall_s"] = round(time.monotonic() - clock, 1)
    record["cpu_s"] = round(usage.ru_utime + usage.ru_stime, 1)
    record["peak_rss_gb"] = round(usage.ru_maxrss / 1024 / 1024, 2)   # Linux reports kilobytes
    record["exit_status"] = proc.returncode
    if proc.returncode == 0 and partial.stat().st_size > 0:
        partial.rename(paf)
        record["paf_bytes"] = paf.stat().st_size
        record["paf_records"] = count_lines(paf)
    else:
        failed = Path("{}.failed-{}".format(paf, start.strftime("%Y%m%d-%H%M%S")))
        partial.rename(failed)
        record.update({"paf": str(failed.relative_to(base)), "paf_bytes": failed.stat().st_size, "paf_records": ""})
    (logs / (run["run"] + ".json")).write_text(json.dumps(record, indent=2) + "\n")
    with manifest_lock:
        manifest = base / "manifest.tsv"
        new = not manifest.exists()
        with open(manifest, "a") as handle:
            if new:
                handle.write("\t".join(MANIFEST_COLUMNS) + "\n")
            handle.write("\t".join(str(record.get(c, "")) for c in MANIFEST_COLUMNS) + "\n")
    return record


# ---------------------------------------------------------------- main

def read_version():
    try:
        return (REPO / "VERSION").read_text().strip()
    except OSError:
        return "unknown version"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    which = parser.add_mutually_exclusive_group(required=True)
    which.add_argument("--assemblies", nargs="+", metavar="ID", help="assembly ids to align")
    which.add_argument("--all", action="store_true", help="every panel assembly except the reference")
    parser.add_argument("--exclude", nargs="+", default=[], metavar="ID",
                        help="leave these out (for example --all --exclude YH_YH1m once it is done)")
    parser.add_argument("--preset", nargs="+", choices=PRESETS, default=["asm5"],
                        help="minimap2 preset(s); several give one run per preset and assembly (default: asm5)")
    parser.add_argument("--reverse", action="store_true", help="the reference as query, the assembly as target")
    parser.add_argument("--version", default="PW_v1", help="the version folder (default: %(default)s)")
    parser.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH,
                        help="the project's scratch folder (default: %(default)s)")
    parser.add_argument("--minimap2", default="minimap2", help="the minimap2 binary (default: on PATH)")
    parser.add_argument("--go", action="store_true", help="run the alignments")
    add_throttle_arguments(parser)
    args = parser.parse_args(argv)
    if not getattr(args, "threads", None):
        parser.error("--threads is required: it becomes minimap2's -t")
    return args


def main(argv=None, deps=None):
    args = parse_args(argv)
    deps = deps or load_deps()
    base, runs = plan_runs(args, deps)
    jobs = max(1, min(getattr(args, "max_jobs", None) or 1, len(runs)))
    blocked = [p for r in runs for p in (r["paf"], Path(str(r["paf"]) + ".partial")) if p.exists()]
    mm2_version = minimap2_version(args.minimap2)
    print("align_pairwise.py ({}): {} runs, preset {}, {} threads each, {} at a time; minimap2 {}".format(
        read_version(), len(runs), " and ".join(args.preset), args.threads, jobs, mm2_version or "NOT FOUND"))
    for run in runs:
        print("  {:<26} {}\n  {:<26} -> {}".format(run["run"], " ".join(command_for(
            args.minimap2, args.threads, run)), "", run["paf"]))
    if blocked:
        print("STOP: these already exist, and this script never overwrites:", file=sys.stderr)
        for path in blocked:
            print("  rm {}".format(path), file=sys.stderr)
        return 2
    if not mm2_version:
        print("STOP: {} --version did not run; activate the ga-alignment environment".format(args.minimap2),
              file=sys.stderr)
        return 2
    if not args.go:
        print("Dry run: nothing has been run or written. Add --go to run these alignments.")
        return 0
    for folder in {r["paf"].parent for r in runs} | {base / "logs"}:
        folder.mkdir(parents=True, exist_ok=True)
    lock = threading.Lock()
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(run_one, args.minimap2, args.threads, run, base, args.version,
                               mm2_version, lock) for run in runs]
        for future in concurrent.futures.as_completed(futures):
            record = future.result()
            results.append(record)
            print("  {} {}: exit {}, {:.0f} s wall, {:.0f} s CPU, peak {} GB, {} records".format(
                "done  " if record["exit_status"] == 0 else "FAILED", record["run"], record["exit_status"],
                record["wall_s"], record["cpu_s"], record["peak_rss_gb"], record["paf_records"]), flush=True)
    failed = [r for r in results if r["exit_status"] != 0]
    print("Finished: {} succeeded, {} failed. Outputs and logs in {}".format(len(results) - len(failed), len(failed), base))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
