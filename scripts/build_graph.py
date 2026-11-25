#!/usr/bin/env python3
"""Build a Minigraph-Cactus graph with the settings of spec section 4 (step 4, agreed 2026-11-24).

Two modes:
  --chromosome LG20   a chromosome run: for each assembly, the sequences our pairwise alignments assign to that
                      reference chromosome (assignment.tsv of the newest PW_v1 summary) are written to a small FASTA,
                      and the reference contributes that chromosome alone
  (no --chromosome)   a whole-genome run: the genome files themselves, with the reference's 22 linkage groups given
                      by --refContigs and its other sequences lumped by --otherContig chrOther
Assemblies: --set subset (the reference and the seven scored assemblies; the default), --set all (the reference and
the 24 others), or --assemblies ID ... (the reference is always added).

The command, run under Apptainer from the Cactus 4.0.0 image:
  taskset -c <cpus> apptainer exec --cleanenv --env TMPDIR=<work> --bind ... <image> cactus-pangenome <jobstore> <seqfile>
      --outDir <out> --outName <name> --reference MZ_NMK2024a --refContigs ... [--otherContig chrOther]
      --gfa full clip (or full alone with --set all) --viz --validate --workDir <work>
      --maxCores N --mgCores N --consCores N --indexCores N --mapCores min(6, N) --maxMemory <memory> [--restart]
--viz adds odgi's 1-D image of each chromosome's full graph (agreed 2026-11-25); it changes no graph.
taskset confines the whole run to the CPUs given (--cpus, default the first --threads CPUs); the per-step options
size each job; --maxMemory is the ceiling --consRetainPages auto sizes against. Run it under nice -n 10 in screen.

Each run has its own folder, <scratch>/graph/<version>/<name>/ (name: the chromosome, or the set): seqfile.txt,
fasta/ (chromosome runs), jobstore/ (Toil's; keep it until the run has succeeded), work/, out/, run.log and
run.json (command, inputs with checksums, the assignment table used, start and end, wall and CPU time, the
largest process's peak memory, exit status). A dry run unless --go. It never overwrites a run folder; after an
interruption, the same command with --restart resumes from the job store.
"""

import argparse
import datetime
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pafkit  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-alignment"
DEFAULT_SIF = Path.home() / "containers" / "cactus_v4.0.0.sif"
SUBSET = ["MZ_NMetZ", "MZ_MZ4f", "AC_SAMEA123396441", "YH_YH1m", "YH_YH2f", "RC_006pf", "DT_SAMEA123396444"]
OTHER_CONTIG = "chrOther"


# ---------------------------------------------------------------- inputs

def reference_linkage_groups(provenance):
    """{label: name} of the reference's linkage groups (assembled molecules typed Linkage Group)."""
    return {e["chr_name"]: n for n, e in provenance.items()
            if e.get("role") == "assembled-molecule" and e.get("molecule_type") == "Linkage Group"}


def newest_assignment(scratch, pairwise_version):
    folders = sorted((Path(scratch) / "pairwise" / pairwise_version / "summary").glob("*/assignment.tsv"))
    return folders[-1] if folders else None


def assigned_sequences(assignment_rows, assembly, preset, label):
    run = "{}.{}".format(assembly, preset)
    return [r["query"] for r in assignment_rows if r["run"] == run and r["chromosome_label"] == label]


def extract(fasta_gz, names, out_path):
    """Write the named sequences (headers cut to their first word) to a plain FASTA; return the names found."""
    wanted, found, keep = set(names), [], False
    opener = gzip.open if str(fasta_gz).endswith(".gz") else open
    with opener(fasta_gz, "rb") as handle, open(out_path, "wb") as out:
        for line in handle:
            if line.startswith(b">"):
                name = line[1:].split(None, 1)[0].decode()
                keep = name in wanted
                if keep:
                    found.append(name)
                    out.write(b">" + name.encode() + b"\n")
            elif keep:
                out.write(line)
    return found


def sha256_file(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def cpu_count(cpus):
    total = 0
    for part in cpus.split(","):
        a, _, b = part.partition("-")
        total += (int(b) - int(a) + 1) if b else 1
    return total


# ---------------------------------------------------------------- the command

def cactus_command(args, folder, name, ref_contigs, other_contig, gfa, binds):
    threads = cpu_count(args.cpus)
    work = folder / "work"
    cmd = [args.taskset, "-c", args.cpus, args.apptainer, "exec", "--cleanenv", "--env", "TMPDIR={}".format(work)]
    for b in sorted(binds):
        cmd += ["--bind", "{0}:{0}".format(b)]
    cmd += [str(args.sif), "cactus-pangenome", str(folder / "jobstore"), str(folder / "seqfile.txt"),
            "--outDir", str(folder / "out"), "--outName", name, "--reference", args.reference_id,
            "--refContigs"] + list(ref_contigs)
    if other_contig:
        cmd += ["--otherContig", other_contig]
    cmd += ["--gfa"] + list(gfa) + ["--viz", "--validate", "--workDir", str(work),
                                    "--maxCores", str(threads), "--mgCores", str(threads), "--consCores", str(threads),
                                    "--indexCores", str(threads), "--mapCores", str(min(6, threads)),
                                    "--maxMemory", args.memory]
    if args.restart:
        cmd.append("--restart")
    return cmd


def plan(args, deps):
    """The run's name, folder, samples [(id, FASTA path)], reference contigs, other contig and graph types."""
    reference = deps.reference_id
    if args.assemblies:
        assemblies = list(args.assemblies)
    elif args.set == "all":
        assemblies = deps.panel_ids()
    else:
        assemblies = list(SUBSET)
    unknown = [a for a in assemblies if a not in deps.accessions or a == reference]
    if unknown:
        raise SystemExit("not panel assemblies (or the reference): {}".format(", ".join(unknown)))
    name = args.chromosome or args.name or ("all" if args.set == "all" and not args.assemblies else "subset")
    folder = (Path(args.scratch).expanduser() / "graph" / args.version / name).resolve()
    gfa = ["full"] if args.set == "all" and not args.assemblies else ["full", "clip"]
    ref_fasta = Path(deps.locate(reference))
    groups = reference_linkage_groups(pafkit.provenance_sequences(ref_fasta))
    samples = [(reference, ref_fasta)] + [(a, Path(deps.locate(a))) for a in assemblies]
    if args.chromosome:
        if args.chromosome not in groups:
            raise SystemExit("{} is not one of the reference's linkage groups".format(args.chromosome))
        return name, folder, samples, [groups[args.chromosome]], None, gfa
    order = sorted(groups, key=lambda label: int(pafkit.label_number(label)))
    return name, folder, samples, [groups[label] for label in order], OTHER_CONTIG, gfa


# ---------------------------------------------------------------- main

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    which = parser.add_mutually_exclusive_group()
    which.add_argument("--set", choices=("subset", "all"), default="subset")
    which.add_argument("--assemblies", nargs="+", metavar="ID")
    parser.add_argument("--chromosome", help="a chromosome run for this reference linkage group, e.g. LG20")
    parser.add_argument("--name", help="the run's folder name (default: the chromosome, subset or all)")
    parser.add_argument("--version", default="GR_v1", help="the graph version folder (default: %(default)s)")
    parser.add_argument("--pairwise-version", default="PW_v1")
    parser.add_argument("--assignment", type=Path, help="assignment.tsv to choose sequences (default: newest summary)")
    parser.add_argument("--preset", default="asm5", help="the pairwise runs whose assignment is used")
    parser.add_argument("--threads", type=int, default=16, help="CPUs for the run, when --cpus is not given")
    parser.add_argument("--cpus", help="the CPU list for taskset, e.g. 0-15 (default: the first --threads CPUs)")
    parser.add_argument("--memory", default="200G", help="--maxMemory (default: %(default)s)")
    parser.add_argument("--sif", type=Path, default=DEFAULT_SIF)
    parser.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH)
    parser.add_argument("--taskset", default="taskset")
    parser.add_argument("--apptainer", default="apptainer")
    parser.add_argument("--restart", action="store_true", help="resume an interrupted run from its job store")
    parser.add_argument("--go", action="store_true", help="run it")
    args = parser.parse_args(argv)
    args.cpus = args.cpus or "0-{}".format(args.threads - 1)
    return args


def main(argv=None, deps=None, now=None):
    args = parse_args(argv)
    deps = deps or pafkit.load_deps()
    args.reference_id = deps.reference_id
    name, folder, samples, ref_contigs, other_contig, gfa = plan(args, deps)
    assignment = args.assignment or newest_assignment(args.scratch, args.pairwise_version)
    selections = {}
    if args.chromosome:
        if not assignment or not Path(assignment).is_file():
            raise SystemExit("no assignment.tsv found; run summarize_pairwise.py first or give --assignment")
        rows = pafkit.read_table(assignment)
        for assembly, _ in samples[1:]:
            selections[assembly] = assigned_sequences(rows, assembly, args.preset, args.chromosome)
        selections[samples[0][0]] = list(ref_contigs)
    binds = {str(folder)} | ({str(p.resolve().parent) for _, p in samples} if not args.chromosome else set())
    cmd = cactus_command(args, folder, name, ref_contigs, other_contig, gfa, binds)
    print("build_graph.py: run {} in {}".format(name, folder))
    for assembly, path in samples:
        extra = " -> {} sequence(s) for {}".format(len(selections[assembly]), args.chromosome) if args.chromosome else ""
        print("  {:<20} {}{}".format(assembly, path, extra))
    print("  " + " ".join(cmd))
    if args.chromosome and any(not selections[a] for a, _ in samples):
        print("STOP: no sequence assigned to {} for: {}".format(args.chromosome, ", ".join(
            a for a, _ in samples if not selections[a])), file=sys.stderr)
        return 2
    if args.restart:
        if not (folder / "jobstore").exists():
            print("STOP: --restart needs the job store {}, which does not exist".format(folder / "jobstore"), file=sys.stderr)
            return 2
    elif folder.exists():
        print("STOP: {} already exists, and this script never overwrites. To resume an interrupted run, add "
              "--restart; to start again, remove the folder by hand.".format(folder), file=sys.stderr)
        return 2
    if not args.go:
        print("Dry run: nothing has been written or run. Add --go to run it.")
        return 0

    inputs = []
    if not args.restart:
        (folder / "work").mkdir(parents=True)
        (folder / "out").mkdir()
        lines = []
        for assembly, path in samples:
            if args.chromosome:
                (folder / "fasta").mkdir(exist_ok=True)
                fasta = folder / "fasta" / "{}.fa".format(assembly)
                found = extract(path, selections[assembly], fasta)
                if sorted(found) != sorted(selections[assembly]):
                    print("STOP: {} lacks {}".format(path, set(selections[assembly]) - set(found)), file=sys.stderr)
                    return 1
                inputs.append({"assembly": assembly, "fasta": str(fasta), "sequences": found, "sha256": sha256_file(fasta),
                               "from": str(path)})
                lines.append("{}\t{}".format(assembly, fasta))
            else:
                real = path.resolve()
                inputs.append({"assembly": assembly, "fasta": str(real), "sha256_recorded": record_sha(path)})
                lines.append("{}\t{}".format(assembly, real))
        (folder / "seqfile.txt").write_text("\n".join(lines) + "\n")
    record = {"name": name, "mode": "chromosome" if args.chromosome else "whole genome", "chromosome": args.chromosome,
              "assemblies": [a for a, _ in samples], "image": str(args.sif), "command": cmd, "restart": args.restart,
              "assignment": str(assignment) if args.chromosome else "",
              "assignment_sha256": sha256_file(assignment) if args.chromosome else "", "inputs": inputs}
    start, clock = datetime.datetime.now(datetime.timezone.utc), time.monotonic()
    with open(folder / "run.log", "ab") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
    code = os.waitstatus_to_exitcode(status) if hasattr(os, "waitstatus_to_exitcode") else os.WEXITSTATUS(status)
    record.update(start_utc=start.strftime("%Y-%m-%dT%H:%M:%SZ"),
                  end_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                  wall_s=round(time.monotonic() - clock, 1), cpu_s=round(usage.ru_utime + usage.ru_stime, 1),
                  largest_process_peak_gb=round(usage.ru_maxrss / 1024 / 1024, 2), exit_status=code)
    runs = folder / "run.json"
    history = json.loads(runs.read_text()) if runs.exists() else []
    history.append(record)
    runs.write_text(json.dumps(history, indent=2) + "\n")
    print("{} {}: exit {}, {:.0f} min wall, {:.0f} CPU-hours, largest process {} GB. Log: {}".format(
        "done" if code == 0 else "FAILED", name, code, record["wall_s"] / 60, record["cpu_s"] / 3600,
        record["largest_process_peak_gb"], folder / "run.log"))
    return 0 if code == 0 else 1


def record_sha(fasta):
    try:
        provenance = json.loads((Path(fasta).parent / "provenance.json").read_text())
        return next(f["sha256"] for f in provenance.get("files", []) if f.get("role") == "genome")
    except (OSError, ValueError, StopIteration, KeyError):
        return ""


if __name__ == "__main__":
    sys.exit(main())
