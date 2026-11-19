#!/usr/bin/env python3
"""Gate 1's input check: is every genome the project will align present, intact and as recorded?

For each assembly of the panel (read with ga_common.panel), plus the supplementary Y haplotig YH_LG10Y:
  - locate genome.fna.gz through ga_common.paths (never a hand-built path) and provenance.json beside it
    (spec, section 2: "provenance.json in the same folder");
  - check that provenance.json belongs to this assembly and accession, and that the file's sha256 and size
    match what provenance.json recorded when the coordinator fetched it;
  - read every sequence name and length from the FASTA (only the first word of a header, as the tools see
    it) and compare them with provenance.json: the same names, the same lengths;
  - write one table per assembly: name, length, NCBI chromosome label, role, molecule type, N bases and
    lower-case (soft-masked) bases.
Across the panel it also writes:
  - summary.tsv: one row per assembly (accessions requested and fetched, NCBI assembly name, sequence
    count, total length, label style, N and lower-case content, check results);
  - reference_refContigs.tsv and reference_otherContigs.tsv: the reference's linkage groups (the candidate
    --refContigs for Cactus) and everything else in it (--otherContig);
  - panel_chromosomes_check.tsv: panel_chromosomes_P1.tsv against the labels in provenance.json;
  - report.md: what passed, what did not, and what needs a decision.

A dry run unless --go: it then only lists the assemblies, where their files are, and the folder it would
write. With --go it writes one new folder, <scratch>/inputs/<UTC timestamp>/, and never overwrites.
It reads the genomes and writes nothing anywhere else. Exit status: 0 if every check passed, 1 if any
failed, 2 if it refused to overwrite.
"""

import argparse
import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-alignment"
# Held beside the panel (spec, section 2); aligned after the 24 panel assemblies (spec, section 6).
SUPPLEMENTARY = [("YH_LG10Y", "h2tg000014l")]
PANEL_LG_COLUMNS = ["LG2", "LG9", "LG10", "LG11", "LG13", "LG20"]
EXPECTED_LINKAGE_GROUPS = 22  # LG1-LG20, LG22, LG23 (spec, section 2)
LOWER = bytes(range(ord("a"), ord("z") + 1))
SEQUENCE_COLUMNS = ["name", "length", "chr_name", "role", "molecule_type", "n_bases", "lowercase_bases"]
SUMMARY_COLUMNS = ["assembly_id", "panel_role", "accession_requested", "accession_fetched", "kind",
                   "assembly_name", "sequences", "total_bp", "chromosome_sequences", "label_style",
                   "n_bp", "lowercase_fraction", "sha256_ok", "names_lengths_ok", "problems"]


# ---------------------------------------------------------------- helpers (pure, tested)

def sha256_file(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta_stats(path):
    """[(name, length, n_bases, lowercase_bases)] for a FASTA, gzip/BGZF or plain.
    The name is the first word of the header, which is all the tools see. Lower-case counts include n."""
    opener = gzip.open if str(path).endswith(".gz") else open
    stats, name, length, n_bases, lower = [], None, 0, 0, 0
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    stats.append((name, length, n_bases, lower))
                words = line[1:].split(None, 1)
                name = words[0].decode() if words else ""
                length = n_bases = lower = 0
            else:
                seq = line.rstrip(b"\r\n")
                length += len(seq)
                n_bases += seq.count(b"N") + seq.count(b"n")
                lower += len(seq) - len(seq.translate(None, LOWER))
    if name is not None:
        stats.append((name, length, n_bases, lower))
    return stats


def label_style(chr_names):
    """'LG' if every chromosome label is LGn, 'num' if every one is a plain number, else 'mixed' or 'none'.
    The mitochondrion (MT) is left out."""
    labels = [c for c in chr_names if c and c.upper() not in ("MT", "MITO")]
    if not labels:
        return "none"
    if all(re.fullmatch(r"LG\d+", c) for c in labels):
        return "LG"
    if all(re.fullmatch(r"\d+", c) for c in labels):
        return "num"
    return "mixed"


def assembly_name(original_name, accession):
    """The NCBI assembly name inside a downloaded file name, e.g. GCF_041146795.1_Mzebra_GT3a_genomic.fna."""
    match = re.fullmatch(re.escape(accession) + r"_(.+)_genomic\.fna(\.gz)?", original_name or "")
    return match.group(1) if match else ""


def is_chromosome(row):
    return row.get("role") == "assembled-molecule" and (row.get("chr_name") or "").upper() not in ("MT", "MITO")


def reference_lists(sequence_rows):
    """(linkage groups, everything else) of the reference: the --refContigs and --otherContig candidates."""
    groups = [r for r in sequence_rows if is_chromosome(r) and r.get("molecule_type") == "Linkage Group"]
    names = {r["name"] for r in groups}
    return groups, [r for r in sequence_rows if r["name"] not in names]


def check_panel_chromosomes(panel_chromosome_rows, results):
    """Each accession that panel_chromosomes_P1.tsv lists must exist in that assembly with that label."""
    rows = []
    for entry in panel_chromosome_rows:
        result = results.get(entry["assembly_id"])
        if result is None:
            continue
        by_name = {r["name"]: r for r in result["sequences"]}
        for column in PANEL_LG_COLUMNS:
            accession = (entry.get(column) or "").strip()
            if not accession:
                continue
            number = column[2:]
            found = by_name.get(accession)
            label = (found or {}).get("chr_name", "")
            ok = found is not None and label in (column, number)
            rows.append({"assembly_id": entry["assembly_id"], "chromosome": column, "accession": accession,
                         "label_in_provenance": label if found else "(not in this assembly)",
                         "label_style_listed": entry.get("label_style", ""),
                         "label_style_found": result["summary"]["label_style"],
                         "status": "ok" if ok else "MISMATCH"})
    return rows


def check_one(assembly_id, accession_version, fasta_path, panel_role):
    """Every check for one assembly. Runs in a worker process; returns plain data."""
    problems, sequences = [], []
    summary = {"assembly_id": assembly_id, "panel_role": panel_role, "accession_requested": accession_version,
               "accession_fetched": "", "kind": "", "assembly_name": "", "sequences": 0, "total_bp": 0,
               "chromosome_sequences": 0, "label_style": "", "n_bp": 0, "lowercase_fraction": "",
               "sha256_ok": "no", "names_lengths_ok": "no", "problems": 0}
    result = {"assembly_id": assembly_id, "fasta": str(fasta_path), "summary": summary,
              "sequences": sequences, "problems": problems}
    fasta = Path(fasta_path)
    provenance_path = fasta.parent / "provenance.json"
    if not fasta.is_file():
        problems.append("genome FASTA not found: {}".format(fasta))
    if not provenance_path.is_file():
        problems.append("provenance.json not found beside the FASTA")
    if problems:
        summary["problems"] = len(problems)
        return result
    provenance = json.loads(provenance_path.read_text())
    if provenance.get("assembly_id") != assembly_id:
        problems.append("provenance.json names assembly {!r}".format(provenance.get("assembly_id")))
    if provenance.get("accession_requested") != accession_version:
        problems.append("provenance.json was requested as {!r}".format(provenance.get("accession_requested")))
    summary["accession_fetched"] = provenance.get("accession_fetched", "")
    summary["kind"] = provenance.get("kind", "")
    genome_files = [f for f in provenance.get("files", []) if f.get("role") == "genome"]
    if len(genome_files) != 1:
        problems.append("provenance.json lists {} genome files, not 1".format(len(genome_files)))
    else:
        recorded = genome_files[0]
        summary["assembly_name"] = assembly_name(recorded.get("original_name", ""), summary["accession_fetched"])
        if fasta.stat().st_size != recorded.get("bytes"):
            problems.append("size {} differs from the recorded {}".format(fasta.stat().st_size, recorded.get("bytes")))
        if sha256_file(fasta) == recorded.get("sha256"):
            summary["sha256_ok"] = "yes"
        else:
            problems.append("sha256 differs from provenance.json")

    stats = read_fasta_stats(fasta)
    recorded_seqs = {s["name"]: s for s in provenance.get("sequences", [])}
    seen = set()
    mismatches = []
    for name, length, n_bases, lower in stats:
        record = recorded_seqs.get(name)
        if name in seen:
            mismatches.append("{} appears twice in the FASTA".format(name))
        seen.add(name)
        if record is None:
            mismatches.append("{} is in the FASTA but not in provenance.json".format(name))
        elif record.get("length") != length:
            mismatches.append("{}: FASTA {} bp, provenance.json {} bp".format(name, length, record.get("length")))
        record = record or {}
        sequences.append({"name": name, "length": length, "chr_name": record.get("chr_name", ""),
                          "role": record.get("role", ""), "molecule_type": record.get("molecule_type", ""),
                          "n_bases": n_bases, "lowercase_bases": lower})
    for name in recorded_seqs:
        if name not in seen:
            mismatches.append("{} is in provenance.json but not in the FASTA".format(name))
    if mismatches:
        problems.append("{} name/length differences, first: {}".format(len(mismatches), "; ".join(mismatches[:3])))
    else:
        summary["names_lengths_ok"] = "yes"

    total = sum(s["length"] for s in sequences)
    summary["sequences"] = len(sequences)
    summary["total_bp"] = total
    summary["chromosome_sequences"] = sum(1 for s in sequences if is_chromosome(s))
    summary["label_style"] = label_style([s["chr_name"] for s in sequences if s["role"] == "assembled-molecule"])
    summary["n_bp"] = sum(s["n_bases"] for s in sequences)
    summary["lowercase_fraction"] = "{:.3f}".format(sum(s["lowercase_bases"] for s in sequences) / total) if total else ""
    summary["problems"] = len(problems)
    return result


# ---------------------------------------------------------------- the coordinator's library

class Deps:
    """What the check needs from ga_common, in a form tests can stand in for."""

    def __init__(self, locate, rows, panel_version, reference_id, panel_chromosomes):
        self.locate, self.rows, self.panel_version = locate, rows, panel_version
        self.reference_id, self.panel_chromosomes = reference_id, panel_chromosomes


def load_deps():
    from ga_common.paths import Paths
    from ga_common import panel
    paths = Paths.load()
    pan = panel.read_panel(str(paths.config_dir / "genome_panel_P1.tsv"))

    def locate(assembly_id, accession_version):
        return Path(paths.local("genome_fasta", assembly_id=assembly_id, accession_version=accession_version))
    return Deps(locate, list(pan.rows), pan.version, pan.reference()["id"],
                Path(paths.config_dir) / "panel_chromosomes_P1.tsv")


def add_throttle_arguments(parser):
    try:
        from ga_common import throttle
        throttle.add_arguments(parser)
    except ImportError:  # the laptop's tests run without ga_common
        parser.add_argument("--threads", type=int, default=None)
        parser.add_argument("--max-jobs", dest="max_jobs", type=int, default=None)


# ---------------------------------------------------------------- outputs

def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_report(outdir, deps, results, chromosome_rows, groups, others, version):
    summaries = [r["summary"] for r in results.values()]
    failed = [r for r in results.values() if r["problems"]]
    lines = ["# Input check ({})".format(version), "",
             "- written: {} UTC; panel {}; reference {}".format(outdir.name, deps.panel_version, deps.reference_id),
             "- assemblies checked: {} ({} with a problem)".format(len(results), len(failed)),
             "- sha256 as recorded: {} of {}".format(sum(s["sha256_ok"] == "yes" for s in summaries), len(summaries)),
             "- names and lengths as recorded: {} of {}".format(
                 sum(s["names_lengths_ok"] == "yes" for s in summaries), len(summaries)), ""]
    lines += ["## Problems", ""] + (["- **{}**: {}".format(r["assembly_id"], p) for r in failed for p in r["problems"]]
                                    or ["- none"])
    styles = {}
    for s in summaries:
        styles.setdefault(s["label_style"], []).append(s["assembly_id"])
    lines += ["", "## Chromosome labels", ""] + ["- {} ({}): {}".format(k, len(v), ", ".join(sorted(v)))
                                                  for k, v in sorted(styles.items())]
    mismatches = [r for r in chromosome_rows if r["status"] != "ok"]
    lines += ["", "panel_chromosomes_P1.tsv against provenance.json: {} entries, {} mismatches.".format(
        len(chromosome_rows), len(mismatches))]
    lines += ["- {assembly_id} {chromosome}: {accession} is labelled {label_in_provenance}".format(**r) for r in mismatches]
    lines += ["", "## The reference ({})".format(deps.reference_id), "",
              "- linkage groups (candidate `--refContigs`): {}, {:,} bp".format(len(groups), sum(r["length"] for r in groups)),
              "- other sequences (`--otherContig`): {}, {:,} bp".format(len(others), sum(r["length"] for r in others))]
    if len(groups) != EXPECTED_LINKAGE_GROUPS:
        lines.append("- **expected {} linkage groups**".format(EXPECTED_LINKAGE_GROUPS))
    lines += ["", "## Per assembly", "",
              "| assembly | NCBI assembly name | sequences | total Mb | chromosomes | labels | N Mb | lower-case |",
              "|---|---|---|---|---|---|---|---|"]
    for s in sorted(summaries, key=lambda x: x["assembly_id"]):
        lines.append("| {} | {} | {} | {:.1f} | {} | {} | {:.2f} | {} |".format(
            s["assembly_id"], s["assembly_name"] or "-", s["sequences"], s["total_bp"] / 1e6,
            s["chromosome_sequences"], s["label_style"], s["n_bp"] / 1e6, s["lowercase_fraction"]))
    (outdir / "report.md").write_text("\n".join(lines) + "\n")


def read_version():
    try:
        return (REPO / "VERSION").read_text().strip()
    except OSError:
        return "unknown version"


# ---------------------------------------------------------------- main

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH,
                        help="the project's scratch folder (default: %(default)s)")
    parser.add_argument("--no-supplementary", action="store_true",
                        help="leave out the supplementary YH_LG10Y")
    parser.add_argument("--go", action="store_true", help="run the checks and write the results")
    add_throttle_arguments(parser)
    return parser.parse_args(argv)


def main(argv=None, deps=None, now=None):
    args = parse_args(argv)
    deps = deps or load_deps()
    version = read_version()
    stamp = (now or datetime.datetime.now(datetime.timezone.utc)).strftime("%Y%m%d-%H%M%S")
    outdir = Path(args.scratch).expanduser() / "inputs" / stamp
    items = [(r["id"], r["accession_genbank"], "reference" if r["id"] == deps.reference_id else "panel")
             for r in deps.rows]
    if not args.no_supplementary:
        items += [(a, acc, "supplementary") for a, acc in SUPPLEMENTARY]
    located = [(a, acc, str(deps.locate(a, acc)), role) for a, acc, role in items]
    if not args.go:
        print("check_inputs.py ({}): dry run. Panel {}, reference {}.".format(version, deps.panel_version, deps.reference_id))
        for a, _, path, role in located:
            print("  {:<20} {:<13} {}".format(a, role, path))
        print("With --go it would check these {} assemblies and write a new folder:\n  {}".format(len(located), outdir))
        print("Nothing has been read or written. Add --go to run the checks.")
        return 0
    if outdir.exists():
        print("STOP: {} already exists, and this script never overwrites.".format(outdir), file=sys.stderr)
        print("To redo it, remove it by hand:  rm -r {}".format(outdir), file=sys.stderr)
        return 2
    jobs = getattr(args, "max_jobs", None) or min(len(located), os.cpu_count() or 1)
    jobs = max(1, min(jobs, len(located)))
    results = {}
    if jobs == 1:
        for item in located:
            result = check_one(*item)
            results[result["assembly_id"]] = result
            print("  checked {}: {} problems".format(result["assembly_id"], len(result["problems"])), flush=True)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
            for result in pool.map(check_one, *zip(*located)):
                results[result["assembly_id"]] = result
                print("  checked {}: {} problems".format(result["assembly_id"], len(result["problems"])), flush=True)

    outdir.mkdir(parents=True)
    (outdir / "sequences").mkdir()
    for assembly_id, result in results.items():
        write_tsv(outdir / "sequences" / (assembly_id + ".tsv"), SEQUENCE_COLUMNS, result["sequences"])
    write_tsv(outdir / "summary.tsv", SUMMARY_COLUMNS, [r["summary"] for r in results.values()])
    with open(deps.panel_chromosomes) as handle:
        chromosome_rows = check_panel_chromosomes(list(csv.DictReader(handle, delimiter="\t")), results)
    write_tsv(outdir / "panel_chromosomes_check.tsv",
              ["assembly_id", "chromosome", "accession", "label_in_provenance", "label_style_listed",
               "label_style_found", "status"], chromosome_rows)
    reference = results.get(deps.reference_id, {"sequences": []})
    groups, others = reference_lists(reference["sequences"])
    write_tsv(outdir / "reference_refContigs.tsv", SEQUENCE_COLUMNS, groups)
    write_tsv(outdir / "reference_otherContigs.tsv", SEQUENCE_COLUMNS, others)
    write_report(outdir, deps, results, chromosome_rows, groups, others, version)

    failed = [a for a, r in results.items() if r["problems"]]
    mismatched = [r for r in chromosome_rows if r["status"] != "ok"]
    print("Results written to {}".format(outdir))
    print("{} assemblies: {} with problems; {} chromosome-label mismatches; reference: {} linkage groups.".format(
        len(results), len(failed), len(mismatched), len(groups)))
    for assembly_id in failed:
        print("  - {}: {}".format(assembly_id, "; ".join(results[assembly_id]["problems"])))
    return 1 if failed or mismatched or len(groups) != EXPECTED_LINKAGE_GROUPS else 0


if __name__ == "__main__":
    sys.exit(main())
