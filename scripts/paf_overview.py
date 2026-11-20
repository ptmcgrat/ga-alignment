#!/usr/bin/env python3
"""A quick look at one PAF: what each target chromosome received. Read-only: it prints and writes nothing.

  paf_overview.py PAF [--labels TSV] [--chromosome NAME_OR_LABEL] [--min-length BP]

Table 1, one row per target sequence that has primary alignments:
  label (from --labels), length, primary records, bases covered by primary blocks (overlaps merged) and
  their fraction, identity (residue matches / alignment block length, primary records), bases in
  minus-strand primary blocks, the query sequence with most aligned bases, and how many query sequences
  have more than 1 Mb aligned there.
Then the record types (the tp tag: P primary, S secondary, ...) and how many records lack a cs tag.

With --chromosome, table 2 lists that target's primary blocks of at least --min-length bp (default 50 kb)
in target order: target start and end, query and its start and end, strand, length and identity.
The labels file is a sequence table from check_inputs.py (columns name and chr_name), for example
<scratch>/inputs/<timestamp>/sequences/MZ_NMK2024a.tsv; without it, sequences are shown by name only.
PAF coordinates are 0-based, half-open; this script prints them as they are.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


def read_labels(path):
    if not path:
        return {}
    with open(path) as handle:
        return {row["name"]: row.get("chr_name", "") for row in csv.DictReader(handle, delimiter="\t")}


def parse_record(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    tags = {}
    for tag in fields[12:]:
        parts = tag.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return {"qname": fields[0], "qlen": int(fields[1]), "qstart": int(fields[2]), "qend": int(fields[3]),
            "strand": fields[4], "tname": fields[5], "tlen": int(fields[6]), "tstart": int(fields[7]),
            "tend": int(fields[8]), "matches": int(fields[9]), "block": int(fields[10]),
            "tp": tags.get("tp", ""), "has_cs": "cs" in tags}


def merged_length(intervals):
    total, end = 0, -1
    for start, stop in sorted(intervals):
        if start > end:
            total += stop - start
            end = stop
        elif stop > end:
            total += stop - end
            end = stop
    return total


def summarize(records):
    """Per-target figures from primary records, plus record-type counts and records without cs."""
    per_target = {}
    tp_counts, missing_cs = defaultdict(int), 0
    for r in records:
        tp_counts[r["tp"] or "(none)"] += 1
        missing_cs += not r["has_cs"]
        if r["tp"] != "P":
            continue
        t = per_target.setdefault(r["tname"], {"tlen": r["tlen"], "records": 0, "intervals": [], "matches": 0,
                                               "block": 0, "minus_bp": 0, "by_query": defaultdict(int)})
        t["records"] += 1
        t["intervals"].append((r["tstart"], r["tend"]))
        t["matches"] += r["matches"]
        t["block"] += r["block"]
        if r["strand"] == "-":
            t["minus_bp"] += r["tend"] - r["tstart"]
        t["by_query"][r["qname"]] += r["tend"] - r["tstart"]
    for t in per_target.values():
        t["covered"] = merged_length(t["intervals"])
    return per_target, dict(tp_counts), missing_cs


def label_order(name, labels):
    label = labels.get(name, "")
    match = re.fullmatch(r"(?:LG)?(\d+)", label)
    return (0, int(match.group(1)), name) if match else (1, 0, name)


def find_target(wanted, per_target, labels):
    if wanted in per_target:
        return wanted
    matches = [n for n in per_target if labels.get(n) == wanted]
    return matches[0] if len(matches) == 1 else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("paf")
    parser.add_argument("--labels", help="a sequence table from check_inputs.py (name, chr_name)")
    parser.add_argument("--chromosome", help="list this target's primary blocks (a name, or a label such as LG9)")
    parser.add_argument("--min-length", type=int, default=50000, help="shortest block listed (default: %(default)s)")
    args = parser.parse_args(argv)
    labels = read_labels(args.labels)
    with open(args.paf) as handle:
        records = [r for r in (parse_record(line) for line in handle) if r]
    per_target, tp_counts, missing_cs = summarize(records)

    print("{}: {} records".format(args.paf, len(records)))
    print("\t".join(["target", "label", "length", "primary_records", "covered_bp", "covered", "identity",
                     "minus_strand_bp", "main_query", "queries_over_1Mb"]))
    for name in sorted(per_target, key=lambda n: label_order(n, labels)):
        t = per_target[name]
        main_query = max(t["by_query"], key=t["by_query"].get)
        print("\t".join(str(x) for x in [
            name, labels.get(name, ""), t["tlen"], t["records"], t["covered"],
            "{:.3f}".format(t["covered"] / t["tlen"]) if t["tlen"] else "",
            "{:.4f}".format(t["matches"] / t["block"]) if t["block"] else "",
            t["minus_bp"], main_query, sum(1 for bp in t["by_query"].values() if bp > 1_000_000)]))
    print("record types (tp): {}".format(", ".join("{} {}".format(k, v) for k, v in sorted(tp_counts.items()))))
    print("records without a cs tag: {}".format(missing_cs))

    if args.chromosome:
        target = find_target(args.chromosome, per_target, labels)
        if target is None:
            print("no single target named or labelled {!r} has primary blocks".format(args.chromosome), file=sys.stderr)
            return 1
        print("\nprimary blocks on {} ({}) of at least {} bp, in target order:".format(
            target, labels.get(target, "no label"), args.min_length))
        print("\t".join(["tstart", "tend", "query", "qstart", "qend", "strand", "length", "identity"]))
        blocks = sorted((r for r in records if r["tp"] == "P" and r["tname"] == target
                         and r["tend"] - r["tstart"] >= args.min_length), key=lambda r: r["tstart"])
        for r in blocks:
            print("\t".join(str(x) for x in [r["tstart"], r["tend"], r["qname"], r["qstart"], r["qend"], r["strand"],
                                             r["tend"] - r["tstart"], "{:.4f}".format(r["matches"] / r["block"])]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
