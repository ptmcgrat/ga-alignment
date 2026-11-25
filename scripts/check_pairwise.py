#!/usr/bin/env python3
"""The exact and accounting checks of spec section 5, for every PAF of a pairwise version folder.

For each run in <scratch>/pairwise/<version>/manifest.tsv that finished (exit status 0), it reads both genomes
(located through ga_common) and every record of the PAF, and checks:
  names and lengths   each query and target name exists in its genome, with the length the PAF gives,
                      and that length matches provenance.json
  coordinates         0 <= start < end <= length, on both sides
  cs counts           the cs string spans exactly the record's target and query intervals; its matched bases, less
                      N facing N inside matched runs, equal column 10 (residue matches); and its columns minus the
                      ambiguous bases minimap2 reports (the nn tag) equal column 11 (block length). Both columns leave
                      out ambiguous bases, such as the 100-N runs that stand for assembly gaps. The N columns are also
                      counted from the genomes and reported beside the nn total
  cs sequence         the cs string rebuilds the aligned stretch base for base: matched runs and substitutions
                      against both genomes, deletions against the target, insertions against the query
                      (the query reverse-complemented for '-' records); compared in upper case
  cs and CIGAR        the CIGAR (cg tag) is exactly what the cs string implies
Accounting: every base of both genomes is put in one class, in this order: aligned by a primary record (tp P or I),
aligned only by secondary records (tp S or i), N, or unaligned. The classes add up to the genome, per sequence.
It also checks that every panel assembly has a forward PAF with the given preset.

A dry run unless --go. With --go it writes one new folder, <version folder>/checks/<UTC timestamp>/, with
exact.tsv, accounting.tsv, accounting_sequences.tsv and report.md, and never overwrites. Exit status: 0 if every
exact check passed and every panel assembly has a PAF, 1 otherwise, 2 if it refused to start.
"""

import argparse
import collections
import concurrent.futures
import datetime
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pafkit  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-alignment"
KINDS = ["names_lengths", "coordinates", "cs_missing", "cs_counts", "cs_sequence", "cs_cigar"]
EXACT_COLUMNS = (["run", "direction", "preset", "query_id", "target_id", "records", "ambiguous_nn", "n_columns"]
                 + KINDS + ["status", "examples"])
ACCOUNT_COLUMNS = ["run", "side", "assembly", "sequences", "length", "primary_bp", "secondary_only_bp", "n_bp",
                   "unaligned_bp", "primary_fraction"]
SEQ_COLUMNS = ["run", "side", "assembly", "sequence", "label", "length", "primary_bp", "secondary_only_bp", "n_bp",
               "unaligned_bp"]


def check_record(rec, query, target, query_prov, target_prov):
    """[(kind, detail)] for one record; an empty list when it passes every exact check."""
    qseq, tseq = query.get(rec["qname"]), target.get(rec["tname"])
    for name, seq, length, prov in ((rec["qname"], qseq, rec["qlen"], query_prov),
                                    (rec["tname"], tseq, rec["tlen"], target_prov)):
        if seq is None:
            return [("names_lengths", "{} is not in its genome".format(name))]
        if len(seq) != length or prov.get(name, {}).get("length") != length:
            return [("names_lengths", "{}: PAF {} bp, FASTA {} bp, provenance {} bp".format(
                name, length, len(seq), prov.get(name, {}).get("length")))]
    if not (0 <= rec["qs"] < rec["qe"] <= rec["qlen"] and 0 <= rec["ts"] < rec["te"] <= rec["tlen"]):
        return [("coordinates", "query {}-{}, target {}-{}".format(rec["qs"], rec["qe"], rec["ts"], rec["te"]))]
    if not rec["cs"]:
        return [("cs_missing", "no cs tag")]
    problems = []
    ops = pafkit.cs_ops(rec["cs"])
    segment = qseq[rec["qs"]:rec["qe"]]
    if rec["strand"] == "-":
        segment = pafkit.revcomp(segment)
    t, q, matched, columns, n_columns, n_matched, bad = rec["ts"], 0, 0, 0, 0, 0, None
    for kind, value in ops:
        if kind == ":":
            if bad is None and tseq[t:t + value] != segment[q:q + value]:
                bad = "matched run differs at target {}".format(t)
            n_here = tseq.count(b"N", t, t + value)      # N facing N: equal bases, but not a residue match
            n_columns, n_matched = n_columns + n_here, n_matched + n_here
            t, q, matched, columns = t + value, q + value, matched + value, columns + value
        elif kind == "*":
            ref, alt = value
            if bad is None and (tseq[t:t + 1] != ref or segment[q:q + 1] != alt):
                bad = "substitution differs at target {}".format(t)
            n_columns += ref == b"N" or alt == b"N"
            t, q, columns = t + 1, q + 1, columns + 1
        elif kind == "+":
            if bad is None and segment[q:q + len(value)] != value:
                bad = "insertion differs at target {}".format(t)
            n_columns += value.count(b"N")
            q, columns = q + len(value), columns + len(value)
        elif kind == "-":
            if bad is None and tseq[t:t + len(value)] != value:
                bad = "deletion differs at target {}".format(t)
            n_columns += value.count(b"N")
            t, columns = t + len(value), columns + len(value)
        else:
            t += value
    rec["n_columns"] = n_columns
    if t != rec["te"] or q != rec["qe"] - rec["qs"]:
        problems.append(("cs_counts", "cs spans target {} and query {} bp; the record {} and {}".format(
            t - rec["ts"], q, rec["te"] - rec["ts"], rec["qe"] - rec["qs"])))
    elif matched - n_matched != rec["matches"] or columns - rec["nn"] != rec["block"]:
        problems.append(("cs_counts", "cs gives {} residue matches ({} in matched runs, {} of them N facing N) in {} columns, "
                         "{} ambiguous (nn); columns 10 and 11 say {} and {}".format(
                             matched - n_matched, matched, n_matched, columns, rec["nn"], rec["matches"], rec["block"])))
    if bad:
        problems.append(("cs_sequence", bad))
    if not rec["cg"]:
        problems.append(("cs_cigar", "no cg tag"))
    elif pafkit.cigar_from_cs(ops) != pafkit.parse_cigar(rec["cg"]):
        problems.append(("cs_cigar", "the CIGAR is not what the cs string implies"))
    return problems


def account(side, assembly, genome, provenance, primary, secondary, run):
    rows = []
    for name, seq in genome.items():
        length = len(seq)
        prim = pafkit.merge(primary.get(name, []))
        sec_only = pafkit.subtract(pafkit.merge(secondary.get(name, [])), prim)
        unaligned = pafkit.subtract(pafkit.subtract([(0, length)], prim), pafkit.merge(sec_only))
        n_bp = sum(seq[a:b].count(b"N") for a, b in unaligned)
        rows.append({"run": run, "side": side, "assembly": assembly, "sequence": name,
                     "label": provenance.get(name, {}).get("chr_name", ""), "length": length,
                     "primary_bp": pafkit.total(prim), "secondary_only_bp": pafkit.total(sec_only), "n_bp": n_bp,
                     "unaligned_bp": pafkit.total(unaligned) - n_bp})
    return rows


def check_run(row, paf, query_fasta, target_fasta):
    """Every exact check and the accounting for one run. Runs in a worker process; returns plain data."""
    query, target = pafkit.read_fasta(query_fasta), pafkit.read_fasta(target_fasta)
    query_prov, target_prov = pafkit.provenance_sequences(query_fasta), pafkit.provenance_sequences(target_fasta)
    counts, examples, records, ambiguous, n_columns = collections.Counter(), [], 0, 0, 0
    spans = {"q": (collections.defaultdict(list), collections.defaultdict(list)),
             "t": (collections.defaultdict(list), collections.defaultdict(list))}
    for rec in pafkit.read_paf(paf, keep_strings=True):
        records += 1
        problems = check_record(rec, query, target, query_prov, target_prov)
        ambiguous += rec["nn"]
        n_columns += rec.get("n_columns", 0)
        for kind in {k for k, _ in problems}:
            counts[kind] += 1
        for kind, detail in problems:
            if len(examples) < 5:
                examples.append("record {} ({}:{}-{}): {}".format(records, rec["qname"], rec["qs"], rec["qe"], detail))
        which = 0 if rec["tp"] in pafkit.PRIMARY else 1 if rec["tp"] in pafkit.SECONDARY else None
        if which is not None:
            spans["q"][which][rec["qname"]].append((rec["qs"], rec["qe"]))
            spans["t"][which][rec["tname"]].append((rec["ts"], rec["te"]))
    exact = {"run": row["run"], "direction": row["direction"], "preset": row["preset"], "query_id": row["query_id"],
             "target_id": row["target_id"], "records": records, "ambiguous_nn": ambiguous, "n_columns": n_columns,
             "examples": " | ".join(examples)}
    exact.update({k: counts.get(k, 0) for k in KINDS})
    exact["status"] = "pass" if records and not any(counts.values()) else "FAIL"
    sequences = (account("query", row["query_id"], query, query_prov, *spans["q"], row["run"])
                 + account("target", row["target_id"], target, target_prov, *spans["t"], row["run"]))
    return exact, sequences


def summarize_accounting(sequences):
    totals = collections.OrderedDict()
    for s in sequences:
        key = (s["run"], s["side"], s["assembly"])
        t = totals.setdefault(key, {"run": s["run"], "side": s["side"], "assembly": s["assembly"], "sequences": 0,
                                    "length": 0, "primary_bp": 0, "secondary_only_bp": 0, "n_bp": 0, "unaligned_bp": 0})
        t["sequences"] += 1
        for k in ("length", "primary_bp", "secondary_only_bp", "n_bp", "unaligned_bp"):
            t[k] += s[k]
    for t in totals.values():
        t["primary_fraction"] = "{:.4f}".format(t["primary_bp"] / t["length"]) if t["length"] else ""
    return list(totals.values())


def write_report(path, version, exact, accounting, missing, preset):
    failed = [e for e in exact if e["status"] != "pass"]
    lines = ["# Exact and accounting checks ({})".format(version), "",
             "- runs checked: {}; passed every exact check: {}; failed: {}".format(len(exact), len(exact) - len(failed), len(failed)),
             "- records checked: {:,}".format(sum(e["records"] for e in exact)),
             "- ambiguous bases in alignments: {:,} by minimap2's nn tags; {:,} N columns counted from the genomes".format(
                 sum(e["ambiguous_nn"] for e in exact), sum(e["n_columns"] for e in exact)),
             "- panel assemblies without a forward {} PAF: {}".format(preset, ", ".join(missing) or "none"), ""]
    if failed:
        lines += ["## Failed runs", ""] + ["- **{}**: {}".format(e["run"], "; ".join(
            "{} {}".format(e[k], k) for k in KINDS if e[k])) + (" — " + e["examples"] if e["examples"] else "") for e in failed] + [""]
    lines += ["## Accounting (every base in one class; the classes add up to each genome)", "",
              "| run | side | assembly | length Mb | primary | secondary only | N | unaligned |", "|---|---|---|---|---|---|---|---|"]
    for t in accounting:
        lines.append("| {} | {} | {} | {:.1f} | {:.1%} | {:.2%} | {:.2%} | {:.2%} |".format(
            t["run"], t["side"], t["assembly"], t["length"] / 1e6, t["primary_bp"] / t["length"],
            t["secondary_only_bp"] / t["length"], t["n_bp"] / t["length"], t["unaligned_bp"] / t["length"]))
    Path(path).write_text("\n".join(lines) + "\n")


def read_version():
    try:
        return (REPO / "VERSION").read_text().strip()
    except OSError:
        return "unknown version"


def main(argv=None, deps=None, now=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--version", default="PW_v1", help="the pairwise version folder (default: %(default)s)")
    parser.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH, help="the project's scratch folder")
    parser.add_argument("--preset", default="asm5", help="the preset whose forward PAFs every assembly needs (default: %(default)s)")
    parser.add_argument("--runs", nargs="+", metavar="RUN", help="check only these runs (names as in the manifest)")
    parser.add_argument("--go", action="store_true", help="run the checks and write the results")
    pafkit.add_throttle_arguments(parser)
    args = parser.parse_args(argv)
    deps = deps or pafkit.load_deps()
    base = Path(args.scratch).expanduser() / "pairwise" / args.version
    rows = [r for r in pafkit.read_manifest(base) if r["exit_status"] == "0" and (not args.runs or r["run"] in args.runs)]
    forward = {r["query_id"] for r in rows if r["direction"] == "forward" and r["preset"] == args.preset}
    missing = [a for a in deps.panel_ids() if a not in forward] if not args.runs else []
    stamp = (now or datetime.datetime.now(datetime.timezone.utc)).strftime("%Y%m%d-%H%M%S")
    outdir = base / "checks" / stamp
    jobs = max(1, min(getattr(args, "max_jobs", None) or 1, len(rows) or 1))
    print("check_pairwise.py ({}): {} runs in {}, {} at a time".format(read_version(), len(rows), base, jobs))
    for r in rows:
        print("  {:<32} {}  ({} -> {})".format(r["run"], r["paf"], r["query_id"], r["target_id"]))
    if missing:
        print("  panel assemblies without a forward {} PAF: {}".format(args.preset, ", ".join(missing)))
    if not args.go:
        print("Dry run: nothing has been read or written. Add --go to run the checks into\n  {}".format(outdir))
        return 0
    if outdir.exists():
        print("STOP: {} already exists, and this script never overwrites.".format(outdir), file=sys.stderr)
        return 2
    tasks = [(r, str(base / r["paf"]), str(deps.locate(r["query_id"])), str(deps.locate(r["target_id"]))) for r in rows]
    results = []
    if jobs == 1:
        results = [check_run(*t) for t in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
            results = list(pool.map(check_run, *zip(*tasks)))
    exact = [e for e, _ in results]
    for e in exact:
        print("  {:<6} {}: {:,} records".format(e["status"], e["run"], e["records"]), flush=True)
    sequences = [s for _, seqs in results for s in seqs]
    accounting = summarize_accounting(sequences)
    outdir.mkdir(parents=True)
    pafkit.write_tsv(outdir / "exact.tsv", EXACT_COLUMNS, exact)
    pafkit.write_tsv(outdir / "accounting.tsv", ACCOUNT_COLUMNS, accounting)
    pafkit.write_tsv(outdir / "accounting_sequences.tsv", SEQ_COLUMNS, sequences)
    write_report(outdir / "report.md", read_version(), exact, accounting, missing, args.preset)
    failed = [e["run"] for e in exact if e["status"] != "pass"]
    print("Results written to {}\n{} runs: {} failed an exact check; {} panel assemblies lack a forward {} PAF.".format(
        outdir, len(exact), len(failed), len(missing), args.preset))
    return 1 if failed or missing else 0


if __name__ == "__main__":
    sys.exit(main())
