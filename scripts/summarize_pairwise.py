#!/usr/bin/env python3
"""The known-feature and measured checks of spec section 5, for every PAF of a pairwise version folder.

Inversions (forward runs, the reference as target), for each unit of Supplementary file 4 (2, 9, 10, 11a, 11b, 13,
20a, 20b; data/kumar2025_suppfile4_breakpoints.tsv):
  blocks      Kumar et al.'s four filters: primary records only; identity (column 10 / column 11) >= 0.30;
              block length (column 11) >= 10,000 bp; and no block from a query sequence with under 2,000,000 bp
              of such blocks in total
  call        the unit's interior runs from the end of its first breakpoint interval to the start of its second.
              The query sequence with most interior bases is the unit's sequence; its background strand is its
              majority strand outside every unit on that chromosome. Inverted if >= 80% of its interior bases lie
              on the other strand, normal if <= 20%, otherwise unresolved; unresolved too when that sequence has
              under 1 Mb outside the units (interior and flanks on different sequences) or under 100 kb inside
  breakpoints for an inverted unit: collinear segments of that sequence (split where the strand changes or the
              query position jumps by more than 1 Mb; segments of >= 100 kb kept). The observed interval runs from
              the last other segment before the inverted segments to the first of them, and from the last of them
              to the next other segment. Scored by overlap with Supplementary file 4's interval; the distance of
              the inverted block ends from that interval is reported
  expected    data/expected_inversions.tsv (from Supplementary file 2). Scored rows decide the known-feature
              check: it fails only on a disagreement; unresolved units are reported, not failed
Measured, per forward run and reference sequence: coverage by primary records, BLAST-like identity (column 10 /
column 11, which counts gap bases), gap-compressed divergence (the de tag, weighted by aligned length), bases on
the minus strand, and the main query sequence with its label. Assignment: each query sequence's chromosome (most
aligned bases) against its own label. Reciprocity: coverage of the reference and the assembly in the forward and
reverse runs. Outliers: chromosomes whose coverage falls more than 5 points below the subset's median, or whose
divergence is more than twice it. Runs: time, memory and disk from the manifest.

A dry run unless --go. With --go it writes one new folder, <version folder>/summary/<UTC timestamp>/, and never
overwrites. Exit status: 0 if the known-feature check passed, 1 if a scored unit disagreed or nothing was scored,
2 if it refused to start.
"""

import argparse
import collections
import concurrent.futures
import datetime
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pafkit  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-alignment"
CHROM_COLUMNS = ["run", "target", "label", "length", "primary_records", "covered_bp", "covered", "blast_identity",
                 "divergence", "minus_strand_bp", "main_query", "main_query_label", "main_query_share", "label_agrees"]
ASSIGN_COLUMNS = ["run", "query", "query_label", "length", "aligned_bp", "chromosome", "chromosome_label", "share",
                  "second", "second_share", "label_agrees"]
SEG_COLUMNS = ["run", "target", "label", "query", "strand", "ts", "te", "qs", "qe", "blocks", "aligned_bp"]
INV_COLUMNS = ["run", "assembly", "preset", "unit", "chromosome", "main_query", "interior_bp", "interior_other_bp",
               "background_bp", "background_strand", "opposite_fraction", "call", "reason", "inverted_start",
               "inverted_end", "bp1_observed", "bp1_published", "bp1_overlap", "bp1_distance", "bp2_observed",
               "bp2_published", "bp2_overlap", "bp2_distance", "expected", "basis", "scored", "verdict"]


class Params:
    def __init__(self, args):
        for key in ("min_identity", "min_block", "min_contig", "min_segment", "max_jump", "call_fraction",
                    "min_background", "min_interior"):
            setattr(self, key, getattr(args, key))


# ---------------------------------------------------------------- inversions

def kumar_blocks(records, p):
    keep = [r for r in records if r["tp"] in pafkit.PRIMARY and r["block"] >= p.min_block
            and r["matches"] >= p.min_identity * r["block"]]
    per_query = collections.Counter()
    for r in keep:
        per_query[r["qname"]] += r["block"]
    return [r for r in keep if per_query[r["qname"]] >= p.min_contig]


def segments(blocks, max_jump):
    """Collinear segments, built per query sequence: a new one starts where the strand changes or the query
    position jumps by more than max_jump."""
    out = []
    by_query = collections.defaultdict(list)
    for b in blocks:
        by_query[b["qname"]].append(b)
    for qname, group in by_query.items():
        current = None
        for b in sorted(group, key=lambda r: (r["ts"], r["te"])):
            if current and current["strand"] == b["strand"]:
                gap = b["qs"] - current["last"][1] if b["strand"] == "+" else current["last"][0] - b["qe"]
                if -max_jump <= gap <= max_jump:
                    current.update(te=max(current["te"], b["te"]), qs=min(current["qs"], b["qs"]),
                                   qe=max(current["qe"], b["qe"]), last=(b["qs"], b["qe"]))
                    current["blocks"] += 1
                    current["aligned_bp"] += b["te"] - b["ts"]
                    continue
            current = {"query": qname, "strand": b["strand"], "ts": b["ts"], "te": b["te"], "qs": b["qs"],
                       "qe": b["qe"], "blocks": 1, "aligned_bp": b["te"] - b["ts"], "last": (b["qs"], b["qe"])}
            out.append(current)
    return sorted(out, key=lambda s: (s["ts"], s["te"]))


def distance(point, lo, hi):
    return 0 if lo <= point <= hi else min(abs(point - lo), abs(point - hi))


def call_unit(unit, chrom_units, blocks, p):
    outer = (min(u["bp1_start"] for u in chrom_units), max(u["bp2_end"] for u in chrom_units))
    lo, hi = unit["bp1_end"], unit["bp2_start"]
    res = {"unit": unit["unit"], "chromosome": unit["chromosome"], "main_query": "", "interior_bp": 0,
           "interior_other_bp": 0, "background_bp": 0, "background_strand": "", "opposite_fraction": "",
           "call": "unresolved", "reason": "", "inverted_start": "", "inverted_end": "",
           "bp1_published": "{}-{}".format(unit["bp1_start"], unit["bp1_end"]),
           "bp2_published": "{}-{}".format(unit["bp2_start"], unit["bp2_end"]), "bp1_observed": "", "bp2_observed": "",
           "bp1_overlap": "", "bp2_overlap": "", "bp1_distance": "", "bp2_distance": ""}
    interior, background = collections.Counter(), collections.defaultdict(collections.Counter)
    for b in blocks:
        inside = pafkit.overlap(b["ts"], b["te"], lo, hi)
        if inside:
            interior[(b["qname"], b["strand"])] += inside
        outside = (b["te"] - b["ts"]) - pafkit.overlap(b["ts"], b["te"], *outer)
        if outside > 0:
            background[b["qname"]][b["strand"]] += outside
    per_query = collections.Counter()
    for (qname, _), bp in interior.items():
        per_query[qname] += bp
    if not per_query:
        res["reason"] = "nothing aligned in the interior"
        return res
    main = per_query.most_common(1)[0][0]
    bg = background.get(main, collections.Counter())
    res.update(main_query=main, interior_bp=per_query[main], interior_other_bp=sum(per_query.values()) - per_query[main],
               background_bp=sum(bg.values()))
    if per_query[main] < p.min_interior:
        res["reason"] = "under {:,} bp aligned in the interior".format(p.min_interior)
        return res
    if sum(bg.values()) < p.min_background:
        res["reason"] = "the interior's sequence has under {:,} bp outside the units (flanks on other sequences)".format(p.min_background)
        return res
    strand = max(bg, key=bg.get)
    other = "-" if strand == "+" else "+"
    fraction = interior[(main, other)] / per_query[main]
    res.update(background_strand=strand, opposite_fraction="{:.3f}".format(fraction))
    if fraction >= p.call_fraction:
        res["call"] = "inverted"
    elif fraction <= 1 - p.call_fraction:
        res["call"] = "normal"
        return res
    else:
        res["reason"] = "mixed strands inside"
        return res
    segs = [s for s in segments([b for b in blocks if b["qname"] == main], p.max_jump) if s["aligned_bp"] >= p.min_segment]
    span = [s for s in segs if s["strand"] == other and pafkit.overlap(s["ts"], s["te"], lo, hi) > 0]
    if not span:
        res["reason"] = "no inverted segment of {:,} bp or more".format(p.min_segment)
        res["bp1_overlap"] = res["bp2_overlap"] = "unresolved"
        return res
    start, end = min(s["ts"] for s in span), max(s["te"] for s in span)
    others = [s for s in segs if s not in span]
    before = max((s["te"] for s in others if s["te"] <= start), default=None)
    after = min((s["ts"] for s in others if s["ts"] >= end), default=None)
    res.update(inverted_start=start, inverted_end=end, bp1_distance=distance(start, unit["bp1_start"], unit["bp1_end"]),
               bp2_distance=distance(end, unit["bp2_start"], unit["bp2_end"]))
    if before is None:
        res["bp1_observed"], res["bp1_overlap"] = "-{}".format(start), "unresolved"
    else:
        res["bp1_observed"] = "{}-{}".format(before, start)
        res["bp1_overlap"] = "yes" if before <= unit["bp1_end"] and start >= unit["bp1_start"] else "no"
    if after is None:
        res["bp2_observed"], res["bp2_overlap"] = "{}-".format(end), "unresolved"
    else:
        res["bp2_observed"] = "{}-{}".format(end, after)
        res["bp2_overlap"] = "yes" if end <= unit["bp2_end"] and after >= unit["bp2_start"] else "no"
    return res


def verdict(expected, call, bp1, bp2):
    if not expected:
        return "no expectation"
    if expected not in ("inverted", "normal"):
        return "not compared ({})".format(expected)
    if call == "unresolved":
        return "unresolved"
    if call != expected:
        return "disagree"
    if call == "inverted" and "no" in (bp1, bp2):
        return "disagree: breakpoint"
    if call == "inverted" and "unresolved" in (bp1, bp2):
        return "agree; a breakpoint could not be placed"
    return "agree"


# ---------------------------------------------------------------- one run

def summarize_run(row, paf, ref_labels, asm_labels, units, p):
    """Everything for one run. Runs in a worker process; returns plain data."""
    records = list(pafkit.read_paf(paf))
    primary = [r for r in records if r["tp"] in pafkit.PRIMARY]
    out = {"row": row, "chromosomes": [], "assignment": [], "segments": [], "inversions": [], "recip": {}}
    forward = row["direction"] == "forward"
    ref_side, asm_side = ("t", "q") if forward else ("q", "t")
    spans = {"q": collections.defaultdict(list), "t": collections.defaultdict(list)}
    for r in primary:
        spans["q"][r["qname"]].append((r["qs"], r["qe"]))
        spans["t"][r["tname"]].append((r["ts"], r["te"]))
    out["recip"] = {"reference": {k: pafkit.total(pafkit.merge(v)) for k, v in spans[ref_side].items()},
                    "assembly_bp": sum(pafkit.total(pafkit.merge(v)) for v in spans[asm_side].values())}
    if not forward:
        return out
    by_target = collections.defaultdict(list)
    for r in primary:
        by_target[r["tname"]].append(r)
    for tname, recs in by_target.items():
        aligned = sum(r["te"] - r["ts"] for r in recs)
        with_de = [r for r in recs if r["de"] is not None]
        by_query = collections.Counter()
        for r in recs:
            by_query[r["qname"]] += r["te"] - r["ts"]
        main, main_bp = by_query.most_common(1)[0]
        tlabel, qlabel = ref_labels.get(tname, ""), asm_labels.get(main, "")
        covered = pafkit.total(pafkit.merge((r["ts"], r["te"]) for r in recs))
        out["chromosomes"].append({
            "run": row["run"], "target": tname, "label": tlabel, "length": recs[0]["tlen"], "primary_records": len(recs),
            "covered_bp": covered, "covered": "{:.4f}".format(covered / recs[0]["tlen"]),
            "blast_identity": "{:.4f}".format(sum(r["matches"] for r in recs) / sum(r["block"] for r in recs)),
            "divergence": "{:.5f}".format(sum(r["de"] * (r["te"] - r["ts"]) for r in with_de) /
                                          sum(r["te"] - r["ts"] for r in with_de)) if with_de else "",
            "minus_strand_bp": sum(r["te"] - r["ts"] for r in recs if r["strand"] == "-"), "main_query": main,
            "main_query_label": qlabel, "main_query_share": "{:.3f}".format(main_bp / aligned),
            "label_agrees": agrees(qlabel, tlabel)})
    by_query = collections.defaultdict(collections.Counter)
    qlen = {}
    for r in primary:
        by_query[r["qname"]][r["tname"]] += r["qe"] - r["qs"]
        qlen[r["qname"]] = r["qlen"]
    for qname, targets in by_query.items():
        ranked = targets.most_common(2)
        aligned = sum(targets.values())
        out["assignment"].append({
            "run": row["run"], "query": qname, "query_label": asm_labels.get(qname, ""), "length": qlen[qname],
            "aligned_bp": aligned, "chromosome": ranked[0][0], "chromosome_label": ref_labels.get(ranked[0][0], ""),
            "share": "{:.3f}".format(ranked[0][1] / aligned), "second": ranked[1][0] if len(ranked) > 1 else "",
            "second_share": "{:.3f}".format(ranked[1][1] / aligned) if len(ranked) > 1 else "",
            "label_agrees": agrees(asm_labels.get(qname, ""), ref_labels.get(ranked[0][0], ""))})
    blocks = kumar_blocks(records, p)
    by_label = {label: name for name, label in ref_labels.items()}
    for chromosome in sorted({u["chromosome"] for u in units}, key=lambda c: int(pafkit.label_number(c))):
        tname = by_label.get(chromosome)
        chrom_blocks = [b for b in blocks if b["tname"] == tname]
        for s in segments(chrom_blocks, p.max_jump):
            if s["aligned_bp"] >= p.min_segment:
                out["segments"].append(dict(s, run=row["run"], target=tname, label=chromosome))
        chrom_units = [u for u in units if u["chromosome"] == chromosome]
        for unit in chrom_units:
            res = call_unit(unit, chrom_units, chrom_blocks, p)
            res.update(run=row["run"], assembly=row["query_id"], preset=row["preset"])
            out["inversions"].append(res)
    return out


def agrees(query_label, target_label):
    a, b = pafkit.label_number(query_label), pafkit.label_number(target_label)
    if not (a.isdigit() and b.isdigit()):
        return "-"
    return "yes" if a == b else "no"


# ---------------------------------------------------------------- across runs

def apply_expectations(inversions, expected, preset):
    for res in inversions:
        exp = expected.get((res["assembly"], res["unit"]), {})
        scored = exp.get("scored") == "yes" and res["preset"] == preset
        res.update(expected=exp.get("expected", ""), basis=exp.get("basis", ""), scored="yes" if scored else "no",
                   verdict=verdict(exp.get("expected", ""), res["call"], res["bp1_overlap"], res["bp2_overlap"]))


def outliers(chromosomes, subset_runs, runs_of_preset, drop, ratio):
    by_label = collections.defaultdict(list)
    for c in chromosomes:
        if c["run"] in runs_of_preset and c["label"].startswith("LG"):
            by_label[c["label"]].append(c)
    flags = []
    for label, rows in by_label.items():
        subset = [c for c in rows if c["run"] in subset_runs]
        if not subset:
            continue
        cov = statistics.median(float(c["covered"]) for c in subset)
        divs = [float(c["divergence"]) for c in subset if c["divergence"]]
        div = statistics.median(divs) if divs else None
        for c in rows:
            why = []
            if float(c["covered"]) < cov - drop:
                why.append("coverage {:.3f} vs subset median {:.3f}".format(float(c["covered"]), cov))
            if div is not None and c["divergence"] and float(c["divergence"]) > ratio * div:
                why.append("divergence {:.4f} vs subset median {:.4f}".format(float(c["divergence"]), div))
            if why:
                flags.append({"run": c["run"], "label": label, "why": "; ".join(why)})
    return sorted(flags, key=lambda f: (f["run"], int(pafkit.label_number(f["label"]))))


def reciprocity(results, ref_labels):
    rows = []
    for r in results:
        if r["row"]["direction"] != "reverse":
            continue
        assembly = r["row"]["target_id"]
        f = next((x for x in results if x["row"]["direction"] == "forward" and x["row"]["query_id"] == assembly
                  and x["row"]["preset"] == r["row"]["preset"]), None)
        if f is None:
            continue
        for name, label in sorted(ref_labels.items(), key=lambda kv: (not kv[1].startswith("LG"), int(pafkit.label_number(kv[1])) if pafkit.label_number(kv[1]).isdigit() else 0)):
            if not label.startswith("LG"):
                continue
            a, b = f["recip"]["reference"].get(name, 0), r["recip"]["reference"].get(name, 0)
            rows.append({"assembly": assembly, "sequence": name, "label": label, "forward_bp": a, "reverse_bp": b, "difference_bp": b - a})
        a, b = f["recip"]["assembly_bp"], r["recip"]["assembly_bp"]
        rows.append({"assembly": assembly, "sequence": "(assembly, all)", "label": "", "forward_bp": a, "reverse_bp": b, "difference_bp": b - a})
    return rows


# ---------------------------------------------------------------- report

def status(scored, failed):
    return "FAIL" if failed else "PASS" if scored else "NOTHING SCORED (check --preset and --expected)"


def write_report(path, version, p, inversions, chromosomes, flags, recip, manifest_rows, preset):
    scored = [r for r in inversions if r["scored"] == "yes"]
    failed = [r for r in scored if r["verdict"].startswith("disagree")]
    unresolved = [r for r in scored if r["verdict"] == "unresolved"]
    lines = ["# Known features and measured checks ({})".format(version), "",
             "Blocks: primary, identity >= {}, length >= {:,} bp, from sequences with >= {:,} bp of such blocks. Calls: "
             ">= {:.0%} of the interior on the other strand. Segments >= {:,} bp, split at jumps > {:,} bp.".format(
                 p.min_identity, p.min_block, p.min_contig, p.call_fraction, p.min_segment, p.max_jump), "",
             "## Known features (scored on the subset): {}".format(status(scored, failed)), "",
             "{} scored units: {} agree, {} disagree, {} unresolved.".format(
                 len(scored), sum(r["verdict"].startswith("agree") for r in scored), len(failed), len(unresolved)), "",
             "| assembly | unit | expected | call | other strand | breakpoint 1 (observed; published) | breakpoint 2 | verdict |",
             "|---|---|---|---|---|---|---|---|"]
    for r in scored:
        lines.append("| {assembly} | {unit} | {expected} | {call} | {opposite_fraction} | {bp1_observed} ({bp1_overlap}); "
                     "{bp1_published} | {bp2_observed} ({bp2_overlap}); {bp2_published} | {verdict} |".format(**r))
    others = [r for r in inversions if r["scored"] != "yes" and r["preset"] == preset]
    shown = [r for r in others if not r["verdict"].startswith("agree") and r["verdict"] != "no expectation"]
    lines += ["", "## Reported, not scored", "",
              "{} other units with an expectation: {} agree; listed below are the rest.".format(
                  sum(r["verdict"] != "no expectation" for r in others), sum(r["verdict"].startswith("agree") for r in others)), "",
              "| assembly | unit | expected | call | other strand | reason | verdict |", "|---|---|---|---|---|---|---|"]
    lines += ["| {assembly} | {unit} | {expected} | {call} | {opposite_fraction} | {reason} | {verdict} |".format(**r) for r in shown]
    runs_of_preset = {m["run"] for m in manifest_rows if m["direction"] == "forward" and m["preset"] == preset}
    lines += ["", "## Coverage and divergence of the 22 linkage groups ({} forward runs)".format(preset), "",
              "| run | covered | blast identity | divergence | label disagreements |", "|---|---|---|---|---|"]
    by_run = collections.defaultdict(list)
    for c in chromosomes:
        if c["label"].startswith("LG"):
            by_run[c["run"]].append(c)
    for run in sorted(by_run):
        rows = by_run[run]
        length = sum(c["length"] for c in rows)
        aligned_w = sum(c["covered_bp"] for c in rows)
        lines.append("| {} | {:.3f} | {:.4f} | {:.4f} | {} |".format(
            run, aligned_w / length,
            sum(float(c["blast_identity"]) * c["covered_bp"] for c in rows) / aligned_w,
            sum(float(c["divergence"] or 0) * c["covered_bp"] for c in rows) / aligned_w,
            ", ".join(c["label"] for c in rows if c["label_agrees"] == "no") or "none"))
    lines += ["", "## Outliers (reported, not scored)", ""] + (
        ["- {run} {label}: {why}".format(**f) for f in flags] or ["- none"])
    lines += ["", "## Reciprocity (reference bases covered, forward against reverse)", "",
              "| assembly | forward Mb | reverse Mb | difference Mb |", "|---|---|---|---|"]
    for assembly in sorted({r["assembly"] for r in recip}):
        rows = [r for r in recip if r["assembly"] == assembly and r["label"].startswith("LG")]
        a, b = sum(r["forward_bp"] for r in rows), sum(r["reverse_bp"] for r in rows)
        lines.append("| {} | {:.1f} | {:.1f} | {:+.1f} |".format(assembly, a / 1e6, b / 1e6, (b - a) / 1e6))
    done = [m for m in manifest_rows if m["exit_status"] == "0"]
    lines += ["", "## Runs", "", "- {} runs; wall time {:.1f} h in total; peak memory up to {:.1f} GB; PAFs {:.2f} GB".format(
        len(done), sum(float(m["wall_s"]) for m in done) / 3600, max(float(m["peak_rss_gb"]) for m in done),
        sum(int(m["paf_bytes"]) for m in done) / 1e9)]
    Path(path).write_text("\n".join(lines) + "\n")
    return failed


# ---------------------------------------------------------------- main

def read_units(path):
    units = []
    for r in pafkit.read_table(path):
        units.append({"unit": r["unit"], "chromosome": r["chromosome"], "bp1_start": int(r["bp1_start"]),
                      "bp1_end": int(r["bp1_end"]), "bp2_start": int(r["bp2_start"]), "bp2_end": int(r["bp2_end"])})
    return units


def read_version():
    try:
        return (REPO / "VERSION").read_text().strip()
    except OSError:
        return "unknown version"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--version", default="PW_v1", help="the pairwise version folder (default: %(default)s)")
    parser.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH, help="the project's scratch folder")
    parser.add_argument("--preset", default="asm5", help="the preset whose forward runs are scored (default: %(default)s)")
    parser.add_argument("--units", type=Path, default=REPO / "data" / "kumar2025_suppfile4_breakpoints.tsv")
    parser.add_argument("--expected", type=Path, default=REPO / "data" / "expected_inversions.tsv")
    parser.add_argument("--min-identity", type=float, default=0.30)
    parser.add_argument("--min-block", type=int, default=10000)
    parser.add_argument("--min-contig", type=int, default=2000000)
    parser.add_argument("--min-segment", type=int, default=100000)
    parser.add_argument("--max-jump", type=int, default=1000000)
    parser.add_argument("--call-fraction", type=float, default=0.80)
    parser.add_argument("--min-background", type=int, default=1000000)
    parser.add_argument("--min-interior", type=int, default=100000)
    parser.add_argument("--coverage-drop", type=float, default=0.05)
    parser.add_argument("--divergence-ratio", type=float, default=2.0)
    parser.add_argument("--go", action="store_true", help="run the summary and write the results")
    pafkit.add_throttle_arguments(parser)
    return parser.parse_args(argv)


def main(argv=None, deps=None, now=None):
    args = parse_args(argv)
    p = Params(args)
    deps = deps or pafkit.load_deps()
    base = Path(args.scratch).expanduser() / "pairwise" / args.version
    manifest_rows = pafkit.read_manifest(base)
    rows = [r for r in manifest_rows if r["exit_status"] == "0"]
    units = read_units(args.units)
    expected = {(r["assembly_id"], r["unit"]): r for r in pafkit.read_table(args.expected)}
    stamp = (now or datetime.datetime.now(datetime.timezone.utc)).strftime("%Y%m%d-%H%M%S")
    outdir = base / "summary" / stamp
    jobs = max(1, min(getattr(args, "max_jobs", None) or 1, len(rows) or 1))
    print("summarize_pairwise.py ({}): {} runs in {}, {} at a time; {} units, {} expectations ({} scored)".format(
        read_version(), len(rows), base, jobs, len(units), len(expected), sum(e["scored"] == "yes" for e in expected.values())))
    if not args.go:
        for r in rows:
            print("  {:<32} {} {}".format(r["run"], r["direction"], r["paf"]))
        print("Dry run: nothing has been read or written. Add --go to write\n  {}".format(outdir))
        return 0
    if outdir.exists():
        print("STOP: {} already exists, and this script never overwrites.".format(outdir), file=sys.stderr)
        return 2
    ref_labels = {n: e.get("chr_name", "") for n, e in pafkit.provenance_sequences(deps.locate(deps.reference_id)).items()}
    tasks = []
    for r in rows:
        assembly = pafkit.run_assembly(r)
        asm_labels = {n: e.get("chr_name", "") for n, e in pafkit.provenance_sequences(deps.locate(assembly)).items()}
        tasks.append((r, str(base / r["paf"]), ref_labels, asm_labels, units, p))
    if jobs == 1:
        results = [summarize_run(*t) for t in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
            results = list(pool.map(summarize_run, *zip(*tasks)))
    chromosomes = [c for r in results for c in r["chromosomes"]]
    inversions = [i for r in results for i in r["inversions"]]
    apply_expectations(inversions, expected, args.preset)
    subset = {a for (a, _), e in expected.items() if e["scored"] == "yes"}
    runs_of_preset = {r["run"] for r in rows if r["direction"] == "forward" and r["preset"] == args.preset}
    subset_runs = {r["run"] for r in rows if r["run"] in runs_of_preset and r["query_id"] in subset}
    flags = outliers(chromosomes, subset_runs, runs_of_preset, args.coverage_drop, args.divergence_ratio)
    recip = reciprocity(results, ref_labels)
    outdir.mkdir(parents=True)
    pafkit.write_tsv(outdir / "chromosomes.tsv", CHROM_COLUMNS, chromosomes)
    pafkit.write_tsv(outdir / "assignment.tsv", ASSIGN_COLUMNS, [a for r in results for a in r["assignment"]])
    pafkit.write_tsv(outdir / "segments.tsv", SEG_COLUMNS, [s for r in results for s in r["segments"]])
    pafkit.write_tsv(outdir / "inversions.tsv", INV_COLUMNS, inversions)
    pafkit.write_tsv(outdir / "outliers.tsv", ["run", "label", "why"], flags)
    pafkit.write_tsv(outdir / "reciprocity.tsv", ["assembly", "sequence", "label", "forward_bp", "reverse_bp",
                                                 "difference_bp"], recip)
    failed = write_report(outdir / "report.md", read_version(), p, inversions, chromosomes, flags, recip, manifest_rows, args.preset)
    scored = [i for i in inversions if i["scored"] == "yes"]
    print("Results written to {}".format(outdir))
    print("Known features: {} ({} scored units: {} agree, {} disagree, {} unresolved)".format(
        status(scored, failed), len(scored), sum(i["verdict"].startswith("agree") for i in scored),
        len(failed), sum(i["verdict"] == "unresolved" for i in scored)))
    for f in failed:
        print("  - {} {}: expected {}, called {} ({})".format(f["assembly"], f["unit"], f["expected"], f["call"], f["verdict"]))
    return 1 if failed or not scored else 0


if __name__ == "__main__":
    sys.exit(main())
