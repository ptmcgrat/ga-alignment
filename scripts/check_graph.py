#!/usr/bin/env python3
"""The graph checks of spec section 5 (agreed 2026-11-25), for runs in <scratch>/graph/<version>/.

For each run (a folder written by build_graph.py whose last run ended with exit status 0), from its full graph
(<name>.full.gfa.gz) and clip graph (<name>.gfa.gz), its inputs (seqfile.txt) and Cactus's stats tables:
  exact           every walk (W line) in both graphs spells its input contig over its start-end range, compared in
                  upper case, with the first differing position; the reference (the H line's RS tag) has one walk
                  per reference contig, complete, forward only, and visiting no node twice
  accounting      per input contig (its walks merged into contiguous pieces, as Cactus counts fragments), its
                  bases in each graph plus Cactus's absent intervals for that graph (the
                  clipped*.beds archives) equal its length; our counts per contig and per reference contig and genome
                  equal Cactus's clipped-by-input-contig.tsv and clipped-by-reference-contig.tsv exactly
  known features  each walk projected onto the reference walk: steps on reference nodes, merged into blocks where
                  the strand holds and both coordinates continue within --gap bp; then the PAF check's rules
                  (summarize_pairwise.py): Kumar et al.'s filters, flank windows, the 80% rule, segments, and YH's
                  breakpoints against Supplementary file 4. Scored on the subset; a disagreement fails the check
  graph vs PAF    the same call per unit as the newest PW_v1 summary; for units inverted in both, the observed
                  breakpoint intervals overlap. Reported
  reference-only  reference intervals of at least --refonly-min bp that no other walk visits in the full graph,
                  against those no subset PAF covers: the share of each inside the other (flagged below 90%), and
                  Cactus's refgaps.bed alongside. Reported
  assignment      whole-genome runs: each contig's graph according to Cactus against our PAF assignment (the newest
                  PW_v1 summary), contigs outside --refContigs named by --otherContig as Cactus names them. Reported
  measured        nodes, edges and length of each graph; each genome's bases in each graph; what clipping removed by
                  reason; Cactus's inversion and graph stats; wall time, CPU and memory from run.json; disk. Reported
A dry run unless --go; it writes one new folder, <scratch>/graph/<version>/checks/<UTC timestamp>/, and never
overwrites. Exit status: 0 if the exact, accounting and known-feature checks pass, 1 if not, 2 if it refused to start.
"""

import argparse
import collections
import concurrent.futures
import datetime
import gzip
import io
import json
import re
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pafkit  # noqa: E402
import summarize_pairwise as sp  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCRATCH = Path.home() / "scratch" / "ga-alignment"
STEP = re.compile(r"([<>])([^<>]+)")
ACCOUNT_COLUMNS = ["contigs", "input_bp", "full_bp", "full_frags", "clip_bp", "clip_frags"]


# ---------------------------------------------------------------- reading

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_gfa(path):
    """{'reference': RS sample, 'nodes': {id: upper-case bytes}, 'walks': [...], 'edges': n}."""
    graph = {"reference": None, "nodes": {}, "walks": [], "edges": 0}
    with open_text(path) as handle:
        for line in handle:
            tag = line[:1]
            if tag == "S":
                fields = line.rstrip("\n").split("\t")
                graph["nodes"][fields[1]] = fields[2].upper().encode()
            elif tag == "W":
                f = line.rstrip("\n").split("\t")
                graph["walks"].append({"sample": f[1], "hap": f[2], "contig": f[3], "start": int(f[4]),
                                       "end": int(f[5]), "walk": f[6]})
            elif tag == "L":
                graph["edges"] += 1
            elif tag == "H":
                for field in line.rstrip("\n").split("\t")[1:]:
                    if field.startswith("RS:Z:"):
                        graph["reference"] = field[5:]
    return graph


def steps(walk):
    return [(m.group(1), m.group(2)) for m in STEP.finditer(walk)]


def read_cactus_table(path):
    """Rows of a Cactus stats table: its header is the last comment line that names the columns."""
    header, rows = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                if "\t" in line:
                    header = line[1:].rstrip("\n").split("\t")
            elif line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def read_beds(path):
    """{(contig, reason): bp} from a clipped*.beds.tar.gz archive."""
    totals = collections.Counter()
    if not Path(path).exists():
        return totals
    with tarfile.open(path, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            for line in io.TextIOWrapper(tar.extractfile(member)):
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4 and f[1].isdigit():
                    totals[(f[0], f[3])] += int(f[2]) - int(f[1])
    return totals


def read_seqfile(folder):
    out = []
    for line in (Path(folder) / "seqfile.txt").read_text().splitlines():
        if line.strip():
            sample, path = line.split("\t")
            out.append((sample, Path(path)))
    return out


def option_values(command, name):
    """The values given to an option in a recorded command (a list), up to the next option."""
    if name not in command:
        return []
    out, i = [], command.index(name) + 1
    while i < len(command) and not command[i].startswith("--"):
        out.append(command[i])
        i += 1
    return out


def graph_namer(command):
    """How Cactus names the graph a reference contig's sequence lands in: contigs left out of --refContigs are lumped
    under --otherContig in a whole-genome run."""
    ref_contigs, other = set(option_values(command, "--refContigs")), (option_values(command, "--otherContig") or [None])[0]
    return lambda contig: other if other and contig not in ref_contigs else contig


def classify_assignment(ours, cactus):
    if not ours:
        return "no PAF alignment"
    if cactus in ("", "_NONE_"):
        return "Cactus: none"
    return "same" if cactus == ours else "different"


# ---------------------------------------------------------------- exact

def spell(walk_steps, nodes):
    return b"".join(nodes[n] if o == ">" else pafkit.revcomp(nodes[n]) for o, n in walk_steps)


def first_difference(a, b):
    for i in range(0, max(len(a), len(b)), 1 << 20):
        if a[i:i + (1 << 20)] != b[i:i + (1 << 20)]:
            x, y = a[i:i + (1 << 20)], b[i:i + (1 << 20)]
            return i + next((k for k in range(min(len(x), len(y))) if x[k] != y[k]), min(len(x), len(y)))
    return None


def exact_checks(graph, phase, genomes, reference):
    """Problems for one graph: spelling of every walk against its input, and the reference walk's form."""
    problems, spelled = [], 0
    for w in graph["walks"]:
        seq = genomes.get(w["sample"], {}).get(w["contig"])
        if seq is None:
            problems.append("{}: walk {}#{} names a contig not in its input".format(phase, w["sample"], w["contig"]))
            continue
        text = spell(steps(w["walk"]), graph["nodes"])
        want = seq[w["start"]:w["end"]]
        if text != want:
            at = first_difference(text, want)
            problems.append("{}: {}#{} [{}-{}] differs from its input at offset {} ({} vs {} bp)".format(
                phase, w["sample"], w["contig"], w["start"], w["end"], at, len(text), len(want)))
        spelled += len(text)
    if graph["reference"] != reference:
        problems.append("{}: the graph's reference (RS) is {!r}, not {}".format(phase, graph["reference"], reference))
    ref_walks = collections.defaultdict(list)
    for w in graph["walks"]:
        if w["sample"] == reference:
            ref_walks[w["contig"]].append(w)
    for contig, seq in genomes.get(reference, {}).items():
        walks = ref_walks.get(contig, [])
        if len(walks) != 1:
            problems.append("{}: reference contig {} has {} walks, not 1".format(phase, contig, len(walks)))
            continue
        w = walks[0]
        walk_steps = steps(w["walk"])
        if (w["start"], w["end"]) != (0, len(seq)):
            problems.append("{}: reference {} walk covers {}-{} of {}".format(phase, contig, w["start"], w["end"], len(seq)))
        if any(o == "<" for o, _ in walk_steps):
            problems.append("{}: reference {} walk is not forward only".format(phase, contig))
        if len({n for _, n in walk_steps}) != len(walk_steps):
            problems.append("{}: reference {} walk visits a node twice".format(phase, contig))
    return problems, spelled


# ---------------------------------------------------------------- reference coordinates and projection

def reference_positions(graph, reference):
    """{node: (reference contig, start, length)} along the reference walks."""
    pos = {}
    for w in graph["walks"]:
        if w["sample"] != reference:
            continue
        offset = w["start"]
        for _, n in steps(w["walk"]):
            length = len(graph["nodes"][n])
            pos[n] = (w["contig"], offset, length)
            offset += length
    return pos


def project(walk, nodes, pos, gap):
    """Blocks of one walk on the reference, shaped like PAF records for summarize_pairwise's calls."""
    blocks, current, q = [], None, walk["start"]
    for o, n in steps(walk["walk"]):
        length = len(nodes[n])
        hit = pos.get(n)
        if hit:
            contig, ts, _ = hit
            strand = "+" if o == ">" else "-"
            te = ts + length
            if current and current["tname"] == contig and current["strand"] == strand:
                ref_gap = ts - current["te"] if strand == "+" else current["ts"] - te
                if 0 <= ref_gap <= gap and 0 <= q - current["qe"] <= gap:
                    current.update(qe=q + length, ts=min(current["ts"], ts), te=max(current["te"], te))
                    current["matches"] += length
                    q += length
                    continue
            current = {"qname": walk["contig"], "qlen": 0, "qs": q, "qe": q + length, "strand": strand,
                       "tname": contig, "ts": ts, "te": te, "matches": length, "tp": "P"}
            blocks.append(current)
        q += length
    for b in blocks:
        b["block"] = b["te"] - b["ts"]
    return blocks


# ---------------------------------------------------------------- one run

def check_run(folder, reference, ref_labels, units, expected, paf_calls, paf_runs, args, paf_assign=None):
    folder = Path(folder)
    name = folder.name
    out = folder / "out"
    result = {"run": name, "exact": [], "accounting": [], "inversions": [], "refonly": [], "measured": [], "assignment": []}
    record = json.loads((folder / "run.json").read_text())[-1]
    named = graph_namer(record.get("command", []))
    genomes = {sample: pafkit.read_fasta(path) for sample, path in read_seqfile(folder)}
    graphs, sizes_by_phase = {}, {}
    for phase, filename in (("clip", "{}.gfa.gz"), ("full", "{}.full.gfa.gz")):   # one graph in memory at a time
        graph = read_gfa(out / filename.format(name))
        problems, spelled = exact_checks(graph, phase, genomes, reference)
        result["exact"].append({"run": name, "graph": phase, "walks": len(graph["walks"]), "bp_spelled": spelled,
                                "problems": len(problems), "status": "pass" if not problems else "FAIL",
                                "details": " | ".join(problems[:5])})
        sizes_by_phase[phase] = (len(graph["nodes"]), graph["edges"], sum(len(x) for x in graph["nodes"].values()))
        if phase == "clip":
            graphs["clip"] = {"walks": [dict(w, walk="") for w in graph["walks"]]}
            del graph
        else:
            graphs["full"] = graph

    # accounting
    pos_full = reference_positions(graphs["full"], reference)
    sizes = collections.Counter()
    for r in read_cactus_table(out / "{}.input-contig-sizes.tsv.gz".format(name)):
        sizes[(r["event"], r["base_contig"])] += int(r["length"])
    mine, spans, chroms = {}, collections.defaultdict(list), collections.defaultdict(collections.Counter)
    for phase, graph in graphs.items():
        for w in graph["walks"]:
            spans[(phase, w["sample"], w["contig"])].append((w["start"], w["end"]))
            if phase == "full":
                for _, n in steps(w["walk"]):
                    if n in pos_full:
                        chroms[(w["sample"], w["contig"])][pos_full[n][0]] += pos_full[n][2]
    for (phase, sample, contig), intervals in spans.items():
        row = mine.setdefault((sample, contig), {"full_bp": 0, "full_frags": 0, "clip_bp": 0, "clip_frags": 0, "ref_chrom": "_NONE_"})
        pieces = pafkit.merge(intervals)               # contiguous pieces: walks that abut or overlap are one
        row[phase + "_bp"], row[phase + "_frags"] = pafkit.total(pieces), len(pieces)
    for key, row in mine.items():
        if row["full_frags"]:
            row["ref_chrom"] = named(chroms[key].most_common(1)[0][0]) if chroms[key] else "_UNKNOWN_"
    absent = {phase: read_beds(out / "{}.stats".format(name) / ("clipped.full.beds.tar.gz" if phase == "full" else "clipped.beds.tar.gz"))
              for phase in graphs}
    problems = []
    for (genome, contig), length in sorted(sizes.items()):
        row = mine.get((genome, contig), {"full_bp": 0, "full_frags": 0, "clip_bp": 0, "clip_frags": 0, "ref_chrom": "_NONE_"})
        full_missing = sum(bp for (c, _), bp in absent["full"].items() if c == contig)
        for phase in graphs:
            missing = sum(bp for (c, _), bp in absent[phase].items() if c == contig)
            # the clip archive may list everything absent from the clip graph, or only what clipping removed
            closes = row[phase + "_bp"] + missing == length or (
                phase == "clip" and row["clip_bp"] + missing + full_missing == length)
            if not closes:
                problems.append("{} {} {}: {} in graph + {} absent != {}".format(phase, genome, contig, row[phase + "_bp"], missing, length))
    cactus_contigs = {(r["genome"], r["contig"]): r for r in read_cactus_table(out / "{}.stats".format(name) / "clipped-by-input-contig.tsv.gz")}
    for (genome, contig), length in sorted(sizes.items()):
        row = mine.get((genome, contig), {"full_bp": 0, "full_frags": 0, "clip_bp": 0, "clip_frags": 0, "ref_chrom": "_NONE_"})
        c = cactus_contigs.get((genome, contig))
        ours = {"ref_chrom": row["ref_chrom"], "input_bp": length, **{k: row[k] for k in ("full_bp", "full_frags", "clip_bp", "clip_frags")}}
        if c is None or any(str(ours[k]) != c.get(k) for k in ours):
            problems.append("contig {} {}: ours {} vs Cactus {}".format(genome, contig, ours, c))
    per_chrom = collections.defaultdict(lambda: collections.Counter())
    for (genome, contig), length in sizes.items():
        row = mine.get((genome, contig), {"full_bp": 0, "full_frags": 0, "clip_bp": 0, "clip_frags": 0, "ref_chrom": "_NONE_"})
        t = per_chrom[(row["ref_chrom"], genome)]
        t["contigs"] += 1
        t["input_bp"] += length
        for k in ("full_bp", "full_frags", "clip_bp", "clip_frags"):
            t[k] += row[k]
    cactus_chroms = {(r["ref_chrom"], r["genome"]): r for r in read_cactus_table(out / "{}.stats".format(name) / "clipped-by-reference-contig.tsv")}
    for key, t in sorted(per_chrom.items()):
        c = cactus_chroms.get(key)
        if c is None or any(str(t[k]) != c.get(k) for k in ACCOUNT_COLUMNS):
            problems.append("reference contig {} {}: ours {} vs Cactus {}".format(key[0], key[1], dict(t), c))
    if set(cactus_chroms) != set(per_chrom):
        problems.append("rows only in Cactus's table: {}".format(sorted(set(cactus_chroms) - set(per_chrom))))
    result["accounting"].append({"run": name, "contigs": len(sizes), "problems": len(problems),
                                 "status": "pass" if not problems else "FAIL", "details": " | ".join(problems[:5])})
    if record.get("mode") == "whole genome" and paf_assign is not None:
        reasons = {c: reason for (c, reason) in absent["full"]}
        for (genome, contig), length in sorted(sizes.items()):
            if genome == reference:
                continue
            ours = paf_assign.get(genome, {}).get(contig)
            cactus = cactus_contigs.get((genome, contig), {}).get("ref_chrom", "")
            our_chrom = named(ours[0]) if ours else ""
            result["assignment"].append({"run": name, "genome": genome, "contig": contig, "length": length,
                                         "cactus": cactus, "cactus_reason": reasons.get(contig, "") if cactus in ("", "_NONE_") else "",
                                         "paf": our_chrom, "paf_label": ref_labels.get(ours[0], "") if ours else "",
                                         "paf_share": ours[1] if ours else "", "status": classify_assignment(our_chrom, cactus)})

    # known features, graph against PAF
    blocks = []
    for w in graphs["full"]["walks"]:
        if w["sample"] != reference:
            for b in project(w, graphs["full"]["nodes"], pos_full, args.gap):
                b["sample"] = w["sample"]
                blocks.append(b)
    by_label = {label: n for n, label in ref_labels.items()}
    for sample in sorted({w["sample"] for w in graphs["full"]["walks"]} - {reference}):
        sample_blocks = sp.kumar_blocks([b for b in blocks if b["sample"] == sample], args)
        for chromosome in sorted({u["chromosome"] for u in units}):
            tname = by_label.get(chromosome)
            if tname not in {c for c, _, _ in pos_full.values()}:
                continue
            chrom_units = [u for u in units if u["chromosome"] == chromosome]
            for unit in chrom_units:
                res = sp.call_unit(unit, chrom_units, [b for b in sample_blocks if b["tname"] == tname], args)
                exp = expected.get((sample, unit["unit"]), {})
                on_source = sample.startswith(args.breakpoint_species + "_")
                res.update(run=name, assembly=sample, expected=exp.get("expected", ""),
                           scored="yes" if exp.get("scored") == "yes" else "no", breakpoints_scored="yes" if on_source else "no",
                           verdict=sp.verdict(exp.get("expected", ""), res["call"], res["bp1_overlap"], res["bp2_overlap"], on_source))
                paf = paf_calls.get((sample, unit["unit"]))
                res["paf_call"] = paf["call"] if paf else ""
                res["graph_vs_paf"] = graph_vs_paf(res, paf)
                result["inversions"].append(res)

    # reference-only stretches
    visited = set()
    for w in graphs["full"]["walks"]:
        if w["sample"] != reference:
            visited.update(n for _, n in steps(w["walk"]))
    graph_only = collections.defaultdict(list)
    for n, (contig, start, length) in pos_full.items():
        if n not in visited:
            graph_only[contig].append((start, start + length))
    refgaps = collections.defaultdict(list)
    gap_file = out / "{}.stats".format(name) / "refgaps.bed.gz"
    if gap_file.exists():
        with open_text(gap_file) as handle:
            for line in handle:
                f = line.split("\t")
                if len(f) >= 3 and f[1].isdigit():
                    refgaps[f[0]].append((int(f[1]), int(f[2])))
    lengths = {w["contig"]: w["end"] for w in graphs["full"]["walks"] if w["sample"] == reference}
    for contig, length in sorted(lengths.items()):
        g = [iv for iv in pafkit.merge(graph_only[contig]) if iv[1] - iv[0] >= args.refonly_min]
        covered = pafkit.merge(iv for run in paf_runs for iv in run.get(contig, []))
        p = [iv for iv in pafkit.subtract([(0, length)], covered) if iv[1] - iv[0] >= args.refonly_min]
        cg = pafkit.merge(refgaps[contig])
        share = lambda a, b: (pafkit.total(pafkit.subtract(pafkit.merge(a), pafkit.subtract(pafkit.merge(a), pafkit.merge(b)))) / pafkit.total(a)) if pafkit.total(a) else 1.0
        row = {"run": name, "contig": contig, "label": ref_labels.get(contig, ""), "graph_only_bp": pafkit.total(g),
               "paf_only_bp": pafkit.total(p), "cactus_refgap_bp": pafkit.total(cg),
               "graph_in_paf": "{:.3f}".format(share(g, p)), "paf_in_graph": "{:.3f}".format(share(p, g)),
               "graph_in_cactus": "{:.3f}".format(share(g, cg))}
        row["flag"] = "below 90%" if min(share(g, p), share(p, g)) < 0.9 else ""
        result["refonly"].append(row)

    # measured
    stats_dir = out / "{}.stats".format(name)
    disk = sum(f.stat().st_size for f in folder.rglob("*") if f.is_file())
    for phase in ("full", "clip"):
        nodes, edges, length_bp = sizes_by_phase[phase]
        result["measured"].append({"run": name, "graph": phase, "nodes": nodes, "edges": edges, "length_bp": length_bp,
                                   "wall_min": round(record.get("wall_s", 0) / 60, 1), "cpu_h": round(record.get("cpu_s", 0) / 3600, 2),
                                   "largest_process_gb": record.get("largest_process_peak_gb", ""), "disk_gb": round(disk / 1e9, 2)})
    result["clip_by_genome"] = [dict(r, run=name) for r in read_cactus_table(stats_dir / "clipped-by-genome.tsv")] \
        if (stats_dir / "clipped-by-genome.tsv").exists() else []
    result["cactus_notes"] = [(stats_dir / f).read_text().strip() for f in ("inversion-stats.tsv", "graph-stats.tsv")
                              if (stats_dir / f).exists()]
    return result


def interval(text):
    a, _, b = (text or "").partition("-")
    return (int(a), int(b)) if a.isdigit() and b.isdigit() else None


def graph_vs_paf(res, paf):
    if not paf:
        return "no PAF call"
    if res["call"] != paf["call"]:
        return "different call (graph {}, PAF {})".format(res["call"], paf["call"])
    if res["call"] != "inverted":
        return "same call"
    out = []
    for k in ("bp1_observed", "bp2_observed"):
        g, q = interval(res[k]), interval(paf[k])
        out.append("?" if not (g and q) else "yes" if g[0] <= q[1] and q[0] <= g[1] else "no")
    return "same call; breakpoints overlap: {}".format(", ".join(out))


# ---------------------------------------------------------------- main

def write_report(path, results, version):
    lines = ["# Graph checks ({})".format(version), ""]
    exact = [e for r in results for e in r["exact"]]
    acc = [a for r in results for a in r["accounting"]]
    inv = [i for r in results for i in r["inversions"]]
    scored = [i for i in inv if i["scored"] == "yes"]
    failed = [i for i in scored if i["verdict"].startswith("disagree")]
    lines += ["- exact: {} of {} graphs pass ({:,} bp spelled)".format(sum(e["status"] == "pass" for e in exact), len(exact),
                                                                        sum(e["bp_spelled"] for e in exact)),
              "- accounting: {} of {} runs match Cactus's tables exactly".format(sum(a["status"] == "pass" for a in acc), len(acc)),
              "- known features: {} ({} scored units: {} agree, {} disagree, {} unresolved)".format(
                  "FAIL" if failed else "PASS" if scored else "NOTHING SCORED", len(scored),
                  sum(i["verdict"].startswith("agree") for i in scored), len(failed), sum(i["verdict"] == "unresolved" for i in scored)), ""]
    for e in exact + acc:
        if e["status"] != "pass":
            lines.append("- **{}**: {}".format(e["run"], e["details"]))
    lines += ["", "## Inversions in the graph, against the PAFs", "",
              "| run | assembly | unit | expected | graph call | other strand | breakpoints (graph) | PAF call | graph vs PAF | verdict |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for i in inv:
        lines.append("| {run} | {assembly} | {unit} | {expected} | {call} | {opposite_fraction} | {bp1_observed}; {bp2_observed} | "
                     "{paf_call} | {graph_vs_paf} | {verdict} |".format(**i))
    lines += ["", "## Reference-only stretches", "", "| run | contig | graph only | PAF only | Cactus refgaps | graph in PAF | PAF in graph | flag |",
              "|---|---|---|---|---|---|---|---|"]
    for r in (x for res in results for x in res["refonly"]):
        lines.append("| {run} | {label} | {graph_only_bp:,} | {paf_only_bp:,} | {cactus_refgap_bp:,} | {graph_in_paf} | {paf_in_graph} | {flag} |".format(**r))
    rows = [x for res in results for x in res.get("assignment", [])]
    if rows:
        lines += ["", "## Chromosome assignment: Cactus against the PAFs (reported)", "",
                  "| run | genome | same | different | Cactus: none | no PAF alignment |", "|---|---|---|---|---|---|"]
        for key in sorted({(r["run"], r["genome"]) for r in rows}):
            group = [r for r in rows if (r["run"], r["genome"]) == key]
            cell = lambda status: "{} ({:.2f} Mb)".format(sum(r["status"] == status for r in group),
                                                         sum(r["length"] for r in group if r["status"] == status) / 1e6)
            lines.append("| {} | {} | {} | {} | {} | {} |".format(key[0], key[1], cell("same"), cell("different"),
                                                                    cell("Cactus: none"), cell("no PAF alignment")))
        worst = sorted((r for r in rows if r["status"] in ("different", "Cactus: none")), key=lambda r: -r["length"])[:12]
        lines += ["", "Largest contigs not placed as the PAFs place them:", ""] + [
            "- {genome} {contig} ({length:,} bp): Cactus {cactus} {cactus_reason}; PAF {paf} {paf_label} (share {paf_share})".format(**r)
            for r in worst]
    lines += ["", "## Measured", "", "| run | graph | nodes | edges | length | wall min | CPU h | largest GB | disk GB |", "|---|---|---|---|---|---|---|---|---|"]
    for m in (x for res in results for x in res["measured"]):
        lines.append("| {run} | {graph} | {nodes:,} | {edges:,} | {length_bp:,} | {wall_min} | {cpu_h} | {largest_process_gb} | {disk_gb} |".format(**m))
    lines += ["", "Clipped, by genome and reason (Cactus):", ""]
    for r in (x for res in results for x in res["clip_by_genome"]):
        if r.get("reason") not in ("TOTAL",) and r.get("bp", "0") != "0":
            lines.append("- {run} {genome}: {reason} {bp} bp ({pct_of_input}%)".format(**r))
    lines += ["", "Cactus's notes:", ""] + ["    " + n.replace("\n", "\n    ") for res in results for n in res["cactus_notes"]]
    Path(path).write_text("\n".join(lines) + "\n")
    return bool(failed) or any(e["status"] != "pass" for e in exact + acc)


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--version", default="GR_v1")
    parser.add_argument("--runs", nargs="+", help="run folders to check (default: every finished run)")
    parser.add_argument("--pairwise-version", default="PW_v1")
    parser.add_argument("--preset", default="asm5")
    parser.add_argument("--scratch", type=Path, default=DEFAULT_SCRATCH)
    parser.add_argument("--units", type=Path, default=REPO / "published" / "kumar2025_suppfile4_breakpoints.tsv")
    parser.add_argument("--expected", type=Path, default=REPO / "published" / "expected_inversions.tsv")
    parser.add_argument("--gap", type=int, default=10000, help="largest gap within a projected block (default: %(default)s)")
    parser.add_argument("--refonly-min", type=int, default=10000)
    for key, default in (("min-identity", 0.30), ("min-block", 10000), ("min-contig", 2000000), ("min-segment", 100000),
                         ("max-jump", 1000000), ("call-fraction", 0.80), ("min-background", 500000), ("min-interior", 100000),
                         ("flank-window", 2000000)):
        parser.add_argument("--" + key, type=type(default), default=default)
    parser.add_argument("--breakpoint-species", default="YH")
    parser.add_argument("--go", action="store_true")
    pafkit.add_throttle_arguments(parser)
    return parser.parse_args(argv)


def main(argv=None, deps=None, now=None):
    args = parse_args(argv)
    deps = deps or pafkit.load_deps()
    base = Path(args.scratch).expanduser() / "graph" / args.version
    finished = []
    for folder in sorted(p for p in base.iterdir() if p.is_dir() and (p / "run.json").exists()) if base.exists() else []:
        if json.loads((folder / "run.json").read_text())[-1].get("exit_status") == 0 and (not args.runs or folder.name in args.runs):
            finished.append(folder)
    stamp = (now or datetime.datetime.now(datetime.timezone.utc)).strftime("%Y%m%d-%H%M%S")
    outdir = base / "checks" / stamp
    print("check_graph.py: {} finished run(s) in {}: {}".format(len(finished), base, ", ".join(f.name for f in finished)))
    if not args.go:
        print("Dry run: nothing has been read or written. Add --go to check them into\n  {}".format(outdir))
        return 0
    if outdir.exists():
        print("STOP: {} already exists, and this script never overwrites.".format(outdir), file=sys.stderr)
        return 2
    reference = deps.reference_id
    ref_labels = {n: e.get("chr_name", "") for n, e in pafkit.provenance_sequences(deps.locate(reference)).items()}
    units = sp.read_units(args.units)
    expected = {(r["assembly_id"], r["unit"]): r for r in pafkit.read_table(args.expected)}
    summaries = sorted((Path(args.scratch).expanduser() / "pairwise" / args.pairwise_version / "summary").glob("*/inversions.tsv"))
    paf_calls = {(r["assembly"], r["unit"]): r for r in pafkit.read_table(summaries[-1]) if r["preset"] == args.preset} if summaries else {}
    pw = Path(args.scratch).expanduser() / "pairwise" / args.pairwise_version
    manifest = pafkit.read_manifest(pw) if (pw / "manifest.tsv").exists() else []
    paf_assign = collections.defaultdict(dict)
    assignments = sorted((Path(args.scratch).expanduser() / "pairwise" / args.pairwise_version / "summary").glob("*/assignment.tsv"))
    for r in (pafkit.read_table(assignments[-1]) if assignments else []):
        sample, _, preset = r["run"].rpartition(".")
        if preset == args.preset:
            paf_assign[sample][r["query"]] = (r["chromosome"], r["share"], r["aligned_bp"])
    tasks = []
    for folder in finished:
        samples = [s for s, _ in read_seqfile(folder) if s != reference]
        paf_runs = []
        for m in manifest:
            if m["direction"] == "forward" and m["preset"] == args.preset and m["query_id"] in samples and m["exit_status"] == "0":
                cov = collections.defaultdict(list)
                for rec in pafkit.read_paf(pw / m["paf"]):
                    if rec["tp"] in pafkit.PRIMARY:
                        cov[rec["tname"]].append((rec["ts"], rec["te"]))
                paf_runs.append(dict(cov))
        tasks.append((str(folder), reference, ref_labels, units, expected, paf_calls, paf_runs, args, dict(paf_assign)))
    jobs = max(1, min(getattr(args, "max_jobs", None) or 1, len(tasks) or 1))
    if jobs == 1:
        results = [check_run(*t) for t in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
            results = list(pool.map(check_run, *zip(*tasks)))
    outdir.mkdir(parents=True)
    for key, cols in (("exact", ["run", "graph", "walks", "bp_spelled", "problems", "status", "details"]),
                      ("accounting", ["run", "contigs", "problems", "status", "details"]),
                      ("inversions", ["run", "assembly", "unit", "chromosome", "main_query", "call", "opposite_fraction", "reason",
                                      "inverted_start", "inverted_end", "bp1_observed", "bp1_published", "bp1_overlap",
                                      "bp2_observed", "bp2_published", "bp2_overlap", "expected", "scored", "breakpoints_scored",
                                      "verdict", "paf_call", "graph_vs_paf"]),
                      ("refonly", ["run", "contig", "label", "graph_only_bp", "paf_only_bp", "cactus_refgap_bp", "graph_in_paf",
                                   "paf_in_graph", "graph_in_cactus", "flag"]),
                      ("measured", ["run", "graph", "nodes", "edges", "length_bp", "wall_min", "cpu_h", "largest_process_gb", "disk_gb"]),
                      ("assignment", ["run", "genome", "contig", "length", "cactus", "cactus_reason", "paf", "paf_label", "paf_share", "status"])):
        pafkit.write_tsv(outdir / (key + ".tsv"), cols, [x for r in results for x in r[key]])
    failed = write_report(outdir / "report.md", results, args.version)
    print("Results written to {}".format(outdir))
    print(open(outdir / "report.md").read().split("\n## ")[0])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
