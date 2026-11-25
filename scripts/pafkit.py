"""Shared helpers for project 2a's checks: PAF records, cs and CIGAR operations, intervals, genomes and inputs.

Used by check_pairwise.py and summarize_pairwise.py; it does nothing on its own.
"""

import csv
import gzip
import json
import re
from pathlib import Path

SUPPLEMENTARY = {"YH_LG10Y": "h2tg000014l"}
PRIMARY = ("P", "I")      # minimap2's tp tag: P primary, I inversion (primary level); S and i are secondary
SECONDARY = ("S", "i")
CS_OP = re.compile(r":(\d+)|\*([a-zA-Z])([a-zA-Z])|\+([a-zA-Z]+)|-([a-zA-Z]+)|~([a-zA-Z]{2})(\d+)([a-zA-Z]{2})")
CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
_COMPLEMENT = bytes.maketrans(b"ACGTRYKMSWBDHVN", b"TGCAYRMKSWVHDBN")


# ---------------------------------------------------------------- PAF records

def parse_record(line, keep_strings=False):
    """One PAF line as a dict: the 12 columns, tp, de, and (if keep_strings) the cs and cg strings."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    rec = {"qname": fields[0], "qlen": int(fields[1]), "qs": int(fields[2]), "qe": int(fields[3]),
           "strand": fields[4], "tname": fields[5], "tlen": int(fields[6]), "ts": int(fields[7]),
           "te": int(fields[8]), "matches": int(fields[9]), "block": int(fields[10]), "mapq": int(fields[11]),
           "tp": "", "de": None, "cs": None, "cg": None}
    for tag in fields[12:]:
        key, _, value = tag.partition(":")
        value = value.partition(":")[2]
        if key == "tp":
            rec["tp"] = value
        elif key == "de":
            rec["de"] = float(value)
        elif keep_strings and key in ("cs", "cg"):
            rec[key] = value
    return rec


def read_paf(path, keep_strings=False):
    with open(path) as handle:
        for line in handle:
            rec = parse_record(line, keep_strings)
            if rec is not None:
                yield rec


def cs_ops(cs):
    """[(kind, value)]: (':', n), ('*', (ref, alt)), ('+', bases), ('-', bases), ('~', n); bases upper-case bytes."""
    ops = []
    for m in CS_OP.finditer(cs):
        if m.group(1):
            ops.append((":", int(m.group(1))))
        elif m.group(2):
            ops.append(("*", (m.group(2).upper().encode(), m.group(3).upper().encode())))
        elif m.group(4):
            ops.append(("+", m.group(4).upper().encode()))
        elif m.group(5):
            ops.append(("-", m.group(5).upper().encode()))
        else:
            ops.append(("~", int(m.group(7))))
    return ops


def cigar_from_cs(ops):
    """The CIGAR a cs string implies: matches and substitutions merge into M, insertions are I, deletions D."""
    out = []
    for kind, value in ops:
        if kind == ":":
            op, n = "M", value
        elif kind == "*":
            op, n = "M", 1
        elif kind == "+":
            op, n = "I", len(value)
        elif kind == "-":
            op, n = "D", len(value)
        else:
            op, n = "N", value
        if out and out[-1][0] == op:
            out[-1][1] += n
        else:
            out.append([op, n])
    return [tuple(x) for x in out]


def parse_cigar(cg):
    out = []
    for n, op in CIGAR_OP.findall(cg):
        op = "M" if op in "=X" else op
        if out and out[-1][0] == op:
            out[-1][1] += int(n)
        else:
            out.append([op, int(n)])
    return [tuple(x) for x in out]


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------- intervals (half-open)

def merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(x) for x in merged]


def total(merged):
    return sum(end - start for start, end in merged)


def subtract(a, b):
    """Merged intervals a minus merged intervals b."""
    out, j = [], 0
    for start, end in a:
        cur = start
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < end:
            if b[k][0] > cur:
                out.append((cur, b[k][0]))
            cur = max(cur, b[k][1])
            k += 1
        if cur < end:
            out.append((cur, end))
    return out


def overlap(start, end, lo, hi):
    return max(0, min(end, hi) - max(start, lo))


# ---------------------------------------------------------------- genomes and provenance

def read_fasta(path, upper=True):
    """{name: sequence bytes} for a FASTA, gzip/BGZF or plain; the name is the header's first word."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, name, chunks = {}, None, []
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    seqs[name] = b"".join(chunks)
                name, chunks = line[1:].split(None, 1)[0].decode(), []
            else:
                chunks.append(line.rstrip(b"\r\n"))
    if name is not None:
        seqs[name] = b"".join(chunks)
    return {k: v.upper() for k, v in seqs.items()} if upper else seqs


def provenance_sequences(fasta):
    """{name: entry} from the provenance.json beside a genome file (entries have length, chr_name, role ...)."""
    provenance = json.loads((Path(fasta).parent / "provenance.json").read_text())
    return {s["name"]: s for s in provenance.get("sequences", [])}


def label_number(label):
    """'LG9' and '9' both give '9'; anything else (B, Un, MT, '') is returned as it is."""
    match = re.fullmatch(r"(?:LG)?(\d+)", label or "")
    return match.group(1) if match else (label or "")


# ---------------------------------------------------------------- the version folder

def read_manifest(base):
    with open(Path(base) / "manifest.tsv") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_assembly(row):
    """The panel assembly a run is about: the query of a forward run, the target of a reverse run."""
    return row["query_id"] if row["direction"] == "forward" else row["target_id"]


def read_table(path):
    """Rows of a tab-separated file whose comment lines start with '#'."""
    with open(path) as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------- the coordinator's library

class Deps:
    """What the checks need from ga_common, in a form tests can stand in for."""

    def __init__(self, locate, accessions, reference_id):
        self.locate, self.accessions, self.reference_id = locate, accessions, reference_id

    def panel_ids(self):
        return [a for a in self.accessions if a != self.reference_id and a not in SUPPLEMENTARY]


def load_deps():
    from ga_common.paths import Paths
    from ga_common import panel
    paths = Paths.load()
    pan = panel.read_panel(str(paths.config_dir / "genome_panel_P1.tsv"))
    accessions = {r["id"]: r["accession_genbank"] for r in pan.rows}
    accessions.update(SUPPLEMENTARY)

    def locate(assembly_id):
        return Path(paths.local("genome_fasta", assembly_id=assembly_id, accession_version=accessions[assembly_id]))
    return Deps(locate, accessions, pan.reference()["id"])


def add_throttle_arguments(parser):
    try:
        from ga_common import throttle
        throttle.add_arguments(parser)
    except ImportError:  # the laptop's tests run without ga_common
        parser.add_argument("--threads", type=int, default=None)
        parser.add_argument("--max-jobs", dest="max_jobs", type=int, default=None)
