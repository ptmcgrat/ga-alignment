"""Tests for scripts/check_pairwise.py, on simulated genomes and alignments whose cs and CIGAR are true by
construction. No real genome, PAF or ga_common is needed."""

import contextlib
import datetime
import gzip
import io
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_pairwise as cp  # noqa: E402
import pafkit  # noqa: E402

random.seed(7)
DNA = lambda n: bytes(random.choice(b"ACGT") for _ in range(n))


def simulate(target, start, ops):
    """Apply edits to target[start:...]; return the query segment, cs, cg, target end, matches and block length."""
    t, seg, cs, cigar, matches, block = start, bytearray(), [], [], 0, 0

    def add(op, n):
        if cigar and cigar[-1][0] == op:
            cigar[-1][1] += n
        else:
            cigar.append([op, n])
    for kind, value in ops:
        if kind == "=":
            seg += target[t:t + value]; cs.append(":%d" % value); add("M", value); t += value; matches += value; block += value
        elif kind == "X":
            seg += value; cs.append("*%s%s" % (chr(target[t]).lower(), value.decode().lower())); add("M", 1); t += 1; block += 1
        elif kind == "I":
            seg += value; cs.append("+" + value.decode().lower()); add("I", len(value)); block += len(value)
        else:
            cs.append("-" + target[t:t + value].decode().lower()); add("D", value); t += value; block += value
    return bytes(seg), "".join(cs), "".join("%d%s" % (n, op) for op, n in cigar), t, matches, block


OPS = [("=", 300), ("X", b"G"), ("=", 200), ("I", b"TTA"), ("=", 150), ("D", 7), ("=", 340)]


class StandIn(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        target = {"T1": DNA(3000), "T2": DNA(2000)}
        target["T2"] = target["T2"][:1700] + b"N" * 300
        seg1, cs1, cg1, te1, m1, b1 = simulate(target["T1"], 100, OPS)
        seg2, cs2, cg2, te2, m2, b2 = simulate(target["T2"], 500, OPS)
        query = {"Q1": DNA(50) + seg1 + DNA(70), "Q2": DNA(30) + pafkit.revcomp(seg2) + DNA(40),
                 "Q3": target["T1"][2000:2100] + DNA(400)}
        self.lines = [
            self.paf_line("Q1", query, 50, 50 + len(seg1), "+", "T1", target, 100, te1, m1, b1, "P", cg1, cs1),
            self.paf_line("Q2", query, 30, 30 + len(seg2), "-", "T2", target, 500, te2, m2, b2, "P", cg2, cs2),
            self.paf_line("Q3", query, 0, 100, "+", "T1", target, 2000, 2100, 100, 100, "S", "100M", ":100"),
        ]
        self.target, self.query = target, query
        self.write_genome("REF", "GCA_1.1", target)
        self.write_genome("ASM", "GCA_2.1", query)
        self.base = self.root / "scratch" / "pairwise" / "PW_v1"
        self.base.mkdir(parents=True)
        (self.base / "ASM.asm5.paf").write_text("".join(self.lines))
        cols = ["version", "run", "direction", "preset", "query_id", "target_id", "paf", "exit_status"]
        (self.base / "manifest.tsv").write_text("\t".join(cols) + "\n" + "\t".join(
            ["PW_v1", "ASM.asm5", "forward", "asm5", "ASM", "REF", "ASM.asm5.paf", "0"]) + "\n")

    @staticmethod
    def paf_line(q, query, qs, qe, strand, t, target, ts, te, matches, block, tp, cg, cs):
        return "\t".join(map(str, [q, len(query[q]), qs, qe, strand, t, len(target[t]), ts, te, matches, block, 60,
                                   "tp:A:" + tp, "de:f:0.01", "cg:Z:" + cg, "cs:Z:" + cs])) + "\n"

    def write_genome(self, assembly, accession, seqs):
        folder = self.root / "Genomes" / assembly / accession
        folder.mkdir(parents=True)
        (folder / "genome.fna.gz").write_bytes(gzip.compress(b"".join(b">%s x\n%s\n" % (k.encode(), v) for k, v in seqs.items())))
        (folder / "provenance.json").write_text(json.dumps(
            {"sequences": [{"name": k, "length": len(v), "chr_name": k[-1]} for k, v in seqs.items()]}))

    def deps(self, accessions=None):
        accessions = accessions or {"REF": "GCA_1.1", "ASM": "GCA_2.1"}
        return pafkit.Deps(lambda a: self.root / "Genomes" / a / accessions[a] / "genome.fna.gz", accessions, "REF")

    def check(self, line):
        prov_q = {k: {"length": len(v)} for k, v in self.query.items()}
        prov_t = {k: {"length": len(v)} for k, v in self.target.items()}
        return cp.check_record(pafkit.parse_record(line, keep_strings=True), self.query, self.target, prov_q, prov_t)

    def run_main(self, *argv, deps=None):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = cp.main(["--scratch", str(self.root / "scratch")] + list(argv), deps=deps or self.deps(),
                           now=datetime.datetime(2026, 1, 2, 3, 4, 5, tzinfo=datetime.timezone.utc))
        return code, out.getvalue(), err.getvalue()


class Records(StandIn):
    def test_true_records_pass_on_both_strands(self):
        self.assertEqual(self.check(self.lines[0]), [])
        self.assertEqual(self.check(self.lines[1]), [])

    def test_a_wrong_substitution_is_caught(self):
        head, cs = self.lines[0].split("cs:Z:")
        i = cs.index("*") + 2                      # the query base of the substitution
        flipped = cs[:i] + ("a" if cs[i] != "a" else "c") + cs[i + 1:]
        self.assertEqual([k for k, _ in self.check(head + "cs:Z:" + flipped)], ["cs_sequence"])

    def test_a_wrong_cigar_is_caught(self):
        head, rest = self.lines[0].split("cg:Z:")
        first = rest.index("M")                    # the first M one base short
        broken = head + "cg:Z:" + str(int(rest[:first]) - 1) + rest[first:]
        self.assertNotEqual(broken, self.lines[0])
        self.assertEqual([k for k, _ in self.check(broken)], ["cs_cigar"])

    def test_wrong_lengths_and_counts_are_caught(self):
        f = self.lines[0].split("\t")
        self.assertEqual([k for k, _ in self.check("\t".join([f[0], str(int(f[1]) + 1)] + f[2:]))], ["names_lengths"])
        self.assertEqual([k for k, _ in self.check("\t".join(f[:8] + [str(int(f[8]) + 1)] + f[9:]))], ["cs_counts"])
        self.assertEqual([k for k, _ in self.check("\t".join(f[:9] + [str(int(f[9]) - 1)] + f[10:]))], ["cs_counts"])

    def test_accounting_adds_up(self):
        exact, sequences = cp.check_run({"run": "r", "direction": "forward", "preset": "asm5", "query_id": "ASM",
                                         "target_id": "REF"}, str(self.base / "ASM.asm5.paf"),
                                        str(self.deps().locate("ASM")), str(self.deps().locate("REF")))
        self.assertEqual(exact["records"], 3)
        for s in sequences:
            self.assertEqual(s["primary_bp"] + s["secondary_only_bp"] + s["n_bp"] + s["unaligned_bp"], s["length"])
        t2 = next(s for s in sequences if s["side"] == "target" and s["sequence"] == "T2")
        self.assertEqual(t2["n_bp"], 300)
        t1 = next(s for s in sequences if s["side"] == "target" and s["sequence"] == "T1")
        self.assertEqual(t1["secondary_only_bp"], 100)


class CommandLine(StandIn):
    def test_dry_run_then_go_then_never_overwrite(self):
        code, out, _ = self.run_main()
        self.assertEqual(code, 0)
        self.assertIn("Dry run", out)
        self.assertFalse((self.base / "checks").exists())
        code, out, _ = self.run_main("--go")
        folder = self.base / "checks" / "20260102-030405"
        for name in ("exact.tsv", "accounting.tsv", "accounting_sequences.tsv", "report.md"):
            self.assertTrue((folder / name).is_file(), name)
        self.assertIn("pass", (folder / "exact.tsv").read_text())
        self.assertEqual(code, 0, out)
        code, _, err = self.run_main("--go")
        self.assertEqual(code, 2)
        self.assertIn("never overwrites", err)

    def test_a_missing_panel_assembly_fails(self):
        deps = self.deps({"REF": "GCA_1.1", "ASM": "GCA_2.1", "OTHER": "GCA_3.1"})
        code, out, _ = self.run_main("--go", deps=deps)
        self.assertEqual(code, 1)
        self.assertIn("OTHER", out)


if __name__ == "__main__":
    unittest.main()
