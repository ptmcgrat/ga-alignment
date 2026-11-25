"""Tests for scripts/summarize_pairwise.py, on synthetic chromosomes built from 1 Mb blocks.
No real PAF, genome or ga_common is needed."""

import contextlib
import datetime
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pafkit  # noqa: E402
import summarize_pairwise as sp  # noqa: E402

MB = 1_000_000
L = 40 * MB
P = SimpleNamespace(min_identity=0.30, min_block=10_000, min_contig=2 * MB, min_segment=100_000, max_jump=MB,
                    call_fraction=0.80, min_background=500_000, min_interior=100_000, flank_window=2 * MB)
UNIT20A = {"unit": "20a", "chromosome": "LG20", "bp1_start": 19_600_000, "bp1_end": 19_620_000, "bp2_start": 29_680_000, "bp2_end": 29_700_000}
UNIT20B = {"unit": "20b", "chromosome": "LG20", "bp1_start": 29_680_000, "bp1_end": 29_700_000, "bp2_start": 33_170_000, "bp2_end": 33_250_000}
UNIT9 = {"unit": "9", "chromosome": "LG9", "bp1_start": 14 * MB, "bp1_end": 14_100_000, "bp2_start": 30 * MB, "bp2_end": 30_200_000}
UNIT11A = {"unit": "11a", "chromosome": "LG11", "bp1_start": 8 * MB, "bp1_end": 8_100_000, "bp2_start": 18 * MB, "bp2_end": 18 * MB}
UNIT11B = {"unit": "11b", "chromosome": "LG11", "bp1_start": 18 * MB, "bp1_end": 18 * MB, "bp2_start": 30 * MB, "bp2_end": 30_100_000}


def block(q, qs, qe, strand, t, ts, te, ident=0.99, tp="P"):
    n = te - ts
    return {"qname": q, "qlen": L, "qs": qs, "qe": qe, "strand": strand, "tname": t, "tlen": L, "ts": ts, "te": te,
            "matches": int(n * ident), "block": n, "mapq": 60, "tp": tp, "de": 0.005, "cs": None, "cg": None}


def pieces(spec, q="Q", t="R", step=MB):
    """spec: [(ref_lo, ref_hi, strand, query_of)] where query_of(x) maps a reference position to the query."""
    out = []
    for lo, hi, strand, query_of in spec:
        for x in range(lo, hi, step):
            y = min(x + step, hi)
            a, b = sorted((query_of(x), query_of(y)))
            out.append(block(q, a, b, strand, t, x, y))
    return out


fwd = lambda x: x
inv = lambda lo, hi: (lambda x: lo + hi - x)


class Calls(unittest.TestCase):
    def call(self, blocks, unit=UNIT9, units=(UNIT9,)):
        return sp.call_unit(unit, list(units), blocks, P)

    def test_a_plain_inversion_and_its_breakpoints(self):
        r = self.call(pieces([(0, 14 * MB, "+", fwd), (14_200_000, 29_900_000, "-", inv(14_200_000, 29_900_000)),
                              (30_300_000, L, "+", fwd)]))
        self.assertEqual((r["call"], r["opposite_fraction"]), ("inverted", "1.000"))
        self.assertEqual((r["bp1_observed"], r["bp1_overlap"]), ("14000000-14200000", "yes"))
        self.assertEqual((r["bp2_observed"], r["bp2_overlap"]), ("29900000-30300000", "yes"))
        self.assertEqual(r["bp1_distance"], 100_000)

    def test_no_inversion(self):
        self.assertEqual(self.call(pieces([(0, L, "+", fwd)]))["call"], "normal")

    def test_mixed_strands_are_unresolved(self):
        r = self.call(pieces([(0, 22 * MB, "+", fwd), (22 * MB, 30_200_000, "-", inv(22 * MB, 30_200_000)),
                              (30_200_000, L, "+", fwd)]))      # flanks agree; half the interior on each strand
        self.assertEqual((r["call"], r["reason"]), ("unresolved", "mixed strands inside"))

    def test_flanks_split_evenly_are_unresolved(self):
        r = self.call(pieces([(0, 22 * MB, "+", fwd), (22 * MB, L, "-", inv(22 * MB, L))]))
        self.assertEqual((r["call"], r["reason"]), ("unresolved", "the flanks are on both strands equally"))

    def test_interior_on_another_sequence_is_unresolved(self):
        blocks = pieces([(0, 14 * MB, "+", fwd), (30_300_000, L, "+", fwd)]) + \
            pieces([(14_200_000, 29_900_000, "-", inv(14_200_000, 29_900_000))], q="Qb")
        r = self.call(blocks)
        self.assertEqual((r["call"], r["main_query"]), ("unresolved", "Qb"))
        self.assertIn("flanks on other sequences", r["reason"])

    def test_a_chromosome_stored_reversed(self):
        rev = lambda x: L - x
        r = self.call(pieces([(0, 14 * MB, "-", rev), (14_200_000, 29_900_000, "+", lambda x: L - (14_200_000 + 29_900_000 - x)),
                              (30_300_000, L, "-", rev)]))
        self.assertEqual((r["call"], r["background_strand"], r["bp1_overlap"], r["bp2_overlap"]), ("inverted", "-", "yes", "yes"))

    def test_tandem_inversions_split_at_the_query_jump(self):
        blocks = pieces([(0, 8 * MB, "+", fwd), (8_200_000, 17_900_000, "-", inv(8_200_000, 17_900_000)),
                         (18_100_000, 29_900_000, "-", inv(18_100_000, 29_900_000)), (30_200_000, L, "+", fwd)])
        a = self.call(blocks, UNIT11A, (UNIT11A, UNIT11B))
        b = self.call(blocks, UNIT11B, (UNIT11A, UNIT11B))
        self.assertEqual((a["call"], a["bp2_observed"], a["bp2_overlap"]), ("inverted", "17900000-18100000", "yes"))
        self.assertEqual((b["call"], b["bp1_observed"], b["bp1_overlap"]), ("inverted", "17900000-18100000", "yes"))

    def test_a_misplaced_flank_does_not_outvote_the_neighbours(self):
        # YH's LG20: reference 0-17.4 Mb sits at the start on "+"; 20a and 20b are each "+" between "-" neighbours
        blocks = pieces([(0, 17_400_000, "+", fwd),
                         (33_270_000, 38_800_000, "-", lambda x: 18 * MB + (38_800_000 - x)),
                         (29_700_000, 33_250_000, "+", lambda x: 23_600_000 + (x - 29_700_000)),
                         (19_620_000, 29_680_000, "+", lambda x: 27_200_000 + (x - 19_620_000)),
                         (18_650_000, 19_620_000, "-", lambda x: 37_300_000 + (19_620_000 - x))])
        units = (UNIT20A, UNIT20B)
        a, b = self.call(blocks, UNIT20A, units), self.call(blocks, UNIT20B, units)
        self.assertEqual((a["call"], a["background_strand"], a["bp1_overlap"], a["bp2_overlap"]), ("inverted", "-", "yes", "yes"))
        self.assertEqual((b["call"], b["bp1_overlap"], b["bp2_overlap"]), ("inverted", "yes", "yes"))

    def test_a_segment_reaching_slightly_into_a_unit_is_not_absorbed(self):
        unit_a = dict(UNIT11A, bp2_start=17_970_000, bp2_end=18_010_000)
        unit_b = dict(UNIT11B, bp1_start=17_970_000, bp1_end=18_010_000)
        blocks = pieces([(0, 8 * MB, "+", fwd), (8_200_000, 17_900_000, "-", inv(8_200_000, 17_900_000)),
                         (17_965_000, 29_900_000, "-", inv(17_965_000, 29_900_000)), (30_200_000, L, "+", fwd)])
        a = self.call(blocks, unit_a, (unit_a, unit_b))      # 11b's segment starts 5 kb inside 11a
        self.assertEqual((a["call"], a["inverted_end"], a["bp2_observed"]), ("inverted", 17_900_000, "17900000-17965000"))

    def test_breakpoints_are_scored_only_for_the_source_species(self):
        rows = [{"assembly": "YH_A", "unit": "9", "preset": "asm5", "call": "inverted", "bp1_overlap": "no", "bp2_overlap": "yes"},
                {"assembly": "DT_B", "unit": "9", "preset": "asm5", "call": "inverted", "bp1_overlap": "no", "bp2_overlap": "yes"}]
        expected = {(a, "9"): {"expected": "inverted", "scored": "yes", "basis": ""} for a in ("YH_A", "DT_B")}
        sp.apply_expectations(rows, expected, "asm5", "YH")
        self.assertEqual([(r["breakpoints_scored"], r["verdict"]) for r in rows],
                         [("yes", "disagree: breakpoint"), ("no", "agree")])

    def test_the_haplotig_is_left_out_of_outliers(self):
        rows = [{"run": "A.asm5", "direction": "forward", "preset": "asm5", "query_id": "A"},
                {"run": "YH_LG10Y.asm5", "direction": "forward", "preset": "asm5", "query_id": "YH_LG10Y"},
                {"run": "A.asm5.reverse", "direction": "reverse", "preset": "asm5", "query_id": "REF"}]
        self.assertEqual(sp.outlier_runs(rows, "asm5"), {"A.asm5"})

    def test_kumar_filters(self):
        keep = pieces([(0, 3 * MB, "+", fwd)])
        low = [block("Q", 0, MB, "+", "R", 3 * MB, 4 * MB, ident=0.2)]
        short = [block("Q", 0, 5000, "+", "R", 5 * MB, 5 * MB + 5000)]
        secondary = [block("Q", 0, MB, "+", "R", 6 * MB, 7 * MB, tp="S")]
        small_contig = pieces([(8 * MB, 9 * MB, "+", fwd)], q="tiny")
        self.assertEqual(sp.kumar_blocks(keep + low + short + secondary + small_contig, P), keep)

    def test_verdicts(self):
        self.assertEqual(sp.verdict("inverted", "inverted", "yes", "yes"), "agree")
        self.assertEqual(sp.verdict("inverted", "inverted", "yes", "no"), "disagree: breakpoint")
        self.assertEqual(sp.verdict("inverted", "inverted", "unresolved", "yes"), "agree; a breakpoint could not be placed")
        self.assertEqual(sp.verdict("inverted", "normal", "", ""), "disagree")
        self.assertEqual(sp.verdict("normal", "unresolved", "", ""), "unresolved")
        self.assertEqual(sp.verdict("heterozygous", "inverted", "", ""), "not compared (heterozygous)")
        self.assertEqual(sp.verdict("", "normal", "", ""), "no expectation")
        self.assertEqual(sp.verdict("inverted", "inverted", "no", "yes", breakpoints_scored=False), "agree")


class CommandLine(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.acc = {"REF": "GCA_1.1", "A": "GCA_2.1", "B": "GCA_3.1"}
        labels = {"REF": {"R": "LG9", "U": "Un"}, "A": {"Q": "9"}, "B": {"Q": "9"}}
        for assembly, seqs in labels.items():
            folder = self.root / "Genomes" / assembly / self.acc[assembly]
            folder.mkdir(parents=True)
            (folder / "provenance.json").write_text(json.dumps(
                {"sequences": [{"name": n, "length": L, "chr_name": c} for n, c in seqs.items()]}))
        self.base = self.root / "scratch" / "pairwise" / "PW_v1"
        self.base.mkdir(parents=True)
        inverted = pieces([(0, 14 * MB, "+", fwd), (14_200_000, 29_900_000, "-", inv(14_200_000, 29_900_000)), (30_300_000, L, "+", fwd)])
        self.write_paf("A.asm5.paf", inverted)
        self.write_paf("B.asm5.paf", pieces([(0, L, "+", fwd)]))
        self.write_paf("reverse/A.asm5.paf", [dict(b, qname=b["tname"], tname=b["qname"], qs=b["ts"], qe=b["te"], ts=b["qs"], te=b["qe"]) for b in inverted])
        cols = ["version", "run", "direction", "preset", "query_id", "target_id", "paf", "paf_bytes", "paf_records",
                "wall_s", "cpu_s", "peak_rss_gb", "exit_status"]
        runs = [["A.asm5", "forward", "asm5", "A", "REF", "A.asm5.paf"], ["B.asm5", "forward", "asm5", "B", "REF", "B.asm5.paf"],
                ["A.asm5.reverse", "reverse", "asm5", "REF", "A", "reverse/A.asm5.paf"]]
        (self.base / "manifest.tsv").write_text("\t".join(cols) + "\n" + "".join(
            "\t".join(["PW_v1"] + r + ["1000", "40", "60", "120", "2.5", "0"]) + "\n" for r in runs))
        (self.root / "units.tsv").write_text("# test units\nunit\tchromosome\tbp1_start\tbp1_end\tbp2_start\tbp2_end\n"
                                             "9\tLG9\t14000000\t14100000\t30000000\t30200000\n")
        self.write_expected("normal")

    def write_paf(self, name, blocks):
        path = self.base / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join("\t".join(map(str, [b["qname"], b["qlen"], b["qs"], b["qe"], b["strand"], b["tname"], b["tlen"],
                                                    b["ts"], b["te"], b["matches"], b["block"], 60, "tp:A:P", "de:f:0.005"])) + "\n"
                                for b in blocks))

    def write_expected(self, b_state):
        (self.root / "expected.tsv").write_text("# test\nassembly_id\tunit\texpected\tscored\tbasis\n"
                                                "A\t9\tinverted\tyes\ttest\nB\t9\t%s\tyes\ttest\n" % b_state)

    def run_main(self, *argv, stamp=5):
        deps = pafkit.Deps(lambda a: self.root / "Genomes" / a / self.acc[a] / "genome.fna.gz", self.acc, "REF")
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = sp.main(["--scratch", str(self.root / "scratch"), "--units", str(self.root / "units.tsv"),
                            "--expected", str(self.root / "expected.tsv")] + list(argv), deps=deps,
                           now=datetime.datetime(2026, 1, 2, 3, 4, stamp, tzinfo=datetime.timezone.utc))
        return code, out.getvalue(), err.getvalue()

    def test_dry_run_writes_nothing(self):
        code, out, _ = self.run_main()
        self.assertEqual(code, 0)
        self.assertFalse((self.base / "summary").exists())

    def test_go_passes_and_writes_everything(self):
        code, out, _ = self.run_main("--go")
        self.assertEqual(code, 0, out)
        folder = self.base / "summary" / "20260102-030405"
        for name in ("chromosomes.tsv", "assignment.tsv", "segments.tsv", "inversions.tsv", "outliers.tsv",
                     "reciprocity.tsv", "report.md"):
            self.assertTrue((folder / name).is_file(), name)
        report = (folder / "report.md").read_text()
        self.assertIn("PASS", report)
        self.assertIn("Tables used: units.tsv (sha256 ", report)
        assign = pafkit.read_table(folder / "assignment.tsv")
        self.assertEqual({(a["run"], a["label_agrees"]) for a in assign}, {("A.asm5", "yes"), ("B.asm5", "yes")})
        recip = pafkit.read_table(folder / "reciprocity.tsv")
        self.assertEqual(recip[0]["difference_bp"], "0")

    def test_nothing_scored_is_not_a_pass(self):
        (self.root / "expected.tsv").write_text("assembly_id\tunit\texpected\tscored\tbasis\n")
        code, out, _ = self.run_main("--go", stamp=7)
        self.assertEqual(code, 1)
        self.assertIn("NOTHING SCORED", out)

    def test_a_disagreement_fails(self):
        self.write_expected("inverted")
        code, out, _ = self.run_main("--go", stamp=6)
        self.assertEqual(code, 1)
        self.assertIn("B 9: expected inverted, called normal", out)


if __name__ == "__main__":
    unittest.main()
