"""Tests for scripts/pafkit.py."""

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pafkit  # noqa: E402

LINE = "q1\t1000\t10\t110\t-\tt1\t5000\t200\t300\t98\t101\t60\ttp:A:P\tde:f:0.0123\tcg:Z:50M1I50M\tcs:Z::50+a:48*ag\n"


class Records(unittest.TestCase):
    def test_parse_record(self):
        r = pafkit.parse_record(LINE, keep_strings=True)
        self.assertEqual((r["qname"], r["qs"], r["qe"], r["strand"], r["tname"], r["ts"], r["te"]),
                         ("q1", 10, 110, "-", "t1", 200, 300))
        self.assertEqual((r["matches"], r["block"], r["tp"], r["de"]), (98, 101, "P", 0.0123))
        self.assertEqual((r["cg"], r["cs"]), ("50M1I50M", ":50+a:48*ag"))
        self.assertIsNone(pafkit.parse_record(LINE)["cs"])

    def test_cs_ops_and_cigar(self):
        ops = pafkit.cs_ops(":50+a:48*ag-tt:2")
        self.assertEqual(ops, [(":", 50), ("+", b"A"), (":", 48), ("*", (b"A", b"G")), ("-", b"TT"), (":", 2)])
        self.assertEqual(pafkit.cigar_from_cs(ops), [("M", 50), ("I", 1), ("M", 49), ("D", 2), ("M", 2)])
        self.assertEqual(pafkit.parse_cigar("50M1I49M2D2M"), [("M", 50), ("I", 1), ("M", 49), ("D", 2), ("M", 2)])
        self.assertEqual(pafkit.parse_cigar("10=1X5="), [("M", 16)])

    def test_revcomp(self):
        self.assertEqual(pafkit.revcomp(b"AACGTN"), b"NACGTT")


class Intervals(unittest.TestCase):
    def test_merge_total_subtract(self):
        merged = pafkit.merge([(5, 10), (0, 3), (8, 12)])
        self.assertEqual(merged, [(0, 3), (5, 12)])
        self.assertEqual(pafkit.total(merged), 10)
        self.assertEqual(pafkit.subtract([(0, 20)], merged), [(3, 5), (12, 20)])
        self.assertEqual(pafkit.subtract([(0, 4), (6, 8)], [(2, 7)]), [(0, 2), (7, 8)])
        self.assertEqual(pafkit.overlap(0, 10, 5, 20), 5)

    def test_label_number(self):
        self.assertEqual([pafkit.label_number(x) for x in ("LG9", "9", "B", "")], ["9", "9", "B", ""])


class Fasta(unittest.TestCase):
    def test_read_fasta_upper(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "g.fna.gz"
            path.write_bytes(gzip.compress(b">s1 words\nacgt\nNNAC\n>s2\nGG\n"))
            self.assertEqual(pafkit.read_fasta(path), {"s1": b"ACGTNNAC", "s2": b"GG"})


if __name__ == "__main__":
    unittest.main()
