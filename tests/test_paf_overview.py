"""Tests for scripts/paf_overview.py, on a small hand-written PAF."""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import paf_overview as po  # noqa: E402


def record(q, qs, qe, strand, t, tlen, ts, te, matches, block, tp, cs=True):
    fields = [q, "10000", str(qs), str(qe), strand, t, str(tlen), str(ts), str(te), str(matches), str(block), "60",
              "tp:A:" + tp]
    if cs:
        fields.append("cs:Z::{}".format(matches))
    return "\t".join(fields) + "\n"


PAF = "".join([
    record("q1", 0, 500, "+", "NC_1.1", 1000, 0, 500, 495, 500, "P"),
    record("q1", 500, 800, "-", "NC_1.1", 1000, 450, 750, 297, 300, "P"),     # overlaps the first by 50
    record("q2", 0, 100, "+", "NC_1.1", 1000, 900, 1000, 100, 100, "S"),      # secondary: not counted
    record("q3", 0, 2000, "+", "NC_2.1", 4000, 0, 2000, 1990, 2000, "P", cs=False),
])


class Overview(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.paf = Path(tmp.name) / "x.paf"
        self.paf.write_text(PAF)
        self.labels = Path(tmp.name) / "labels.tsv"
        self.labels.write_text("name\tlength\tchr_name\nNC_1.1\t1000\tLG1\nNC_2.1\t4000\tLG2\n")

    def test_merged_length(self):
        self.assertEqual(po.merged_length([(0, 500), (450, 750)]), 750)
        self.assertEqual(po.merged_length([(10, 20), (0, 5), (12, 15)]), 15)

    def test_summary(self):
        records = [po.parse_record(line) for line in PAF.splitlines(True)]
        per_target, tp_counts, missing_cs = po.summarize(records)
        t = per_target["NC_1.1"]
        self.assertEqual((t["records"], t["covered"], t["minus_bp"], t["matches"], t["block"]), (2, 750, 300, 792, 800))
        self.assertEqual(tp_counts, {"P": 3, "S": 1})
        self.assertEqual(missing_cs, 1)

    def test_printed_tables(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = po.main([str(self.paf), "--labels", str(self.labels), "--chromosome", "LG1", "--min-length", "100"])
        self.assertEqual(code, 0)
        text = out.getvalue()
        self.assertIn("NC_1.1\tLG1\t1000\t2\t750\t0.750\t0.9900\t300\tq1\t0", text)
        self.assertIn("records without a cs tag: 1", text)
        blocks = text.split("in target order:")[1].strip().splitlines()[1:]
        self.assertEqual([b.split("\t")[:2] for b in blocks], [["0", "500"], ["450", "750"]])

    def test_unknown_chromosome(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(po.main([str(self.paf), "--chromosome", "LG9"]), 1)


if __name__ == "__main__":
    unittest.main()
