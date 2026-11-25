"""Tests for scripts/check_graph.py, on a synthetic 40 kb graph with an inversion (scaled down from LG9)."""

import contextlib
import datetime
import gzip
import io
import json
import random
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_graph as cg  # noqa: E402
import pafkit  # noqa: E402

random.seed(11)
N, L = 400, 100                                   # 400 nodes of 100 bp: a 40 kb reference
NODES = {str(i): bytes(random.choice(b"ACGT") for _ in range(L)) for i in range(1, N + 1)}
SKIP = set(range(351, 361))                       # reference 35,000-36,000 bp: visited by no other walk
REF_WALK = "".join(">%d" % i for i in range(1, N + 1))
A_WALK = ("".join(">%d" % i for i in range(1, 142)) + "".join("<%d" % i for i in range(300, 141, -1))
          + "".join(">%d" % i for i in range(301, N + 1) if i not in SKIP))      # inverted over 14,100-30,000
B_WALK = "".join(">%d" % i for i in range(1, N + 1) if i not in SKIP)
SCALED = ["--gap", "1000", "--min-block", "500", "--min-contig", "2000", "--min-segment", "1000", "--max-jump", "10000",
          "--flank-window", "2000", "--min-background", "500", "--min-interior", "1000", "--refonly-min", "300"]


def gz(path, text):
    Path(path).write_bytes(gzip.compress(text.encode()))


def beds(path, folder, rows):
    with tarfile.open(path, "w:gz") as tar:
        data = "".join("\t".join(map(str, r)) + "\n" for r in rows).encode()
        info = tarfile.TarInfo("{}/MZ_B.bed".format(folder))
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))


class Fixture(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = self.root = Path(tmp.name)
        graph = {"nodes": NODES, "walks": [], "reference": "MZ_REF", "edges": 0}
        spelled = {w: cg.spell(cg.steps(w), NODES) for w in (REF_WALK, A_WALK, B_WALK)}
        self.genomes = {"MZ_REF": {"R9": spelled[REF_WALK]}, "YH_A": {"a9": spelled[A_WALK]},
                        "MZ_B": {"b9": spelled[B_WALK], "bX": b"ACGT" * 125}}
        # the reference genome, with provenance labels
        ref_dir = root / "Genomes" / "MZ_REF" / "GCA_1.1"
        ref_dir.mkdir(parents=True)
        (ref_dir / "provenance.json").write_text(json.dumps({"sequences": [{"name": "R9", "length": N * L, "chr_name": "LG9"}]}))
        # the run folder
        run = self.run = root / "scratch" / "graph" / "GR_v1" / "LG9"
        out, stats = run / "out", run / "out" / "LG9.stats"
        stats.mkdir(parents=True)
        (run / "fasta").mkdir()
        lines = []
        for sample, seqs in self.genomes.items():
            fasta = run / "fasta" / (sample + ".fa")
            fasta.write_bytes(b"".join(b">%s\n%s\n" % (n.encode(), s) for n, s in seqs.items()))
            lines.append("{}\t{}".format(sample, fasta))
        (run / "seqfile.txt").write_text("\n".join(lines) + "\n")
        (run / "run.json").write_text(json.dumps([{"exit_status": 0, "wall_s": 60, "cpu_s": 7200, "largest_process_peak_gb": 1.0}]))
        walks = [("MZ_REF", "R9", REF_WALK), ("YH_A", "a9", A_WALK), ("MZ_B", "b9", B_WALK)]
        gfa = "H\tVN:Z:1.1\tRS:Z:MZ_REF\n" + "".join("S\t{}\t{}\n".format(n, s.decode()) for n, s in NODES.items()) + \
            "L\t1\t+\t2\t+\t0M\n" + "".join("W\t{}\t0\t{}\t0\t{}\t{}\n".format(s, c, len(cg.spell(cg.steps(w), NODES)), w) for s, c, w in walks)
        gz(out / "LG9.full.gfa.gz", gfa)
        gz(out / "LG9.gfa.gz", gfa)
        sizes = {(s, c): len(seq) for s, seqs in self.genomes.items() for c, seq in seqs.items()}
        gz(out / "LG9.input-contig-sizes.tsv.gz", "#event\tpansn_prefix\tbase_contig\toffset\tlength\n" + "".join(
            "{}\t{}#0\t{}\t0\t{}\n".format(s, s, c, n) for (s, c), n in sizes.items()))
        rows = [(s, c, "R9" if c != "bX" else "_NONE_", n, n if c != "bX" else 0, 1 if c != "bX" else 0) for (s, c), n in sizes.items()]
        gz(stats / "clipped-by-input-contig.tsv.gz", "# test\n#genome\tcontig\tref_chrom\tinput_bp\tfull_bp\tfull_frags\tclip_bp\tclip_frags\n"
           + "".join("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(s, c, r, n, b, f, b, f) for s, c, r, n, b, f in rows))
        (stats / "clipped-by-reference-contig.tsv").write_text(
            "# test\n#ref_chrom\tgenome\tcontigs\tinput_bp\tfull_bp\tfull_frags\tclip_bp\tclip_frags\n" + "".join(
                "{}\t{}\t1\t{}\t{}\t{}\t{}\t{}\n".format(r, s, n, b, f, b, f) for s, c, r, n, b, f in rows))
        beds(stats / "clipped.full.beds.tar.gz", "clipped.full.beds", [("bX", 0, 500, "unassigned")])
        beds(stats / "clipped.beds.tar.gz", "clipped.beds", [("bX", 0, 500, "unassigned")])
        (stats / "clipped-by-genome.tsv").write_text("#genome\treason\tintervals\tbp\tpct_of_input\nMZ_B\tunassigned\t1\t500\t1.3\n")
        gz(stats / "refgaps.bed.gz", "R9\t35000\t36000\trefgap\n")
        (stats / "inversion-stats.tsv").write_text("#ref_chrom\tsevered\tfound\tmessage\nR9\t0\t1\tnone\n")
        # the pairwise side: PAFs covering all but the reference-only stretch, and the PAF calls
        pw = root / "scratch" / "pairwise" / "PW_v1"
        (pw / "summary" / "20260101-000000").mkdir(parents=True)
        (pw / "manifest.tsv").write_text("run\tdirection\tpreset\tquery_id\tpaf\texit_status\n"
                                         "YH_A.asm5\tforward\tasm5\tYH_A\tA.paf\t0\nMZ_B.asm5\tforward\tasm5\tMZ_B\tB.paf\t0\n")
        for name in ("A.paf", "B.paf"):
            (pw / name).write_text("".join("q\t1\t0\t1\t+\tR9\t{}\t{}\t{}\t1\t1\t60\ttp:A:P\n".format(N * L, a, b)
                                           for a, b in ((0, 35000), (36000, 40000))))
        (pw / "summary" / "20260101-000000" / "inversions.tsv").write_text(
            "assembly\tunit\tpreset\tcall\tbp1_observed\tbp2_observed\n"
            "YH_A\t9\tasm5\tinverted\t14100-14100\t30000-30000\nMZ_B\t9\tasm5\tnormal\t\t\n")
        (root / "units.tsv").write_text("unit\tchromosome\tbp1_start\tbp1_end\tbp2_start\tbp2_end\n9\tLG9\t14000\t14100\t30000\t30200\n")
        (root / "expected.tsv").write_text("assembly_id\tunit\texpected\tscored\tbasis\nYH_A\t9\tinverted\tyes\tt\nMZ_B\t9\tnormal\tyes\tt\n")
        self.graph = graph
        graph["walks"] = [{"sample": s, "hap": "0", "contig": c, "start": 0, "end": len(cg.spell(cg.steps(w), NODES)), "walk": w}
                          for s, c, w in walks]

    def run_main(self, *extra, stamp=5):
        deps = pafkit.Deps(lambda a: self.root / "Genomes" / a / "GCA_1.1" / "genome.fna.gz", {"MZ_REF": "GCA_1.1"}, "MZ_REF")
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            code = cg.main(["--scratch", str(self.root / "scratch"), "--units", str(self.root / "units.tsv"),
                            "--expected", str(self.root / "expected.tsv")] + SCALED + list(extra), deps=deps,
                           now=datetime.datetime(2026, 1, 2, 3, 4, stamp, tzinfo=datetime.timezone.utc))
        return code, out.getvalue()


class Pieces(Fixture):
    def test_projection_gives_the_inverted_block(self):
        pos = cg.reference_positions(self.graph, "MZ_REF")
        a = next(w for w in self.graph["walks"] if w["sample"] == "YH_A")
        blocks = cg.project(a, NODES, pos, 1000)
        self.assertEqual([(b["strand"], b["ts"], b["te"]) for b in blocks],
                         [("+", 0, 14100), ("-", 14100, 30000), ("+", 30000, 40000)])

    def test_exact_passes_and_finds_a_wrong_base(self):
        problems, spelled = cg.exact_checks(self.graph, "full", self.genomes, "MZ_REF")
        self.assertEqual(problems, [])
        bad = dict(self.genomes, YH_A={"a9": self.genomes["YH_A"]["a9"][:5000] + b"N" + self.genomes["YH_A"]["a9"][5001:]})
        problems, _ = cg.exact_checks(self.graph, "full", bad, "MZ_REF")
        self.assertEqual(len(problems), 1)
        self.assertIn("at offset 5000", problems[0])

    def test_the_reference_must_be_forward_and_visit_no_node_twice(self):
        self.graph["walks"][0]["walk"] = REF_WALK.replace(">5>", "<5>").replace(">7>", ">5>", 1)
        problems, _ = cg.exact_checks(self.graph, "full", self.genomes, "MZ_REF")
        self.assertTrue(any("not forward only" in p for p in problems))
        self.assertTrue(any("visits a node twice" in p for p in problems))


class CommandLine(Fixture):
    def test_a_good_run_passes_every_check(self):
        code, out = self.run_main("--go")
        self.assertEqual(code, 0, out)
        folder = self.root / "scratch" / "graph" / "GR_v1" / "checks" / "20260102-030405"
        inv = {r["assembly"]: r for r in pafkit.read_table(folder / "inversions.tsv")}
        self.assertEqual((inv["YH_A"]["call"], inv["YH_A"]["verdict"]), ("inverted", "agree"))
        self.assertEqual(inv["YH_A"]["graph_vs_paf"], "same call; breakpoints overlap: yes, yes")
        self.assertEqual((inv["MZ_B"]["call"], inv["MZ_B"]["graph_vs_paf"]), ("normal", "same call"))
        ref = pafkit.read_table(folder / "refonly.tsv")[0]
        self.assertEqual((ref["graph_only_bp"], ref["paf_only_bp"], ref["graph_in_paf"], ref["graph_in_cactus"]), ("1000", "1000", "1.000", "1.000"))
        self.assertEqual([r["status"] for r in pafkit.read_table(folder / "accounting.tsv")], ["pass"])
        self.assertIn("known features: PASS", out)

    def test_an_accounting_difference_fails(self):
        table = self.run / "out" / "LG9.stats" / "clipped-by-reference-contig.tsv"
        table.write_text(table.read_text().replace("\t1\t39000\t39000\t1\t39000\t1\n", "\t1\t39000\t38999\t1\t39000\t1\n", 1))
        code, out = self.run_main("--go", stamp=6)
        self.assertEqual(code, 1)
        self.assertIn("accounting: 0 of 1 runs", out)

    def test_dry_run_writes_nothing(self):
        code, out = self.run_main()
        self.assertEqual(code, 0)
        self.assertIn("1 finished run(s)", out)
        self.assertFalse((self.root / "scratch" / "graph" / "GR_v1" / "checks").exists())


if __name__ == "__main__":
    unittest.main()
