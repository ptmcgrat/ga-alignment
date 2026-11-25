"""Tests for scripts/build_graph.py, with stand-in taskset and apptainer executables and stand-in genomes.
No Cactus, image, genome or ga_common is needed."""

import contextlib
import gzip
import io
import json
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_graph as bg  # noqa: E402
import pafkit  # noqa: E402

TASKSET = '''#!{python}
import os, sys
os.execvp(sys.argv[3], sys.argv[3:])          # drop "-c <cpus>"
'''
APPTAINER = '''#!{python}
import gzip, os, sys
args = sys.argv[1:]
i = args.index("cactus-pangenome")
os.makedirs(args[i + 1], exist_ok=True)       # the job store
out, name = args[args.index("--outDir") + 1], args[args.index("--outName") + 1]
open(os.path.join(out, "argv.txt"), "a").write(" ".join(args) + "\\n")
gzip.open(os.path.join(out, name + ".full.gfa.gz"), "wb").close()
sys.exit(3 if os.environ.get("STAND_IN_FAIL") else 0)
'''


class StandIn(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        for name, text in (("taskset", TASKSET), ("apptainer", APPTAINER)):
            path = self.bin / name
            path.write_text(text.format(python=sys.executable))
            path.chmod(path.stat().st_mode | stat.S_IXUSR)
        self.acc = {"REF": "GCA_1.1", "A": "GCA_2.1", "B": "GCA_3.1"}
        self.genome("REF", {"R1": ("LG1", "Linkage Group"), "R20": ("LG20", "Linkage Group"), "U1": ("Un", "")})
        self.genome("A", {"a1": ("1", "Chromosome"), "a20": ("20", "Chromosome"), "u": ("Un", "")})
        self.genome("B", {"b20": ("LG20", "Linkage Group"), "bX": ("LG1", "Linkage Group")})
        self.scratch = self.root / "scratch"
        assign = self.scratch / "pairwise" / "PW_v1" / "summary" / "20260101-000000"
        assign.mkdir(parents=True)
        rows = [("A.asm5", "a1", "LG1"), ("A.asm5", "a20", "LG20"), ("A.asm5", "u", "LG20"),
                ("B.asm5", "b20", "LG20"), ("B.asm5", "bX", "LG1")]
        (assign / "assignment.tsv").write_text("run\tquery\tchromosome_label\n" + "".join("\t".join(r) + "\n" for r in rows))

    def genome(self, assembly, seqs):
        folder = self.root / "Genomes" / assembly / self.acc[assembly]
        folder.mkdir(parents=True)
        (folder / "genome.fna.gz").write_bytes(gzip.compress(b"".join(
            b">%s description\nACGTACGT\nNNAC\n" % n.encode() for n in seqs)))
        (folder / "provenance.json").write_text(json.dumps({
            "files": [{"role": "genome", "sha256": "0" * 64}],
            "sequences": [{"name": n, "length": 12, "chr_name": c, "role": "assembled-molecule" if m else "unplaced-scaffold",
                           "molecule_type": m} for n, (c, m) in seqs.items()]}))

    def deps(self):
        return pafkit.Deps(lambda a: self.root / "Genomes" / a / self.acc[a] / "genome.fna.gz", dict(self.acc), "REF")

    def run_main(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        argv = ["--scratch", str(self.scratch), "--sif", str(self.root / "cactus.sif"), "--taskset", str(self.bin / "taskset"),
                "--apptainer", str(self.bin / "apptainer")] + list(argv)
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = bg.main(argv, deps=self.deps())
        return code, out.getvalue(), err.getvalue()

    def folder(self, name):
        return (self.scratch / "graph" / "GR_v1" / name).resolve()


class Planning(StandIn):
    def test_the_subset_is_the_scored_seven(self):
        rows = pafkit.read_table(Path(__file__).resolve().parents[1] / "published" / "expected_inversions.tsv")
        self.assertEqual(set(bg.SUBSET), {r["assembly_id"] for r in rows if r["scored"] == "yes"})

    def test_dry_run_of_a_chromosome_run(self):
        code, out, _ = self.run_main("--assemblies", "A", "B", "--chromosome", "LG20", "--threads", "8")
        self.assertEqual(code, 0, out)
        for piece in ("taskset -c 0-7", "--refContigs R20 --gfa full clip --viz --validate", "--maxCores 8", "--mapCores 6",
                      "--maxMemory 200G", "-> 2 sequence(s) for LG20"):
            self.assertIn(piece, out)
        self.assertNotIn("--otherContig", out)
        self.assertFalse(self.folder("LG20").exists())

    def test_whole_genome_commands(self):
        code, out, _ = self.run_main("--assemblies", "A", "B", "--cpus", "16-31")
        self.assertEqual(code, 0)
        self.assertIn("taskset -c 16-31", out)
        self.assertIn("--refContigs R1 R20 --otherContig chrOther --gfa full clip", out)
        self.assertIn("--maxCores 16", out)
        code, out, _ = self.run_main("--set", "all")
        self.assertIn("--gfa full --viz --validate", out)

    def test_an_assembly_without_sequences_for_the_chromosome_stops(self):
        code, _, err = self.run_main("--assemblies", "A", "B", "--chromosome", "LG1", "--name", "x")
        self.assertEqual(code, 0)                       # both have LG1 sequences
        rows = (self.scratch / "pairwise" / "PW_v1" / "summary" / "20260101-000000" / "assignment.tsv")
        rows.write_text("run\tquery\tchromosome_label\nA.asm5\ta20\tLG20\n")
        code, _, err = self.run_main("--assemblies", "A", "B", "--chromosome", "LG20")
        self.assertEqual(code, 2)
        self.assertIn("no sequence assigned to LG20 for: B", err)


class Running(StandIn):
    def test_a_chromosome_run(self):
        code, out, _ = self.run_main("--assemblies", "A", "B", "--chromosome", "LG20", "--threads", "8", "--go")
        self.assertEqual(code, 0, out)
        folder = self.folder("LG20")
        names = lambda f: [l[1:].strip() for l in (folder / "fasta" / f).read_text().splitlines() if l.startswith(">")]
        self.assertEqual((names("REF.fa"), names("A.fa"), names("B.fa")), (["R20"], ["a20", "u"], ["b20"]))
        self.assertEqual([l.split("\t")[0] for l in (folder / "seqfile.txt").read_text().splitlines()], ["REF", "A", "B"])
        record = json.loads((folder / "run.json").read_text())[-1]
        self.assertEqual((record["exit_status"], record["mode"], len(record["inputs"])), (0, "chromosome", 3))
        self.assertTrue((folder / "out" / "LG20.full.gfa.gz").is_file())

    def test_never_overwrites_and_resumes_with_restart(self):
        self.run_main("--assemblies", "A", "B", "--chromosome", "LG20", "--go")
        code, _, err = self.run_main("--assemblies", "A", "B", "--chromosome", "LG20", "--go")
        self.assertEqual(code, 2)
        self.assertIn("never overwrites", err)
        code, _, _ = self.run_main("--assemblies", "A", "B", "--chromosome", "LG20", "--restart", "--go")
        self.assertEqual(code, 0)
        folder = self.folder("LG20")
        self.assertTrue((folder / "out" / "argv.txt").read_text().splitlines()[-1].endswith("--restart"))
        self.assertEqual(len(json.loads((folder / "run.json").read_text())), 2)

    def test_restart_needs_a_job_store(self):
        code, _, err = self.run_main("--assemblies", "A", "B", "--chromosome", "LG20", "--restart", "--go")
        self.assertEqual(code, 2)
        self.assertIn("needs the job store", err)

    def test_a_failed_run_is_recorded(self):
        with mock.patch.dict("os.environ", {"STAND_IN_FAIL": "1"}):
            code, out, _ = self.run_main("--assemblies", "A", "B", "--chromosome", "LG20", "--go")
        self.assertEqual(code, 1)
        self.assertIn("FAILED", out)
        self.assertEqual(json.loads((self.folder("LG20") / "run.json").read_text())[-1]["exit_status"], 3)


if __name__ == "__main__":
    unittest.main()
