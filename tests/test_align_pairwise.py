"""Tests for scripts/align_pairwise.py, with a stand-in minimap2 (a small executable written here) and
stand-in genomes. No real minimap2, genome or ga_common is needed, so they run on the laptop and on mzebra."""

import contextlib
import hashlib
import io
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import align_pairwise as ap  # noqa: E402

STAND_IN_MINIMAP2 = '''#!{python}
import os, sys
args = sys.argv[1:]
if args == ["--version"]:
    print("2.31-stand-in"); sys.exit(0)
if os.environ.get("STAND_IN_FAIL"):
    sys.stderr.write("[ERROR] stand-in failure\\n"); sys.exit(3)
target, query = args[-2], args[-1]
sys.stderr.write("[M::main] stand-in run\\n")
for i in range(2):
    print("\\t".join(["q{{}}".format(i), "100", "0", "100", "+", "t", "100", "0", "100", "99", "100", "60",
                     "tp:A:P", "cs:Z::99", "target=" + target, "query=" + query]))
'''


class StandIn(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.minimap2 = self.root / "minimap2"
        self.minimap2.write_text(STAND_IN_MINIMAP2.format(python=sys.executable))
        self.minimap2.chmod(self.minimap2.stat().st_mode | stat.S_IXUSR)
        self.accessions = {"MZ_REF": "GCA_1.1", "YH_Q": "GCA_2.1", "AC_Q": "GCA_3.1", "YH_LG10Y": "h2tg000014l"}
        for assembly, accession in self.accessions.items():
            folder = self.root / "Genomes" / assembly / accession
            folder.mkdir(parents=True)
            fasta = folder / "genome.fna.gz"
            fasta.write_bytes(assembly.encode())
            (folder / "provenance.json").write_text(json.dumps(
                {"files": [{"role": "genome", "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest()}]}))
        self.scratch = self.root / "scratch"
        self.base = self.scratch / "pairwise" / "PW_v1"

    def deps(self):
        return ap.Deps(lambda a: self.root / "Genomes" / a / self.accessions[a] / "genome.fna.gz",
                       dict(self.accessions), "MZ_REF")

    def run_main(self, *argv, minimap2=None):
        out, err = io.StringIO(), io.StringIO()
        argv = ["--scratch", str(self.scratch), "--minimap2", str(minimap2 or self.minimap2)] + list(argv)
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = ap.main(argv, deps=self.deps())
        return code, out.getvalue(), err.getvalue()

    def manifest(self):
        lines = (self.base / "manifest.tsv").read_text().splitlines()
        header = lines[0].split("\t")
        return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


class Planning(StandIn):
    def test_threads_are_required(self):
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            self.run_main("--assemblies", "YH_Q")

    def test_dry_run_writes_nothing(self):
        code, out, _ = self.run_main("--assemblies", "YH_Q", "--threads", "16")
        self.assertEqual(code, 0)
        self.assertIn("-c -x asm5 --cs -t 16", out)
        self.assertIn("Dry run", out)
        self.assertFalse(self.scratch.exists())

    def test_all_leaves_out_the_reference_and_the_haplotig(self):
        code, out, _ = self.run_main("--all", "--threads", "4")
        self.assertEqual(code, 0)
        self.assertIn("YH_Q.asm5", out)
        self.assertIn("AC_Q.asm5", out)
        self.assertNotIn("MZ_REF.asm5", out)
        self.assertNotIn("YH_LG10Y", out)

    def test_exclude_leaves_assemblies_out(self):
        code, out, _ = self.run_main("--all", "--exclude", "YH_Q", "--threads", "4")
        self.assertEqual(code, 0)
        self.assertNotIn("YH_Q.asm5", out)
        self.assertIn("AC_Q.asm5", out)
        with self.assertRaises(SystemExit):
            self.run_main("--all", "--exclude", "NOT_AN_ID", "--threads", "4")

    def test_unknown_or_reference_ids_are_refused(self):
        for bad in ("NOT_AN_ID", "MZ_REF"):
            with self.assertRaises(SystemExit):
                self.run_main("--assemblies", bad, "--threads", "4")

    def test_missing_minimap2_stops(self):
        code, _, err = self.run_main("--assemblies", "YH_Q", "--threads", "4", minimap2=self.root / "absent")
        self.assertEqual(code, 2)
        self.assertIn("ga-alignment", err)


class Running(StandIn):
    def test_forward_run_writes_paf_log_and_manifest(self):
        code, out, _ = self.run_main("--assemblies", "YH_Q", "--threads", "16", "--go")
        self.assertEqual(code, 0, out)
        paf = self.base / "YH_Q.asm5.paf"
        self.assertEqual(len(paf.read_text().splitlines()), 2)
        self.assertFalse(Path(str(paf) + ".partial").exists())
        log = json.loads((self.base / "logs" / "YH_Q.asm5.json").read_text())
        self.assertEqual((log["exit_status"], log["minimap2_version"], log["paf_records"]), (0, "2.31-stand-in", 2))
        reference = str(self.deps().locate("MZ_REF"))
        query = str(self.deps().locate("YH_Q"))
        self.assertTrue(log["command"].endswith(reference + " " + query))      # target first, then query
        self.assertEqual(log["query_sha256_recorded"], hashlib.sha256(b"YH_Q").hexdigest())
        self.assertIn("stand-in run", (self.base / "logs" / "YH_Q.asm5.stderr").read_text())
        rows = self.manifest()
        self.assertEqual([(r["run"], r["direction"], r["exit_status"]) for r in rows], [("YH_Q.asm5", "forward", "0")])

    def test_reverse_run_swaps_query_and_target(self):
        code, _, _ = self.run_main("--assemblies", "AC_Q", "--reverse", "--threads", "8", "--go")
        self.assertEqual(code, 0)
        self.assertTrue((self.base / "reverse" / "AC_Q.asm5.paf").is_file())
        log = json.loads((self.base / "logs" / "AC_Q.asm5.reverse.json").read_text())
        self.assertEqual((log["query_id"], log["target_id"]), ("MZ_REF", "AC_Q"))
        self.assertTrue(log["command"].endswith(str(self.deps().locate("AC_Q")) + " " + str(self.deps().locate("MZ_REF"))))

    def test_two_at_a_time_and_a_second_preset(self):
        code, _, _ = self.run_main("--assemblies", "YH_Q", "AC_Q", "--threads", "4", "--max-jobs", "2", "--go")
        self.assertEqual(code, 0)
        code, _, _ = self.run_main("--assemblies", "YH_Q", "--preset", "asm10", "--threads", "4", "--go")
        self.assertEqual(code, 0)
        self.assertEqual(sorted(r["run"] for r in self.manifest()), ["AC_Q.asm5", "YH_Q.asm10", "YH_Q.asm5"])

    def test_two_presets_in_one_call(self):
        code, out, _ = self.run_main("--assemblies", "YH_Q", "--preset", "asm5", "asm10", "--threads", "4",
                                     "--max-jobs", "2", "--go")
        self.assertEqual(code, 0, out)
        self.assertEqual(sorted(r["run"] for r in self.manifest()), ["YH_Q.asm10", "YH_Q.asm5"])
        self.assertEqual(len((self.base / "manifest.tsv").read_text().splitlines()), 3)   # one header

    def test_never_overwrites(self):
        self.run_main("--assemblies", "YH_Q", "--threads", "4", "--go")
        code, _, err = self.run_main("--assemblies", "YH_Q", "--threads", "4", "--go")
        self.assertEqual(code, 2)
        self.assertIn("never overwrites", err)
        self.assertIn("rm ", err)

    def test_a_failed_run_keeps_its_output_aside(self):
        with mock.patch.dict(os.environ, {"STAND_IN_FAIL": "1"}):
            code, out, _ = self.run_main("--assemblies", "YH_Q", "--threads", "4", "--go")
        self.assertEqual(code, 1)
        self.assertIn("FAILED", out)
        self.assertFalse((self.base / "YH_Q.asm5.paf").exists())
        self.assertEqual(len(list(self.base.glob("YH_Q.asm5.paf.failed-*"))), 1)
        self.assertEqual(self.manifest()[0]["exit_status"], "3")


if __name__ == "__main__":
    unittest.main()
