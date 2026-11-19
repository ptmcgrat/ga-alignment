"""Tests for scripts/check_inputs.py, on small stand-in genomes with their own provenance.json.
No real genome, no data root and no ga_common is needed, so they run on the laptop and on mzebra."""

import contextlib
import datetime
import gzip
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_inputs as ci  # noqa: E402

# name, sequence, chr_name, role, molecule_type
REFERENCE = [("NC_1.1", "ACGTACGTNNacgt", "LG1", "assembled-molecule", "Linkage Group"),
             ("NC_2.1", "ACGTACGTACGTACGTAC", "LG2", "assembled-molecule", "Linkage Group"),
             ("NW_9.1", "acgtNNNN", "Un", "unplaced-scaffold", ""),
             ("NC_MT.1", "ACGTAC", "MT", "assembled-molecule", "Mitochondrion")]
QUERY = [("OZ_1.1", "ACGTACGTAC", "1", "assembled-molecule", "Chromosome"),
         ("OZ_2.1", "ACGTACGTACGTAC", "2", "assembled-molecule", "Chromosome")]


def write_assembly(root, assembly_id, accession, fetched, records, kind="genome", original=None):
    folder = root / "Genomes" / assembly_id / accession
    folder.mkdir(parents=True)
    fasta = folder / "genome.fna.gz"
    text = "".join(">{} description words\n{}\n".format(name, seq) for name, seq, *_ in records)
    fasta.write_bytes(gzip.compress(text.encode()))
    sequences = []
    for name, seq, chr_name, role, molecule_type in records:
        entry = {"name": name, "length": len(seq)}
        if chr_name is not None:
            entry.update({"chr_name": chr_name, "role": role, "molecule_type": molecule_type})
        sequences.append(entry)
    provenance = {"schema": "ga-provenance-1", "assembly_id": assembly_id, "accession_requested": accession,
                  "accession_fetched": fetched, "kind": kind, "sequences": sequences,
                  "files": [{"role": "genome", "name": "genome.fna.gz", "bytes": fasta.stat().st_size,
                             "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
                             "original_name": original or "{}_Stand_in_v1_genomic.fna".format(fetched)}]}
    (folder / "provenance.json").write_text(json.dumps(provenance))
    return fasta


class StandIn(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.ref = write_assembly(self.root, "MZ_REF", "GCA_000000001.1", "GCF_000000001.1", REFERENCE)
        self.query = write_assembly(self.root, "AC_Q1", "GCA_000000002.1", "GCA_000000002.1", QUERY)
        write_assembly(self.root, "YH_LG10Y", "h2tg000014l", "h2tg000014l",
                       [("h2tg000014l", "ACGTNNACGT", None, None, None)], kind="supplementary",
                       original="YH_Y.fasta")
        self.panel_chromosomes = self.root / "panel_chromosomes_P1.tsv"
        self.write_panel_chromosomes("NC_2.1", "OZ_2.1")
        self.scratch = self.root / "scratch"
        self.now = datetime.datetime(2026, 1, 2, 3, 4, 5, tzinfo=datetime.timezone.utc)

    def write_panel_chromosomes(self, ref_lg2, query_lg2):
        header = "panel_version\tassembly_id\tlabel_style\tLG2\tLG9\tLG10\tLG11\tLG13\tLG20\n"
        rows = "PX\tMZ_REF\tLG\t{}\t\t\t\t\t\nPX\tAC_Q1\tnum\t{}\t\t\t\t\t\n".format(ref_lg2, query_lg2)
        self.panel_chromosomes.write_text(header + rows)

    def deps(self):
        rows = [{"id": "MZ_REF", "accession_genbank": "GCA_000000001.1"},
                {"id": "AC_Q1", "accession_genbank": "GCA_000000002.1"}]
        return ci.Deps(lambda a, acc: self.root / "Genomes" / a / acc / "genome.fna.gz",
                       rows, "PX", "MZ_REF", self.panel_chromosomes)

    def run_main(self, *extra):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err), \
                mock.patch.object(ci, "EXPECTED_LINKAGE_GROUPS", 2):
            code = ci.main(["--scratch", str(self.scratch), "--max-jobs", "1"] + list(extra),
                           deps=self.deps(), now=self.now)
        return code, out.getvalue(), err.getvalue()


class Helpers(StandIn):
    def test_read_fasta_stats(self):
        stats = ci.read_fasta_stats(self.ref)
        self.assertEqual(stats[0], ("NC_1.1", 14, 2, 4))      # ACGTACGTNNacgt: 2 N, 4 lower-case
        self.assertEqual(stats[2], ("NW_9.1", 8, 4, 4))
        self.assertEqual([s[0] for s in stats], ["NC_1.1", "NC_2.1", "NW_9.1", "NC_MT.1"])

    def test_label_style(self):
        self.assertEqual(ci.label_style(["LG1", "LG2", "MT"]), "LG")
        self.assertEqual(ci.label_style(["1", "2"]), "num")
        self.assertEqual(ci.label_style(["LG1", "2"]), "mixed")
        self.assertEqual(ci.label_style([]), "none")

    def test_assembly_name(self):
        self.assertEqual(ci.assembly_name("GCF_041146795.1_Mzebra_GT3a_genomic.fna", "GCF_041146795.1"), "Mzebra_GT3a")
        self.assertEqual(ci.assembly_name("YH_Y.fasta", "h2tg000014l"), "")

    def test_reference_lists(self):
        result = ci.check_one("MZ_REF", "GCA_000000001.1", self.ref, "reference")
        groups, others = ci.reference_lists(result["sequences"])
        self.assertEqual([r["name"] for r in groups], ["NC_1.1", "NC_2.1"])
        self.assertEqual([r["name"] for r in others], ["NW_9.1", "NC_MT.1"])


class Checks(StandIn):
    def test_a_good_assembly_passes(self):
        result = ci.check_one("MZ_REF", "GCA_000000001.1", self.ref, "reference")
        self.assertEqual(result["problems"], [])
        summary = result["summary"]
        self.assertEqual((summary["sha256_ok"], summary["names_lengths_ok"]), ("yes", "yes"))
        self.assertEqual((summary["label_style"], summary["chromosome_sequences"]), ("LG", 2))
        self.assertEqual(summary["assembly_name"], "Stand_in_v1")
        self.assertEqual(summary["accession_fetched"], "GCF_000000001.1")

    def test_a_changed_file_is_caught(self):
        self.query.write_bytes(gzip.compress(b">OZ_1.1\nACGT\n>OZ_2.1\nACGTACGTACGTAC\n"))
        result = ci.check_one("AC_Q1", "GCA_000000002.1", self.query, "panel")
        self.assertEqual(result["summary"]["sha256_ok"], "no")
        self.assertEqual(result["summary"]["names_lengths_ok"], "no")
        self.assertTrue(any("OZ_1.1: FASTA 4 bp" in p for p in result["problems"]))

    def test_a_wrong_accession_is_caught(self):
        result = ci.check_one("AC_Q1", "GCA_999999999.1", self.query, "panel")
        self.assertTrue(any("requested as" in p for p in result["problems"]))

    def test_the_supplementary_haplotig_has_no_labels(self):
        path = self.root / "Genomes" / "YH_LG10Y" / "h2tg000014l" / "genome.fna.gz"
        result = ci.check_one("YH_LG10Y", "h2tg000014l", path, "supplementary")
        self.assertEqual(result["problems"], [])
        self.assertEqual(result["summary"]["label_style"], "none")

    def test_a_label_mismatch_is_flagged(self):
        self.write_panel_chromosomes("NC_2.1", "OZ_1.1")       # OZ_1.1 is labelled 1, not 2
        results = {a: ci.check_one(a, acc, self.root / "Genomes" / a / acc / "genome.fna.gz", "panel")
                   for a, acc in (("MZ_REF", "GCA_000000001.1"), ("AC_Q1", "GCA_000000002.1"))}
        import csv
        with open(self.panel_chromosomes) as handle:
            rows = ci.check_panel_chromosomes(list(csv.DictReader(handle, delimiter="\t")), results)
        self.assertEqual([r["status"] for r in rows], ["ok", "MISMATCH"])


class Runs(StandIn):
    def test_dry_run_writes_nothing(self):
        code, out, _ = self.run_main()
        self.assertEqual(code, 0)
        self.assertIn("dry run", out)
        self.assertIn("YH_LG10Y", out)
        self.assertFalse(self.scratch.exists())

    def test_go_writes_every_output(self):
        code, out, _ = self.run_main("--go")
        self.assertEqual(code, 0, out)
        folder = self.scratch / "inputs" / "20260102-030405"
        for name in ("summary.tsv", "report.md", "panel_chromosomes_check.tsv",
                     "reference_refContigs.tsv", "reference_otherContigs.tsv",
                     "sequences/MZ_REF.tsv", "sequences/AC_Q1.tsv", "sequences/YH_LG10Y.tsv"):
            self.assertTrue((folder / name).is_file(), name)
        self.assertEqual(len((folder / "reference_refContigs.tsv").read_text().splitlines()), 3)
        self.assertIn("num (1): AC_Q1", (folder / "report.md").read_text())

    def test_never_overwrites(self):
        self.run_main("--go")
        code, _, err = self.run_main("--go")
        self.assertEqual(code, 2)
        self.assertIn("never overwrites", err)

    def test_without_the_haplotig(self):
        code, out, _ = self.run_main("--go", "--no-supplementary")
        self.assertEqual(code, 0, out)
        self.assertFalse((self.scratch / "inputs" / "20260102-030405" / "sequences" / "YH_LG10Y.tsv").exists())

    def test_a_problem_gives_exit_status_1(self):
        self.query.write_bytes(gzip.compress(b">OZ_1.1\nACGT\n"))
        code, out, _ = self.run_main("--go")
        self.assertEqual(code, 1)
        self.assertIn("AC_Q1", out)


if __name__ == "__main__":
    unittest.main()
