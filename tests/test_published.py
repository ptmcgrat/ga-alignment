"""The published tables the checks read must be in the repository, unchanged.

They live in published/, not data/: the repository's .gitignore keeps data/ out (data lives under the data root),
which is how alignment_007's copies were left out of git. These tests fail if either table is missing or differs
from what was recorded, so a pull that lacks them fails the test suite instead of a later run.
"""

import hashlib
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import pafkit  # noqa: E402

BREAKPOINTS = REPO / "published" / "kumar2025_suppfile4_breakpoints.tsv"
EXPECTED = REPO / "published" / "expected_inversions.tsv"
SHA256 = {BREAKPOINTS.name: "edc9a2cddd0e34446ac896dd870a89fb7271d76445341ddb2361f4522280cff2",
          EXPECTED.name: "ba046c393445b1c19a94d9ae4d6b5c359137b465f2810f9907b7b9f2fc8a81d4"}
SUBSET = {"YH_YH1m", "YH_YH2f", "AC_SAMEA123396441", "MZ_NMetZ", "MZ_MZ4f", "RC_006pf", "DT_SAMEA123396444"}


class PublishedTables(unittest.TestCase):
    def test_both_tables_are_present_and_unchanged(self):
        for path in (BREAKPOINTS, EXPECTED):
            self.assertTrue(path.is_file(), "missing from the repository: {}".format(path))
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), SHA256[path.name], path.name)

    def test_the_eight_units(self):
        units = pafkit.read_table(BREAKPOINTS)
        self.assertEqual([u["unit"] for u in units], ["2", "9", "10", "11a", "11b", "13", "20a", "20b"])
        for u in units:
            self.assertLessEqual(int(u["bp1_start"]), int(u["bp1_end"]))
            self.assertLessEqual(int(u["bp2_start"]), int(u["bp2_end"]))

    def test_the_expectations(self):
        rows = pafkit.read_table(EXPECTED)
        self.assertEqual(len(rows), 193)
        scored = [r for r in rows if r["scored"] == "yes"]
        self.assertEqual(len(scored), 55)
        self.assertEqual({r["assembly_id"] for r in scored}, SUBSET)
        self.assertTrue(all(r["expected"] in ("inverted", "normal") for r in scored))


if __name__ == "__main__":
    unittest.main()
