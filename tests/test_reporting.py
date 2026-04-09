from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import pandas as pd

from aou_workbench.reporting import load_table_if_exists


class ReportingTests(unittest.TestCase):
    def test_load_table_if_exists_reads_compressed_tsv(self) -> None:
        root = Path(tempfile.mkdtemp(prefix="aou-workbench-reporting-"))
        path = root / "results.tsv.gz"
        pd.DataFrame([{"variant_id": "1-100-A-G", "fisher_p": 0.05}]).to_csv(path, sep="\t", index=False)

        loaded = load_table_if_exists(str(path))

        self.assertEqual(loaded.shape, (1, 2))
        self.assertEqual(loaded.loc[0, "variant_id"], "1-100-A-G")


if __name__ == "__main__":
    unittest.main()
