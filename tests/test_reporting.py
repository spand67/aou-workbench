from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import pandas as pd

from aou_workbench.reporting import dataframe_markdown, load_table_if_exists


class ReportingTests(unittest.TestCase):
    def test_dataframe_markdown_falls_back_without_tabulate(self) -> None:
        frame = pd.DataFrame([{"variant_id": "1-100-A-G", "fisher_p": 0.05}])

        with unittest.mock.patch.object(pd.DataFrame, "to_markdown", side_effect=ImportError):
            rendered = dataframe_markdown(frame)

        self.assertIn("1-100-A-G", rendered)
        self.assertIn("fisher_p", rendered)

    def test_load_table_if_exists_reads_compressed_tsv(self) -> None:
        root = Path(tempfile.mkdtemp(prefix="aou-workbench-reporting-"))
        path = root / "results.tsv.gz"
        pd.DataFrame([{"variant_id": "1-100-A-G", "fisher_p": 0.05}]).to_csv(path, sep="\t", index=False)

        loaded = load_table_if_exists(str(path))

        self.assertEqual(loaded.shape, (1, 2))
        self.assertEqual(loaded.loc[0, "variant_id"], "1-100-A-G")


if __name__ == "__main__":
    unittest.main()
