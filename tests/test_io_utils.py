from __future__ import annotations

import unittest

from aou_workbench.io_utils import stable_hash


class IoUtilsTests(unittest.TestCase):
    def test_stable_hash_changes_when_payload_changes_beyond_prefix(self) -> None:
        first = stable_hash({"vat_path": "same", "max_af": 0.001})
        second = stable_hash({"vat_path": "same", "max_af": 0.01})

        self.assertEqual(len(first), 16)
        self.assertNotEqual(first, second)
        self.assertEqual(first, stable_hash({"max_af": 0.001, "vat_path": "same"}))


if __name__ == "__main__":
    unittest.main()
