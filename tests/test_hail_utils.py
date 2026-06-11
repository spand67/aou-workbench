from __future__ import annotations

from types import SimpleNamespace
import unittest
from unittest import mock

from aou_workbench import hail_utils


class HailUtilsTests(unittest.TestCase):
    def test_resolve_requester_pays_project_prefers_explicit_value(self) -> None:
        with mock.patch.dict("os.environ", {"GCS_REQUESTER_PAYS_PROJECT": "env-project"}, clear=True):
            self.assertEqual(hail_utils.resolve_requester_pays_project("explicit-project"), "explicit-project")

    def test_resolve_requester_pays_project_uses_gcloud_fallback(self) -> None:
        with (
            mock.patch.dict("os.environ", {}, clear=True),
            mock.patch(
                "aou_workbench.hail_utils.subprocess.run",
                return_value=SimpleNamespace(returncode=0, stdout="workspace-project\n"),
            ) as run,
        ):
            self.assertEqual(hail_utils.resolve_requester_pays_project(), "workspace-project")

        run.assert_called_once_with(
            ["gcloud", "config", "get-value", "project"],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_configure_hail_bootstrap_adds_requester_pays_spark_conf(self) -> None:
        calls: list[tuple[list[str], dict[str, object]]] = []

        def fake_run(cmd: list[str], **kwargs: object) -> SimpleNamespace:
            calls.append((cmd, kwargs))
            if cmd[:4] == ["gcloud", "config", "get-value", "project"]:
                return SimpleNamespace(returncode=0, stdout="workspace-project\n")
            return SimpleNamespace(returncode=0, stdout="")

        with (
            mock.patch.dict("os.environ", {}, clear=True),
            mock.patch("aou_workbench.hail_utils.subprocess.run", side_effect=fake_run),
        ):
            hail_utils.configure_aou_hail_bootstrap(None, ["gs://vwb-aou-datasets-controlled"])
            args = hail_utils.os.environ["PYSPARK_SUBMIT_ARGS"]
            self.assertIn("spark.serializer=org.apache.spark.serializer.KryoSerializer", args)
            self.assertIn("spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator", args)
            self.assertIn("spark.hadoop.fs.gs.requester.pays.project.id=workspace-project", args)
            self.assertIn("spark.hadoop.fs.gs.requester.pays.mode=CUSTOM", args)
            self.assertIn("spark.hadoop.fs.gs.requester.pays.buckets=vwb-aou-datasets-controlled", args)
            self.assertEqual(hail_utils.os.environ["GCS_REQUESTER_PAYS_PROJECT"], "workspace-project")

        self.assertIn(
            (["hailctl", "config", "set", "gcs_requester_pays/project", "workspace-project"], mock.ANY),
            calls,
        )
        self.assertIn(
            (
                [
                    "hailctl",
                    "config",
                    "set",
                    "gcs_requester_pays/buckets",
                    "vwb-aou-datasets-controlled",
                ],
                mock.ANY,
            ),
            calls,
        )

    def test_active_spark_requester_pays_conf_uses_hadoop_keys(self) -> None:
        seen: dict[str, str] = {}

        class FakeHadoopConf:
            def set(self, key: str, value: str) -> None:
                seen[key] = value

        class FakeJsc:
            def hadoopConfiguration(self) -> FakeHadoopConf:
                return FakeHadoopConf()

        class FakeSc:
            _jsc = FakeJsc()

        fake_hail = SimpleNamespace(current_backend=lambda: SimpleNamespace(sc=FakeSc()))

        hail_utils._set_active_spark_requester_pays_conf(  # noqa: SLF001
            fake_hail,
            "workspace-project",
            ["gs://vwb-aou-datasets-controlled"],
        )

        self.assertEqual(seen["fs.gs.requester.pays.project.id"], "workspace-project")
        self.assertEqual(seen["fs.gs.requester.pays.mode"], "CUSTOM")
        self.assertEqual(seen["fs.gs.requester.pays.buckets"], "vwb-aou-datasets-controlled")


if __name__ == "__main__":
    unittest.main()
