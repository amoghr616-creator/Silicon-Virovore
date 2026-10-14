import tempfile
import unittest
from pathlib import Path

from prototype2.run_p2_rapidock import (
    DEFAULT_RECEPTOR,
    _backend_preflight,
    _build_dock_command,
    _preflight_report_path,
)


class P2PreflightTests(unittest.TestCase):
    def test_missing_backend_and_long_checkpoint_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            result = _backend_preflight(
                "MKQAVNALLVFFAGSSDAIRR",
                52.0,
                20,
                receptor=DEFAULT_RECEPTOR,
                output_root=Path(directory) / "new_run",
            )
        self.assertFalse(result["ready"])
        self.assertTrue(any("not available on PATH" in error for error in result["errors"]))
        self.assertTrue(any("longer_local.pt" in error for error in result["errors"]))
        self.assertTrue(result["receptor_mapping_valid"])
        self.assertTrue(result["coordinate_frame"]["coordinate_frame_comparable"])

    def test_nonempty_output_directory_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "existing"
            output.mkdir()
            (output / "old-result.json").write_text("preserve")
            result = _backend_preflight(
                "MKQAVNALLVFFAGSSDAIRR",
                52.0,
                20,
                receptor=DEFAULT_RECEPTOR,
                output_root=output,
            )
        self.assertTrue(any("non-empty" in error for error in result["errors"]))

    def test_seed_is_passed_to_backend_and_reports_stay_outside_run_dir(self):
        command = _build_dock_command(
            "P2-01",
            "MKQAVNALLVFFAGSSDAIRR",
            DEFAULT_RECEPTOR,
            Path("/tmp/p2-pilot/P2-01"),
            52.0,
            20,
            616,
        )
        self.assertEqual(command[-2:], ["--seed", "616"])
        self.assertEqual(
            _preflight_report_path(Path("/tmp/p2-pilot")),
            Path("/tmp/p2-pilot_preflight.json"),
        )


if __name__ == "__main__":
    unittest.main()
