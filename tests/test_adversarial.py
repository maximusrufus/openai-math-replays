"""Adversarial tests: the replay must FAIL when an input is perturbed, and must match the expected output when not.

  python -m unittest discover -s tests -v
"""

import contextlib
import importlib.util
import io
import pathlib
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DIR = ROOT / "replays" / "wall_d2"
SCRIPT = DIR / "replay_wall.py"


def run(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True
    )


def load():
    spec = importlib.util.spec_from_file_location("replay_wall", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class Replay(unittest.TestCase):
    def test_baseline_matches_expected_output(self):
        p = run()
        self.assertEqual(p.returncode, 0, p.stdout[-500:])
        expected = (DIR / "expected_output.txt").read_text(encoding="utf-8")
        self.assertEqual(p.stdout.replace("\r\n", "\n"), expected)
        self.assertIn("21/21 claims passed", p.stdout)

    def test_perturbed_matrix_entry_fails(self):
        p = run("--perturb-matrix")
        self.assertEqual(p.returncode, 1)
        self.assertIn("[FAIL] C5:", p.stdout)
        self.assertIn("RESULT: FAILURES", p.stdout)

    def test_perturbed_relator_fails(self):
        mod = load()
        mod.RELS["r3"] = mod.word("x1 a2^-12")  # was a2^-13
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            rc = mod.main()
        self.assertEqual(rc, 1)
        self.assertIn("[FAIL] C1:", buf.getvalue())

    def test_perturbed_character_fails(self):
        mod = load()
        orig = mod.Qchar
        mod.Qchar = lambda: {**orig(), "a2": mod.Q(1)}  # was -1
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            rc = mod.main()
        self.assertEqual(rc, 1)
        self.assertIn("[FAIL]", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
