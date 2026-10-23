"""Adversarial tests for the eight additional replays: each replay has a --perturb mode that changes exactly one input of the
computation; the replay must then exit 1 and print the corresponding FAIL line.

  python -m unittest discover -s tests -v                 # all
  python -m unittest tests.test_replays_adversarial -v    # only these
"""

import pathlib
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
REPLAYS = ROOT / "replays"

# slug, script, text of the claim that must FAIL, what is perturbed
CASES = [
    (
        "e059_zariski_multiplicity",
        "replay_zariski.py",
        "mu(d0) = 51 for the perturbed curve",
        "q1: w^8 -> w^7",
    ),
    (
        "e067_campana_peternell",
        "replay_campana_peternell.py",
        "N (B_i H_j) = 0 for the perturbed witness",
        "witness entry Ha[4] + 1",
    ),
    (
        "e101_simplex_isotropic",
        "replay_simplex.py",
        "(2.5)-(2.6) step bound m0 on U1 raised",
        "m0(U1) .68 -> .75",
    ),
    (
        "e138_subset_sum",
        "replay_subset_sum.py",
        "Table 1 (nodes 0..32) reproduced",
        "sign of r6 flipped",
    ),
    (
        "e153_bernoulli_non_pisot",
        "replay_bernoulli.py",
        "root-disk certificate",
        "one polynomial coefficient + 1",
    ),
    (
        "e268_haldane_gap",
        "replay_haldane.py",
        "(C.8) at N = 72",
        "MPS diagonal parameter D_1 + 1",
    ),
    (
        "e318_chromatic_splitting",
        "replay_chromatic.py",
        "det d3[R3,C3] = 1024",
        "row index 61 -> 62",
    ),
    (
        "e352_calabi_flow",
        "replay_calabi.py",
        "(A.1) |q| <= .25 R at the centre",
        "k shifted by 1e-9",
    ),
]


def run(slug, script, *args):
    return subprocess.run(
        [sys.executable, str(REPLAYS / slug / script), *args],
        capture_output=True,
        text=True,
        cwd=REPLAYS / slug,
    )


class Perturbed(unittest.TestCase):
    pass


def make(slug, script, line, what):
    def test(self):
        p = run(slug, script, "--perturb")
        self.assertEqual(p.returncode, 1, p.stdout[-600:] + p.stderr[-600:])
        self.assertIn(
            line,
            p.stdout,
            f"perturbation ({what}) did not produce the expected FAIL line",
        )
        self.assertIn("RESULT: FAILURES", p.stdout)

    return test


for slug, script, line, what in CASES:
    setattr(
        Perturbed, f"test_{slug}_perturbed_input_fails", make(slug, script, line, what)
    )


if __name__ == "__main__":
    unittest.main()
