#!/usr/bin/env python3
"""Independent replay of the finite computations in 'The periodic spin-one Haldane gap' (OpenAI, 2026-09-24).

Written from the manuscript text only; no code from the authors is used.  Three parts:
  haldane_scalar.py  exact rational certificates of Sections 4-5 (the Lemma 5.1 tables are inputs there);
  haldane_mps.py     Appendix C: exact integer matrix-product trial vectors at N = 72 and 120 (E0(N) < -N a);
  haldane_traces.py  Lemma 5.1: the thermal tables Z_n(b,g) re-computed by exact diagonalisation in float64 and compared with the
                     printed centres within the printed errors (10e-6 at b = 21/2, 30e-6 at b = 49/4).  NOT a rigorous enclosure.

  python replay_haldane.py             # prints every claim with PASS/FAIL; exit 0 only if all pass (about 6-8 min)
  python replay_haldane.py --perturb   # adversarial: the MPS diagonal parameter D_1 changed by 1; the (C.8) interval at N = 72 must FAIL
  python replay_haldane.py --selftest  # also runs the built-in wrong-input controls (each must be rejected)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import importlib.util
import pathlib
import re
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
RES = []
_FLT = re.compile(r"-?\d+\.\d{7,}(?:e[+-]?\d+)?|-?\d+(?:\.\d+)?e[+-]?\d+")


def _norm(s):
    return _FLT.sub(lambda m: f"{float(m.group(0)):.6g}", s)


def chk(name, ok):
    if name.startswith("WRONG INPUT") and not SELFTEST:
        return
    RES.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] H{len(RES):02d}: {_norm(name)}", flush=True)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def finish():
    npass = sum(1 for _, o in RES if o)
    ok = npass == len(RES)
    print(
        f"\n{npass}/{len(RES)} claims passed; RESULT: {'ALL PASS' if ok else 'FAILURES'}"
    )
    return 0 if ok else 1


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.chk = chk
    return mod


def traces_part(T):
    thetas = {"1": 0.0, "P": np.pi, "C": 2 * np.pi / 3}
    # (a) the gauge/momentum reduction reproduces literal diagonalisation of (2.5) on small rings
    worst = 0.0
    for n in (4, 5, 6):
        for th in thetas.values():
            dsp = np.sort(T.direct_spectrum(n, th))
            red = T.reduced_spectrum(n, th, verbose=False)
            allr = np.sort(
                np.concatenate(
                    [np.concatenate([v] * (1 if m == 0 else 2)) for m, v in red.items()]
                )
            )
            assert len(dsp) == len(allr)
            worst = max(worst, float(np.abs(dsp - allr).max()))
    chk(
        f"reduced (gauge + momentum) spectra equal literal diagonalisation of (2.5) for n = 4, 5, 6 and theta = 0, pi, 2 pi/3 (max difference < 1e-10: {worst < 1e-10})",
        worst < 1e-10,
    )
    res = {}
    CUT = 2.5
    # discarded-eigenvalue bound for the Lanczos runs (n >= 10): at most 3^n eigenvalues above the cut, each contributing < exp(-b CUT)
    worst_bound = max(3**n * float(np.exp(-b * CUT)) for n, b in [(10, 21 / 2), (10, 49 / 4), (11, 49 / 4), (12, 49 / 4)])
    chk(f"Lanczos cut E + a n < {CUT}: the discarded eigenvalues change Z by less than {worst_bound:.1e} (3^n e^(-b cut) for the Lanczos cases n = 10, b = 21/2 and n = 10..12, b = 49/4), below 1e-6 (units of the tables)", worst_bound < 1e-6)
    todo = [(n, g) for n in range(4, 13) for g in "1P"] + [(n, "C") for n in (6, 8, 10)]
    for n, g in todo:
        t0 = time.time()
        if n <= 9:
            out = T.reduced_spectrum(n, thetas[g], verbose=False)
        else:
            out = T.reduced_spectrum(
                n, thetas[g], cutoff=CUT - T.A_SHIFT * n, verbose=False, lanczos=True
            )
        res[(n, g)] = (T.partition(out, n, 21 / 2), T.partition(out, n, 49 / 4))
        log(f"   n={n} g={g} ({time.time() - t0:.0f}s)")
    for n, row in T.TABLE_21_2.items():
        for g, c in zip("1PC", row):
            z = res[(n, g)][0] * 1e6
            chk(
                f"Lemma 5.1, b = 21/2, n = {n}, g = {g}: replayed {round(z)} vs table {c} (|difference| < 10 in units 1e-6)",
                abs(z - c) < 10,
            )
    for n, row in T.TABLE_49_4.items():
        for g, c in zip("1P", row):
            z = res[(n, g)][1] * 1e6
            chk(
                f"Lemma 5.1, b = 49/4, n = {n}, g = {g}: replayed {round(z)} vs table {c} (|difference| < 30 in units 1e-6)",
                abs(z - c) < 30,
            )
    z = res[(10, "1")][0] * 1e6
    chk(
        "WRONG INPUT: table entry n = 10, g = I (b = 21/2) shifted by +50 fails the 10e-6 tolerance",
        not (abs(z - (2327547 + 50)) < 10),
    )
    z = res[(8, "P")][1] * 1e6
    chk(
        "WRONG INPUT: table entry n = 8, g = P (b = 49/4) shifted by +100 fails the 30e-6 tolerance",
        not (abs(z - (510140 + 100)) < 30),
    )


def run():
    t0 = time.time()
    S = load("haldane_scalar")
    S.run(SELFTEST)
    M = load("haldane_mps")
    if PERTURB:  # perturb the INPUT (MPS diagonal parameter D_1 changed by 1), then run the unchanged main checks
        M.DDIAG[0] += 1
    M.run(SELFTEST)
    T = load("haldane_traces")
    traces_part(T)
    log(f"elapsed {time.time() - t0:.1f}s")
    return finish()


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    sys.exit(run())
