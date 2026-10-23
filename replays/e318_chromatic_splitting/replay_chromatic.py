"""Independent replay of the finite Chevalley-Eilenberg certificate in 'A rational obstruction to strong chromatic splitting at height three'
(OpenAI, 2026-09-25), Appendix A.1.

Written from the manuscript text only; no code from the authors (including the accompanying certificate script) is used.
Exact integer arithmetic (python-flint).  Claims: the matrix of d1 for sl2 (determinant 4) and d2 = 0; dims of C^0..C^4 for sl3 = 1, 8, 28, 56, 70;
the three nonsingular minors det d1[R1,C1] = -1, det d2[R2,C2] = 2916, det d3[R3,C3] = 1024; ranks 8, 20, 35 of the full matrices;
H1 = H2 = 0, H3 = Q.[omega]; the trace cocycle is closed, omega(E12,E21,E11-E22) = 2, and it is not a coboundary.

  python replay_chromatic.py             # prints every claim with PASS/FAIL; exit 0 only if all pass
  python replay_chromatic.py --perturb   # adversarial: one row index of R3 changed (61 -> 62); the printed minor must FAIL (exit 1)
  python replay_chromatic.py --selftest  # also runs the built-in wrong-input controls (each must be rejected)

Deterministic stdout (no timing).
"""

import itertools
import re
import sys

from flint import fmpz_mat

RES = []


SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
PFX = "CH"
_FLT = re.compile(r"-?\d+\.\d{7,}(?:e[+-]?\d+)?|-?\d+(?:\.\d+)?e[+-]?\d+")


def _norm(s):
    return _FLT.sub(lambda m: f"{float(m.group(0)):.6g}", s)


def chk(name, ok):
    if name.startswith("WRONG INPUT") and not SELFTEST:
        return
    RES.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {PFX}{len(RES):02d}: {_norm(name)}", flush=True)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def finish():
    npass = sum(1 for _, o in RES if o)
    ok = npass == len(RES)
    print(f"\n{npass}/{len(RES)} claims passed; RESULT: {'ALL PASS' if ok else 'FAILURES'}")
    return 0 if ok else 1


def sl_basis(n):
    """matrix units E_ij (i != j) lexicographic, then H_i = E_ii - E_nn for 1 <= i < n.  Matrices as dict {(r,c): v} 0-indexed."""
    B = []
    for i in range(n):
        for j in range(n):
            if i != j:
                B.append({(i, j): 1})
    for i in range(n - 1):
        B.append({(i, i): 1, (n - 1, n - 1): -1})
    return B


def mat_mul(A, B):
    out = {}
    for (i, k), a in A.items():
        for (k2, j), b in B.items():
            if k == k2:
                out[(i, j)] = out.get((i, j), 0) + a * b
    return {k: v for k, v in out.items() if v}


def bracket(A, B):
    AB, BA = mat_mul(A, B), mat_mul(B, A)
    out = dict(AB)
    for k, v in BA.items():
        out[k] = out.get(k, 0) - v
    return {k: v for k, v in out.items() if v}


def coords(M, n):
    """expand a traceless matrix in the basis (E_ij, H_i): off-diagonal entries give E_ij; diagonal d_i gives sum_i d_i H_i (H_i = E_ii - E_nn)"""
    c = {}
    idx = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                if (i, j) in M:
                    c[idx] = M[(i, j)]
                idx += 1
    for i in range(n - 1):
        if (i, i) in M:
            c[idx + i] = M[(i, i)]
    # trace check
    tr = sum(M.get((i, i), 0) for i in range(n))
    assert tr == 0
    return c


def perm_sign(seq):
    seq = list(seq)
    sgn = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                sgn = -sgn
    return sgn


def ce_matrices(n):
    B = sl_basis(n)
    dim = len(B)
    br = {}
    for a in range(dim):
        for b in range(dim):
            br[(a, b)] = coords(bracket(B[a], B[b]), n)
    tuples = {k: list(itertools.combinations(range(dim), k)) for k in range(dim + 1)}
    tidx = {k: {t: i for i, t in enumerate(tuples[k])} for k in tuples}
    mats = {}
    for k in range(5 if n == 3 else 3):
        if k + 1 > dim:
            continue
        rows = tuples[k + 1]
        cols = tuples[k]
        M = [[0] * len(cols) for _ in rows]
        for r, T in enumerate(rows):
            # (d omega)(X_T) = sum_{i<j} (-1)^{i+j} omega([X_i, X_j], X_T without i, j)
            for i in range(k + 1):
                for j in range(i + 1, k + 1):
                    rest = [T[x] for x in range(k + 1) if x not in (i, j)]
                    for m, cm in br[(T[i], T[j])].items():
                        if m in rest:
                            continue
                        tup = [m] + rest
                        S = tuple(sorted(tup))
                        sgn = perm_sign(tup)
                        M[r][tidx[k][S]] += (-1) ** (i + j) * cm * sgn
        mats[k] = M
    return mats, dim, tuples


def minor(M, R, C):
    return fmpz_mat([[M[r][c] for c in C] for r in R]).det()


def rank(M):
    return fmpz_mat(M).rank()


def run():
    selftest = SELFTEST
    # ---------------- sl2
    m2, dim2, tup2 = ce_matrices(2)
    d1 = m2[1]
    chk(
        f"sl2: basis (E12, E21, H1), dim {dim2}; d1 is {len(d1)} x {len(d1[0])}",
        dim2 == 3 and len(d1) == 3 and len(d1[0]) == 3,
    )
    detd1 = fmpz_mat(d1).det()
    chk(
        f"sl2: det d1 = {detd1} (stated 4); d1 = {d1} vs printed [[0,0,-1],[2,0,0],[0,-2,0]]",
        abs(int(detd1)) == 4
        and [[int(x) for x in row] for row in d1] == [[0, 0, -1], [2, 0, 0], [0, -2, 0]],  # entry by entry, signs and order
    )
    chk(
        "sl2: d2 = 0 (C^2 -> C^3 = Lambda^3 of a 3-dimensional space, d2 has all zero entries)",
        all(x == 0 for row in m2[2] for x in row),
    )
    # ---------------- sl3
    mats, dim, tuples = ce_matrices(3)
    dims = [len(tuples[k]) for k in range(5)]
    chk(
        f"sl3: dims of C^0..C^4 = {dims} (stated 1, 8, 28, 56, 70)",
        dims == [1, 8, 28, 56, 70],
    )
    # d o d = 0 (sanity of the construction)
    ok = True
    for k in range(3):
        P = fmpz_mat(mats[k + 1]) * fmpz_mat(mats[k])
        ok &= all(P[i, j] == 0 for i in range(P.nrows()) for j in range(P.ncols()))
    chk("sanity: d_{k+1} d_k = 0 for k = 0..2 on the constructed matrices", ok)
    R1 = (1, 2, 3, 5, 7, 9, 15, 16)
    C1 = (0, 1, 2, 3, 4, 5, 6, 7)
    R2 = (0, 2, 4, 6, 7, 8, 9, 11, 12, 13, 19, 22, 24, 28, 30, 37, 39, 41, 43, 52)
    C2 = (0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 14, 15, 16, 17, 20, 22, 23, 25, 27)
    R3 = (
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        9,
        10,
        11,
        12,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        25,
        26,
        27,
        31,
        35,
        36,
        37,
        39,
        40,
        41,
        45,
        55,
        56,
        61,
    )
    C3 = (
        0,
        1,
        2,
        3,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        20,
        21,
        22,
        24,
        25,
        28,
        31,
        32,
        36,
        37,
        38,
        40,
        42,
        43,
        44,
        45,
        52,
    )
    if PERTURB:  # perturb the INPUT (last row index of R3: 61 -> 62), then run the unchanged main checks
        R3 = tuple(list(R3[:-1]) + [62])
    chk(
        f"index sets: |R1|=|C1|=8, |R2|=|C2|=20, |R3|=|C3|={len(R3)}=35",
        len(R1) == len(C1) == 8
        and len(R2) == len(C2) == 20
        and len(R3) == len(C3) == 35,
    )
    m1 = int(minor(mats[1], R1, C1))
    m2_ = int(minor(mats[2], R2, C2))
    m3_ = int(minor(mats[3], R3, C3))
    log(f"   minors: d1 {m1}, d2 {m2_}, d3 {m3_}")
    chk(f"det d1[R1,C1] = {m1} (stated -1)", m1 == -1)
    chk(
        f"|det d2[R2,C2]| = {abs(m2_)} (stated 2916), signed value {m2_}",
        abs(m2_) == 2916,
    )
    chk(
        f"|det d3[R3,C3]| = {abs(m3_)} (stated 1024), signed value {m3_}",
        abs(m3_) == 1024,
    )
    chk(
        f"sign conventions: signed minors ({m1}, {m2_}, {m3_}) vs stated (-1, 2916, 1024)",
        (m1, m2_, m3_) == (-1, 2916, 1024),
    )
    r1, r2, r3 = rank(mats[1]), rank(mats[2]), rank(mats[3])
    chk(
        f"full ranks: rank d1 = {r1} (8), rank d2 = {r2} (20), rank d3 = {r3} (35)",
        (r1, r2, r3) == (8, 20, 35),
    )
    h1 = dims[1] - r1 - 0
    h2 = dims[2] - r2 - r1
    h3 = dims[3] - r3 - r2
    chk(
        f"H^1 = {h1}, H^2 = {h2}, H^3 = {h3} (stated 0, 0, 1)",
        (h1, h2, h3) == (0, 0, 1),
    )
    # trace cocycle omega(X,Y,Z) = tr(X [Y,Z]) as a cochain in C^3 (coordinates on the 56 basis triples)
    B = sl_basis(3)

    def trace_prod(A, Bm):
        P = mat_mul(A, Bm)
        return sum(P.get((i, i), 0) for i in range(3))

    omega = []
    for T in tuples[3]:
        X, Y, Z = (B[t] for t in T)
        omega.append(trace_prod(X, bracket(Y, Z)))
    # closedness: d3 omega = 0
    d3 = fmpz_mat(mats[3])
    ov = fmpz_mat([[w] for w in omega])
    prod = d3 * ov
    chk(
        "trace cocycle omega(X,Y,Z) = tr(X[Y,Z]) is closed: d3 omega = 0",
        all(prod[i, 0] == 0 for i in range(prod.nrows())),
    )
    # not a coboundary: rank [im d2 | omega] = 21
    aug = [
        row + [w]
        for row, w in zip(
            [list(map(int, [mats[2][i][j] for j in range(28)])) for i in range(56)],
            omega,
        )
    ]
    chk(
        f"omega is not in im d2: rank [d2 | omega] = {rank(aug)} (= 21 = 20 + 1)",
        rank(aug) == 21,
    )
    # sl2 block: omega(E12, E21, E11 - E22) = 2
    E12, E21, H = {(0, 1): 1}, {(1, 0): 1}, {(0, 0): 1, (1, 1): -1}
    val = trace_prod(E12, bracket(E21, H))
    chk(f"omega(E12, E21, E11-E22) = {val} (stated 2)", val == 2)
    if selftest:
        # deliberately wrong claims must be detected
        chk(
            "WRONG INPUT: claimed det d2 minor 2917 is rejected", not (abs(m2_) == 2917)
        )
        chk("WRONG INPUT: claimed rank d3 = 34 is rejected", r3 != 34)
        chk("WRONG INPUT: claimed H^3 dimension 2 is rejected", h3 != 2)
        R3bad = tuple(list(R3[:-1]) + [62])
        mb = int(minor(mats[3], R3bad, C3))
        chk(
            f"WRONG INPUT: replacing the last row index 61 by 62 changes the minor to {mb} != 1024",
            abs(mb) != 1024,
        )
        d1_bad = [row[:] for row in mats[1]]
        d1_bad[1][1] += 1
        P = fmpz_mat(mats[2]) * fmpz_mat(d1_bad)
        chk(
            "WRONG INPUT: perturbed d1 breaks d2 d1 = 0 (the checker's sanity test catches it)",
            any(P[i, j] != 0 for i in range(P.nrows()) for j in range(P.ncols())),
        )
    return finish()


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    sys.exit(run())
