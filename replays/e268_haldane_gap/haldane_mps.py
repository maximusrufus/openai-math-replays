"""Exact integer matrix-product trial vectors of Appendix C of 'The periodic spin-one Haldane gap' (module of replay_haldane.py).
From the manuscript text only (Lemma C.1, Tables 1-2, entry rules (C.4)-(C.5), contractions (C.2)-(C.3), result (C.8)); exact big-integer
arithmetic (python-flint fmpz_mat); Lemma C.1 is also cross-checked by brute force on a 4-site ring."""

import itertools
import sys
import time
from fractions import Fraction as Fr

from flint import fmpz_mat

RES = []


def chk(name, ok):
    raise RuntimeError("chk is provided by replay_haldane.py")


LELL = [1, 1, 1, 1, 3, 3, 3, 5]  # l_alpha, alpha = 1..8
DDIAG = [-289, -84, -9, 33, -22, 27, 71, -11]
PTAB = {
    (1, 5): -11,
    (1, 6): -35,
    (1, 7): -28,
    (2, 5): -23,
    (2, 6): 16,
    (2, 7): -49,
    (3, 5): -6,
    (3, 6): -25,
    (3, 7): -20,
    (4, 5): -34,
    (4, 6): 9,
    (4, 7): -32,
    (5, 8): 1,
    (6, 8): 2,
    (7, 8): -16,
}


def labels():
    out = []
    for a in range(1, 9):
        l = LELL[a - 1]
        for d in range(-l, l + 1, 2):
            out.append((a, d))
    return out


def fvals(l, d, kind, s):
    r = (d + l) // 2
    assert (d + l) % 2 == 0
    if kind == "up":  # k = l + 2
        return {
            -1: 2 * (l + 1 - r) * (l + 2 - r),
            0: 2 * (r + 1) * (l - r + 1),
            1: (r + 1) * (r + 2),
        }[s]
    if kind == "same":
        return {-1: 2 * (l - r + 1), 0: d, 1: -(r + 1)}[s]
    if kind == "down":
        return {-1: 2, 0: -2, 1: 1}[s]


def build_A(diag=DDIAG, ptab=PTAB):
    labs = labels()
    idx = {x: i for i, x in enumerate(labs)}
    D = len(labs)
    assert D == 26
    A = {s: [[0] * D for _ in range(D)] for s in (-1, 0, 1)}
    for a, d in labs:
        l = LELL[a - 1]
        for b, e in labs:
            k = LELL[b - 1]
            for s in (-1, 0, 1):
                if e != d + 2 * s:
                    continue
                if abs(k - l) > 2:
                    continue
                if k == l and a != b:
                    continue
                if k == l:
                    p, kind = diag[a - 1], "same"
                elif k == l + 2:
                    p, kind = ptab[(a, b)], "up"
                else:
                    p, kind = -ptab[(b, a)], "down"
                A[s][idx[(a, d)]][idx[(b, e)]] = p * fvals(l, d, kind, s)
    return A


def kron(X, Y):
    n, m = len(X), len(Y)
    return [
        [X[i // m][j // m] * Y[i % m][j % m] for j in range(n * m)]
        for i in range(n * m)
    ]


def mm(X, Y):
    return fmpz_mat(X) * fmpz_mat(Y)


def tolist(M):
    return [[int(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]


def build_TO(A):
    D = 26
    # T = sum_s A_s (x) A_s
    T = None
    for s in (-1, 0, 1):
        K = fmpz_mat(kron(A[s], A[s]))
        T = K if T is None else T + K
    B = {
        (s, t): fmpz_mat(A[s]) * fmpz_mat(A[t]) for s in (-1, 0, 1) for t in (-1, 0, 1)
    }
    O = None

    def add(X):
        nonlocal O
        O = X if O is None else O + X

    for s in (-1, 0, 1):
        for t in (-1, 0, 1):
            Kt = fmpz_mat(kron(tolist(B[(s, t)]), tolist(B[(s, t)])))
            add(Kt * (s * t) if s * t != 0 else Kt * 0)
    for s in (-1, 0):
        for t in (0, 1):
            add(fmpz_mat(kron(tolist(B[(s, t)]), tolist(B[(s + 1, t - 1)]))))
            add(fmpz_mat(kron(tolist(B[(s + 1, t - 1)]), tolist(B[(s, t)]))))
    return T, O


def matpow(M, n):
    R = None
    P = M
    while n:
        if n & 1:
            R = P if R is None else R * P
        n >>= 1
        if n:
            P = P * P
    return R


def trace(M):
    return sum(M[i, i] for i in range(M.nrows()))


def contraction(A, N):
    T, O = build_TO(A)
    TN2 = matpow(T, N - 2)
    TN = TN2 * T * T
    V = trace(TN)
    # U_N = Tr(T^{N-2} O)
    U = trace(TN2 * O)
    return U, V


def lemma_c1_bruteforce(A, N):
    """<psi|H|psi> / (N) and |psi|^2 by direct enumeration for small N (D=26 is too big for 3^N amplitudes only when N large; N=4 gives 81 amplitudes)."""
    import numpy as np

    Am = {s: np.array(A[s], dtype=object) for s in (-1, 0, 1)}
    psi = {}
    for word in itertools.product((-1, 0, 1), repeat=N):
        M = Am[word[0]]
        for s in word[1:]:
            M = M.dot(Am[s])
        psi[word] = sum(M[i, i] for i in range(26))
    norm = sum(v * v for v in psi.values())
    # H = sum_j h_{j,j+1}, h = Sz Sz + (S+S- + S-S+)/2 with ladder entries sqrt2 -> off-diagonal coefficient 1 between (s,t) and (s+1,t-1)
    Hpsi = {w: 0 for w in psi}
    for w in psi:
        for j in range(N):
            j2 = (j + 1) % N
            s, t = w[j], w[j2]
            Hpsi[w] += s * t * psi[w]
            # <w| h |w'> = 1 when w' differs by (s+1,t-1) or (s-1,t+1) from w
            for ds in (1, -1):
                s2, t2 = s + ds, t - ds
                if -1 <= s2 <= 1 and -1 <= t2 <= 1:
                    w2 = list(w)
                    w2[j], w2[j2] = s2, t2
                    Hpsi[w] += psi[tuple(w2)]
    expect = sum(psi[w] * Hpsi[w] for w in psi)
    return norm, expect


def run(selftest):
    A = build_A()
    nnz = sum(1 for s in A for row in A[s] for x in row if x != 0)
    chk(
        f"matrices A_s built from Tables 1-2 and rules (C.4)-(C.5): D=26, {nnz} nonzero entries",
        len(A[0]) == 26,
    )
    # block dimensions of the d-e decomposition (C.6)
    labs = labels()
    from collections import Counter

    cnt = Counter()
    for a, d in labs:
        for b, e in labs:
            cnt[d - e] += 1
    dims = [cnt[x] for x in range(-10, 11, 2)]
    chk(
        f"(C.6) difference-block dimensions {dims} = [1, 8, 32, 80, 136, 162, 136, 80, 32, 8, 1], total {sum(dims)} = 676",
        dims == [1, 8, 32, 80, 136, 162, 136, 80, 32, 8, 1] and sum(dims) == 676,
    )
    # Lemma C.1 cross-check on a small ring (N=4): <psi|H|psi> = N U_N and |psi|^2 = V_N
    T, O = build_TO(A)
    U4, V4 = trace(matpow(T, 2) * O), trace(matpow(T, 4))
    nrm, ex = lemma_c1_bruteforce(A, 4)
    chk(
        f"Lemma C.1 on N=4 by brute force: |psi|^2 = {nrm} equals V_4 = {V4}; <psi|H|psi> = {ex} equals 4*U_4 = {4 * U4}",
        nrm == V4 and ex == 4 * U4,
    )
    for N, cN in ((72, 1401482112601), (120, 1401482105174)):
        t1 = time.time()
        U, V = contraction(A, N)
        U, V = int(U), int(V)
        ok_pos = V > 0
        ok1 = -cN * V <= 10**12 * U < (-cN + 1) * V
        ratio = Fr(U, V)
        a = Fr(700741, 500000)
        chk(
            f"N={N}: V_N > 0 ({len(str(V))} digits); -c_N V <= 10^12 U < (-c_N+1) V with c_N = {cN}; U/V = {float(ratio):.12f}",
            ok_pos and ok1,
        )
        chk(
            f"N={N}: U_N/V_N < -a = -700741/500000 (variational bound E0 < -N a; margin {float(-a - ratio):.3e})",
            ratio < -a,
        )
        print(f"   ({time.time() - t1:.0f}s)", file=sys.stderr, flush=True)
        if selftest and N == 72:
            Abad = build_A(diag=[DDIAG[0] + 1] + DDIAG[1:])
            Ub, Vb = contraction(Abad, N)
            Ub, Vb = int(Ub), int(Vb)
            chk(
                f"WRONG INPUT: D_1 -289 -> -288 changes U/V to {float(Fr(Ub, Vb)):.9f}; the stated (C.8) interval for N=72 must FAIL",
                not (Vb > 0 and -cN * Vb <= 10**12 * Ub < (-cN + 1) * Vb),
            )
            pt = dict(PTAB)
            pt[(1, 5)] = -12
            Abad2 = build_A(ptab=pt)
            Ub, Vb = contraction(Abad2, N)
            Ub, Vb = int(Ub), int(Vb)
            chk(
                f"WRONG INPUT: P_15 -11 -> -12: (C.8) interval must FAIL (U/V = {float(Fr(Ub, Vb)) if Vb else 'n/a'})",
                not (Vb > 0 and -cN * Vb <= 10**12 * Ub < (-cN + 1) * Vb),
            )
            chk(
                "WRONG INPUT: a claimed c_72 off by one fails the interval",
                not (-(cN + 1) * V <= 10**12 * U < (-(cN + 1) + 1) * V),
            )
    return 0
