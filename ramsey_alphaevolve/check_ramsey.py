"""Re-check Ramsey lower-bound graphs: no K_s and no independent set of size t.

A file graphs/R_<s>_<t>_ge_<m>.txt must hold the adjacency matrix of a simple undirected graph on m-1 vertices
(square, symmetric, entries 0/1, zero diagonal). The graph then certifies R(s,t) >= m iff it has no clique of
size s and no independent set of size t (an independent set of size t is a clique of size t in the complement).

Two clique-search methods are run on both sides:
  A  bitset depth-first search (uint64 words, popcount pruning), compiled with numba, parallel over the
     smallest vertex of the clique.
  B  branch and bound with a greedy-colouring bound on Python big-integer bitsets, multiprocessing over root vertices.

Usage:
  python check_ramsey.py graphs/R_3_12_ge_53.txt [more files]   # both methods
  python check_ramsey.py --cheap                                # the graphs listed in CHEAP
  python check_ramsey.py --selftest                             # planted K_s / I_t must be detected
Options: --method A|B|both (default both), --procs N (method B processes).
Exit status is nonzero on any malformed file, any found K_s or I_t, or a failed self-test.
Progress and timing go to stderr; stdout is deterministic.
"""

import argparse
import hashlib
import os
import random
import re
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
CHEAP = [
    "R_3_12_ge_53",
    "R_3_13_ge_61",
    "R_3_15_ge_74",
    "R_3_17_ge_93",
    "R_3_18_ge_100",
    "R_4_13_ge_139",
    "R_4_14_ge_148",
    "R_4_15_ge_163",
    "R_4_16_ge_174",
]


def parse_name(path):
    m = re.fullmatch(r"R_(\d+)_(\d+)_ge_(\d+)\.txt", os.path.basename(path))
    if not m:
        raise ValueError("file name must look like R_<s>_<t>_ge_<m>.txt: " + path)
    return tuple(int(x) for x in m.groups())


def parse_matrix(text, expected_n=None):
    tok = text.replace("[", " ").replace("]", " ").split()
    if not all(x in ("0", "1") for x in tok):
        raise ValueError("entry that is not 0 or 1")
    n = int(round(len(tok) ** 0.5))
    if n * n != len(tok) or n == 0:
        raise ValueError("matrix is not square")
    A = np.array(tok, dtype=np.uint8).reshape(n, n)
    if not (A == A.T).all():
        raise ValueError("matrix is not symmetric")
    if not (np.diag(A) == 0).all():
        raise ValueError("nonzero diagonal (self loop)")
    if expected_n is not None and n != expected_n:
        raise ValueError("expected %d vertices (m-1), found %d" % (expected_n, n))
    return A


def load(path):
    s, t, m = parse_name(path)
    with open(path, "r", encoding="ascii") as f:
        return s, t, m, parse_matrix(f.read(), m - 1)


def to_bits(A):
    n = A.shape[0]
    W = (n + 63) // 64
    B = np.zeros((n, W), dtype=np.uint64)
    for i in range(n):
        for j in np.nonzero(A[i])[0]:
            B[i, j >> 6] |= np.uint64(1) << np.uint64(j & 63)
    return B


@njit(cache=False)
def popc(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + (
        (x >> np.uint64(2)) & np.uint64(0x3333333333333333)
    )
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return (x * np.uint64(0x0101010101010101)) >> np.uint64(56)


@njit(cache=False)
def ctz(x):
    c = 0
    while (x & np.uint64(1)) == np.uint64(0):
        x >>= np.uint64(1)
        c += 1
    return c


@njit(parallel=True, cache=False)
def find_clique_a(B, n, k, out):
    """out[v,:k] = a k-clique whose smallest vertex is v; found[v]=1. Parallel over v."""
    W = B.shape[1]
    found = np.zeros(n, dtype=np.int64)
    for v in prange(n):
        cand = np.zeros((k + 1, W), dtype=np.uint64)
        path = np.zeros(k + 1, dtype=np.int64)
        path[0] = v
        for w in range(W):
            cand[1, w] = B[v, w]
        for u in range(v + 1):
            cand[1, u >> 6] &= ~(np.uint64(1) << np.uint64(u & 63))
        d = 1
        if k == 1:
            found[v] = 1
            out[v, 0] = v
            continue
        while d >= 1:
            cnt = 0
            for w in range(W):
                cnt += popc(cand[d, w])
            if cnt < k - d:
                d -= 1
                continue
            u = -1
            for w in range(W):
                if cand[d, w] != np.uint64(0):
                    b = ctz(cand[d, w])
                    u = w * 64 + b
                    cand[d, w] &= ~(np.uint64(1) << np.uint64(b))
                    break
            if u < 0:
                d -= 1
                continue
            path[d] = u
            if d + 1 == k:
                for i in range(k):
                    out[v, i] = path[i]
                found[v] = 1
                break
            for w in range(W):
                cand[d + 1, w] = cand[d, w] & B[u, w]
            d += 1
    return found


def is_clique(A, w, k):
    return len(set(w)) == k and all(A[a, b] for a in w for b in w if a != b)


def check_a(A, k):
    """Return a k-clique (list of vertices) or None."""
    n = A.shape[0]
    out = np.zeros((n, k), dtype=np.int64)
    f = find_clique_a(to_bits(A), n, k, out)
    idx = np.nonzero(f)[0]
    if len(idx) == 0:
        return None
    w = [int(x) for x in out[idx[0]]]
    assert is_clique(A, w, k)
    return w


_ADJ = None


def _init(adj):
    global _ADJ
    _ADJ = adj


def _root(args):
    k, v = args
    P = _ADJ[v] & ~((1 << (v + 1)) - 1)
    return _expand(_ADJ, P, [v], k)


def _expand(adj, P, R, k):
    if len(R) == k:
        return list(R)
    if P.bit_count() < k - len(R):
        return None
    order, colors = [], []
    U, c = P, 0
    while U:
        c += 1
        Q = U
        while Q:
            low = Q & -Q
            v = low.bit_length() - 1
            order.append(v)
            colors.append(c)
            U &= ~low
            Q &= ~low
            Q &= ~adj[v]
    for i in range(len(order) - 1, -1, -1):
        if len(R) + colors[i] < k:
            return None
        v = order[i]
        r = _expand(adj, P & adj[v], R + [v], k)
        if r:
            return r
        P &= ~(1 << v)
    return None


def check_b(A, k, procs=1):
    """Return a k-clique (list of vertices) or None."""
    import multiprocessing as mp

    n = A.shape[0]
    adj = [
        int("".join("1" if A[i, j] else "0" for j in range(n - 1, -1, -1)), 2)
        for i in range(n)
    ]
    tasks = [(k, v) for v in range(n)]
    pool = None
    if procs <= 1:
        _init(adj)
        it = map(_root, tasks)
    else:
        pool = mp.Pool(procs, initializer=_init, initargs=(adj,))
        it = pool.imap_unordered(_root, tasks, chunksize=2)
    res = None
    for r in it:
        if r:
            res = r
            break
    if pool:
        pool.terminate()
        pool.join()
    if res:
        assert is_clique(A, res, k)
    return res


def complement(A):
    C = (1 - A).astype(np.uint8)
    np.fill_diagonal(C, 0)
    return C


def run_methods(A, k, methods, procs):
    """dict method -> witness or None"""
    r = {}
    if "A" in methods:
        r["A"] = check_a(A, k)
    if "B" in methods:
        r["B"] = check_b(A, k, procs)
    return r


def fmt(r):
    return "/".join("FOUND" if r[m] else "free" for m in sorted(r))


def check_file(path, methods, procs):
    s, t, m, A = load(path)
    n = A.shape[0]
    with open(path, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    t0 = time.time()
    rk = run_methods(A, s, methods, procs)
    ri = run_methods(complement(A), t, methods, procs)
    dt = time.time() - t0
    ok = not any(rk.values()) and not any(ri.values())
    verdict = "VERIFIED by %s" % ("+".join(sorted(methods))) if ok else "FAIL"
    if not ok:
        for name, rr in (("K%d" % s, rk), ("I%d" % t, ri)):
            for meth, w in rr.items():
                if w:
                    print("  witness %s (method %s): %s" % (name, meth, w))
    print(
        "R(%d,%d)>=%d  n=%d  sha256=%s  K%d-free(%s)  I%d-free(%s)  %s"
        % (s, t, m, n, sha, s, fmt(rk), t, fmt(ri), verdict)
    )
    print("  %s: %.1f s" % (os.path.basename(path), dt), file=sys.stderr)
    return ok


def plant_clique(A, k, rng):
    B = A.copy()
    vs = rng.sample(range(A.shape[0]), k)
    for a in vs:
        for b in vs:
            if a != b:
                B[a, b] = 1
    return B, vs


def plant_indep(A, k, rng):
    B = A.copy()
    vs = rng.sample(range(A.shape[0]), k)
    for a in vs:
        for b in vs:
            B[a, b] = 0
    return B, vs


def selftest(procs):
    """Planted K_s and I_t must be found by both methods; the unplanted graph must pass; bad input must be rejected."""
    import contextlib
    import io
    import tempfile

    path = os.path.join(HERE, "graphs", "R_3_12_ge_53.txt")
    s, t, m, A = load(path)
    rng = random.Random(20261022)
    good = True

    def expect(label, cond):
        nonlocal good
        print("selftest %-58s %s" % (label, "ok" if cond else "FAILED"))
        good = good and bool(cond)

    expect(
        "unplanted: no K%d (A,B)" % s, not any(run_methods(A, s, "AB", procs).values())
    )
    expect(
        "unplanted: no I%d (A,B)" % t,
        not any(run_methods(complement(A), t, "AB", procs).values()),
    )
    for k in (s, s + 1):
        P, vs = plant_clique(A, k, rng)
        r = run_methods(P, k, "AB", procs)
        expect(
            "planted K%d detected by A and B, valid witness" % k,
            all(r[x] is not None and is_clique(P, r[x], k) for x in "AB"),
        )
    P, vs = plant_indep(A, t, rng)
    C = complement(P)
    r = run_methods(C, t, "AB", procs)
    expect(
        "planted I%d detected by A and B, valid witness" % t,
        all(r[x] is not None and is_clique(C, r[x], t) for x in "AB"),
    )
    P, vs = plant_clique(A, s, rng)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "R_%d_%d_ge_%d.txt" % (s, t, m))
        np.savetxt(p, P, fmt="%d")
        with (
            contextlib.redirect_stdout(io.StringIO()),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            res = check_file(p, "AB", procs)
        expect("check_file rejects a file with a planted K%d" % s, res is False)

    def rejects(text, n=None):
        try:
            parse_matrix(text, n)
            return False
        except ValueError:
            return True

    expect("rejects non-symmetric matrix", rejects("0 1 0\n0 0 0\n0 0 0"))
    expect("rejects nonzero diagonal", rejects("1 0\n0 0"))
    expect("rejects entry 2", rejects("0 2\n2 0"))
    expect("rejects non-square", rejects("0 1 0\n1 0 0"))
    expect("rejects wrong vertex count", rejects("0 1\n1 0", 3))
    print("selftest %s" % ("PASSED" if good else "FAILED"))
    return good


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("files", nargs="*")
    ap.add_argument(
        "--cheap", action="store_true", help="check the graphs listed in CHEAP"
    )
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--method", choices=["A", "B", "both"], default="both")
    ap.add_argument("--procs", type=int, default=min(8, os.cpu_count() or 1))
    a = ap.parse_args()
    methods = "AB" if a.method == "both" else a.method
    if a.selftest:
        return 0 if selftest(a.procs) else 1
    files = list(a.files)
    if a.cheap:
        files += [os.path.join(HERE, "graphs", n + ".txt") for n in CHEAP]
    if not files:
        ap.error("give graph files, --cheap or --selftest")
    ok = True
    for p in files:
        try:
            ok = check_file(p, methods, a.procs) and ok
        except (ValueError, OSError) as e:
            print("%s: INVALID: %s" % (os.path.basename(p), e))
            ok = False
    print("ALL PASSED" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
