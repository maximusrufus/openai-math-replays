"""Independent replay of the finite computations in 'The Campana-Peternell conjecture in dimension six' (OpenAI, 2026-10; Sections 3-5).

Written from the manuscript text only (Sections 3, 4, 5); no code from the authors is used.  Exact arithmetic (python-flint
fmpq_mpoly / fmpz_mat / nmod_mat, sympy) or modular (a rank mod a prime is a rigorous lower bound for the rational rank when no
denominator is divisible by the prime; this is checked).

  python replay_campana_peternell.py             # every claim with PASS/FAIL; exit 0 only if all pass (about 3-6 min)
  python replay_campana_peternell.py --perturb   # adversarial: one witness entry changed; the slice identity must FAIL (exit 1)
  python replay_campana_peternell.py --selftest  # also runs the built-in wrong-input controls (each must be rejected)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import itertools
import re
import sys
import time
from math import comb, factorial

from flint import fmpq, fmpq_mat, fmpq_mpoly_ctx, fmpz_mat, nmod_mat

RES = []
SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
PFX = "CP"
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


W = (1, 2, 3, 4, 5, 6, 1, 1, 2, 3)  # weights of h,x2..x6,l,r,s,t
NAMES = ("h", "x2", "x3", "x4", "x5", "x6", "l", "r", "s", "t")
ctx = fmpq_mpoly_ctx.get(NAMES, "lex")
h, x2, x3, x4, x5, x6, l, r, s, t = ctx.gens()
ONE = ctx.from_dict({(0,) * 10: 1}) if hasattr(ctx, "from_dict") else h**0
c = 2 * h - l


# --- monomial enumeration (weights) ------------------------------------------------
def monos(idx, k):
    """exponent tuples (length 10) of weight k in the variable indices idx."""
    if k < 0:
        return []
    out = []

    def rec(pos, rem, cur):
        if pos == len(idx):
            if rem == 0:
                e = [0] * 10
                for ii, ex in zip(idx, cur):
                    e[ii] = ex
                out.append(tuple(e))
            return
        w = W[idx[pos]]
        for ex in range(rem // w + 1):
            rec(pos + 1, rem - ex * w, cur + [ex])

    rec(0, k, [])
    return out


BASE = [0, 1, 2, 3, 4, 5]
IDEAL = [6, 7, 8, 9]
ALL = list(range(10))


def mono_poly(e):
    return ctx.from_dict({tuple(e): 1})


def B(k):
    return [mono_poly(e) for e in monos(BASE, k)]


def I(k):
    return [mono_poly(e) for e in monos(IDEAL, k)]


def T(k):
    return [mono_poly(e) for e in monos(ALL, k)]


# --- Newton / q / Todd / psi / D --------------------------------------------------
NMAX = 10
p_ = [ctx.from_dict({(0,) * 10: 3})]
sj = {1: r, 2: s, 3: t}
for i in range(1, NMAX + 1):
    acc = ctx.from_dict({})
    for j in range(1, min(3, i) + 1):
        term = sj[j] * (i if j == i else p_[i - j])
        acc = acc + term if (j + 1) % 2 == 0 else acc - term
    p_.append(acc)
# sign: (-1)^(j+1): j=1 -> +, j=2 -> -, j=3 -> +   (fix below, written explicitly)
p_ = [p_[0]]
for i in range(1, NMAX + 1):
    acc = ctx.from_dict({})
    for j in range(1, min(3, i) + 1):
        term = sj[j] * (i if j == i else p_[i - j])
        acc = acc + term if (j % 2 == 1) else acc - term
    p_.append(acc)

q_ = []
for i in range(NMAX + 1):
    acc = ctx.from_dict({})
    for j in range(i + 1):
        acc = acc + comb(i, j) * p_[j] * (-h) ** (i - j)
    q_.append(acc)

beta = {
    1: fmpq(1, 2),
    2: fmpq(-1, 12),
    4: fmpq(1, 720),
    6: fmpq(-1, 30240),
    8: fmpq(1, 1209600),
}
theta = [ONE]
for i in range(1, NMAX + 1):
    acc = ctx.from_dict({})
    for j in range(1, i + 1):
        if j in beta:
            acc = acc + beta[j] * q_[j] * theta[i - j]
    theta.append(acc * fmpq(1, i))

psi = []
for i in range(NMAX + 1):
    acc = ctx.from_dict({})
    for j in range(i + 1):
        acc = acc + fmpq((-1) ** j, factorial(j)) * q_[j] * theta[i - j]
    psi.append(acc)


def Dfun(i):
    first = (1 - (-1) ** i) * c ** (i - 1)
    acc = ctx.from_dict({})
    for j in range(i):
        inner = ctx.from_dict({})
        for k in range(i - j):
            inner = inner + h**k * (l - h) ** (i - j - 1 - k)
        acc = acc + comb(i, j) * p_[j] * ((-1) ** (i - j)) * inner
    return first - acc


Dl = {i: Dfun(i) for i in range(1, NMAX + 1)}

# --- Delta recursion (4.8) ----------------------------------------------------------
_memo = {}
ZERO = ctx.from_dict({})


def Delta(e):
    """e: exponent tuple in the six base generators only (others 0)."""
    e = tuple(e)
    if e in _memo:
        return _memo[e]
    if sum(e) == 0:
        res = ZERO
    else:
        g = next(
            i for i in range(6) if e[i] > 0
        )  # first generator dividing m in displayed order
        n = list(e)
        n[g] -= 1
        n = tuple(n)
        gp = mono_poly([1 if i == g else 0 for i in range(10)])
        dg = (
            ctx.from_dict({(0,) * 10: 1}) if g == 0 else Dl[g + 1]
        )  # Delta h = 1, Delta x_i = D_i (i = g+1)
        res = dg * mono_poly(n) + (gp - c * dg) * Delta(n)
    _memo[e] = res
    return res


def Delta_poly(m):
    ((e, co),) = m.to_dict().items()
    return Delta(e) * co


# --- build A ------------------------------------------------------------------------
def build_rows(drop=()):
    rows = []
    for i in range(7):
        for m in B(i):
            for w in I(9 - i):
                rows.append(("4.11", w * (m - h * Delta_poly(m))))
    for i in range(7):
        for m in B(i):
            for w in I(10 - i):
                rows.append(("4.12", w * Delta_poly(m)))
    for i in range(7, 11):
        for m in B(i):
            for w in T(9 - i):
                rows.append(("4.13", w * m))
    for i in range(7, 11):
        for m in B(i):
            for w in T(10 - i):
                rows.append(("4.14", w * Delta_poly(m)))
    for i in range(7, 11):
        for w in T(10 - i):
            rows.append(("4.15", w * Dl[i]))
    for i in range(4, 10):
        for m in B(9 - i):
            rows.append(("4.16t", m * theta[i]))
            rows.append(("4.16p", m * psi[i]))
    return [(f, p) for f, p in rows if f not in drop]


def dense_rows(rows, cols):
    idx = {e: k for k, e in enumerate(cols)}
    out = []
    for f, p in rows:
        d = p.to_dict()
        if not d:
            continue
        v = [fmpq(0)] * len(cols)
        for e, co in d.items():
            assert sum(a * b for a, b in zip(e, W)) == 9, "row of wrong weight"
            v[idx[e]] = co
        out.append((f, v))
    return out


def to_fmpz_rows(rows_dense):
    M = []
    for f, v in rows_dense:
        den = 1
        for x in v:
            den = den * int(x.q) // __import__("math").gcd(den, int(x.q))
        M.append([int(x * den) for x in v])
    return fmpz_mat(M)


def iszero(M):
    return all(M[i, j] == 0 for i in range(M.nrows()) for j in range(M.ncols()))


def apply_pair_check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    return cond


def sec3():
    """Section 3: moment identities and lambda = 9/4 (exact, sympy)."""
    import sympy as sp
    global m
    m = sp.symbols("m0:8")




    def E(i, sign=-1):
        """sum_{j<=i} (-1)^j C(i,j) m_j - sign*m_i, with m0 = m5 = m6 = m7 = 0 substituted"""
        e = sum((-1) ** j * comb(i, j) * m[j] for j in range(i + 1)) - sign * m[i]
        return sp.expand(e.subs({m[0]: 0, m[5]: 0, m[6]: 0, m[7]: 0}))


    chk("i=2 identity reads 2(m2 - m1) = 0, i.e. m2 = m1", sp.expand(E(2) - 2 * (m[2] - m[1])) == 0)
    e4 = E(4).subs(m[2], m[1])
    chk("i=4 identity with m2 = m1 equals 2*(m1 - 2 m3 + m4)", sp.expand(e4 - 2 * (m[1] - 2 * m[3] + m[4])) == 0)
    e6 = E(6).subs(m[2], m[1])
    chk("i=6 identity with m2 = m1 equals 9 m1 - 20 m3 + 15 m4", sp.expand(e6 - (9 * m[1] - 20 * m[3] + 15 * m[4])) == 0)
    s = sp.solve([m[1] - 2 * m[3] + m[4], 9 * m[1] - 20 * m[3] + 15 * m[4]], [m[1], m[3]], dict=True)[0]
    chk(f"solution m1 = {s[m[1]]}, m3 = {s[m[3]]} (stated m1 = m2 = 5 m4, m3 = 3 m4)", sp.simplify(s[m[1]] - 5 * m[4]) == 0 and sp.simplify(s[m[3]] - 3 * m[4]) == 0)
    lam = sp.symbols("lam")
    ls = sp.solve(sp.Eq(3**4 * s[m[3]] - 4 * 3**3 * lam * m[4], 0), lam)[0]
    chk(f"0 = 3^4 m3 - 4*3^3*lam*m4 gives lam = {ls} = 9/4, not an integer", ls == sp.Rational(9, 4) and ls.q != 1)
    bad = sp.solve([m[1] - 2 * m[3] + m[4], 9 * m[1] - 21 * m[3] + 15 * m[4]], [m[1], m[3]], dict=True)[0]
    chk(f"WRONG INPUT: coefficient 20 -> 21 changes m3 to {bad[m[3]]} != 3 m4", sp.simplify(bad[m[3]] - 3 * m[4]) != 0)
    chk("WRONG INPUT: a sign error in the i=4 identity (sign=+1) no longer reproduces the stated equation", sp.expand(E(4, sign=+1).subs(m[2], m[1]) - 2 * (m[1] - 2 * m[3] + m[4])) != 0)
    chk("WRONG INPUT: claimed lambda = 2 is rejected", ls != 2)


def run():
    selftest = SELFTEST
    t0 = time.time()
    res = []
    cols = monos(ALL, 9)
    res.append(("T9 has 581 monomials", len(cols) == 581))
    rows = build_rows()
    dr = dense_rows(rows, cols)
    res.append((f"A has 1785 nonzero rows (got {len(dr)})", len(dr) == 1785))
    Az = to_fmpz_rows(dr)
    rk = Az.rank()
    res.append((f"rank_Q(A) = 516 (got {rk}), nullity {581 - rk}", rk == 516))
    X, nul = Az.nullspace()
    res.append((f"dim ker A = 65 (got {nul})", nul == 65))
    K = [[X[i, j] for j in range(nul)] for i in range(581)]  # 581 x 65 integer basis
    # AK = 0 exactly
    from flint import fmpz_mat as ZM

    Km = ZM([[K[i][j] for j in range(nul)] for i in range(581)])
    res.append(
        (
            "A K = 0 over Q",
            iszero(Az * Km),
        )
    )
    res.append(("rank(K) = 65", Km.rank() == 65))
    # selected rows
    Mexp = [
        (6, 0, 0, 0, 0, 0),
        (4, 1, 0, 0, 0, 0),
        (2, 2, 0, 0, 0, 0),
        (0, 3, 0, 0, 0, 0),
        (3, 0, 1, 0, 0, 0),
        (1, 1, 1, 0, 0, 0),
        (0, 0, 2, 0, 0, 0),
        (2, 0, 0, 1, 0, 0),
        (0, 1, 0, 1, 0, 0),
        (1, 0, 0, 0, 1, 0),
        (0, 0, 0, 0, 0, 1),
    ]
    Gexp = [  # exponents of (l, r, s, t)
        (3, 0, 0, 0),
        (2, 1, 0, 0),
        (1, 2, 0, 0),
        (0, 3, 0, 0),
        (1, 0, 1, 0),
        (0, 1, 1, 0),
        (0, 0, 0, 1),
        (4, 0, 0, 0),
        (3, 1, 0, 0),
        (2, 2, 0, 0),
        (1, 3, 0, 0),
        (0, 4, 0, 0),
        (2, 0, 1, 0),
        (1, 1, 1, 0),
        (0, 2, 1, 0),
        (0, 0, 2, 0),
        (1, 0, 0, 1),
        (0, 1, 0, 1),
    ]
    wG = [sum(a * b for a, b in zip(g, (1, 1, 2, 3))) for g in Gexp]
    J = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17]
    # NB: J indexes G.  The text's M index set (i) runs over 0..10, G index j over J.
    C1 = []
    for i in range(11):
        hdiv = Mexp[i][0] > 0
        for j in J:
            if wG[j] == 3 or hdiv:
                C1.append((i, j))
    res.append((f"|C1| = 129 (got {len(C1)})", len(C1) == 129))
    cidx = {e: k for k, e in enumerate(cols)}

    def sel_monomial(i, j):
        e = list(Mexp[i]) + list(Gexp[j])
        e[0] -= wG[j] - 3
        assert e[0] >= 0
        return tuple(e)

    P = [[K[cidx[sel_monomial(i, j)]][k] for k in range(nul)] for (i, j) in C1]
    Pm = fmpq_mat(len(C1), nul, [x for rw in P for x in rw])
    rP = Pm.rank()
    res.append(
        (
            f"rank_Q(P) = 60 (got {rP}), P is {len(C1)} x {nul}",
            rP == 60 and len(C1) == 129,
        )
    )
    # N: rref of a basis of ker(P^T)
    PT = fmpq_mat(nul, len(C1), [Pm[i, j] for j in range(nul) for i in range(len(C1))])
    PTz = fmpz_mat([[int(PT[i, j] * 1) for j in range(len(C1))] for i in range(nul)])
    Y, ny = PTz.nullspace()
    Nm = fmpq_mat(
        ny, len(C1), [fmpq(int(Y[j, i])) for i in range(ny) for j in range(len(C1))]
    )
    Nr, rankN = Nm.rref()
    Nrows = [[Nr[i, j] for j in range(len(C1))] for i in range(rankN)]
    res.append(
        (
            f"N has 69 rows (got {rankN}) and N P = 0",
            rankN == 69 and iszero(fmpq_mat(rankN, len(C1), [x for rw in Nrows for x in rw]) * Pm),
        )
    )
    Nfm = fmpq_mat(rankN, len(C1), [x for rw in Nrows for x in rw])
    Nden = 1
    for rw in Nrows:
        for x in rw:
            Nden = Nden * int(x.q) // __import__("math").gcd(Nden, int(x.q))
    res.append((f"all N denominators prime to 10007 (lcm {Nden})", Nden % 10007 != 0))

    if PERTURB:
        for name, ok in res:
            chk(name, ok)
        Ha_p = [fmpq(x) for x in (1, 4, 16, 64, 6, 24, 4, 3, 12, 48, 192, 768, 17, 68, 272, 96, 9, 36)]
        Ha_p[4] += 1
        B0 = [fmpq(x) for x in (1, 7, 45, 275, 11, 65, 85, 17, 95, 25, 35)]
        vec = fmpq_mat(len(C1), 1, [B0[i] * Ha_p[j] for (i, j) in C1])
        chk("N (B_i H_j) = 0 for the perturbed witness (Ha[4] increased by 1)", iszero(Nfm * vec))
        return finish()

    # --- K2 / K3 ------------------------------------------------------------------
    jpos = {j: k for k, j in enumerate(J)}
    c1pos = {e: k for k, e in enumerate(C1)}

    def Cn(n):
        out = []
        for z in itertools.combinations_with_replacement(J, n):
            has3 = any(wG[j] == 3 for j in z)
            for i in range(11):
                if has3 or Mexp[i][0] > 0:
                    out.append((i, z))
        return out

    C2, C3 = Cn(2), Cn(3)
    res.append(
        (
            f"|C2| = 1140, |C3| = 6820 (got {len(C2)}, {len(C3)})",
            len(C2) == 1140 and len(C3) == 6820,
        )
    )

    def kernel_matrix_mod(Ccols, n, prime, drop_pair=None):
        pos = {e: k for k, e in enumerate(Ccols)}
        rowsM = []
        for combo in itertools.combinations_with_replacement(J, n - 1):
            if drop_pair is not None and combo == drop_pair:
                continue
            for rr in range(rankN):
                row = {}
                for (i, j), kk in c1pos.items():
                    co = Nrows[rr][kk]
                    if co == 0:
                        continue
                    z = tuple(sorted((j,) + combo))
                    col = pos[(i, z)]
                    row[col] = (
                        row.get(col, 0) + int(co.p) * pow(int(co.q), -1, prime)
                    ) % prime
                rowsM.append(row)
        M = nmod_mat(len(rowsM), len(Ccols), prime)
        for ri, row in enumerate(rowsM):
            for cj, v in row.items():
                M[ri, cj] = v
        return M

    for prime in (10007, 1000003):
        # K2
        M2 = kernel_matrix_mod(C2, 2, prime)
        d2 = len(C2) - M2.rank()
        # K3
        M3 = kernel_matrix_mod(C3, 3, prime)
        rk3 = M3.rank()
        d3 = len(C3) - rk3
        if prime == 10007:
            res.append(
                (f"mod 10007: dim K1 = 60 (129-69), dim K2 = 163 (got {d2})", d2 == 163)
            )
            res.append(
                (
                    f"mod 10007: K3 matrix {M3.nrows()} x {M3.ncols()} (8280 x 6820), rank {rk3} (6814), dim K3 = {d3} (6)",
                    M3.nrows() == 8280 and rk3 == 6814 and d3 == 6,
                )
            )
        else:
            res.append((f"mod {prime} (second prime): dim K3 = {d3} <= 6", d3 == 6))
        if prime == 10007 and selftest:
            Mbad = kernel_matrix_mod(C3, 3, prime, drop_pair=(0, 0))
            res.append(
                (
                    "WRONG INPUT: drop one pair of slice conditions -> dim K3 must exceed 6",
                    len(C3) - Mbad.rank() > 6,
                )
            )

    # --- six explicit exact columns ------------------------------------------------
    Bv = lambda v: [
        fmpq(1),
        7 - fmpq(2) * v / 5,
        45 - 4 * v,
        275 - 30 * v,
        fmpq(11),
        65 - 2 * v,
        fmpq(85),
        fmpq(17),
        95 - 2 * v,
        fmpq(25),
        fmpq(35),
    ]
    Bp = [
        fmpq(1),
        fmpq(11, 5),
        fmpq(5),
        fmpq(11),
        fmpq(-1),
        fmpq(-1),
        fmpq(10),
        fmpq(-3),
        fmpq(-1),
        fmpq(0),
        fmpq(14),
    ]
    Ha = [
        fmpq(x)
        for x in (1, 4, 16, 64, 6, 24, 4, 3, 12, 48, 192, 768, 17, 68, 272, 96, 9, 36)
    ]
    Hb = [
        fmpq(x, 3)
        for x in (3, 8, 21, 54, 9, 24, 6, 14, 37, 97, 252, 648, 37, 97, 252, 96, 16, 42)
    ]
    res.append(
        (
            "sanity: |B|=11,|B'|=11,|Ha|=18,|Hb|=18",
            len(Bv(fmpq(0))) == 11
            and len(Bp) == 11
            and len(Ha) == 18
            and len(Hb) == 18,
        )
    )
    # (N B_i H_j)_{C1} = 0 for the three product-form witnesses (the paper's "direct rational substitution")
    cols_prod = []
    for Bb, Hh in ((Bv(fmpq(0)), Ha), (Bv(fmpq(1)), Ha), (Bp, Hb)):
        vec = fmpq_mat(len(C1), 1, [Bb[i] * Hh[j] for (i, j) in C1])
        res.append(("N (B_i H_j) = 0 (product witness)", iszero(Nfm * vec)))
        cols_prod.append(
            {e: Bb[e[0]] * Hh[e[1][0]] * Hh[e[1][1]] * Hh[e[1][2]] for e in C3}
        )
    # sparse columns
    s1 = {
        (5, (15, 15, 15)): fmpq(3),
        (7, (15, 15, 15)): fmpq(1),
        (9, (15, 15, 15)): fmpq(-5),
    }
    s2 = {
        (7, (10, 10, 15)): fmpq(2),
        (7, (10, 10, 10)): fmpq(4),
        (9, (10, 10, 10)): fmpq(25),
    }
    s3 = {(7, (10, 10, 10)): fmpq(1)}

    def slice_check(vdict):
        """exact check: N applied to every double slice (120 pairs) of the column v on C3."""
        bad = 0
        for a_, b_ in itertools.combinations_with_replacement(J, 2):
            vec = fmpq_mat(
                len(C1),
                1,
                [vdict.get((i, tuple(sorted((j, a_, b_)))), fmpq(0)) for (i, j) in C1],
            )
            if not iszero(Nfm * vec):
                bad += 1
        return bad

    names6 = ["B(0)Ha", "B(1)Ha", "B'Hb", "s1", "s2", "s3"]
    dicts6 = cols_prod + [s1, s2, s3]
    for nm, dct in zip(names6, dicts6):
        bad = slice_check(dct)
        res.append(
            (
                f"column {nm}: N kills all 120 double slices exactly (bad slices {bad})",
                bad == 0,
            )
        )
    # independence minor
    rows6 = [
        (0, (0, 0, 0)),
        (0, (0, 0, 1)),
        (1, (0, 0, 0)),
        (5, (15, 15, 15)),
        (7, (10, 10, 15)),
        (7, (10, 10, 10)),
    ]
    Mn = fmpq_mat(
        6,
        6,
        [
            dct.get(rw, fmpq(0)) if isinstance(dct, dict) else 0
            for rw in rows6
            for dct in dicts6
        ],
    )
    det = Mn.det()
    res.append(
        (f"6x6 minor on the displayed rows = -16/5 (got {det})", det == fmpq(-16, 5))
    )
    ul = fmpq_mat(3, 3, [Mn[i, j] for i in range(3) for j in range(3)]).det()
    res.append((f"upper-left det = -8/15 (got {ul})", ul == fmpq(-8, 15)))

    # functional (5.2)
    def func(dct):
        return (
            32 * dct.get((0, (0, 0, 0)), fmpq(0))
            - 20 * dct.get((0, (0, 0, 1)), fmpq(0))
            + 3 * dct.get((0, (0, 1, 1)), fmpq(0))
        )

    vals = [func(d) for d in dicts6]
    res.append(
        (
            f"functional (5.2) vanishes on all six columns (values {[str(v) for v in vals]})",
            all(v == 0 for v in vals),
        )
    )
    # u values
    res.append(
        (
            "32-20u+3u^2 vanishes at u=4 and u=8/3 (and (u-4)(3u-8) = 3u^2-20u+32)",
            all(32 - 20 * u + 3 * u * u == 0 for u in (fmpq(4), fmpq(8, 3))),
        )
    )
    if selftest:
        # deliberately wrong inputs: each must be detected
        wrong = []
        Hbad = list(Ha)
        Hbad[4] = Hbad[4] + 1
        vec = fmpq_mat(len(C1), 1, [Bv(fmpq(0))[i] * Hbad[j] for (i, j) in C1])
        wrong.append(
            (
                "WRONG INPUT: perturbed Ha[4] -> N kill must FAIL",
                not iszero(Nfm * vec),
            )
        )
        s1b = dict(s1)
        s1b[(5, (15, 15, 15))] = fmpq(4)
        wrong.append(
            (
                "WRONG INPUT: s1 entry 3->4 -> slice check must FAIL",
                slice_check(s1b) > 0,
            )
        )
        wrong.append(
            (
                "WRONG INPUT: functional with 33 instead of 32 must NOT vanish on B(0)Ha",
                (
                    33 * cols_prod[0][(0, (0, 0, 0))]
                    - 20 * cols_prod[0][(0, (0, 0, 1))]
                    + 3 * cols_prod[0][(0, (0, 1, 1))]
                )
                != 0,
            )
        )
        wrong.append(
            (
                "WRONG INPUT: minor with s3 entry 1->0 must differ from -16/5",
                fmpq_mat(
                    6,
                    6,
                    [
                        (
                            fmpq(0)
                            if (rw == (7, (10, 10, 10)) and dct is s3)
                            else (dct.get(rw, fmpq(0)))
                        )
                        for rw in rows6
                        for dct in dicts6
                    ],
                ).det()
                != fmpq(-16, 5),
            )
        )
        # remove the (4.16) rows: nullity of A must exceed 65
        dr2 = dense_rows(build_rows(drop=("4.16t", "4.16p")), cols)
        n2 = 581 - to_fmpz_rows(dr2).rank()
        wrong.append(
            (
                f"WRONG INPUT: drop row family (4.16) -> dim ker A must exceed 65 (got {n2})",
                n2 > 65,
            )
        )
        res.extend(wrong)
    for name, ok in res:
        chk(name, ok)
    sec3()
    log(f"elapsed {time.time() - t0:.1f}s")
    return finish()


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    sys.exit(run())
