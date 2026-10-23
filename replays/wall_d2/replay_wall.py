#!/usr/bin/env python3
"""Independent replay of the computable claims in 'A Counterexample to Wall's D(2) Problem' (OpenAI, 2026-10-06).

Written from the paper text only (pages 7-9 were also read as rendered images to confirm the transcription of presentation (10),
character (14) and matrix (15)). It does not import or reuse any code from the paper's authors.
Exact arithmetic throughout: (1) a self-contained implementation of Q(w), w^2+w+1=0, with Fractions; (2) a self-contained
free-group word engine and Fox calculus; (3) sympy for Smith normal form and the Laurent-polynomial minor as a cross-check.

  python replay_wall.py                  # prints every claim with PASS/FAIL; exit 0 only if all pass
  python replay_wall.py --perturb-matrix # adversarial self-test: changes one entry of matrix (15); must exit 1 (C5 FAIL)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import itertools
import sys
import time
from fractions import Fraction as Fr
from math import gcd


# ----------------------------------------------------------------------------- Q(w), w^2 = -1 - w
class Q:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a, self.b = Fr(a), Fr(b)

    def __add__(s, o):
        o = o if isinstance(o, Q) else Q(o)
        return Q(s.a + o.a, s.b + o.b)

    __radd__ = __add__

    def __neg__(s):
        return Q(-s.a, -s.b)

    def __sub__(s, o):
        return s + (-(o if isinstance(o, Q) else Q(o)))

    def __rsub__(s, o):
        return Q(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, Q) else Q(o)
        return Q(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a - s.b * o.b)

    __rmul__ = __mul__

    def inv(s):
        n = s.a * s.a - s.a * s.b + s.b * s.b  # norm
        if n == 0:
            raise ZeroDivisionError
        return Q((s.a - s.b) / n, -s.b / n)

    def __truediv__(s, o):
        return s * (o if isinstance(o, Q) else Q(o)).inv()

    def __pow__(s, n):
        if n < 0:
            return s.inv() ** (-n)
        r = Q(1)
        for _ in range(n):
            r = r * s
        return r

    def __eq__(s, o):
        o = o if isinstance(o, Q) else Q(o)
        return s.a == o.a and s.b == o.b

    def iszero(s):
        return s.a == 0 and s.b == 0

    def __repr__(s):
        return f"({s.a}+{s.b}w)"


W = Q(0, 1)  # primitive cube root of unity
assert W * W * W == Q(1) and W * W + W + 1 == Q(0)


def rank_det(M):
    """Gaussian elimination over Q(w). -> (rank, det if square else None)."""
    M = [row[:] for row in M]
    n, m = len(M), len(M[0])
    r, det = 0, Q(1)
    for c in range(m):
        piv = next((i for i in range(r, n) if not M[i][c].iszero()), None)
        if piv is None:
            continue
        if piv != r:
            M[r], M[piv] = M[piv], M[r]
            det = -det
        det = det * M[r][c]
        inv = M[r][c].inv()
        M[r] = [x * inv for x in M[r]]
        for i in range(r + 1, n):
            f = M[i][c]
            if not f.iszero():
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
        if r == n:
            break
    return r, (det if n == m and r == n else (Q(0) if n == m else None))


# ----------------------------------------------------------------------------- free group words
def inv_word(w):
    return [(g, -e) for g, e in reversed(w)]


def reduce_word(w):
    out = []
    for g, e in w:
        if out and out[-1][0] == g and out[-1][1] == -e:
            out.pop()
        else:
            out.append((g, e))
    return out


def word(s):
    """'x1 a1 x1^-1 a1^-4' -> list of (gen, +-1); gens are tokens like x1,a1,x2,a2,s."""
    out = []
    for tok in s.split():
        if "^" in tok:
            g, e = tok.split("^")
            e = int(e)
        else:
            g, e = tok, 1
        out += [(g, 1 if e > 0 else -1)] * abs(e)
    return out


def conj(c, r):
    return c + r + inv_word(c)


GENS = ["x1", "a1", "x2", "a2", "s"]
RELS = {
    "r1": word("x1 a1 x1^-1 a1^-4"),
    "r2": word("x2 a2 x2^-1 a2^-3"),
    "r3": word("x1 a2^-13"),
    "r4": word("s x2 s^-1 a1^-5"),
}


# ----------------------------------------------------------------------------- Fox calculus, evaluated under a character
def fox_row(w, rho):
    """d w / d g evaluated at rho (a dict gen -> element of a commutative ring with inverses). Free-differential rules:
    d(uv)=du + u dv, dg/dg = 1, d(g^-1)/dg = -g^-1. Scans left to right keeping the value of the prefix."""
    zero = rho[GENS[0]] * 0
    one = rho[GENS[0]] ** 0 if hasattr(rho[GENS[0]], "__pow__") else 1
    row = {g: zero for g in GENS}
    p = one
    for g, e in w:
        if e == 1:
            row[g] = row[g] + p
            p = p * rho[g]
        else:
            ginv = rho[g] ** -1
            row[g] = row[g] - p * ginv
            p = p * ginv
    return [row[g] for g in GENS], p


def Qchar():
    return {"x1": Q(-1), "a1": W, "x2": W * W, "a2": Q(-1), "s": Q(1)}


# ----------------------------------------------------------------------------- certificates in the free group
def conj_chain(x, a, n, k, R):
    """Factors f (as words) with  x^k a x^-k a^-(n^k)  ==  product of conjugates of R^{+-1},  where R = x a x^-1 a^-n.
    Uses x a^m x^-1 = (R a^n)^m = prod_{j<m} (a^{nj} R a^{-nj}) a^{nm}."""
    ax = lambda m: [(a, 1 if m > 0 else -1)] * abs(m)
    if k == 1:
        return [R]
    prev = conj_chain(
        x, a, n, k - 1, R
    )  # x^{k-1} a x^{-(k-1)} a^-(n^{k-1}) = prod prev
    m = n ** (k - 1)
    out = [conj([(x, 1)], f) for f in prev]  # conjugate by x
    out += [conj(ax(n * j), R) for j in range(m)]  # x a^m x^-1 a^-(nm)
    return out


def product(ws):
    out = []
    for w in ws:
        out += w
    return reduce_word(out)


def check_identity(lhs, rhs_extra, factors):
    """True iff  lhs  ==  (product of factors) * rhs_extra  in the free group (so lhs = rhs_extra modulo the normal closure of the factors)."""
    return reduce_word(lhs) == reduce_word(product(factors) + rhs_extra)


A = lambda g, m: [(g, 1 if m > 0 else -1)] * abs(m)

# ----------------------------------------------------------------------------- the claims
RESULTS = []


def claim(cid, text, ok, detail=""):
    RESULTS.append((cid, text, bool(ok), detail))
    print(
        f"[{'PASS' if ok else 'FAIL'}] {cid}: {text}"
        + (f"  ({detail})" if detail else ""),
        flush=True,
    )


def main(perturb_matrix=False):
    t0 = time.time()
    del RESULTS[:]
    rho = Qchar()

    # C1: rho respects the four relations (paper (14)) and the paper's displayed exponent checks.
    vals = {k: fox_row(w, rho)[1] for k, w in RELS.items()}
    claim(
        "C1",
        "rho(x1,a1,x2,a2,s)=(-1,w,w^2,-1,1) sends every relator r1..r4 to 1",
        all(v == Q(1) for v in vals.values()),
        str(vals),
    )
    claim(
        "C1b",
        "the four displayed identities w=w^4, -1=(-1)^3, -1=(-1)^13, w^2=w^5",
        W == W**4 and Q(-1) == Q(-1) ** 3 and Q(-1) == Q(-1) ** 13 and W * W == W**5,
    )
    claim(
        "C2",
        "rho(t1)=rho(x1^2)=1 and rho(t2)=rho(x2^3)=1 (so rho factors through G)",
        rho["x1"] ** 2 == Q(1) and rho["x2"] ** 3 == Q(1),
    )
    claim("C3", "rho(x1) = -1 != 1 (x1 nontrivial in G)", not (rho["x1"] == Q(1)))
    # the epimorphism to Z: s -> 1, others -> 0 respects relators and kills t1,t2
    phi = {"x1": 0, "a1": 0, "x2": 0, "a2": 0, "s": 1}
    claim(
        "C4",
        "s->1, others->0 respects all four relators and t1,t2 (G surjects onto Z)",
        all(sum(phi[g] * e for g, e in w) == 0 for w in RELS.values())
        and 2 * phi["x1"] == 0
        and 3 * phi["x2"] == 0,
    )

    # C5: Fox matrix, computed by the left-to-right engine, equals the paper's (15)
    D = [fox_row(RELS[k], rho)[0] for k in ("r1", "r2", "r3", "r4")]
    Dp = [
        [Q(1) - W, Q(-3 if perturb_matrix else -2), Q(0), Q(0), Q(0)],
        [Q(0), Q(0), Q(2), W * W - 1, Q(0)],
        [Q(1), Q(0), Q(0), Q(-1), Q(0)],
        [Q(0), -(Q(1) + W), Q(1), Q(0), Q(1) - W * W],
    ]
    claim(
        "C5",
        "Fox matrix of (r1..r4) at rho (computed) equals matrix (15) of the paper entry by entry",
        all(D[i][j] == Dp[i][j] for i in range(4) for j in range(5)),
    )

    # the paper's closed-form entries for a relator x a x^-1 a^-n
    def closed(xv, av, n):
        return Q(1) - av**n, xv - sum((av**j for j in range(n)), Q(0))

    claim(
        "C5b",
        "closed-form entries (1-rho(a)^n, rho(x)-sum rho(a)^j) match the computed rows for r1,r2 (and the quoted sums 1, 1)",
        closed(rho["x1"], rho["a1"], 4) == (D[0][0], D[0][1])
        and closed(rho["x2"], rho["a2"], 3) == (D[1][2], D[1][3])
        and sum((rho["a1"] ** j for j in range(4)), Q(0)) == Q(1)
        and sum((rho["a2"] ** j for j in range(3)), Q(0)) == Q(1),
    )
    claim(
        "C5c",
        "r3 row: a2-entry = -(sum of 13 alternating terms) = -1; r4 row: a1-entry = -(1+w), x2-entry = rho(s) = 1, s-entry = 1-rho(x2)",
        sum((rho["a2"] ** j for j in range(13)), Q(0)) == Q(1)
        and D[2][3] == Q(-1)
        and D[3][1] == -(Q(1) + W)
        and D[3][2] == Q(1)
        and D[3][4] == Q(1) - rho["x2"],
    )

    # C6: rank 4, the minor on columns (x1,a1,x2,s), value -4(1-w^2)
    r, _ = rank_det(D)
    minor = [[D[i][j] for j in (0, 1, 2, 4)] for i in range(4)]
    _, dm = rank_det(minor)
    claim(
        "C6",
        "rank of the 4x5 twisted Fox matrix over Q(w) is 4 (d2: C^4 -> C^5 injective)",
        r == 4,
        f"rank={r}",
    )
    claim(
        "C6b",
        "minor on columns (x1,a1,x2,s) = -4(1-w^2) exactly, and it is nonzero",
        dm == Q(-4) * (Q(1) - W * W) and not dm.iszero(),
        f"det={dm}",
    )

    # C7: twisted Betti numbers of B: d1 row (1-rho(g)), d1 d2 = 0, h0,h1,h2 = 0
    d1 = [Q(1) - rho[g] for g in GENS]
    d1d2 = [sum((D[i][j] * d1[j] for j in range(5)), Q(0)) for i in range(4)]
    claim(
        "C7a",
        "d1 . d2 = 0 (the twisted chain complex of B is a complex)",
        all(x.iszero() for x in d1d2),
    )
    rk1, _ = rank_det([d1])
    h0 = 1 - rk1
    h1 = 5 - rk1 - r
    h2 = 4 - r
    claim(
        "C7b",
        "twisted Betti numbers of B: h0=h1=h2=0 (hence chi(B)=0 consistent: 1-5+4)",
        (h0, h1, h2) == (0, 0, 0),
        f"h0,h1,h2={h0},{h1},{h2}",
    )
    claim(
        "C7c",
        "cell count: X has 1+5+(4+2)+2 = 14 cells and chi(X) = 1-5+6-2 = 0",
        1 + 5 + 6 + 2 == 14 and 1 - 5 + 6 - 2 == 0,
    )

    # C8: Lemma 4.1 word identities in P, by explicit certificates in the free group F(x1,a1,x2,a2,s)
    R1, R2, R3, R4 = (RELS[k] for k in ("r1", "r2", "r3", "r4"))
    # t1 = x1^2 = a2^26 : x1^2 a2^-26 = R3 (a2^13 R3 a2^-13)
    c_t1 = [R3, conj(A("a2", 13), R3)]
    claim(
        "C8a",
        "t1 = x1^2 = a2^26 in P (certificate: R3 * a2^13 R3 a2^-13)",
        check_identity(A("x1", 2), A("a2", 26), c_t1),
    )
    # conjugation by t2 = x2^3 sends a2 to a2^27; [t2,a2] = a2^26
    ch2 = conj_chain("x2", "a2", 3, 3, R2)
    claim(
        "C8b",
        "x2^3 a2 x2^-3 = a2^27 in P (27 = 3^3), hence [t2,a2] = t2 a2 t2^-1 a2^-1 = a2^26",
        check_identity(A("x2", 3) + A("a2", 1) + A("x2", -3), A("a2", 27), ch2)
        and check_identity(
            A("x2", 3) + A("a2", 1) + A("x2", -3) + A("a2", -1), A("a2", 26), ch2
        ),
    )
    # conjugation by t1 = x1^2 sends a1 to a1^16
    ch1 = conj_chain("x1", "a1", 4, 2, R1)
    claim(
        "C8c",
        "x1^2 a1 x1^-2 = a1^16 in P (16 = 4^2), hence [t1,a1] = a1^15",
        check_identity(A("x1", 2) + A("a1", 1) + A("x1", -2), A("a1", 16), ch1)
        and check_identity(
            A("x1", 2) + A("a1", 1) + A("x1", -2) + A("a1", -1), A("a1", 15), ch1
        ),
    )
    # s t2 s^-1 = (s x2 s^-1)^3 = a1^15 : (R4 a1^5)^3 = prod_{j<3} a1^{5j} R4 a1^{-5j} * a1^15
    c_s = [conj(A("a1", 5 * j), R4) for j in range(3)]
    claim(
        "C8d",
        "s t2 s^-1 = s x2^3 s^-1 = a1^15 in P (certificate (R4 a1^5)^3)",
        check_identity([("s", 1)] + A("x2", 3) + [("s", -1)], A("a1", 15), c_s),
    )
    claim(
        "C8e",
        "arithmetic of Lemma 4.1: 4^2-1=15=3*5, 3^3-1=26=2*13, gcd(15,26)=1",
        4**2 - 1 == 15 == 3 * 5 and 3**3 - 1 == 26 == 2 * 13 and gcd(15, 26) == 1,
    )

    # C9: abelianisation sanity consequence of perfectness: t1,t2 already trivial in P^ab, so G^ab = P^ab
    import sympy as sp
    from sympy.matrices.normalforms import smith_normal_form

    # columns x1,a1,x2,a2,s ; rows = abelianised relators
    Pab = sp.Matrix(
        [[0, 1 - 4, 0, 0, 0], [0, 0, 0, 1 - 3, 0], [1, 0, 0, -13, 0], [0, -5, 1, 0, 0]]
    )
    Gab = sp.Matrix(Pab.tolist() + [[2, 0, 0, 0, 0], [0, 0, 3, 0, 0]])
    sP = [abs(x) for x in smith_normal_form(Pab, domain=sp.ZZ).diagonal() if x != 0]
    sG = [abs(x) for x in smith_normal_form(Gab, domain=sp.ZZ).diagonal() if x != 0]
    claim(
        "C9",
        "P^ab = Z + Z/6 (rank 1; torsion from Smith form) and G^ab = P^ab (adding 2*x1, 3*x2 changes nothing), as K perfect would imply",
        sorted(sP) == sorted(sG) and Pab.rank() == 4,
        f"SNF(P^ab) diag={sP}, SNF(G^ab) diag={sG}",
    )

    # C10: Lemma 2.1 mechanism on B along s -> Z: Fox matrix with monodromy rho(g) z^{phi(g)}; rank 4 over C(z), minor nonzero Laurent poly, equals -4(1-w^2) at z=1
    z, w = sp.symbols("z w")
    rz = {"x1": -1, "a1": w, "x2": w**2, "a2": -1, "s": z}
    rzi = {"x1": -1, "a1": w**2, "x2": w, "a2": -1, "s": 1 / z}  # inverses (w^3 = 1)

    def fox_sym(wd):
        row = {g: 0 for g in GENS}
        p = 1
        for g, e in wd:
            if e == 1:
                row[g] += p
                p = p * rz[g]
            else:
                row[g] -= p * rzi[g]
                p = p * rzi[g]
        return [sp.expand(row[g]) for g in GENS]

    Dz = sp.Matrix([fox_sym(RELS[k]) for k in ("r1", "r2", "r3", "r4")])
    cyc = sp.Poly(w**2 + w + 1, w)

    def red(e):
        e = sp.expand(sp.together(e))
        return sp.expand(sp.rem(sp.Poly(e, w), cyc).as_expr())

    dets = {}
    for cols in itertools.combinations(range(5), 4):
        dets[cols] = red(Dz[:, list(cols)].det())
    nz = [c for c, v in dets.items() if v != 0]
    dm_z = dets[(0, 1, 2, 4)]
    at1 = red(dm_z.subs(z, 1))
    target = red(-4 * (1 - w**2))
    claim(
        "C10",
        "with monodromy rho(g) z^phi(g) the Fox matrix has a nonzero 4x4 minor in Q(w)[z]  (rank 4 over C(z)); minor(x1,a1,x2,s) at z=1 equals -4(1-w^2)",
        len(nz) >= 1 and sp.simplify(at1 - target) == 0,
        f"minor(x1,a1,x2,s)={sp.factor(dm_z)}; nonzero minors: {len(nz)} of 5",
    )

    # C11: bounded solvable-quotient sanity: K perfect => in every solvable image of P, t1 and t2 die. Search affine groups AGL(1,p).
    def agl(p):
        return [(m, c) for m in range(1, p) for c in range(p)]

    def mul(
        g, h, p
    ):  # (m,c)(m',c') acts x -> m x + c ; g*h = g after h : x -> m(m'x+c')+c
        return ((g[0] * h[0]) % p, (g[0] * h[1] + g[1]) % p)

    def inv(g, p):
        mi = pow(g[0], -1, p)
        return (mi, (-mi * g[1]) % p)

    def pw(g, n, p):
        r = (1, 0)
        for _ in range(n):
            r = mul(r, g, p)
        return r

    def conj_(x, a, p):
        return mul(mul(x, a, p), inv(x, p), p)

    checked, bad, nonab = 0, 0, 0
    for p in (2, 3, 5, 7):
        G = agl(p)
        for a1, a2, x2, s in itertools.product(G, repeat=4):
            x1 = pw(a2, 13, p)
            if conj_(x1, a1, p) != pw(a1, 4, p):
                continue
            if conj_(x2, a2, p) != pw(a2, 3, p):
                continue
            if mul(mul(s, x2, p), inv(s, p), p) != pw(a1, 5, p):
                continue
            checked += 1
            gl = [x1, a1, x2, a2, s]
            nonab += any(mul(u, v, p) != mul(v, u, p) for u in gl for v in gl)
            if pw(x1, 2, p) != (1, 0) or pw(x2, 3, p) != (1, 0):
                bad += 1
    claim(
        "C11",
        "sanity (bounded, necessary condition only): in every homomorphism P -> AGL(1,p), p in {2,3,5,7}, t1 and t2 map to 1",
        bad == 0,
        f"{checked} homomorphisms examined ({nonab} with non-abelian image), {bad} violations",
    )

    elapsed = time.time() - t0
    ok = all(r[2] for r in RESULTS)
    print(
        f"\n{sum(r[2] for r in RESULTS)}/{len(RESULTS)} claims passed; RESULT: {'ALL PASS' if ok else 'FAILURES'}"
    )
    print(f"elapsed {elapsed:.1f}s", file=sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    sys.exit(main(perturb_matrix="--perturb-matrix" in sys.argv[1:]))
