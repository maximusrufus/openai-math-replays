"""Independent replay of the finite/explicit computations in 'Ambiently homeomorphic isolated hypersurface germs in C^4 with
multiplicities four and five' (OpenAI, 2026-09-27).

Written from the manuscript text only; no code from the authors is used.  Exact arithmetic: python-flint (fmpz_mpoly, fmpz_mat), Python integers.
Scope: Lemma 3.1 / 9.1 arithmetic on concrete instances, Lemma 2.1 plane data (expansions, Milnor numbers 51 and 6, rank/trace equations),
Proposition 3.2 initial forms, gradient estimate (3.1).  Homological and topological arguments are not finite computations and are not replayed.

  python replay_zariski.py             # prints every claim with PASS/FAIL; exit 0 only if all pass
  python replay_zariski.py --perturb   # adversarial: perturbs q1 (w^8 -> w^7); the Milnor-number claim must FAIL (exit 1)
  python replay_zariski.py --selftest  # also runs the built-in wrong-input controls (each must be rejected)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import re
import sys
import time
from fractions import Fraction as Fr
from math import comb, gcd

from flint import fmpz_mat, fmpz_mpoly_ctx

RES = []


SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
PFX = "Z"
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


# ---------------------------------------------------------------- arithmetic: Lemma 3.1 and 9.1
def lemma31_instance(r, n0):
    """returns dict of derived values, or None if n0 does not satisfy the hypotheses of Lemma 3.1"""
    if not (n0 % 2**r == 1 and n0 % 5 == 2 and n0 > 5 * 2**r):
        return None
    n1 = n0 + 5 * 2**r
    a = (n0 * n0 - 1) // 4
    b = (n1 * n1 - 1) // 4
    assert (n0 * n0 - 1) % 4 == 0 and (n1 * n1 - 1) % 4 == 0
    alpha = Fr(4 * b, a + b)
    beta = Fr(4 * a, a + b)
    return dict(r=r, n0=n0, n1=n1, a=a, b=b, alpha=alpha, beta=beta)


def check31(d):
    a, b, n0, n1, al, be = d["a"], d["b"], d["n0"], d["n1"], d["alpha"], d["beta"]
    ok = True
    ok &= n0 % 4 == 1 and n1 % 4 == 1
    ok &= gcd(n0, n1) == 1
    ok &= a > 0 and b > 0
    ok &= al + be == 4
    ok &= a * al == b * be
    ok &= n0 < n1 < 2 * n0
    ok &= n1 == n0 + 5 * 2 ** d["r"]  # the shift of Lemma 3.1, asserted independently of how the instance was built
    ok &= a == (n0 * n0 - 1) // 4 and b == (n1 * n1 - 1) // 4
    return ok


def strict_bounds(d):
    """the two strict lower bounds of Lemma 3.1; they hold only for n0 'sufficiently large', so a failure means: n0 too small, skip it"""
    a, b, al, be = d["a"], d["b"], d["alpha"], d["beta"]
    return (a - 2) * al / 4 > 4 and (b - 2) * be / 4 > 4 and a > 4 and b > 4


def check91(d, qlist=(1, 4)):
    a, b, r = d["a"], d["b"], d["r"]
    g = gcd(a, b)
    m = 4 * a * b + a + b
    ok = a % 5 == 2 and b % 5 == 2 and g % 5 != 0 and g >= 2 ** (r - 1) > 5
    found = {}
    for q in qlist:
        # brute force over j (no use of the paper's solution recipe): need i>0, j>0, i = j (mod 4), (b i + a j)/m = q/5 exactly
        target = Fr(q, 5) * m  # b i + a j = q m / 5, integer because 5 | m
        if target.denominator != 1:
            ok = False
            continue
        T = int(target)
        sol = None
        # i ranges over positive integers with b i < T, j = (T - b i)/a must be positive integer with j = i (mod 4).
        # iterate over j instead (cheaper when a large): j in 1..: b i = T - a j.
        j = 1
        while a * j < T and sol is None:
            rem = T - a * j
            if rem % b == 0:
                i = rem // b
                if i > 0 and (i - j) % 4 == 0:
                    sol = (i, j)
            j += 1
            if j > 8 * b + 8:
                break
        found[q] = sol
        ok &= sol is not None
        if sol is not None:  # verify the found pair directly against the stated conditions
            ok &= sol[0] > 0 and sol[1] > 0 and (sol[0] - sol[1]) % 4 == 0 and 5 * (b * sol[0] + a * sol[1]) == q * m
    return ok, found


def arithmetic(selftest):
    count, bad = 0, 0
    for r in range(4, 12):
        n0 = None
        # first 3 admissible n0 per r that also give the strict bounds
        cand = 0
        n = 5 * 2**r + 1
        while cand < 3:
            n += 1
            d = lemma31_instance(r, n)
            if d is None:
                continue
            if not strict_bounds(d):
                continue  # the lemma requires n0 sufficiently large for the two strict bounds: skip n0 below that threshold
            cand += 1
            okk, f = check91(d)
            count += 1
            bad += (not okk) or (not check31(d))  # every other property must hold for every admissible n0
    chk(
        f"Lemma 3.1 + 9.1: {count} concrete (r, n0) instances (r=4..11, 3 each; r>=4 needed for 2^(r-1)>5), all properties hold (failures {bad})",
        bad == 0 and count == 24,
    )
    if selftest:
        # deliberately wrong inputs
        d = lemma31_instance(4, 17 + 80 * 20)
        chk(
            "WRONG INPUT: n0 not 2 mod 5 is rejected by hypothesis test",
            lemma31_instance(4, 16 * 70 + 1 if (16 * 70 + 1) % 5 != 2 else 16 * 71 + 1)
            is None,
        )
        d = lemma31_instance(4, 1217)
        dd = dict(d)
        dd["n1"] = d["n0"] + 3 * 2**4  # wrong shift (5 | n1)
        dd["b"] = (dd["n1"] ** 2 - 1) // 4
        dd["alpha"] = Fr(4 * dd["b"], dd["a"] + dd["b"])
        dd["beta"] = Fr(4 * dd["a"], dd["a"] + dd["b"])
        chk(
            "WRONG INPUT: n1 = n0 + 3*2^r gives mod-5 failure in Lemma 9.1 checker (must FAIL)",
            not check91(dd)[0] or not check31(dd),
        )
        d3 = lemma31_instance(4, 1217)
        d3["alpha"] = d3["alpha"] + 1
        chk("WRONG INPUT: alpha+beta != 4 is rejected by check31", not check31(d3))


# ---------------------------------------------------------------- series / polynomial computations
ctx2 = fmpz_mpoly_ctx.get(("x", "w"), "lex")
X, Wv = ctx2.gens()


def trunc(p, N):
    return ctx2.from_dict({e: c for e, c in p.to_dict().items() if sum(e) < N})


def mult_trunc(p, q, N):
    return trunc(p * q, N)


def milnor_number(f, N):
    """dim_Q O/(f_x, f_y, m^N) via exact rank of truncated multiples; equals mu(f) when m^N lies in the Jacobian ideal."""
    fx = trunc(f.derivative(0), N)
    fy = trunc(f.derivative(1), N)
    mons = [(i, j) for d in range(N) for i in range(d + 1) for j in [d - i]]
    idx = {m: k for k, m in enumerate(mons)}
    rows = []
    for g in (fx, fy):
        gd = g.to_dict()
        if not gd:
            continue
        omin = min(sum(e) for e in gd)
        for i, j in mons:
            if i + j + omin >= N:
                continue
            row = [0] * len(mons)
            nz = False
            for e, c in gd.items():
                ee = (e[0] + i, e[1] + j)
                if sum(ee) < N:
                    row[idx[ee]] += int(c)
                    nz = True
            if nz:
                rows.append(row)
    M = fmpz_mat(rows)
    return len(mons) - M.rank(), M


def build_plane(variant="paper"):
    u = (
        X + Wv
    )  # u = x + w   (first graph uses (u,w) as coordinates; we name them (x,w) for ctx reuse)
    return u


def hfun(X_, T_):
    """h(X, T) = T (X - T^2)(X - 2 T^2)"""
    return T_ * (X_ - T_ * T_) * (X_ - 2 * T_ * T_)


def plane_polys(w8=8, shift=0):
    # first graph: coordinates (u, w); we reuse ctx variables (X := u, Wv := w)
    u, w = X, Wv
    q1 = (u - w - w**w8) * (u - w - 2 * w**w8) * (u - w**2)
    q2 = q1 - u * u
    d0 = u * q1 * q2  # = h(phi(u,w), u), phi = u^2 + q1  (verified symbolically below)
    # second graph: coordinates (x, w) with u = x + w
    x = X
    uu = x + w
    q1b = (uu - w - w**w8) * (uu - w - 2 * w**w8) * (uu - w**2)
    q2b = q1b - uu * uu
    psi = uu * (uu - w**2) * q2b
    d1 = psi * (x - psi * psi) * (x - 2 * psi * psi)
    return q1, q2, d0, psi, d1


def plane_checks(selftest):
    N_list = (55, 60)
    u, w = X, Wv
    q1, q2, d0, psi, d1 = plane_polys()
    # h o pi_0 equals u q1 q2 (definition check): h(P,T) = T (P - T^2)(P - 2 T^2), P = phi = u^2 + q1, T = u
    phi = u * u + q1
    h0 = u * (phi - u * u) * (phi - 2 * u * u)
    chk("h(pi_0(u,w)) = u q1 q2 as polynomials", h0 == d0)
    chk("d1 = h(pi_1(x,w)) = psi (x - psi^2)(x - 2 psi^2) with h(X,T) = T (X - T^2)(X - 2 T^2) at (X,T) = (x,psi)", d1 == hfun(X, psi))
    chk("d0 = h(pi_0(u,w)) with h(X,T) = T (X - T^2)(X - 2 T^2) at (X,T) = (phi,u)", d0 == hfun(phi, u))
    # phi(0,w) = -w^4 (1+w^7)(1+2w^7)
    ev = lambda p, xv: (
        ctx2.from_dict({(0, e[1]): c for e, c in p.to_dict().items() if e[0] == 0})
        if xv == 0
        else None
    )
    p0 = ev(phi, 0)
    chk("phi(0,w) = -w^4 (1+w^7)(1+2 w^7)", p0 == -(w**4) * (1 + w**7) * (1 + 2 * w**7))
    # psi(0,w) = -w^4 (1 - w) + O(w^19), exact: -w^4(1-w) + 2 w^19 (1-w)^2
    ps0 = ev(psi, 0)
    chk(
        "psi(0,w) = -w^4(1-w) + 2 w^19 (1-w)^2 exactly (so O(w^19))",
        ps0 == -(w**4) * (1 - w) + 2 * w**19 * (1 - w) ** 2,
    )
    # q2(u,0) = u^2 (u-1)
    q2w0 = ctx2.from_dict({(e[0], 0): c for e, c in q2.to_dict().items() if e[1] == 0})
    chk("q2(u,0) = u^2 (u-1)", q2w0 == u**2 * (u - 1))
    # q2(v w^2, w) = w^4 (v - 1 - v^2) + O(w^5): substitute u = v w^2 with v as ctx variable x
    v = X
    q1s = (v * w**2 - w - w**8) * (v * w**2 - w - 2 * w**8) * (v * w**2 - w**2)
    q2s = q1s - (v * w**2) ** 2
    lowdeg = {e: c for e, c in q2s.to_dict().items() if e[1] <= 4}
    chk(
        "q2(v w^2, w): w-degree<=4 part equals w^4 (v-1-v^2)",
        ctx2.from_dict(lowdeg) == w**4 * (v - 1 - v * v),
    )
    # quadratic roots simple: disc of v^2 - v + 1 = -3 != 0
    chk(
        "v^2 - v + 1 has discriminant -3 (simple roots, f = (1 +- i sqrt3)/2)",
        1 - 4 == -3,
    )
    # delta and mu arithmetic from the stated intersection multiplicities
    delta = 8 + comb(4, 2) * 2 + 2 * 4
    chk(
        f"delta(d_i) = 8 + C(4,2)*2 + 2*4 = {delta} = 28; mu = 2*delta - 6 + 1 = {2 * delta - 6 + 1} = 51",
        delta == 28 and 2 * delta - 6 + 1 == 51,
    )
    # h: three smooth branches, pairwise intersections 1,1,2 -> delta = 4, mu = 2*4-3+1 = 6; verify via Milnor number of h
    hX = (
        X * X
    )  # placeholder (kept to reuse ctx): real check below uses variables (X,T) = (x,w)
    Xv, T = X, Wv
    hh = hfun(Xv, T)
    mu_h, _ = milnor_number(hh, 12)
    chk(f"mu(h) = 6 by exact truncated Jacobian algebra (got {mu_h})", mu_h == 6)
    # gradient formulas (3.1): h_X = T(2X - 3T^2), h_T = X^2 - 9T^2 X + 10 T^4
    hXd, hTd = hh.derivative(0), hh.derivative(1)
    chk(
        "h_X = T(2X - 3T^2) and h_T = X^2 - 9 T^2 X + 10 T^4",
        hXd == T * (2 * Xv - 3 * T * T) and hTd == Xv * Xv - 9 * T * T * Xv + 10 * T**4,
    )
    c0 = Fr(3, 2)
    chk(
        "2c-3 and c^2-9c+10 have no common root (c=3/2 gives c^2-9c+10 = -5/4; and 2c-3 is linear)",
        c0 * c0 - 9 * c0 + 10 == Fr(-5, 4),
    )
    # representation ranks: a + 4 b = 6, a - b = 1 -> (a,b) = (2,1)
    sol = [
        (a_, b_)
        for a_ in range(7)
        for b_ in range(2)
        if a_ + 4 * b_ == 6 and a_ - b_ == 1
    ]
    chk(
        f"Lefschetz rank/trace equations a+4b=6, a-b=1 have the unique solution {sol} -> dim ker(M-I)=2, dim ker Phi5 = 4",
        sol == [(2, 1)],
    )
    # Milnor numbers by exact truncated Jacobian algebra
    t0 = time.time()
    mus = {}
    for nm, f in (("d0", d0), ("d1", d1)):
        for N in N_list:
            mu, _ = milnor_number(f, N)
            mus[(nm, N)] = mu
            log("   ", nm, "N=", N, "dim O/(J+m^N) =", mu, f"({time.time() - t0:.0f}s)")
    chk(
        f"mu(d0) = 51 and mu(d1) = 51 by exact truncated Jacobian algebra, N=55,60 stable ({mus})",
        all(v == 51 for v in mus.values()),
    )
    if selftest:
        chk("WRONG INPUT: claim mu(d0) = 52 is rejected", mus[("d0", 60)] != 52)
        # perturbed germ (w^8 -> w^7 in q1) must NOT have mu = 51
        q1p, q2p, d0p, psip, d1p = plane_polys(w8=7)
        mup, _ = milnor_number(d0p, 60)
        chk(f"WRONG INPUT: perturbed d0 (w^8 -> w^7) has mu = {mup} != 51", mup != 51)
        chk("WRONG INPUT: h Milnor number is not 5", mu_h != 5)
        chk(
            "WRONG INPUT: wrong psi(0,w) = -w^4 (1+w) is rejected",
            ps0 != -(w**4) * (1 + w) + 2 * w**19 * (1 - w) ** 2,
        )


def initial_forms():
    # In weighted terms: f0: X = T^2 + Y + O>=3 with Y = s^-4 yz of degree 2: h = T (X-T^2)(X-2T^2) -> T*Y*(Y-T^2) lowest degree 5
    ctx = fmpz_mpoly_ctx.get(("T", "Y", "E"), "lex")  # E marks a degree-3 error term
    T, Y, E = ctx.gens()
    Xs = T * T + Y + E
    hh = T * (Xs - T * T) * (Xs - 2 * T * T)
    # weights: T:1, Y:2, E:3 (E stands for terms of order >=3 and in particular carries degree 3); lowest weighted-degree part
    wts = (1, 2, 3)
    d = hh.to_dict()
    mind = min(sum(e_ * w_ for e_, w_ in zip(e, wts)) for e in d)
    low = ctx.from_dict(
        {e: c for e, c in d.items() if sum(e_ * w_ for e_, w_ in zip(e, wts)) == mind}
    )
    chk(
        f"in(f0,s): lowest-degree part of h(T^2+Y+O3, T) with deg T=1,deg Y=2 has degree {mind} = 5 (no E term enters) and equals T*Y*(Y-T^2)",
        mind == 5 and low == T * Y * (Y - T * T),
    )
    # f1: T = Y + O>=3 where Y degree 2; h = T (X - T^2)(X - 2T^2), X degree 1: T^2 degree 4 -> X^2 * Y lowest, degree 4
    ctx = fmpz_mpoly_ctx.get(("X", "Y", "E"), "lex")
    Xv, Y, E = ctx.gens()
    Ts = Y + E
    hh = Ts * (Xv - Ts * Ts) * (Xv - 2 * Ts * Ts)
    wts = (1, 2, 3)
    d = hh.to_dict()
    mind = min(sum(e_ * w_ for e_, w_ in zip(e, wts)) for e in d)
    low = ctx.from_dict(
        {e: c for e, c in d.items() if sum(e_ * w_ for e_, w_ in zip(e, wts)) == mind}
    )
    chk(
        f"in(f1,s): lowest-degree part has degree {mind} = 4 and equals Y*X^2",
        mind == 4 and low == Y * Xv * Xv,
    )
    # wrong claim must fail
    chk("WRONG INPUT: claimed order 5 for f1 is rejected", mind != 5)
    chk(
        "WRONG INPUT: wrong initial form Y*X*(X-T) is rejected",
        low != Y * Xv * (Xv - Y),
    )


def run():
    t0 = time.time()
    if PERTURB:  # perturb the INPUT (q1: w^8 -> w^7) and run the unchanged main checks
        global plane_polys
        _orig = plane_polys
        plane_polys = lambda w8=7, shift=0: _orig(w8=7, shift=shift)
    arithmetic(SELFTEST)
    initial_forms()
    plane_checks(SELFTEST)
    log(f"elapsed {time.time() - t0:.1f}s")
    return finish()


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    sys.exit(run())
