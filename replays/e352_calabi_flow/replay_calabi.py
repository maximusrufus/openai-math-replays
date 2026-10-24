"""Independent replay of the finite computations in 'A finite-time singularity of Calabi flow on projective space' (OpenAI, 2026-09-24),
Section 3 and Appendix A.

Written from the manuscript text only; no code from the authors (including the printed certificate) is used.  mpmath at 40 digits (Taylor-series ODE
solver, tolerance 1e-34), sympy, exact fractions.  The printed matrix of (3.6) and the squares in the bounds of Lemma 3.3 were read from the rendered pages.
Replayed: (3.5) coefficients from Q(D); (3.7); Lemma 3.2 constants; Lemma 3.3 root/projection/convolution bounds; Lemma 3.4 tail residual (sampled);
Proposition A.1 shooting estimate (A.1) by high-precision integration at 8 SAMPLE POINTS of the parameter square Z (centre, four corners, three interior points).
NOT replayed: the rigorous integer-interval enclosure over ALL z* in Z.

  python replay_calabi.py             # prints every claim with PASS/FAIL; exit 0 only if all pass (about 2-3 min)
  python replay_calabi.py --perturb   # adversarial: k shifted by 1e-9 at the centre; the bound |q| <= .25 R must FAIL (exit 1)
  python replay_calabi.py --selftest  # also runs the built-in wrong-input controls (each must be rejected)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import re
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
RES = []


SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
PFX = "CF"
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


B_ = Fr(751, 4)
C_ = Fr(-1215, 4)
R_ = Fr(5, 10**10)
S0 = 56
K0 = Fr("-.094617116666041")
D0 = Fr("0.400682721856639")
MAT = [[25190, -63731], [721855, 359157]]


def mpq(x):
    return mp.mpf(x.numerator) / x.denominator


def kd(zs, Rr=1e-8):
    """(k, d) = (k0, d0) + 1e-8 * MAT * z"""
    return (
        K0 + Fr(1, 10**8) * (MAT[0][0] * zs[0] + MAT[0][1] * zs[1]),
        D0 + Fr(1, 10**8) * (MAT[1][0] * zs[0] + MAT[1][1] * zs[1]),
    )


# ------------------------------------------------------------------ origin series (3.8)-(3.10), A.4
def origin_coeffs(k, d, N=61, prec=None):
    mpf = mp.mpf
    z = [mpf(0)] * (N + 1)
    z[0] = 1 + d
    z[1] = k
    bb = [1 / z[0]]
    for i in range(2, N + 1):
        for n in range(len(bb), i - 1):
            bb.append(-(sum(z[t] * bb[n - t] for t in range(1, n + 1))) / z[0])
        z[i] = (bb[i - 2] - (1 if i == 2 else 0)) / (
            mpf(27) / 8 * (i - 1) * (i + 9) * (i + 10)
        )
    return z  # z[0] = 1 + d  (so p = z - 1 at r)


def origin_jet(k, d, s, N=61):
    z = origin_coeffs(k, d, N)
    r = s ** mp.mpf(1.5)
    p = d + sum(z[i] * r**i for i in range(1, N + 1))
    p1 = sum(z[i] * (mp.mpf(3 * i) / 2) * r**i / s for i in range(1, N + 1))
    p2 = sum(
        z[i] * (mp.mpf(3 * i) / 2) * (mp.mpf(3 * i) / 2 - 1) * r**i / s**2
        for i in range(1, N + 1)
    )
    return p, p1, p2, z


def pa_coeffs(d):
    a = mpq(C_) * d
    b = a * a + 567 * a
    d3 = 2 * a * b - a**3 + mp.mpf(2025) / 4 * b
    d4 = b * b + 2 * a * d3 - 3 * a * a * b + a**4 + mp.mpf(567) / 2 * d3
    return [a, b, d3, d4]


def pa_jet(d, s):
    c = pa_coeffs(d)
    p = sum(c[j] * s ** (-3 * (j + 1)) for j in range(4))
    p1 = sum(c[j] * (-3 * (j + 1)) * s ** (-3 * (j + 1) - 1) for j in range(4))
    p2 = sum(
        c[j] * (-3 * (j + 1)) * (-3 * (j + 1) - 1) * s ** (-3 * (j + 1) - 2)
        for j in range(4)
    )
    return p, p1, p2


def Ed_residual(d, s):
    c = pa_coeffs(d)
    p, p1, p2 = pa_jet(d, s)
    p3 = sum(
        c[j]
        * (-3 * (j + 1))
        * (-3 * (j + 1) - 1)
        * (-3 * (j + 1) - 2)
        * s ** (-3 * (j + 1) - 3)
        for j in range(4)
    )
    return (
        p3
        + 30 / s * p2
        + mpq(B_) / s**2 * p1
        + mpq(C_) / s**3 * (p - d)
        + p / (1 + p)
    )


def shoot(k, d, upto=S0, record_min=True):
    p, p1, p2, z = origin_jet(k, d, mp.mpf(1))
    Bm, Cm = mpq(B_), mpq(C_)

    def F(s, y):
        return [
            y[1],
            y[2],
            -30 / s * y[2]
            - Bm / s**2 * y[1]
            - Cm / s**3 * (y[0] - d)
            - y[0] / (1 + y[0]),
        ]

    f = mp.odefun(F, 1, [p, p1, p2], tol=mp.mpf(10) ** -34)
    mn = mp.mpf(10)
    if record_min:
        for s in range(1, upto + 1):
            mn = min(mn, 1 + f(s)[0])
    return f(upto), mn


def D0poly(lam, S=S0):
    return (
        lam**3
        + mp.mpf(15) / S * lam**2
        + mpq(B_) / (2 * S**2) * lam
        + 1
        + mpq(C_) / (2 * S**3)
    )


def run():
    selftest = SELFTEST
    t0 = time.time()
    if PERTURB:  # perturb the INPUT (k0 shifted by 1e-9), then run the unchanged main checks
        global K0
        K0 = K0 + Fr(1, 10**9)
    # --- (3.5) coefficients from Q(D) = (D-1)(D+9)(D+10)
    Dv, s_ = sp.symbols("D s")
    Q = sp.expand((Dv - 1) * (Dv + 9) * (Dv + 10))
    chk(
        f"Q(D) = D^3 + 18 D^2 + 71 D - 90 (got {Q})",
        sp.Poly(Q, Dv).all_coeffs() == [1, 18, 71, -90],
    )
    th = sp.symbols("th")
    # D = (2/3) theta with theta = s d/ds:  Q((2/3) th) = (8/27)(th-3/2)(th+27/2)(th+15); expand in s-derivatives
    P3 = sp.expand((th - sp.Rational(3, 2)) * (th + sp.Rational(27, 2)) * (th + 15))
    c3, c2, c1, c0 = sp.Poly(P3, th).all_coeffs()
    # th^3 = s^3 d^3 + 3 s^2 d^2 + s d,  th^2 = s^2 d^2 + s d,  th = s d
    coef_s3, coef_s2, coef_s1 = c3, 3 * c3 + c2, c3 + c2 + c1
    chk(
        f"(3.5) coefficients: s^3 p''' : {coef_s3}, s^2 p'' : {coef_s2} (=30), s p' : {coef_s1} (=751/4), constant {c0} (=-1215/4)",
        (coef_s3, coef_s2, coef_s1, c0)
        == (1, 30, sp.Rational(751, 4), sp.Rational(-1215, 4)),
    )
    # the nonlinear side: Q(D)(p-d) = (y^2/(2c^2))(1/(1+p) - 1), y^2 = c^2 (16/27) s^3 -> (8/27) s^3 (1/(1+p) - 1); with Q((2/3)th) = (8/27) P3 -> s^3 p/(1+p) sign
    chk(
        "y = c (16/27)^(1/2) s^(3/2) from s = (27/16)^(1/3) (y/c)^(2/3): y^2/(2 c^2) = (8/27) s^3",
        sp.simplify((sp.Rational(16, 27) * s_**3) / 2 - sp.Rational(8, 27) * s_**3)
        == 0,
    )
    # --- (3.7)
    chk(
        "transcription of (3.6): k0 = -.094617116666041, d0 = .400682721856639, matrix [[25190, -63731], [721855, 359157]], R = 5e-10, S0 = 56",
        K0 == Fr(-94617116666041, 10**15)
        and D0 == Fr(400682721856639, 10**15)
        and MAT == [[25190, -63731], [721855, 359157]]
        and R_ == Fr(5, 10**10)
        and S0 == 56,
    )
    kmax = Fr(1, 10**8) * sum(abs(x) for x in MAT[0]) * R_  # row sum of |entries| times R
    dmax = Fr(1, 10**8) * sum(abs(x) for x in MAT[1]) * R_
    chk(
        f"(3.7) |k-k0| <= {float(kmax):.5e} (stated 4.44605e-13) and |d-d0| <= {float(dmax):.5e} (stated 5.40506e-12) < .011 R; exact values lie within half a unit of the last printed digit",
        abs(kmax - Fr("4.44605e-13")) <= Fr(1, 2 * 10**18)  # exact decimal interval [v - 1/2 ulp, v + 1/2 ulp], ulp = 1e-18
        and abs(dmax - Fr("5.40506e-12")) <= Fr(1, 2 * 10**17)  # ulp = 1e-17
        and dmax < Fr(11, 1000) * R_,
    )
    # --- Lemma 3.3
    al = Fr(15, 56)
    be = B_ / (2 * 56**2)
    H0 = 1 + C_ / (2 * 56**3)
    chk(
        f"alpha^2 - 3 beta = {al**2 - 3 * be} = -453/25088 < 0 (D0 increasing on R)",
        al**2 - 3 * be == Fr(-453, 25088) < 0,
    )

    def D0f(x):
        return x**3 + al * x**2 + be * x + H0

    v1, v2 = D0f(Fr(-108659, 100000)), D0f(Fr(-108655, 100000))
    chk(
        f"D0(-1.08659) = {v1} = -17566501197647/343000000000000000 < 0 and D0(-1.08655) = {v2} = 187622583137/2744000000000000 > 0",
        v1 == Fr(-17566501197647, 343000000000000000)
        and v2 == Fr(187622583137, 2744000000000000)
        and v1 < 0 < v2,
    )
    mp.mp.dps = 40
    roots = mp.polyroots(
        [
            1,
            mp.mpf(al.numerator) / al.denominator,
            mp.mpf(be.numerator) / be.denominator,
            mp.mpf(H0.numerator) / H0.denominator,
        ],
        maxsteps=200,
        extraprec=100,
    )
    lam0 = [r for r in roots if abs(mp.im(r)) < 1e-30][0]
    lam1 = [r for r in roots if mp.im(r) > 1e-3][0]
    lam2 = mp.conj(lam1)
    chk(
        f"(3.14) lambda0 = {mp.nstr(mp.re(lam0), 8)} in (-1.08659, -1.08655); Re lambda1 = {mp.nstr(mp.re(lam1), 8)} in (.40933, .40938); Im lambda1 = {mp.nstr(mp.im(lam1), 8)} in (.8671, .86721)",
        -1.08659 < mp.re(lam0) < -1.08655
        and 0.40933 < mp.re(lam1) < 0.40938
        and 0.8671 < mp.im(lam1) < 0.86721,
    )
    um, up = Fr(114617, 280000), Fr(573113, 1400000)
    chk(
        "u- = 114617/280000 > .40933 and u+ = 573113/1400000 < .40938 (rational bounds on Re lambda1)",
        um > Fr(40933, 100000) and up < Fr(40938, 100000),
    )
    chk(
        "H0/1.08655 < .9592 (|lambda1| < .959)",
        H0 / Fr(108655, 100000) < Fr(9592, 10000),
    )
    # projection matrix Gamma: rows from (3.17): complex row i: (lj lk, -lj - lk, 1)/D0'(li)
    lam = [lam0, lam1, lam2]

    def dD0(x):
        return (
            3 * x**2
            + 2 * mp.mpf(al.numerator) / al.denominator * x
            + mp.mpf(be.numerator) / be.denominator
        )

    def row(i):
        j, k = [x for x in (0, 1, 2) if x != i]
        dd = dD0(lam[i])
        return [lam[j] * lam[k] / dd, (-lam[j] - lam[k]) / dd, 1 / dd]

    r0, r1 = row(0), row(1)
    Gam = [[mp.re(x) for x in r0], [mp.re(x) for x in r1], [mp.im(x) for x in r1]]
    Lam = [
        [mp.mpf(x) * mp.mpf(10) ** -5 for x in rr]
        for rr in (
            [30756, -27384, 33447],
            [34622, 13692, -16724],
            [-2925, -34040, -28850],
        )
    ]
    log("   Gamma =", [[mp.nstr(x, 6) for x in rr] for rr in Gam])
    maxdiff = max(abs(Gam[i][j] - Lam[i][j]) for i in range(3) for j in range(3))
    chk(
        f"every entry of Gamma - Lambda has modulus < .001 (max {mp.nstr(maxdiff, 5)})",
        maxdiff < 0.001,
    )
    rs = [sum(abs(x) for x in rr) for rr in Gam]
    chk(
        f"absolute row sums of Gamma = {[mp.nstr(x, 5) for x in rs]} < (.92, .68, .68)",
        rs[0] < 0.92 and rs[1] < 0.68 and rs[2] < 0.68,
    )
    Ms = [
        sum(
            abs(lam[j]) ** i / (abs(mp.re(lam[j])) * abs(dD0(lam[j]))) for j in range(3)
        )
        for i in range(3)
    ]
    chk(
        f"M_i = {[mp.nstr(x, 6) for x in Ms]} < (1.96, 1.92, 1.88)",
        Ms[0] < 1.96 and Ms[1] < 1.92 and Ms[2] < 1.88,
    )
    inv = 1 / (mp.re(lam1) * abs(dD0(lam1)))
    chk(f"1/(Re lambda1 |D0'(lambda1)|) = {mp.nstr(inv, 6)} < .82", inv < 0.82)
    # printed bounds: .8671^2 < ... < .86721^2 and E > (1.08655 + u-)^2 + .8671^2
    lo_v2 = H0 / Fr(108659, 100000) - up**2
    hi_v2 = H0 / Fr(108655, 100000) - um**2
    chk(f".8671^2 = {float(Fr(8671, 10000) ** 2):.6f} < H0/1.08659 - u+^2 = {float(lo_v2):.6f} <= v^2 <= H0/1.08655 - u-^2 = {float(hi_v2):.6f} < .86721^2 = {float(Fr(86721, 100000) ** 2):.6f}", Fr(8671, 10000) ** 2 < lo_v2 <= hi_v2 < Fr(86721, 100000) ** 2)
    vsq_num = mp.im(lam1) ** 2
    chk(f"numerical v^2 = {mp.nstr(vsq_num, 8)} lies in [{float(lo_v2):.6f}, {float(hi_v2):.6f}]", lo_v2 <= Fr(mp.nstr(vsq_num, 20)) <= hi_v2)
    E_exact = (Fr(108655, 100000) + um) ** 2 + Fr(8671, 10000) ** 2
    chk(f"E > (1.08655 + u-)^2 + .8671^2 = {E_exact} = 46876434629/15680000000 > 2.9895 > 1.729^2 = {float(Fr(1729, 1000) ** 2):.6f}", E_exact == Fr(46876434629, 15680000000) > Fr(29895, 10000) > Fr(1729, 1000) ** 2)
    E_num = (-lam0 + mp.re(lam1)) ** 2 + mp.im(lam1) ** 2
    chk(f"numerical E = (r0+u)^2 + v^2 = {mp.nstr(E_num, 8)} > 2.9895, and D0'(-r0) = {mp.nstr(dD0(lam0), 8)} equals E", E_num > 2.9895 and abs(dD0(lam0) - E_num) < 1e-25)
    chk(
        "1/(2*.40933*.8671*1.729) = 500000000000/613674044347 < .82",
        Fr(1, 1) / (2 * Fr(40933, 100000) * Fr(8671, 10000) * Fr(1729, 1000))
        < Fr(82, 100),
    )
    # --- Lemma 3.2 constants at the centre and corners
    ok32 = True
    for zs in [(0, 0)] + [(sx * R_, sy * R_) for sx in (-1, 1) for sy in (-1, 1)]:
        k, d = kd(zs)
        k, d = mp.mpf(k.numerator) / k.denominator, mp.mpf(d.numerator) / d.denominator
        z = origin_coeffs(k, d)
        nrm = 3 * abs(k) + sum(abs(z[i]) * 3**i for i in range(2, 62))
        ok32 &= nrm < 0.8
    chk(
        f"(3.9) 3|k| + sum_(i>=2) |a_i| 3^i < .8 at centre and 4 corners (centre: {mp.nstr(3 * abs(mp.mpf(K0.numerator) / K0.denominator) + sum(abs(origin_coeffs(mp.mpf(K0.numerator) / K0.denominator, mp.mpf(D0.numerator) / D0.denominator)[i]) * 3**i for i in range(2, 62)), 6)})",
        ok32,
    )
    kc, dc = (
        mp.mpf(K0.numerator) / K0.denominator,
        mp.mpf(D0.numerator) / D0.denominator,
    )
    p, p1, p2, z = origin_jet(kc, dc, mp.mpf(1))
    chk(
        f"jet at s=1: |p|,|p'|,|p''| = {mp.nstr(abs(p), 5)}, {mp.nstr(abs(p1), 5)}, {mp.nstr(abs(p2), 5)} all <= 2 (and < .768, .4, .534 as in A.2)",
        abs(p) < 0.768 and abs(p1) <= 0.4 and abs(p2) < 0.534,
    )
    pos = min(1 + origin_jet(kc, dc, mp.mpf(i) / 50)[0] for i in range(1, 51))
    chk(f"1 + p > 0 on 0 < s <= 1 (grid min {mp.nstr(pos, 6)})", pos > 0)
    # the origin series solves (3.5): residual of the ODE at s0 = 1/2, with the derivatives of the summed series taken by numerical
    # differentiation (independent of the coefficient recurrence and of the analytic jets)
    s0 = mp.mpf(1) / 2

    def ps(x):
        return origin_jet(kc, dc, x)[0]

    resid = abs(
        mp.diff(ps, s0, 3)
        + 30 / s0 * mp.diff(ps, s0, 2)
        + mpq(B_) / s0**2 * mp.diff(ps, s0, 1)
        + mpq(C_) / s0**3 * (ps(s0) - dc)
        + ps(s0) / (1 + ps(s0))
    )
    chk(f"the origin series solves (3.5): residual of p(3) + 30 p''/s + B p'/s^2 + C (p - d)/s^3 + p/(1+p) at s = 1/2 is {mp.nstr(resid, 3)} < 1e-20 (numerical derivatives of the summed series)", resid < mp.mpf(10) ** -20)
    # --- Lemma 3.4
    sup = max(
        abs(Ed_residual(dd, mp.mpf(s)))
        for dd in (
            dc - mp.mpf(dmax.numerator) / dmax.denominator,
            dc,
            dc + mp.mpf(dmax.numerator) / dmax.denominator,
        )
        for s in (56, 56.5, 57, 60, 70, 100, 200, 1000)
    )
    chk(
        f"Lemma 3.4: sup_(s>=56) |E_d(p_a)| sampled = {mp.nstr(sup, 4)} < 3e-13",
        sup < 3e-13,
    )
    Pa0 = pa_jet(dc, mp.mpf(56))
    dev = 0
    for dd in (
        dc - mp.mpf(dmax.numerator) / dmax.denominator,
        dc + mp.mpf(dmax.numerator) / dmax.denominator,
    ):
        Pa = pa_jet(dd, mp.mpf(56))
        dev = max(dev, max(abs(Pa[i] - Pa0[i]) for i in range(3)) / abs(dd - dc))
    chk(
        f"(3.19) ||Pa(d) - Pa(d0)|| / |d - d0| = {mp.nstr(dev, 5)} <= .01 (at the extreme d of the square)",
        dev <= 0.01,
    )
    # --- shooting (A.1)
    samples = (
        [(0, 0)]
        + [(sx * R_, sy * R_) for sx in (-1, 1) for sy in (-1, 1)]
        + [
            (R_ / 2, -R_ / 3),
            (-R_ / 5, R_ * Fr(7, 10)),
            (Fr(3, 10) * R_, Fr(1, 10) * R_),
        ]
    )
    Gm = Gam
    worst_q, worst_xi, minpos = 0, 0, 10
    for zs in samples:
        k, d = kd(zs)
        kk, dd = (
            mp.mpf(k.numerator) / k.denominator,
            mp.mpf(d.numerator) / d.denominator,
        )
        (pj, p1j, p2j), mn = shoot(kk, dd)
        pa = pa_jet(dd, mp.mpf(S0))
        delta = [pj - pa[0], p1j - pa[1], p2j - pa[2]]
        q = sum(Gm[0][j] * delta[j] for j in range(3))
        xr = sum(Gm[1][j] * delta[j] for j in range(3))
        xi = sum(Gm[2][j] * delta[j] for j in range(3))
        zr, zi = (
            mp.mpf(zs[0].numerator) / zs[0].denominator
            if hasattr(zs[0], "numerator")
            else mp.mpf(zs[0]),
            mp.mpf(zs[1].numerator) / zs[1].denominator
            if hasattr(zs[1], "numerator")
            else mp.mpf(zs[1]),
        )
        e1 = max(abs(xr - zr), abs(xi - zi))
        worst_q = max(worst_q, abs(q))
        worst_xi = max(worst_xi, e1)
        minpos = min(minpos, mn)
    Rf = mp.mpf(R_.numerator) / R_.denominator
    chk(
        f"(A.1) |q| <= .25 R at all {len(samples)} sampled z* (max {mp.nstr(worst_q / Rf, 5)} R)",
        worst_q <= 0.25 * Rf,
    )
    chk(
        f"(A.1) ||(Re xi, Im xi) - z*||_inf <= .13 R at all sampled z* (max {mp.nstr(worst_xi / Rf, 5)} R)",
        worst_xi <= 0.13 * Rf,
    )
    chk(
        f"1 + p > 0 on [1, 56] at all sampled z* (min over integer s: {mp.nstr(minpos, 5)})",
        minpos > 0,
    )
    if selftest:
        # wrong inputs: perturb the central parameter / use a wrong matrix / wrong eigen-coefficient
        k, d = kd((0, 0))
        kk, dd = (
            mp.mpf(k.numerator) / k.denominator,
            mp.mpf(d.numerator) / d.denominator,
        )
        (pj, p1j, p2j), mn = shoot(kk + mp.mpf(10) ** -9, dd)
        pa = pa_jet(dd, mp.mpf(S0))
        delta = [pj - pa[0], p1j - pa[1], p2j - pa[2]]
        q = sum(Gm[0][j] * delta[j] for j in range(3))
        chk(
            f"WRONG INPUT: k shifted by 1e-9 (about 2000x the allowed spread) violates |q| <= .25 R (|q|/R = {float(abs(q) / Rf):.2f})",
            abs(q) > 0.25 * Rf,
        )
        chk("WRONG INPUT: claimed bound M0 < 1.9 is rejected", not (Ms[0] < 1.9))
        chk(
            "WRONG INPUT: claimed C = -1215/5 would change the (3.5) coefficient check",
            not (
                (coef_s3, coef_s2, coef_s1, c0)
                == (1, 30, sp.Rational(751, 4), sp.Rational(-1215, 5))
            ),
        )
        chk(
            "WRONG INPUT: wrong projection entry (Lambda11 + .01) is rejected by the .001 entry bound",
            not (abs(Gam[0][0] - (Lam[0][0] + 0.01)) < 0.001),
        )
    log(f"elapsed {time.time() - t0:.1f}s")
    return finish()


def _checker_sha256():
    """sha256 of this script's source (CRLF normalised to LF); several files are concatenated in the order listed."""
    import hashlib
    import os

    here = os.path.dirname(os.path.abspath(__file__))
    names = [os.path.basename(__file__)]
    h = hashlib.sha256()
    for n in names:
        with open(os.path.join(here, n), "rb") as f:
            h.update(f.read().replace(b"\r\n", b"\n"))
    return h.hexdigest()


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    print("checker-sha256: " + _checker_sha256(), flush=True)
    sys.exit(run())
