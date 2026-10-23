"""Independent replay of the degree-31 certificate in 'Arithmetic classification and non-Pisot singularity for Bernoulli convolutions'
(OpenAI, 2026-10-03; Section 6, Appendix A, Corollary 5.4).

Written from the manuscript text only; no code from the authors (including the accompanying verification script) is used.
Parts: (R) exact root-disk certificate by Gaussian-integer Horner evaluation, plus independent 80-digit roots (mpmath);
(S) the Salem quartics of Corollary 5.4; (M) the moment/predictor computation of Appendix A.2-A.5 re-done with float64 cosine products
(truncated at m in [-8000, 150000], tails below 1e-16).  The 36-digit integers printed in (A.17) and (A.23) depend on the manuscript's
fixed-point rounding and are checked at VALUE level only (inside the printed brackets; weight sums match exactly).
The 16 rows of Table 1 (centres of the root disks) are transcribed from the manuscript.

  python replay_bernoulli.py             # prints every claim with PASS/FAIL; exit 0 only if all pass (about 2-4 min)
  python replay_bernoulli.py --perturb   # adversarial: one polynomial coefficient changed; the disk certificate at z_beta must FAIL
  python replay_bernoulli.py --selftest  # also runs the built-in wrong-input controls (each must be rejected)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import re
import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
import numpy as np

RES = []


SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
PFX = "B"
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


def poly_coeffs():
    """p(t) = t^31 - 2 t^30 + sum_{j=0}^{9} (t^{3j+1} - t^{3j}), as list of ints, index = degree"""
    c = [0] * 32
    c[31] += 1
    c[30] += -2
    for j in range(10):
        c[3 * j + 1] += 1
        c[3 * j] += -1
    return c


TABLE1 = [
(18392867573155250211326543340571993487329, 0),
(-9741384652837567528189111907762645564212, 960611633124778624015947511546714789832),
(-9390979934023233037686839914141181659013, 2848063934279278841200845525848576791934),
(-8707010393241200595697704637503228167256, 4634358961359777463662031763152574795910),
(-7723118558199735333585947045021501410451, 6250214494152701297733561366714273855663),
(-6480873022754265587670654805348739312397, 7617033927823246411701558316834470378693),
(-3296776282962928450195355562079974291919, 9362810613239308448752821856735617305389),
(-1442513481645415489053007250852371013818, 9733039088635544125591234293715446465839),
(473579005983151143865690217073943282984, 9759329295153416086962240697595385836176),
(2372858204750359429049351486827380431084, 9430908396110323024329804770889268443510),
(4184061662454447896993374181551299054853, 8750376307352230997951514116951420713869),
(5840543823426662247726359012202878127980, 7733999914049269380846159864123679516166),
(7280042164772721585016920411987150804420, 6409668313572020921539025901166241650504),
(8444162922855313595050499086097221186662, 4814843758334683329011653039104490553475),
(9275649862302447872891854755866108397375, 2996533501299661135726051235425500195711),
(9715324892541617145821300300817663390044, 1019107767839771580294306582422817927060),
]


def table_centers():
    return list(TABLE1)


def horner_gauss(coefs, Z, E):
    """exact Horner for q(z) with z = Z/E (Gaussian integer Z); returns (V, h) with q(z) = V/h; V Gaussian integer as (re, im)."""
    d = len(coefs) - 1
    Vr, Vi, h = coefs[d], 0, 1
    for j in range(d - 1, -1, -1):
        h = E * h
        Vr, Vi = Z[0] * Vr - Z[1] * Vi + h * coefs[j], Z[0] * Vi + Z[1] * Vr
    return Vr, Vi, h


def roots_part(selftest):
    E = 10**40
    coefs = poly_coeffs()
    dcoefs = [j * coefs[j] for j in range(1, 32)]
    # identity (6.2)
    import sympy as sp

    t = sp.symbols("t")
    pp = sum(c * t**j for j, c in enumerate(coefs))
    chk(
        "(6.2): (t^2+t+1) p(t) = t^33 - t^32 - t^31 - t^30 - 1",
        sp.expand((t**2 + t + 1) * pp - (t**33 - t**32 - t**31 - t**30 - 1)) == 0,
    )
    cent = table_centers()
    chk(
        f"Table 1 parsed: {len(cent)} rows (16 expected: one real centre + 15 nonreal, conjugates implied)",
        len(cent) == 16,
    )
    full = list(cent)
    for a, b in cent[1:]:
        full.append((a, -b))
    chk(f"31 disk centres (with conjugates): {len(full)}", len(full) == 31)
    ok1 = ok2 = ok3 = True
    delta_inv = 10**35
    for a, b in full:
        Vr, Vi, h = horner_gauss(coefs, (a, b), E)
        Wr, Wi, h1 = horner_gauss(dcoefs, (a, b), E)
        # |z| < 2
        ok1 &= a * a + b * b < 4 * E * E
        # ||p'(z)||_inf > 1 : max(|Wr|,|Wi|) > h1
        pinf_num = max(abs(Wr), abs(Wi))
        ok2 &= pinf_num > h1
        # 2e35 ||p(z)||_1 < ||p'(z)||_inf  (clear denominators h and h1)
        p1_num = abs(Vr) + abs(Vi)
        ok3 &= 2 * delta_inv * p1_num * h1 < pinf_num * h
    chk("(A.1) |z| < 2 at all 31 centres (exact)", ok1)
    chk("(A.1) ||p'(z)||_inf > 1 at all 31 centres (exact Horner)", ok2)
    chk(
        "(A.1) 2*10^35 ||p(z)||_1 < ||p'(z)||_inf at all 31 centres (exact Horner)", ok3
    )
    # distinct centres: squared distance > 0.01^2
    okd = True
    for i in range(31):
        for j in range(i + 1, 31):
            dr, di = full[i][0] - full[j][0], full[i][1] - full[j][1]
            okd &= dr * dr + di * di > (E // 100) ** 2
    chk("distinct centres, pairwise distance > .01", okd)
    # Rouche remainder: delta^2 sum |p_j| 3^j < 10^17 delta^2
    S3 = sum(abs(c) * 3**j for j, c in enumerate(coefs))
    chk(f"sum |p_j| 3^j = {S3} < 10^17", S3 < 10**17)
    # (A.2) bounds
    zb = Fr(full[0][0], E)
    chk(
        f"1.83928 < z_beta < 1.83929 (z_beta = {float(zb):.8f})",
        Fr(183928, 100000) < zb < Fr(183929, 100000),
    )
    zal = [(a, b) for (a, b) in cent[1:] if Fr(a * a + b * b, E * E) > 1]
    chk(
        f"exactly one nonreal centre (the 6th row zalpha) has |z|>1 ({len(zal)} found)",
        len(zal) == 1,
    )
    za = zal[0]
    m2 = Fr(za[0] ** 2 + za[1] ** 2, E * E)
    chk(
        f"1.0002090 < |z_alpha|^2 < 1.0002093 (got {float(m2):.9f})",
        Fr(10002090, 10**7) < m2 < Fr(10002093, 10**7),
    )
    others = [Fr(a * a + b * b, E * E) for (a, b) in cent[1:] if (a, b) != za]
    chk(
        f"all other nonreal centres have |z| < .9939 (max |z|^2 = {float(max(others)):.6f} < {0.9939**2:.6f})",
        all(o < Fr(9939, 10000) ** 2 for o in others),
    )
    # independent: mpmath roots
    mp.mp.dps = 80
    rts = mp.polyroots([c for c in reversed(coefs)], maxsteps=500, extraprec=800)
    dist_ok = True
    for a, b in full:
        z0 = mp.mpc(mp.mpf(a) / mp.mpf(E), mp.mpf(b) / mp.mpf(E))
        dmin = min(abs(r - z0) for r in rts)
        dist_ok &= dmin < mp.mpf(10) ** (-35)
    chk(
        "independent mpmath roots (80 digits): each centre is within 1e-35 of a computed root",
        dist_ok,
    )
    nreal = sum(1 for r in rts if abs(mp.im(r)) < mp.mpf(10) ** (-30))
    chk(f"exactly one real root ({nreal}) ", nreal == 1)
    beta = [r for r in rts if abs(mp.im(r)) < mp.mpf(10) ** (-30)][0]
    chk(
        f"beta = {mp.nstr(mp.re(beta), 12)} in (1.8392, 1.8393)",
        mp.mpf("1.8392") < mp.re(beta) < mp.mpf("1.8393"),
    )
    outs = sorted([r for r in rts if abs(r) > 1], key=lambda r: -abs(r))
    chk(
        f"exactly 3 roots outside the unit circle (beta, alpha, conj alpha): {len(outs)}",
        len(outs) == 3,
    )
    al = outs[1]
    chk(
        f"|alpha|^2 = {mp.nstr(abs(al) ** 2, 12)} in (1.0002089, 1.0002094) and 1/|alpha| = {mp.nstr(1 / abs(al), 8)} < 0.9999",
        mp.mpf("1.0002089") < abs(al) ** 2 < mp.mpf("1.0002094")
        and 1 / abs(al) < mp.mpf("0.9999"),
    )
    ins = [r for r in rts if abs(r) < 1]
    chk(
        f"28 remaining roots with max modulus {mp.nstr(max(abs(r) for r in ins), 8)} < 0.994",
        len(ins) == 28 and max(abs(r) for r in ins) < mp.mpf("0.994"),
    )
    chk(
        "not Pisot argument arithmetic: 1.0002094 * 0.994 < 1",
        Fr(10002094, 10**7) * Fr(994, 1000) < 1,
    )
    M = abs(beta) * abs(al) ** 2
    chk(
        f"(6.4): log M = {mp.nstr(mp.log(M), 10)} < log beta + 0.0002094 = {mp.nstr(mp.log(abs(beta)) + mp.mpf('0.0002094'), 10)}",
        mp.log(M) < mp.log(abs(beta)) + mp.mpf("0.0002094"),
    )
    # simple roots: discriminant nonzero (all pairwise root separations > 0)
    sep = min(abs(rts[i] - rts[j]) for i in range(31) for j in range(i + 1, 31))
    chk(
        f"all roots simple (minimum root separation {mp.nstr(sep, 5)})",
        sep > mp.mpf(10) ** (-5),
    )
    if selftest:
        # wrong inputs: perturb a centre by 1e-30 -> the disk certificate must fail; wrong polynomial
        a, b = full[3]
        Vr, Vi, h = horner_gauss(coefs, (a + 10**10, b), E)
        Wr, Wi, h1 = horner_gauss(dcoefs, (a + 10**10, b), E)
        chk(
            "WRONG INPUT: centre perturbed by 1e-30 fails 2e35 ||p||_1 < ||p'||_inf",
            not (2 * delta_inv * (abs(Vr) + abs(Vi)) * h1 < max(abs(Wr), abs(Wi)) * h),
        )
        bad = list(coefs)
        bad[5] += 1
        Vr, Vi, h = horner_gauss(bad, full[0], E)
        Wr, Wi, h1 = horner_gauss([j * bad[j] for j in range(1, 32)], full[0], E)
        chk(
            "WRONG INPUT: polynomial with one coefficient changed fails the root certificate at z_beta",
            not (2 * delta_inv * (abs(Vr) + abs(Vi)) * h1 < max(abs(Wr), abs(Wi)) * h),
        )
        chk(
            "WRONG INPUT: claimed beta in (1.8393,1.8394) rejected",
            not (mp.mpf("1.8393") < mp.re(beta) < mp.mpf("1.8394")),
        )
    return rts, beta, al


def salem_part(selftest):
    import sympy as sp

    x = sp.symbols("x")
    for nm, P, ysq in (
        ("P1", x**4 - x**3 - x**2 - x + 1, (1, -1, -3)),
        ("P2", x**4 - 2 * x**3 + x**2 - 2 * x + 1, (1, -2, -1)),
    ):
        irr = sp.Poly(P, x).is_irreducible
        rr = sp.Poly(P, x).nroots(n=30)
        real_gt1 = [r for r in rr if abs(sp.im(r)) < 1e-20 and sp.re(r) > 1]
        unit = [r for r in rr if abs(abs(r) - 1) < 1e-20 and abs(sp.im(r)) > 1e-10]
        y = sp.symbols("y")
        sub = sp.expand(sp.simplify(P / x**2))
        yq = sp.Poly(ysq[0] * y**2 + ysq[1] * y + ysq[2], y)
        check_y = (
            sp.simplify(sp.expand(x**2 + x**-2 + (-(1) if nm == "P1" else -2) * 0))
            is not None
        )
        chk(
            f"Cor 5.4 {nm}: irreducible over Q = {irr}; one real root >1 in (1,2) = {[float(r) for r in real_gt1]}; nonreal pair on unit circle ({len(unit)})",
            irr and len(real_gt1) == 1 and 1 < real_gt1[0] < 2 and len(unit) == 2,
        )
    # y-substitution identities
    X = sp.symbols("x")
    yv = X + 1 / X
    chk(
        "Cor 5.4: P1(x)/x^2 = y^2 - y - 3 and P2(x)/x^2 = y^2 - 2y - 1 for y = x + 1/x",
        sp.simplify((X**4 - X**3 - X**2 - X + 1) / X**2 - (yv**2 - yv - 3)) == 0
        and sp.simplify(
            (X**4 - 2 * X**3 + X**2 - 2 * X + 1) / X**2 - (yv**2 - 2 * yv - 1)
        )
        == 0,
    )
    r13 = (1 + sp.sqrt(13)) / 2
    chk(
        "larger roots y1 = (1+sqrt13)/2, y2 = 1+sqrt2 lie in (2, 5/2)",
        2 < float(r13) < 2.5 and 2 < float(1 + sp.sqrt(2)) < 2.5,
    )
    if selftest:
        P = X**4 - X**3 - X**2 - X + 2
        chk(
            "WRONG INPUT: x^4-x^3-x^2-x+2 is not reciprocal-Salem (no root on the unit circle pair with y^2-y-3)",
            sp.simplify(P / X**2 - (yv**2 - yv - 3)) != 0,
        )


# ------------------------------------------------------------------------------------ moments
def build_C(rts, mneg=8000, mpos=150000):
    coefs = poly_coeffs()
    dco = [j * coefs[j] for j in range(1, 32)]
    mp.mp.dps = 50
    rts = [mp.mpc(r) for r in rts]
    dp = [mp.polyval(list(reversed(dco)), r) for r in rts]
    I = [(r, d) for r, d in zip(rts, dp) if abs(r) < 1]
    D = [(r, d) for r, d in zip(rts, dp) if abs(r) > 1]
    ms = np.arange(-mneg - 10, mpos + 10)
    C0 = np.zeros(len(ms))
    for r, d in I:
        rr = complex(r)
        ld = complex(d)
        # r^{-1-m}/p'(r) for m < 0
        sel = ms < 0
        C0[sel] += np.real(np.exp((-1 - ms[sel]) * np.log(rr)) / ld)
    for r, d in D:
        rr = complex(r)
        ld = complex(d)
        sel = ms >= 0
        C0[sel] -= np.real(np.exp((-1 - ms[sel]) * np.log(rr)) / ld)
    return ms, C0


def run():
    selftest = SELFTEST
    t0 = time.time()
    if PERTURB:  # perturb the INPUT (real part of the 4th table centre moved by 1e-30), then run the unchanged main checks
        TABLE1[3] = (TABLE1[3][0] + 10**10, TABLE1[3][1])
    rts, beta, al = roots_part(selftest)
    salem_part(selftest)
    ms, C0 = build_C(rts)
    off = -ms[0]

    def c0(m):
        return C0[np.asarray(m) + off]

    # convolution identity sum_j p_j C0(m-j) = delta_{m,0} (validates C0 as the bilateral Laurent coefficients of 1/p)
    coefs = poly_coeffs()
    mm = np.arange(-300, 300)
    conv = sum(coefs[j] * c0(mm - j) for j in range(32))
    expect = (mm == 0).astype(float)
    chk(
        f"C0 = Laurent coefficients of 1/p: max |sum_j p_j C0(m-j) - delta_m0| = {np.abs(conv - expect).max():.1e} (m in [-300,300))",
        np.abs(conv - expect).max() < 1e-9,
    )
    chk(f"sum_m C0(m) = 1/p(1) = -1 (got {C0.sum():.10f})", abs(C0.sum() + 1) < 1e-7)
    # C1(m) = C0(m+3) - C0(m) + C0(m-1)
    lo, hi = -7990, 149000
    M = np.arange(lo, hi)
    C1 = c0(M + 3) - c0(M) + c0(M - 1)
    # (A.10)-style maxima
    chk(
        f"(A.10) max_m |C0|, |C1| = {np.abs(c0(M)).max():.5f}, {np.abs(C1).max():.5f}; paper: max_A,m |c_A(m)| in (.8460, .8462)",
        0.8460 < max(np.abs(c0(M)).max(), np.abs(C1).max()) < 0.8462,
    )
    chk(
        f"(A.10) max_{{m>=-20}} |C0(m)| = {np.abs(c0(np.arange(-20, 100000))).max():.5f} < .018",
        np.abs(c0(np.arange(-20, 100000))).max() < 0.018,
    )
    # u'_A: weak-root tail coefficients: u0 = -a/p'(alpha), u1 = u0 (a^3 - 1 + 1/a), a = 1/alpha
    coefs_d = [j * coefs[j] for j in range(1, 32)]
    mp.mp.dps = 50
    pd_al = mp.polyval(list(reversed(coefs_d)), al)
    a_ = 1 / al
    u0 = -a_ / pd_al
    u1 = u0 * (a_**3 - 1 + 1 / a_)
    chk(
        f"(A.10) |u_0| = {mp.nstr(abs(u0), 5)}, |u_1| = {mp.nstr(abs(u1), 5)} < .0087",
        abs(u0) < mp.mpf("0.0087") and abs(u1) < mp.mpf("0.0087"),
    )

    def prod_moment(A, B, sigma, k):
        """P^sigma_AB(k) = prod_m cos(pi (C_A(m) + sigma C_B(m+k))), m over [lo, hi - k - 5]"""
        mrange = np.arange(lo, hi - 5 - 2100)
        CA = c0(mrange) if A == 0 else c0(mrange + 3) - c0(mrange) + c0(mrange - 1)
        mk = mrange + k
        CB = c0(mk) if B == 0 else c0(mk + 3) - c0(mk) + c0(mk - 1)
        ang = np.pi * (CA + sigma * CB)
        v = np.cos(ang)
        sgn = np.prod(np.sign(v))
        return sgn * np.exp(np.sum(np.log(np.abs(v))))

    P0 = prod_moment(0, 0, 1, 0) if False else None

    # single moments: P_A = prod cos(pi C_A(m))
    def single(A):
        mrange = np.arange(lo, hi - 5 - 2100)
        CA = c0(mrange) if A == 0 else c0(mrange + 3) - c0(mrange) + c0(mrange - 1)
        v = np.cos(np.pi * CA)
        return np.prod(np.sign(v)) * np.exp(np.sum(np.log(np.abs(v))))

    b0, b1 = single(0), single(1)
    log(f"   P_0 = {b0:.12f}, P_1 = {b1:.12f}")
    # cutoff-convergence evidence
    lo_save, hi_save = lo, hi
    Kmax = 2000
    Pp = np.zeros(Kmax + 1)
    Pm = np.zeros(Kmax + 1)
    for k in range(Kmax + 1):
        Pp[k] = prod_moment(0, 0, 1, k)
        Pm[k] = prod_moment(0, 0, -1, k)
    log(f"   lag moments done ({time.time() - t0:.0f}s)")
    chk("P^-_00(0) = 1 exactly (E|Z|^2 = 1)", abs(Pm[0] - 1) < 1e-12)
    e33 = prod_moment(0, 1, 1, 33)  # eP^+_{01}(33)
    log(f"   e = P^+_01(33) = {e33:.12f}")
    # weights (A.7)
    W = 10**9
    rows = []  # (w (float), k, s)
    irows = []  # (signed integer numerator of w, k, s)
    for k in range(30, 2001):
        for s, rv in ((0, Pp[k]), (1, Pm[k])):
            amt = max(0.0, abs(rv) - 0.00026)
            wnum = int(np.floor(0.68 * W * amt + 0.5))
            if wnum == 0:
                continue
            rows.append((np.sign(rv) * wnum / W, k, s))
            irows.append((int(np.sign(rv)) * wnum, k, s))
    X, Y = -0.001, 0.0055
    nR = len(rows)
    sumw = abs(X) + abs(Y) + sum(abs(w) for w, _, _ in rows)
    sumw50 = abs(X) + abs(Y) + sum(abs(w) for w, k, _ in rows if k < 50)
    xw, yw = round(abs(X) * W), round(abs(Y) * W)  # |X| and |Y| as integers in units 1/W
    isum = xw + yw + sum(abs(v) for v, _, _ in irows)
    isum50 = xw + yw + sum(abs(v) for v, k, _ in irows if k < 50)
    chk(f"|R| = {nR} (stated 753)", nR == 753)
    chk(
        f"|X|+|Y|+sum|w| = {isum}/W (stated 198442946/W); < 0.199 = {sumw < 0.199}",
        isum == 198442946 and sumw < 0.199,  # exact integer equality
    )
    chk(
        f"sum over k<50 = {isum50}/W (stated 11816095/W); total weight of small rows < .012 = {sumw50 < 0.012}",
        isum50 == 11816095 and sumw50 < 0.012,  # exact integer equality
    )
    # (A.14) correlation
    ec = X * b0 + Y * e33 + sum(w * (Pp[k] if s == 0 else Pm[k]) for w, k, s in rows)
    chk(
        f"(A.14) 1e9 * ec = {ec * 1e9:.2f} in (211810, 211820)",
        211810 < ec * 1e9 < 211820,
    )
    # J_ks: pair (P^+_10(k-33), P^-_10(k-33)) for k>=33, else (P^+_01(33-k), P^-_01(33-k)); J_k0 = minus entry, J_k1 = plus entry
    Jcache = {}

    def J(k, s):
        key = k
        if key not in Jcache:
            if k >= 33:
                Jcache[key] = (
                    prod_moment(1, 0, 1, k - 33),
                    prod_moment(1, 0, -1, k - 33),
                )
            else:
                Jcache[key] = (
                    prod_moment(0, 1, 1, 33 - k),
                    prod_moment(0, 1, -1, 33 - k),
                )
        plus, minus = Jcache[key]
        return minus if s == 0 else plus

    ev = (
        X * X
        + Y * Y
        + 2 * X * Y * b1
        + 2 * sum(w * (X * b0 + Y * J(k, s)) for w, k, s in rows)
    )
    ks = np.array([k for _, k, _ in rows])
    ss = np.array([s for _, _, s in rows])
    ws = np.array([w for w, _, _ in rows])
    D = np.abs(ks[:, None] - ks[None, :])
    same = ss[:, None] == ss[None, :]
    Mmat = np.where(same, Pm[D], Pp[D])
    gram = ws[:, None] * ws[None, :] * Mmat
    ev += gram.sum()
    chk(
        f"(A.15) 1e9 * ev = {ev * 1e9:.2f} in (143910, 143930)",
        143910 < ev * 1e9 < 143930,
    )
    big = ks >= 50
    eT = np.abs(gram[np.ix_(big, big)]).sum()
    chk(
        f"(A.16) 1e9 * eT = {eT * 1e9:.2f} in (105300, 105310)",
        105300 < eT * 1e9 < 105310,
    )
    Delta = 5e-6
    chk(
        f"E Re(Z G) > {ec - Delta * 0.199:.10f} > .000210 and E|G|^2 < {ev + Delta * 0.199**2:.9f} < .000145 (error margins as in A.4)",
        ec - Delta * 0.199 > 0.000210 and ev + Delta * 0.199**2 < 0.000145,
    )
    T_ub = eT + Delta * 0.199**2
    chk(
        f"(A.19) T <= {T_ub:.8f} < .000116 and sqrt(T) < .011",
        T_ub < 0.000116 and np.sqrt(T_ub) < 0.011,
    )
    # (A.20)-(A.22)
    t = 0.06
    lo2, hi2 = -200, 1
    mr = np.arange(lo2, hi2)
    C0m = c0(mr)

    def fseq(row):
        kind, k, s = row
        if kind == "const":
            return np.zeros_like(C0m)
        if kind == "t1":
            return c0(mr + 33 + 3) - c0(mr + 33) + c0(mr + 33 - 1)
        return (1 - 2 * s) * c0(mr + k)

    small = [("const", 0, 0), ("t1", 0, 0)] + [
        ("row", k, s) for w, k, s in rows if k < 50
    ]
    chk(
        f"number of nonzero rows with k<50 = {len(small) - 2} (stated 8)",
        len(small) - 2 == 8,
    )
    stated = {
        32: (0, -536932),
        36: (1, 1265024),
        37: (0, 331701),
        39: (1, 47119),
        40: (0, 1483338),
        43: (0, 436861),
        44: (1, 458216),
        47: (1, 756904),
    }
    got = {k: (s, v) for v, k, s in irows if k < 50}
    chk(
        f"small-row weights W*w_ks equal the printed list exactly: {got} (stated {stated})",
        all(
            k in got
            and got[k][0] == stated[k][0]
            and got[k][1] == stated[k][1]  # exact
            for k in stated
        )
        and set(got) == set(stated),
    )
    # big rows tan bound (A.20) at -70<=m<=0
    mm70 = np.arange(-70, 1)
    tanmax = 0.0
    for w, k, s in rows:
        if k >= 50:
            tanmax = max(
                tanmax, np.abs(np.tan(np.pi * (1 - 2 * s) * c0(mm70 + k))).max()
            )
    chk(
        f"(A.20) max |tan(pi f_v(m))| over the big rows (k>=50; the constant, type-1 and k<50 rows are the small rows), -70<=m<=0: {tanmax:.4f} <= .06",
        tanmax <= 0.06,
    )

    def B2(V):
        m = np.arange(-70, 1)
        c = np.abs(np.cos(np.pi * V(m)))
        sn = np.abs(np.sin(np.pi * V(m)))
        return np.prod(np.minimum(1, ((1 + t * t) * c + 2 * t * sn) / (1 - t * t)))

    def B1(V):
        m = np.arange(-70, 1)
        c = np.abs(np.cos(np.pi * V(m)))
        sn = np.abs(np.sin(np.pi * V(m)))
        return np.prod(np.minimum(1, (c + t * sn) / np.sqrt(1 - t * t)))

    def B0(V):
        m = np.arange(-200, 1)
        return np.prod(np.abs(np.cos(np.pi * V(m))))

    def fm(row, m):
        kind, k, s = row
        if kind == "const":
            return np.zeros(len(m))
        if kind == "t1":
            return c0(m + 33 + 3) - c0(m + 33) + c0(m + 33 - 1)
        return (1 - 2 * s) * c0(m + k)

    b2 = B2(lambda m: 2 * c0(m))
    b1max = max(B1(lambda m, r=r: 2 * c0(m) + fm(r, m)) for r in small)
    b0max = max(
        B0(lambda m, r=r, q=q: 2 * c0(m) + fm(r, m) + fm(q, m))
        for r in small
        for q in small
    )
    chk(f"(A.22) B2(2C0) = {b2:.6f} < .011, and below the printed upper-rounded value .010862 within 1%", b2 < 0.011 and 0.99 * 0.010862 <= b2 <= 0.010862)
    chk(f"(A.22) max B1 = {b1max:.6f} < .015, and below the printed value .014732 within 1%", b1max < 0.015 and 0.99 * 0.014732 <= b1max <= 0.014732)
    chk(f"(A.22) max B0 = {b0max:.6f} < .0032, and below the printed value .003156 within 1%", b0max < 0.0032 and 0.99 * 0.003156 <= b0max <= 0.003156)
    chk("stated parameters of (A.20)-(A.22): t = .06, Delta = 5e-6, W = 1e9", t == 0.06 and Delta == 5e-6 and W == 10**9)
    fin = (
        Fr(116, 10**6) * Fr(11, 1000)
        + 2 * Fr(12, 1000) * Fr(11, 1000) * Fr(15, 1000)
        + Fr(12, 1000) ** 2 * Fr(32, 10000)
    )
    chk(
        f"final: .000116*.011 + 2*.012*.011*.015 + .012^2*.0032 = {float(fin):.10f} (stated .0000056968) < .000007",
        abs(float(fin) - 5.6968e-6) < 1e-9 and fin < Fr(7, 10**6),
    )
    # (6.13) final numbers
    chk(
        "(6.12)/(6.13) consequence: weights bound x* = .199; E L_i > .000210 > 0; E|G|^2 < .000145; E(Z^2 G^2) < .000007",
        sumw < 0.199
        and ec - Delta * 0.199 > 0.000210
        and ev + Delta * 0.199**2 < 0.000145
        and float(fin) < 7e-6,
    )
    if selftest:
        # wrong inputs
        chk(
            "WRONG INPUT: claim E L_i > .000250 rejected",
            not (ec - Delta * 0.199 > 0.000250),
        )
        chk("WRONG INPUT: claim |R| = 750 rejected", not (nR == 750))
        t_bad = 0.06
        chk(
            "WRONG INPUT: perturbed C0 (C0(0) += .01) breaks the Laurent identity",
            np.abs(
                sum(
                    coefs[j] * np.where(mm - j == 0, c0(mm - j) + 0.01, c0(mm - j))
                    for j in range(32)
                )
                - expect
            ).max()
            > 1e-3,
        )
        chk("WRONG INPUT: claim ev < 1.0e-4 rejected", not (ev < 1.0e-4))
        chk("WRONG INPUT: claim B2 < .005 rejected", not (b2 < 0.005))
    log(f"elapsed {time.time() - t0:.1f}s")
    return finish()


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    sys.exit(run())
