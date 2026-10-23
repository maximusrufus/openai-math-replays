"""Exact scalar certificates of 'The periodic spin-one Haldane gap' (OpenAI, 2026-09-24): Sections 4-6 arithmetic (module of replay_haldane.py).
Exact rational arithmetic (fractions, sympy).  The tables of Lemma 5.1 are INPUTS here; the physical traces behind them are replayed in haldane_traces.py.
In the R(t) comparisons the printed '0.994' and '0.9984' are read as 0.99^4 and 0.998^4."""

import sys
from fractions import Fraction as Fr
from math import ceil, floor, log

import sympy as sp

RES = []


def chk(name, ok):
    raise RuntimeError("chk is provided by replay_haldane.py")


E6 = Fr(1, 10**6)
T21 = {
    6: (8945303, 346734, 939987),
    8: (3741281, 567167, 984108),
    10: (2327547, 721090, 1002613),
}
T49 = {
    4: (124885379, 65757),
    5: (2925, 2960742),
    6: (12870614, 286599),
    7: (54069, 1920414),
    8: (4639283, 510140),
    9: (195455, 1514846),
    10: (2654823, 675419),
    11: (377629, 1314544),
    12: (1900547, 787405),
}


def Ihat(n):
    z, zp, zc = [Fr(x) * E6 for x in T21[n]]
    return (z + 3 * zp + 8 * zc) / 12


def Ohat(n):
    z, zp, zc = [Fr(x) * E6 for x in T21[n]]
    return (z + 3 * zp - 4 * zc) / 12


def Nhat21(n):
    z, zp, zc = [Fr(x) * E6 for x in T21[n]]
    return (z - zp) / 4


# the two parameter rows (n0, b0, u0, C) of (4.18)
ROWS = [dict(n0=60, b0=Fr(105, 4), u0=Fr(1, 8), C=Fr(79, 10)), dict(n0=2304, b0=Fr(784), u0=Fr(1, 100), C=Fr(5))]


def f(u):
    return u * u / (2 * (1 - u) ** 2)


def Gf(u):
    return Fr(1, 2) * ((1 - u) ** -1 + (1 - u) ** -2 + (1 - u) ** -3) ** 2


def calC(m, h, r):
    q = floor(m / h)
    return q * h**r + (m - q * h) ** r


def up7(x):
    return Fr(ceil(x * 10**7), 10**7)


def run(selftest):
    # ------------------------------------------------------------ Section 4
    u = sp.symbols("u")
    fu = u**2 / (2 * (1 - u) ** 2)
    Gu = sp.Rational(1, 2) * ((1 - u) ** -1 + (1 - u) ** -2 + (1 - u) ** -3) ** 2
    chk(
        "(4.9) f(1-(1-u)^3) = G(u) u^2 (symbolic)",
        sp.simplify(fu.subs(u, 1 - (1 - u) ** 3) - Gu * u**2) == 0,
    )
    G18 = Gf(Fr(1, 8))
    chk(
        f"G(1/8) = {G18} = 913952/117649 < 79/10",
        G18 == Fr(913952, 117649) < Fr(79, 10),
    )
    chk(
        "row 1 of (4.18): u0 = 1/8 < 1/6, r = C u0 = 79/80 < 1, C(1-3u0) = 79/16 > 3, 3/(C(1-3u0)) = 48/79 < 1",
        Fr(1, 8) < Fr(1, 6)
        and Fr(79, 10) * Fr(1, 8) == Fr(79, 80) < 1
        and Fr(79, 10) * (1 - Fr(3, 8)) == Fr(79, 16) > 3
        and 3 / (Fr(79, 10) * (1 - Fr(3, 8))) == Fr(48, 79),
    )
    chk(
        "row 2: 9/(2(1-3u)^2) at u=1/100 = 45000/9409 < 5; G(1/100) <= that; r = 1/20; 3/(C(1-3u0)) = 60/97 < 1",
        Fr(9, 2) / (1 - Fr(3, 100)) ** 2 == Fr(45000, 9409) < 5
        and Gf(Fr(1, 100)) <= Fr(45000, 9409)
        and Fr(5) * Fr(1, 100) == Fr(1, 20)
        and 3 / (5 * (1 - Fr(3, 100))) == Fr(60, 97),
    )
    chk(
        "G(u) <= 9/(2(1-3u)^2) on a grid of u in (0, 1/100] (using 1-(1-u)^3 <= 3u)",
        all(
            Gf(Fr(k, 10000)) <= Fr(9, 2) / (1 - 3 * Fr(k, 10000)) ** 2
            for k in range(1, 101)
        ),
    )
    printed_rows = [(60, Fr(105, 4), Fr(1, 8), Fr(79, 10), Fr(79, 80)), (2304, Fr(784), Fr(1, 100), Fr(5), Fr(1, 20))]
    chk(
        "(4.18): the parameter rows (n0, b0, u0, C, r = C u0) are (60, 105/4, 1/8, 79/10, 79/80) and (2304, 784, 1/100, 5, 1/20)",
        [(w["n0"], w["b0"], w["u0"], w["C"], w["C"] * w["u0"]) for w in ROWS] == printed_rows,
    )
    chk(
        "(4.10) for both rows: 0 < u0 < 1/6, G(u0) <= C, r = C u0 < 1, C(1 - 3 u0) > 3",
        all(0 < w["u0"] < Fr(1, 6) and Gf(w["u0"]) <= w["C"] and w["C"] * w["u0"] < 1 and w["C"] * (1 - 3 * w["u0"]) > 3 for w in ROWS),
    )
    gaps = [log(1 / (w["C"] * w["u0"])) / w["b0"] for w in ROWS]  # (4.14): gamma_L > log(1/r)/b0
    chk(
        f"(4.14) gap bounds log(1/r)/b0 = {gaps[0]:.8f} and {gaps[1]:.8f} equal the printed closed forms (4/105) log(80/79) and (log 20)/784",
        abs(gaps[0] - 4 / 105 * log(80 / 79)) < 1e-15 and abs(gaps[1] - log(20) / 784) < 1e-15,
    )
    # ------------------------------------------------------------ Section 5.2
    chk(
        "(5.1) sector formulas: I=(Z+3ZP+8ZC)/12, O=(Z+3ZP-4ZC)/12, N=(Z-ZP)/4 reproduce Z = I + 2 O + 3 N for J = I + 2 O: I+2O = (Z+3ZP)/4",
        all(
            Ihat(n) + 2 * Ohat(n) == (Fr(T21[n][0]) + 3 * Fr(T21[n][1])) * E6 / 4
            for n in T21
        ),
    )
    mI, mO, mN = Fr(1042655, 10**6), Fr(40041, 10**6), Fr(401625, 10**6)
    eps = Fr(1, 10**5)
    marg = (mI - (Ihat(10) + eps), mO - (Ohat(10) + eps), mN - (Nhat21(10) + eps))
    chk(
        f"margins of mI, mO, mN over centre + 1e-5: {[str(x * 10**6) for x in marg]} (stated 19/12, 7/12, 3/4 in units 1e-6)",
        marg == (Fr(19, 12) * E6, Fr(7, 12) * E6, Fr(3, 4) * E6),
    )

    def poly_coeffs_mom(
        c0, c1
    ):  # x^6 (c0 + c1 x^2)^2 = c0^2 x^6 + 2 c0 c1 x^8 + c1^2 x^10
        return {6: c0 * c0, 8: 2 * c0 * c1, 10: c1 * c1}

    cI = poly_coeffs_mom(-2, 10)
    cN = poly_coeffs_mom(-25, 100)
    BI = sum(cI[j] * Ihat(j) for j in (6, 8, 10)) + eps * (abs(-2) + abs(10)) ** 2
    BN = sum(cN[j] * Nhat21(j) for j in (6, 8, 10)) + eps * (abs(-25) + abs(100)) ** 2
    U, V = Fr(100245, 100000), Fr(88643, 100000)
    FI = lambda x: x**6 * (-2 + 10 * x * x) ** 2
    FN = lambda x: x**6 * (-25 + 100 * x * x) ** 2
    chk(
        f"filter FI(U) - B_I = {float(FI(U) - BI):.5f} > 0.0301 (B_I = {float(BI):.6f}); FN(V) - B_N = {float(FN(V) - BN):.5f} > 0.3463 (B_N = {float(BN):.6f})",
        FI(U) - BI > Fr(301, 10000) and FN(V) - BN > Fr(3463, 10000),
    )
    chk(
        "filters strictly increasing beyond the endpoints: c1 > 0 and c0 + c1 x^2 > 0 at x = U and x = V",
        (-2 + 10 * U * U) > 0 and (-25 + 100 * V * V) > 0,
    )
    chk(
        "U^10 < mI < 2 U^10 and V^10 < mN < 2 V^10",
        U**10 < mI < 2 * U**10 and V**10 < mN < 2 * V**10,
    )

    # A_k, B_k, C_k
    def Ak(k):
        return U**k + (mI - U**10) ** (k // 10)

    def Bk(k):
        return float(mO) ** (k / 10)  # mO^{k/10}: k=30,120 -> integer powers 3 and 12

    def Bk_exact(k):
        return mO ** (k // 10)

    def Ck(k):
        return V**k + (mN - V**10) ** (k // 10)

    m_ = Ak(30) + 2 * Bk_exact(30) + 3 * Ck(30)
    p_ = Ak(30) + 2 * Bk_exact(30)
    q_ = Ak(120) + 2 * Bk_exact(120) + 3 * Ck(120)
    Dd, Ee, tau, tau0, Rq = (
        Fr(1450, 1000),
        Fr(1202, 1000),
        Fr(1225, 1000),
        Fr(1264, 1000),
        Fr(684, 10000),
    )
    chk(
        f"(5.5) m^5 < D^2 with margin {float(Dd**2 - m_**5):.6f} > 0.001068 (m = {float(m_):.6f})",
        m_**5 < Dd**2 and Dd**2 - m_**5 > Fr(1068, 10**6),
    )
    chk(
        f"(5.5) p^5 < E^2 with margin {float(Ee**2 - p_**5):.6f} > 0.000433 (p = {float(p_):.6f})",
        p_**5 < Ee**2 and Ee**2 - p_**5 > Fr(433, 10**6),
    )
    chk(
        f"(5.5) q >= 1 (q = {float(q_):.6f}) and (q-1)^5 < Rq^2 with margin {float(Rq**2 - (q_ - 1) ** 5):.7f} > 0.0000472",
        q_ >= 1 and (q_ - 1) ** 5 < Rq**2 and Rq**2 - (q_ - 1) ** 5 > Fr(472, 10**7),
    )
    chk(
        f"(5.5) D + 11E/12 ... : (D+11E)/12 = {float((Dd + 11 * Ee) / 12):.5f} <= tau = 1.225 and (D+3E)/4 = {float((Dd + 3 * Ee) / 4):.5f} <= tau0 = 1.264",
        (Dd + 11 * Ee) / 12 <= tau and (Dd + 3 * Ee) / 4 <= tau0,
    )

    # R(t)
    def Rfun(t):
        verts = [
            (Fr(0), Fr(0)),
            (tau - t, Fr(0)),
            (tau - t, tau0 - tau),
            (Fr(0), tau0 - t),
        ]
        return max(x * x + y * y / 2 + (Dd - t - x - y) ** 2 / 3 for x, y in verts)

    Rv = {t: Rfun(Fr(t)) for t in ("0.88", "0.99", "0.998")}
    pass
    chk(
        f"R(.88), R(.99), R(.998) = {[round(float(v), 6) for v in Rv.values()]} vs printed 0.1359, 0.0721, 0.068404 ",
        Rv["0.88"] == Fr(1359, 10000) and Rv["0.99"] == Fr(721, 10000) and Rv["0.998"] == Fr(68404, 10**6),  # exact
    )
    c1 = 1 - Fr(99, 100) ** 4 - Fr(1359, 10000) ** 2
    c2 = 1 - Fr(998, 1000) ** 4 - Fr(721, 10000) ** 2
    chk(
        f"1 - 0.99^4 - R(.88)^2 = {float(c1):.8f} (printed 0.02093518, using R=.1359); 1 - 0.998^4 - R(.99)^2 = {float(c2):.12f} (printed 0.002777621984, using R=.0721)",
        c1 == Fr(2093518, 10**8) and c2 == Fr(2777621984, 10**12),  # equal to the printed values at their full printed precision
    )
    c2b = 1 - Fr(998, 1000) ** 4 - Rv["0.99"] ** 2
    S60 = (1 + Rv["0.998"] / Fr(998, 1000) ** 2) ** -2
    T120 = (1 + Rq) ** -2
    chk(
        f"S(60,b) >= (1 + R(.998)/.998^2)^-2 = {float(S60):.6f} > .8756 > 7/8; T(120,b) >= (1+Rq)^-2 = {float(T120):.6f} > .8760 > 7/8",
        S60 > Fr(8756, 10000) > Fr(7, 8) and T120 > Fr(8760, 10000) > Fr(7, 8),
    )
    # ------------------------------------------------------------ Section 5.3
    Jn = lambda n: (Fr(T49[n][0]) + 3 * Fr(T49[n][1])) * E6 / 4
    Nn = lambda n: (Fr(T49[n][0]) - Fr(T49[n][1])) * E6 / 4
    x = sp.symbols("x")
    FJ = sp.Poly(
        sp.expand(x**4 * (7 + 17 * x - 280 * x**2 - 185 * x**3 + 1000 * x**4) ** 2), x
    )
    FN = sp.Poly(sp.expand(x**4 * (12 - 394 * x**2 + 1000 * x**4) ** 2), x)
    eps2 = Fr(30, 10**6)

    def Bsum(P, mom):
        tot = Fr(0)
        eabs = Fr(0)
        for (deg,), co in P.terms():
            co = Fr(int(co))
            tot += co * mom(deg)
            eabs += abs(co)
        return tot + eps2 * eabs

    BJ = Bsum(FJ, Jn)
    BN = Bsum(FN, Nn)
    chk(
        f"B_J = {float(BJ):.8f} (stated 318604.54848175), B_N = {float(BN):.6f} (stated 48169.880699)",
        abs(float(BJ) - 318604.54848175) < 1e-6
        and abs(float(BN) - 48169.880699) < 1e-5,
    )
    vJ, vN = Fr(100139, 100000), Fr(872, 1000)
    FJf = lambda xx: sum(Fr(int(co)) * xx**d for (d,), co in FJ.terms())
    FNf = lambda xx: sum(Fr(int(co)) * xx**d for (d,), co in FN.terms())
    chk(
        f"FJ(vJ) - BJ = {float(FJf(vJ) - BJ):.3f} > 180; FJ(-vJ) - BJ = {float(FJf(-vJ) - BJ):.3f} > 496888; FN(+-vN) - BN = {float(FNf(vN) - BN):.3f}, {float(FNf(-vN) - BN):.3f} > 654",
        FJf(vJ) - BJ > 180
        and FJf(-vJ) - BJ > 496888
        and FNf(vN) - BN > 654
        and FNf(-vN) - BN > 654,
    )
    y = sp.symbols("y")
    Qx = 7 + 17 * x - 280 * x**2 - 185 * x**3 + 1000 * x**4
    chk(
        "Q(1+y) = 559 + 2902 y + 5165 y^2 + 3815 y^3 + 1000 y^4 and Q(-1-y) = 895 + 3978 y + 6275 y^2 + 4185 y^3 + 1000 y^4",
        sp.expand(Qx.subs(x, 1 + y))
        == sp.expand(559 + 2902 * y + 5165 * y**2 + 3815 * y**3 + 1000 * y**4)
        and sp.expand(Qx.subs(x, -1 - y))
        == sp.expand(895 + 3978 * y + 6275 * y**2 + 4185 * y**3 + 1000 * y**4),
    )
    chk(
        "12 - 394 r^2 + 1000 r^4 positive at vN = .872 and derivative 4r(1000 r^2 - 197) > 0 for r >= vN",
        (12 - 394 * vN**2 + 1000 * vN**4) > 0 and 1000 * vN**2 - 197 > 0,
    )
    mJ, mN2 = Fr(10657205, 10**7), Fr(2783155, 10**7)
    chk(
        f"twelfth-moment upper bounds: J12 + 3e-5 = {float(Jn(12) + eps2)} = mJ = 1.0657205; N12 + 3e-5 = {float(Nn(12) + eps2)} = mN = 0.2783155",
        Jn(12) + eps2 == mJ and Nn(12) + eps2 == mN2,
    )
    ell = Fr(999, 1000)
    chk(f"mJ > ell^12 = {float(ell**12):.6f}", mJ > ell**12)
    Bk5 = lambda k: 3 * calC(mN2, vN**12, Fr(k, 12))
    RA = lambda k: Bk5(k) + (mJ - ell**12) ** (k // 12)
    # C(m,h,r) with r = k/12 integer for k = 36, 72
    lhs = calC(mJ, ell**12, 6) + Bk5(72)
    chk(
        f"C(mJ, l^12, 6) + B(72) = {float(lhs):.6f} < 1 with margin {float(1 - lhs):.5f} > 0.0693",
        lhs < 1 and 1 - lhs > Fr(693, 10000),
    )
    p0 = 1 - (1 + RA(36) / ell**36) ** -2
    q0 = 1 - (vJ**72 + RA(72)) ** -2
    chk(
        f"p0 = {float(p0):.7f} < 0.047916 and q0 = {float(q0):.7f} < 0.181521",
        p0 < Fr(47916, 10**6) and q0 < Fr(181521, 10**6),
    )
    # rounded updates
    pq = []
    pj, qj = p0, q0
    for j in range(1, 7):
        pj1 = up7(f(1 - (1 - pj) ** 2 * (1 - qj)))
        qj1 = up7(f(1 - (1 - qj) ** 2 * (1 - pj1)))
        pq.append((pj1, qj1))
        pj, qj = pj1, qj1
    printed = [
        (Fr(151249, 2500000), Fr(433453, 2500000)),
        (Fr(343303, 5000000), Fr(1632381, 10000000)),
        (Fr(44601, 625000), Fr(361773, 2500000)),
        (Fr(632939, 10000000), Fr(1055161, 10000000)),
        (Fr(375793, 10000000), Fr(445941, 10000000)),
        (Fr(84513, 10000000), Fr(27493, 5000000)),
    ]
    chk(
        f"six rounded updates reproduce the printed table exactly: {[(str(a), str(b)) for a, b in pq]}",
        pq == printed,
    )
    chk(
        f"1 - S(2304, 784) <= 84513/10^7 = {float(pq[-1][0])} < 1/100 and 1 - T(4608, 784) <= 27493/(5e6) = {float(pq[-1][1])} < 1/100",
        pq[-1][0] < Fr(1, 100) and pq[-1][1] < Fr(1, 100),
    )
    chk(
        "lengths and temperatures: 36 * 2^6 = 2304 and (49/4) * 2^6 = 784",
        36 * 64 == 2304 and Fr(49, 4) * 64 == 784,
    )
    # Corollary 5.6 numbers: ell = .999, U = vJ, V = vN
    chk(
        "Cor 5.6 interface: p1 = 151249/2500000 is the printed bound on 1 - S(72, 49/2)",
        pq[0][0] == Fr(151249, 2500000),
    )
    if selftest:
        # wrong inputs
        bad_tab = dict(T49)
        bad_tab[12] = (1900547 + 5000, 787405)
        Jb = lambda n: (Fr(bad_tab[n][0]) + 3 * Fr(bad_tab[n][1])) * E6 / 4
        chk(
            "WRONG INPUT: Z_12 table entry +0.005 changes m_J away from 1.0657205",
            Jb(12) + eps2 != mJ,
        )
        chk(
            "WRONG INPUT: claim margin of FI(U)-B_I > 0.05 is rejected",
            not (FI(U) - BI > Fr(5, 100)),
        )
        chk("WRONG INPUT: claim p0 < 0.04 is rejected", not (p0 < Fr(4, 100)))
        chk(
            "WRONG INPUT: update rounding ceil -> floor changes the printed table",
            [(Fr(floor(f(1 - (1 - p0) ** 2 * (1 - q0)) * 10**7), 10**7)) for _ in [0]][
                0
            ]
            != printed[0][0],
        )
        chk("WRONG INPUT: claim G(1/8) < 7 is rejected", not (G18 < 7))
    return 0

