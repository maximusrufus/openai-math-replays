"""Independent replay of the numerical certificate in 'Subset Sum in Time O(2^{0.49n})' (OpenAI, 2026-10-04; Proposition 3.1, Appendix A)
and of the exact elementary inequalities of Section 2.2 of 'A Low-Space Algorithm for Worst-Case Subset Sum' (OpenAI, 2026-09-26).

Written from the manuscript text only (Section 3 definitions, Appendix A.1-A.3); no code from the authors is used.
Conventions read from the rendered pages 8-9: F = 2^r, G = 2^s; the integrands are M(x)^2 * 2^{8.5(Delta-u)} and M(x)^2 * 2^{8.3(a v - u)};
all logarithms are base 2.
Method: (1) the 64-point hats of Appendix A.1-A.2 in ball arithmetic compared with Table 1 and (A.3);
(2) the continuous integrals of Proposition 3.1 computed directly from (3.1)-(3.2) (not via (A.2)) by trapezoid sums on
N = 256, 512, 1024, 2048 nodes (analytic periodic integrands: geometric convergence).  This is a high-precision numerical replay,
NOT a rigorous enclosure of the quadrature remainder.

  python replay_subset_sum.py             # prints every claim with PASS/FAIL; exit 0 only if all pass
  python replay_subset_sum.py --perturb   # adversarial: sign of r6 flipped; the Table 1 comparison must FAIL (exit 1)
  python replay_subset_sum.py --selftest  # also runs the built-in wrong-input controls (each must be rejected)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import re
import sys
import time

import numpy as np
from flint import arb, ctx

ctx.prec = 160
R6 = (0.851, 0.466, 0.325, 0.253, 0.224, -0.314)
S6 = (1.281, -0.938, 0.731, -0.559, 0.391, 0.110)
KAPPA = 1.0315
ETA, TAU = 0.51772, 0.48996
DELTA = ETA - TAU
GAMMA = 2 * TAU - ETA
AA = DELTA / GAMMA
RES = []


SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
PFX = "SS"
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


# ---------------------------------------------------------------- float64 continuous replay
def continuous(N, rr=None, ss=S6, kappa=KAPPA, c1=8.5, c2=8.3, extra=None):
    x = 2 * np.pi * np.arange(N) / N
    rr = R6 if rr is None else rr
    r = sum(c * np.cos((j + 1) * x) for j, c in enumerate(rr))
    s = sum(c * np.cos((j + 1) * x) for j, c in enumerate(ss))
    F = 2.0**r
    G = 2.0**s
    Zf, Zg = F.mean(), G.mean()
    Zfk, Zgk = (F**kappa).mean(), (G**kappa).mean()
    f = F / Zf
    g = G / Zg
    idx = (np.arange(N)[:, None] + np.arange(N)[None, :]) % N  # x_i + w_j
    fx = f[idx]  # f(x_i + w_j)
    Kn = fx * g[None, :]
    Q = Kn.mean(axis=1)  # Q(x_i) = int f(x+w) g(w) dw
    K = Kn / Q[:, None]
    logK = np.log2(K)
    logfs = kappa * r[idx] - np.log2(
        Zfk
    )  # log f*(x+w) = kappa*log2 F - log2 Zf,kappa  (derived from f* = f^k/int f^k)
    # independent route: build f* numerically from its definition and take log
    fk = f**kappa
    fstar = fk / fk.mean()
    logfs_direct = np.log2(fstar)[idx]
    gk = g**kappa
    gstar = gk / gk.mean()
    loggs = np.log2(gstar)[None, :]
    D = (K * logK).mean(axis=1)
    A = (K * logfs_direct).mean(axis=1)
    E = (K * loggs).mean(axis=1)
    u = A + E - D
    v = 2 * D - A - E
    M2 = 2 + 2 * np.cos(x)
    I1 = (M2 * 2.0 ** (c1 * (DELTA - u))).mean()
    I2 = (M2 * 2.0 ** (c2 * (AA * v - u))).mean()
    out = dict(
        I1=I1,
        I2=I2,
        D=D,
        A=A,
        E=E,
        u=u,
        v=v,
        Q=Q,
        logfs_gap=np.abs(logfs - logfs_direct).max(),
        Zf=Zf,
        Zg=Zg,
    )
    return out


# ---------------------------------------------------------------- 64-point hats in ball arithmetic
def hats64(prec_digits=None, c1=8.5, c2=8.3, rr=None, ss=S6, kappa=KAPPA):
    rr = R6 if rr is None else rr
    n = 64
    pi = arb.pi()
    xs = [2 * pi * j / n for j in range(n)]
    two = arb(2)
    lg = lambda y: y.log() / arb(2).log()
    r = [
        sum(
            arb(str(c)) * (j1 + 1 > 0 and (xs[i] * (j1 + 1)).cos())
            for j1, c in enumerate(rr)
        )
        for i in range(n)
    ]
    s = [
        sum(arb(str(c)) * (xs[i] * (j1 + 1)).cos() for j1, c in enumerate(ss))
        for i in range(n)
    ]
    F = [two**r_ for r_ in r]
    G = [two**s_ for s_ in s]
    Zf = sum(F) / n
    Zg = sum(G) / n
    kap = arb(str(kappa))
    Zfk = sum(f_**kap for f_ in F) / n
    Zgk = sum(g_**kap for g_ in G) / n
    Ph, Yh = [], []
    for j in range(n):
        Ph.append(sum(F[(i + j) % n] * G[i] for i in range(n)) / n)
        Yh.append(
            sum(F[(i + j) % n] * G[i] * (r[(i + j) % n] + s[i]) for i in range(n)) / n
        )
    Hh = lg(Zfk) + lg(Zgk) - kap * (lg(Zf) + lg(Zg))
    Lh = [lg(Ph[j]) - lg(Zf) - lg(Zg) for j in range(n)]
    Bh = [Yh[j] / Ph[j] - lg(Zf) - lg(Zg) for j in range(n)]
    uh = [Lh[j] + (kap - 1) * Bh[j] - Hh for j in range(n)]
    vh = [(2 - kap) * Bh[j] - 2 * Lh[j] + Hh for j in range(n)]
    Dl = arb(str(DELTA))
    aa = arb(str(AA))
    M2 = [2 + 2 * xs[j].cos() for j in range(n)]
    X1 = [M2[j] * two ** (arb(str(c1)) * (Dl - uh[j])) for j in range(n)]
    X2 = [M2[j] * two ** (arb(str(c2)) * (aa * vh[j] - uh[j])) for j in range(n)]
    return dict(
        Zf=Zf, Zg=Zg, Zfk=Zfk, Zgk=Zgk, P=Ph, Y=Yh, X1=X1, X2=X2, H=Hh, r=r, s=s, xs=xs, L=Lh, B=Bh, u=uh, v=vh
    )


TABLE = """0 1.7255739 2.3814527 1.71479 1.72487
1 1.7267975 2.3812891 1.70072 1.70914
2 1.7298232 2.3798409 1.66451 1.66858
3 1.7329951 2.3749116 1.62068 1.61932
4 1.7343372 2.3645517 1.58491 1.57884
5 1.7323079 2.3482031 1.56813 1.55943
6 1.7263040 2.3267607 1.57405 1.56556
7 1.7167948 2.3018438 1.59915 1.59372
8 1.7051334 2.2751159 1.63351 1.63323
9 1.6931319 2.2480097 1.66273 1.66812
10 1.6824618 2.2214392 1.67219 1.68188
11 1.6740048 2.1950779 1.65371 1.66486
12 1.6673951 2.1665832 1.61070 1.62013
13 1.6610380 2.1317550 1.55729 1.56268
14 1.6527012 2.0862092 1.51133 1.51181
15 1.6404589 2.0278888 1.48634 1.48247
16 1.6235129 1.9586974 1.48675 1.48020
17 1.6024686 1.8838599 1.50690 1.50018
18 1.5789685 1.8092831 1.53309 1.52883
19 1.5549462 1.7387429 1.54784 1.54794
20 1.5318669 1.6725904 1.53630 1.54099
21 1.5101580 1.6080342 1.49268 1.50055
22 1.4888661 1.5397614 1.42405 1.43289
23 1.4655935 1.4601151 1.34881 1.35678
24 1.4368999 1.3596108 1.29024 1.29629
25 1.3993369 1.2292982 1.26827 1.27187
26 1.3509552 1.0651814 1.29131 1.29177
27 1.2926991 0.8725981 1.34397 1.33998
28 1.2289553 0.6674810 1.36465 1.35461
29 1.1668804 0.4731120 1.23032 1.21474
30 1.1147571 0.3140122 0.82842 0.81341
31 1.0800632 0.2102374 0.27604 0.26987
32 1.0678994 0.1742636 0 0"""
TAB = [tuple(map(float, ln.split())) for ln in TABLE.splitlines()]
ERR = (4e-6, 2e-5, 4e-4, 4e-4)


def mid(a):
    return float(a.mid())



def lowspace_exact(selftest):
    """Section 2.2 of the companion 'A Low-Space Algorithm for Worst-Case Subset Sum' (second manuscript): exact inequalities."""
    from fractions import Fraction as Fr
    from math import log2

    chk("3^12/2^19 > 1 + 1/74", Fr(3**12, 2**19) > 1 + Fr(1, 74))
    chk("(1+1/74)^4 > 1 + 4/74 = 39/37 > 256/243", (1 + Fr(1, 74)) ** 4 > 1 + Fr(4, 74) == Fr(39, 37) > Fr(256, 243))
    chk("3^53 > 2^84", 3**53 > 2**84)
    chk("84/53 > 1.5848", Fr(84, 53) > Fr(15848, 10000))
    hstar_ub = 2 - Fr(3, 4) * Fr(84, 53)
    chk(f"h* = 2 - (3/4) log 3 < 0.8114 (using log 3 > 84/53: bound {float(hstar_ub):.6f}); float check {2 - 0.75 * log2(3):.6f}", hstar_ub < Fr(8114, 10000) and 2 - 0.75 * log2(3) < 0.8114)
    chk("3^5 < 2^8 and h(1/3) = log 3 - 2/3 < 8/5 - 2/3 = 14/15", 3**5 < 2**8 and Fr(8, 5) - Fr(2, 3) == Fr(14, 15) and log2(3) - 2 / 3 < 14 / 15)
    chk("6^2 > 2^5, (6/5)^4 > 2 and h(1/6) > 5/8", 6**2 > 2**5 and Fr(6, 5) ** 4 > 2 and Fr(1, 6) * Fr(5, 2) + Fr(5, 6) * Fr(1, 4) == Fr(5, 8) and (1 / 6) * log2(6) + (5 / 6) * log2(6 / 5) > 5 / 8)
    chk("(39) arithmetic: .032 - 1.25 + 1.5 h* < 0 with h* < .8114", Fr(32, 1000) - Fr(125, 100) + Fr(3, 2) * Fr(8114, 10000) < 0)
    chk("(39) arithmetic: 1.032 - h* > .2206 and -.051 + .2206*.24 >= .0015", Fr(1032, 1000) - Fr(8114, 10000) == Fr(2206, 10000) and -Fr(51, 1000) + Fr(2206, 10000) * Fr(24, 100) >= Fr(15, 10000))
    chk("exponent sum .249 + .25 = .499 < .5", Fr(249, 1000) + Fr(1, 4) == Fr(499, 1000) < Fr(1, 2))
    if selftest:
        chk("WRONG INPUT: claim 3^53 > 2^85 is rejected", not (3**53 > 2**85))
        chk("WRONG INPUT: claim h(1/6) > 2/3 is rejected", not ((1 / 6) * log2(6) + (5 / 6) * log2(6 / 5) > 2 / 3))
        chk("WRONG INPUT: claim h* < .81 is rejected", not (2 - 0.75 * log2(3) < 0.81))


def run():
    selftest = SELFTEST
    t0 = time.time()
    if PERTURB:  # perturb the INPUT (sign of r6 flipped), then run the unchanged main checks
        global R6
        R6 = tuple(list(R6[:5]) + [0.314])
    # parameters
    chk(
        "Delta = eta - tau = 0.02776 and Gamma = 2 tau - eta = 0.46220 (exact decimal arithmetic)",
        abs(DELTA - 0.02776) < 1e-12 and abs(GAMMA - 0.46220) < 1e-12,
    )
    chk(
        "|r| <= sum |r_j| = 2.433 and |s| <= 4.010",
        abs(sum(abs(c) for c in R6) - 2.433) < 1e-12
        and abs(sum(abs(c) for c in S6) - 4.010) < 1e-12,
    )
    # 64-point hats
    H64 = hats64()
    cz = {"Zf": 1.1963625, "Zg": 1.2877407, "Zfk": 1.2110790, "Zgk": 1.3038111}
    okz = all(abs(mid(H64[k]) - v) <= 5e-7 for k, v in cz.items())
    chk(
        f"Zhat_f, Zhat_g, Zhat_f,kappa, Zhat_g,kappa within 5e-7 of the stated centres ({[round(mid(H64[k]), 7) for k in cz]})",
        okz,
    )
    # T64(2^{2 r_beta}) and T64(2^{2 s_beta})  

    def T64_exp(coefs, beta=0.19):
        n = 64
        tot = arb(0)
        for i in range(n):
            xi = 2 * arb.pi() * i / n
            val = sum(
                arb(str(c)) * ((arb(str(beta)) * (j + 1)).cosh()) * (xi * (j + 1)).cos()
                for j, c in enumerate(coefs)
            )
            tot += arb(2) ** (2 * val)
        return tot / n

    t_r, t_s = mid(T64_exp(R6)), mid(T64_exp(S6))
    chk(
        f"T64(2^(2 r_beta)) = {t_r:.7f} (stated 2.4901514 +-5e-6), T64(2^(2 s_beta)) = {t_s:.7f} (stated 2.2607200 +-5e-6)",
        abs(t_r - 2.4901514) <= 5e-6 and abs(t_s - 2.2607200) <= 5e-6,
    )
    worst = [0, 0, 0, 0]
    for j, (jj, p, y, x1, x2) in enumerate(TAB):
        got = (mid(H64["P"][j]), mid(H64["Y"][j]), mid(H64["X1"][j]), mid(H64["X2"][j]))
        for c in range(4):
            worst[c] = max(worst[c], abs(got[c] - (p, y, x1, x2)[c]))
    chk(
        f"Table 1 (nodes 0..32): max |ours - stated| = {['%.2e' % w for w in worst]} vs allowed {['%.0e' % e for e in ERR]}",
        all(worst[c] <= ERR[c] for c in range(4)),
    )
    sym = max(abs(mid(H64["X1"][j]) - mid(H64["X1"][64 - j])) for j in range(1, 32))
    chk(
        f"symmetry j -> 64-j of the node columns (max asymmetry {sym:.1e})", sym < 1e-12
    )
    a1 = sum(H64["X1"]) / 64
    a2 = sum(H64["X2"]) / 64
    chk(
        f"(A.3) mean of Xhat_j = {a1.str(12)} < 1.4295 and mean of Xhat'_j = {a2.str(12)} < 1.4295 (ball arithmetic, 64 nodes)",
        a1 < arb("1.4295") and a2 < arb("1.4295"),
    )
    bound = 2**0.5168
    chk(
        f"1.4295 < 2^0.5168 = {bound:.6f}; margin for the outer quadrature remainder = {bound - 1.4295:.6f}",
        1.4295 < bound,
    )
    # continuous integrals
    vals = {}
    for N in (256, 512, 1024, 2048):
        o = continuous(N)
        vals[N] = (o["I1"], o["I2"])
        log(f"   N={N}: I1={o['I1']:.15f} I2={o['I2']:.15f}")
    conv = max(abs(vals[2048][0] - vals[1024][0]), abs(vals[2048][1] - vals[1024][1]))
    chk("continuous integrals converge: |I(2048)-I(1024)| < 1e-12", conv < 1e-12)
    I1, I2 = vals[2048]
    chk(
        f"Prop 3.1 (3.6): I1 = {I1:.10f} < 2^0.5168 = {bound:.10f}; I2 = {I2:.10f} < {bound:.10f}",
        I1 < bound and I2 < bound,
    )
    o = continuous(1024)
    chk(
        f"(3.3): min D = {o['D'].min():.4f} >= 0 (the other two relative entropies are checked below)",
        o["D"].min() >= -1e-12,
    )
    # D-A >= 0 and D-E >= 0 via recomputation
    chk(
        f"(3.3): min (D-A) = {(o['D'] - o['A']).min():.4f} >= 0 and min (D-E) = {(o['D'] - o['E']).min():.4f} >= 0 (direct from the definitions)",
        (o["D"] - o["A"]).min() >= -1e-12 and (o["D"] - o["E"]).min() >= -1e-12,
    )
    # (A.2): D = B - L, u = L + (kappa-1) B - H, v = (2-kappa) B - 2 L + H, with L = log Q at the 64 nodes; the hats (64-point sums)
    # against the continuous quantities of (3.1)-(3.2) computed on 1024 nodes (every 16th node is a 64-point node)
    fl = lambda xs_: np.array([float(mid(x_)) for x_ in xs_])
    Lh, Bh, uh, vh = fl(H64["L"]), fl(H64["B"]), fl(H64["u"]), fl(H64["v"])
    Hh = float(mid(H64["H"]))
    k_ = KAPPA
    sl = slice(None, None, 16)
    eL = np.abs(np.log2(o["Q"][sl]) - Lh).max()
    eD = np.abs(o["D"][sl] - (Bh - Lh)).max()
    eu = np.abs(o["u"][sl] - (Lh + (k_ - 1) * Bh - Hh)).max()
    ev = np.abs(o["v"][sl] - ((2 - k_) * Bh - 2 * Lh + Hh)).max()
    chk("(A.2) L = log Q: 64-point hats agree with the continuous log2 Q at the 64 nodes to better than 1e-9", eL < 1e-9)
    chk("(A.2) D = B - L: continuous D agrees with Bhat - Lhat at the 64 nodes to better than 1e-9", eD < 1e-9)
    chk("(A.2) u = L + (kappa-1) B - H: continuous u agrees with the hat formula to better than 1e-9", eu < 1e-9)
    chk("(A.2) v = (2-kappa) B - 2 L + H: continuous v agrees with the hat formula to better than 1e-9", ev < 1e-9)
    if selftest:
        # deliberately wrong inputs
        w1 = continuous(512, c1=8.5 * 1.2)
        chk(
            f"WRONG INPUT: exponent 8.5 -> 10.2 makes I1 = {w1['I1']:.4f} (must NOT satisfy < 2^0.5168)",
            not (w1["I1"] < bound),
        )
        r_bad = tuple(list(R6[:5]) + [0.314])  # flipped sign of r6
        w2 = continuous(512, rr=r_bad)
        chk(
            f"WRONG INPUT: r6 sign flipped changes I1 to {w2['I1']:.6f} != {I1:.6f}",
            abs(w2["I1"] - I1) > 1e-6,
        )
        w3 = continuous(512, kappa=1.1)
        chk(
            f"WRONG INPUT: kappa 1.0315 -> 1.1 changes I2 to {w3['I2']:.6f} != {I2:.6f}",
            abs(w3["I2"] - I2) > 1e-6,
        )
        chk(
            "WRONG INPUT: a Table-1 entry perturbed by 0.01 exceeds its allowed error 4e-4",
            abs((TAB[10][3] + 0.01) - mid(H64["X1"][10])) > 4e-4,
        )
        chk("WRONG INPUT: claim mean Xhat < 1.40 is rejected", not (a1 < arb("1.40")))
    lowspace_exact(selftest)
    log(f"elapsed {time.time() - t0:.1f}s")
    return finish()


if __name__ == "__main__":
    sys.stdout.reconfigure(newline="\n")  # byte-identical output on Windows and Linux
    sys.exit(run())
