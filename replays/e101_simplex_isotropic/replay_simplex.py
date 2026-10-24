#!/usr/bin/env python3
"""Independent replay of the finite scalar estimates in 'A sharp entropy bound and the simplex inequality for isotropic constants'
(OpenAI, 2026-10; Section 2.1, Lemma 2.2, Lemma 4.1, Appendix A).

Written from the manuscript text only (definitions (2.2), (2.3), (3.8), (3.14), (3.16), (3.24), (4.1), (4.5)); no code from the
authors is used.  Arithmetic: python-flint Arb (ball arithmetic, 150 bits) with Taylor-model enclosures inside a bisection prover:
a claim is REPLAYED only when the prover decides it on every ball of the stated interval.

Scope: the compact range a in [-6, 10] (the manuscript treats a <= -6 and a >= 10 by analytic tail estimates; those are NOT replayed),
except (4.6) for a >= -1, which is proved here on [-1, 6] only (the method does not finish on (6, 10] in reasonable time).

  python replay_simplex.py             # full run (about 10-15 min on one core); exit 0 only if all claims are proved
  python replay_simplex.py --quick     # documented subset (skips the two longest claims: T1 table row and (4.6) a >= -1; about 1 min)
  python replay_simplex.py --perturb   # adversarial: the step bound m0 on U1 raised to .75; must FAIL (exit 1)
  python replay_simplex.py --selftest  # the built-in wrong-claim controls (each must be rejected)

Deterministic stdout (no timing); elapsed time goes to stderr.
"""

import re
import sys
import time

from flint import arb, ctx

RES = []
SELFTEST = "--selftest" in sys.argv
PERTURB = "--perturb" in sys.argv
PFX = "S"
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


QUICK = "--quick" in sys.argv

ctx.prec = 150
SQ2 = arb(2).sqrt()
SQ2PI = (2 * arb.pi()).sqrt()
R = lambda p, q=1: arb(p) / arb(q)

# constants (3.8), (3.14), (3.24), (4.1)
ALPHA = R(3, 40)
S = 1 / (1 - ALPHA)
BETA = R(7, 25)
CC = R(3, 2)
GAMMA = CC / 2
EPS = R(3, 20)
KAPPA = R(49, 100)
Dstar = 1 / (4 * (1 - EPS))
Astar = S + 1 + S / 8
etastar = S * BETA + (1 + S / 8) * GAMMA
nustar = S * BETA**2 + (9 * S / 8) * GAMMA**2
Qm = [[R(43, 100), R(13, 100)], [R(13, 100), R(272, 100)]]
# the three lower/upper bounds of (4.2) that are used as numbers in claim_42
C42 = dict(h0=R(67, 100), p_pos=R(25, 100), p_a4=R(46, 100))


def Phi(x):
    return (-x / SQ2).erfc() / 2


def phi(x):
    return (-(x * x) / 2).exp() / SQ2PI


NTERMS = 8  # Taylor order of the enclosures


def _series_fns(base, N):
    """All scalar functions of the paper as truncated power series in x, where the argument is a = base + x (base an arb point or ball)."""
    from flint import arb_series

    a = arb_series([base, 1], prec=N)
    P = (a / SQ2).erfc() / 2  # Phi(-a)
    t = -P.log()
    ph = (-(a * a) / 2).exp() / SQ2PI
    q = ph / P
    d = q - a
    k = t / q
    b = 1 - k * d
    l = k - q + a * b
    lp = 2 * b - q * d + a * l
    f = k * (1 + b)
    fp = b * (1 + b) + k * l
    h0 = f - CC * k - 2 * (b - BETA) * l / ALPHA
    p = fp - b * b - l * l / ALPHA
    e = fp - 2 * b * b - 2 * l * l / ALPHA
    E = e + GAMMA * b
    NAN = arb_series([arb("nan")], prec=N)

    def sf(fn):
        try:
            return fn()
        except (ValueError, ZeroDivisionError):
            return NAN

    D = sf(lambda: Dstar / p)
    dlogM = sf(lambda: 2 * (1 - q * d) / (d * d) - 2 * b * l / (d * (1 - b * b)))
    z = sf(lambda: lp / l)
    hstar = 1 + 3 * b + k * z
    chi = sf(lambda: l / d)
    return dict(q=q, d=d, k=k, b=b, l=l, lp=lp, p=p, h0=h0, e=e, E=E, D=D, L=l / ALPHA, K0=D * E, J0=D * CC * b, dlogM=dlogM, z=z, hstar=hstar,
                chi=chi)


_CACHE = {}


def fns_tm(lo, hi, N=NTERMS):
    key = (lo, hi, N)
    if key in _CACHE:
        v = _CACHE[key]
        if isinstance(v, Exception):
            raise v
        return v
    try:
        out = _fns_tm(lo, hi, N)
    except (ValueError, ZeroDivisionError) as ex:
        _CACHE[key] = ex
        raise
    _CACHE[key] = out
    return out


def _fns_tm(lo, hi, N=NTERMS):
    """Taylor-model enclosure of every function over [lo, hi]: sum_{k<N} c_k(m) h^k + c_N(ball) h^N, h in [-w, w], where c_k(m) are the
    Taylor coefficients at the midpoint m and c_N(ball) the N-th coefficient computed over the whole interval (Lagrange remainder)."""
    m = arb((lo + hi) / 2)
    w = arb((hi - lo) / 2)
    h = arb(0, w.upper())
    pt = _series_fns(m, N + 1)
    # remainder coefficients: union over M sub-balls of [m-w, m+w] (M doubled until no divisor ball touches 0)
    M = 1
    while True:
        try:
            pieces = []
            wr = w.upper()
            for j in range(M):
                cj = m - wr + wr * (2 * j + 1) / M
                pieces.append(_series_fns(arb(cj, (wr / M).upper()), N + 1))
            break
        except (ValueError, ZeroDivisionError):
            M *= 2
            if M > 32:
                raise
    out = {}
    for name in pt:
        cp = pt[name].coeffs()
        cbN = None
        for pc in pieces:
            cc = pc[name].coeffs()
            cc = cc[N] if len(cc) > N else arb(0)
            cbN = cc if cbN is None else cbN.union(cc)
        cp = list(cp) + [arb(0)] * (N + 1 - len(cp))
        val = arb(0)
        for kk in range(N):
            val += cp[kk] if kk == 0 else cp[kk] * arb(0, w.upper() ** kk)
        val += cbN * arb(0, w.upper() ** N)
        out[name] = val
    out["a"] = arb(m, w.upper())
    return out


def tern(true_cond, false_cond):
    return True if true_cond else (False if false_cond else None)


def prove(lo, hi, pred, depth=0, maxdepth=34):
    if depth == 0:
        t0 = time.time()
        r = _prove(lo, hi, pred, 0, maxdepth)
        log('  .. prove', lo, hi, r[0], 'nodes', r[2], round(time.time() - t0, 1), 's')
        return r
    return _prove(lo, hi, pred, depth, maxdepth)


def _prove(lo, hi, pred, depth=0, maxdepth=34):
    """pred(fns_dict, a_ball) -> (certainly_true, certainly_false) booleans. Returns (ok, worst_counterexample_or_None, nodes)."""
    a = arb((lo + hi) / 2, (hi - lo) / 2)
    try:
        ct, cf = pred(fns_tm(lo, hi), a)
    except (ValueError, ZeroDivisionError):  # a divisor ball contains 0 on this (too wide) interval: split
        ct, cf = False, False
    if ct:
        return True, None, 1
    if cf and (hi - lo) < 1e-3:
        return False, (lo, hi), 1
    if depth >= maxdepth:
        return False, (lo, hi, "undecided"), 1
    mid = (lo + hi) / 2
    ok1, c1, n1 = _prove(lo, mid, pred, depth + 1, maxdepth)
    if not ok1:
        return False, c1, n1
    ok2, c2, n2 = _prove(mid, hi, pred, depth + 1, maxdepth)
    return ok2, c2, n1 + n2


def le(x, y):  # certainly x <= y / certainly x > y
    return (x <= y), (x > y)


def ge(x, y):
    return (x >= y), (x < y)


def all_of(*pairs):
    return all(p[0] for p in pairs), any(p[1] for p in pairs)


# ----- the interval table (4.3) and (2.5): (name, lo, hi)
INTERVALS = {
    "O": (-6.0, -3.0),
    "P1": (-3.0, -2.0),
    "R1": (-2.0, -1.0),
    "U1": (-1.0, 0.5),
    "V1": (0.5, 4.0),
    "T1": (4.0, 10.0),
}
M0 = {
    "O": R(0),
    "P1": R(19, 100),
    "R1": R(38, 100),
    "U1": R(68, 100),
    "V1": R(68, 100),
    "T1": R(68, 100),
}
# (4.3): b range, L<=, D<= (O: .295(7+a^2)), K0 range, J0<=, ell1
TAB = {
    "O": dict(
        b=(R(0), R(9, 100)),
        L=R(60, 100),
        D=None,
        K0=(R(37, 100), R(56, 100)),
        J0=R(72, 100),
        l1=R(13, 100),
    ),
    "P1": dict(
        b=(R(83, 1000), R(15, 100)),
        L=R(101, 100),
        D=R(43, 10),
        K0=(R(29, 100), R(46, 100)),
        J0=R(72, 100),
        l1=R(13, 100),
    ),
    "R1": dict(
        b=(R(14, 100), R(23, 100)),
        L=R(119, 100),
        D=R(3),
        K0=(R(29, 100), R(371, 1000)),
        J0=R(68, 100),
        l1=R(13, 100),
    ),
    "U1": dict(
        b=(R(22, 100), R(342, 1000)),
        L=R(119, 100),
        D=R(172, 100),
        K0=(R(29, 100), R(398, 1000)),
        J0=R(61, 100),
        l1=R(123, 1000),
    ),
    "V1": dict(
        b=(R(337, 1000), R(45, 100)),
        L=R(83, 100),
        D=R(87, 100),
        K0=(R(31, 100), R(4, 10)),
        J0=R(5, 10),
        l1=R(60, 1000),
    ),
    "T1": dict(
        b=(R(445, 1000), R(5, 10)),
        L=R(2, 10),
        D=R(64, 100),
        K0=(R(32, 100), R(4, 10)),
        J0=R(48, 100),
        l1=R(22, 1000),
    ),
}
# (4.6): z in range, h* range, chi <= ; a <= -1 and a >= -1
TAB46 = {
    "a<=-1": (
        (-6.0, -1.0),
        (R(-8, 100), R(70, 100)),
        (R(8, 10), R(177, 100)),
        R(70, 1000),
    ),
    "a>=-1": (
        (-1.0, 6.0),
        (R(-475, 1000), R(0)),
        (R(8, 10), R(177, 100)),
        R(97, 1000),
    ),
}


def claim_table(name, tab, ell_scale=1, b_hi_shift=0):
    lo, hi = INTERVALS[name]
    T = tab[name]

    def pred(F, a):
        D_ok = (
            le(F["D"], T["D"])
            if T["D"] is not None
            else le(F["D"], R(295, 1000) * (7 + a * a))
        )
        checks = [
            ge(F["b"], T["b"][0]),
            le(F["b"], T["b"][1] + b_hi_shift),
            le(F["L"], T["L"]),
            D_ok,
            ge(F["K0"], T["K0"][0]),
            le(F["K0"], T["K0"][1]),
            le(F["J0"], T["J0"]),
        ]
        # (4.4): matrix inequality, 2x2 PSD test of N - diag(ell1,.75)
        l2a = F["l"] ** 2 / ALPHA
        Dv, Ev, cb = F["D"], F["E"], CC * F["b"]
        m11 = (
            Qm[0][0]
            - (Astar * F["b"] ** 2 - etastar * F["b"] + nustar / 4 + l2a)
            - Dv * Ev * Ev
            - T["l1"] * ell_scale
        )
        m12 = Qm[0][1] - (nustar / 2 - etastar * F["b"]) - Dv * Ev * cb
        m22 = Qm[1][1] - nustar - Dv * cb * cb - R(75, 100)
        det = m11 * m22 - m12 * m12
        checks += [ge(m11, R(0)), ge(m22, R(0)), ge(det, R(0))]
        return all_of(*checks)

    return prove(lo, hi, pred)


def claim_m0(name, scale=1):
    lo, hi = INTERVALS[name]
    return prove(lo, hi, lambda F, a: ge(F["dlogM"], M0[name] * scale))


def claim_scalar(lo, hi, expr):
    return prove(lo, hi, expr)


def claim_42(shift=0):
    """(4.2) and Lemma 2.2: p>0, |h0| <= .67 sqrt(p), p >= 1/(7+a^2) (a<=0), p >= .25 (a>=0), p >= .46 (a>=4), 0<b<1/2, l>0."""
    res = {}
    res["0<b<1/2, l>0"] = prove(
        -6.0,
        10.0,
        lambda F, a: all_of(ge(F["b"], R(0)), le(F["b"], R(1, 2)), ge(F["l"], R(0))),
    )
    res["p>0"] = prove(-6.0, 10.0, lambda F, a: ge(F["p"], R(0)))
    res["|h0|<=.67 sqrt p"] = prove(
        -6.0, 10.0, lambda F, a: le(F["h0"] ** 2, C42["h0"] ** 2 * F["p"] + shift)
    )
    res["p>=1/(7+a^2), a<=0"] = prove(
        -6.0, 0.0, lambda F, a: ge(F["p"] * (7 + a * a), R(1) + shift)
    )
    res["p>=.25, a>=0"] = prove(0.0, 10.0, lambda F, a: ge(F["p"], C42["p_pos"] + shift))
    res["p>=.46, a>=4"] = prove(4.0, 10.0, lambda F, a: ge(F["p"], C42["p_a4"] + shift))
    return res


def claim_46(key, shrink=0):
    (lo, hi), hr, hs, chimax = TAB46[key]
    return prove(
        lo,
        hi,
        lambda F, a: all_of(
            ge(F["z"], hr[0]),
            le(F["z"], hr[1]),
            ge(F["hstar"], hs[0]),
            le(F["hstar"], hs[1]),
            le(F["chi"], chimax - shrink),
        ),
    )


def exact_checks():
    """small exact identities of the text (Fractions)."""
    from fractions import Fraction as Fr

    out = {}
    # (3.29): Q(1,u-1/2) + 7u/50 <= 1-u on [0,1/2], difference = (1-2u)(1+68u)/50
    import sympy as sp

    u = sp.symbols("u")
    Q = lambda x, y: (
        sp.Rational(43, 100) * x**2
        + sp.Rational(26, 100) * x * y
        + sp.Rational(68, 25) * y**2
    )
    diff = sp.expand(1 - u - (Q(1, u - sp.Rational(1, 2)) + sp.Rational(7, 50) * u))
    out["(3.29) factorisation"] = (
        sp.simplify(diff - sp.Rational(1, 50) * (1 - 2 * u) * (1 + 68 * u)) == 0
    )
    # (4.7): sqrt( s ( (b-beta)^2 + alpha h0^2 / (4 (1-eps-kappa) p) ) ) with 0 < b < 1/2 and |h0| <= .67 sqrt(p):
    # the manuscript's middle expression is s (.28^2 + .075 * .67^2 / (4 * .36)) < .111
    s = Fr(1) / (1 - Fr(3, 40))
    beta = Fr(7, 25)
    bmax2 = max(beta**2, (Fr(1, 2) - beta) ** 2)  # max of (b - beta)^2 over [0, 1/2]
    out["(4.7) max_b (b-beta)^2 over [0,1/2] equals .28^2"] = bmax2 == Fr(28, 100) ** 2
    val = s * (
        Fr(28, 100) ** 2
        + Fr(3, 40) * Fr(67, 100) ** 2 / (4 * (1 - Fr(3, 20) - Fr(49, 100)))
    )
    out["(4.7) s(.28^2 + alpha .67^2/(4(1-eps-kappa))) = .11003 < .111; and .075/(4*.36) form: 4(1-eps-kappa) = 1.44"] = val < Fr(111, 1000) and 4 * (1 - Fr(3, 20) - Fr(49, 100)) == Fr(144, 100)
    # (3.8)-(3.24): the stated rational constants
    out["0 < alpha, 1-eps-kappa > 0"] = Fr(3, 20) + Fr(49, 100) < 1
    # the constants used by the claims equal the printed ones (second, independent transcription)
    def eq(x, num, den):
        return abs(x - R(num, den)) < arb("1e-40")

    out["stated constants: alpha = 3/40, beta = 7/25, c = 3/2, eps = 3/20, kappa = 49/100"] = (
        eq(ALPHA, 3, 40) and eq(BETA, 7, 25) and eq(CC, 3, 2) and eq(EPS, 3, 20) and eq(KAPPA, 49, 100)
    )
    out["stated constants: Q = [[.43, .13], [.13, 2.72]] (Proposition 3.6)"] = (
        eq(Qm[0][0], 43, 100) and eq(Qm[0][1], 13, 100) and eq(Qm[1][0], 13, 100) and eq(Qm[1][1], 272, 100)
    )
    out["stated constants of (4.2): .67, .25, .46"] = eq(C42["h0"], 67, 100) and eq(C42["p_pos"], 25, 100) and eq(C42["p_a4"], 46, 100)
    out["the six intervals O, P1, R1, U1, V1, T1 tile [-6, 10] without gap or overlap"] = (
        INTERVALS["O"][0] == -6.0
        and INTERVALS["T1"][1] == 10.0
        and all(INTERVALS[x][1] == INTERVALS[y][0] for x, y in (("O", "P1"), ("P1", "R1"), ("R1", "U1"), ("U1", "V1"), ("V1", "T1")))
    )
    out["stated (2.5)-(2.6) step table m0 = 0, .19, .38, .68, .68, .68"] = all(
        eq(M0[k], n, d) for k, n, d in (("O", 0, 1), ("P1", 19, 100), ("R1", 38, 100), ("U1", 68, 100), ("V1", 68, 100), ("T1", 68, 100))
    )
    return out


def report(name, res):
    ok, w, n = res
    chk(f"{name}" + ("" if ok or w is None else f"  witness {w}") + (f" [nodes {n}]" if n else ""), ok)


def run():
    t0 = time.time()
    if PERTURB:  # perturb the INPUT (step bound m0 on U1 raised from .68 to .75), then run the unchanged main checks
        M0["U1"] = R(75, 100)
    if SELFTEST:
        wrong = [
            ("m0 on U1 raised to .75", claim_m0("U1", scale=R(75, 68))),
            ("m0 on P1 scaled x2", claim_m0("P1", scale=2)),
            ("(4.4) on V1 with ell1 x3", claim_table("V1", TAB, ell_scale=3)),
            ("(4.3) b upper bound on P1 lowered by .05", claim_table("P1", TAB, b_hi_shift=R(-5, 100))),
            ("(4.2) p >= .46 replaced by p >= .6 on a>=4", prove(4.0, 10.0, lambda F, a: ge(F["p"], R(6, 10)))),
        ]
        if not QUICK:
            wrong.append(("(4.6) a>=-1 chi bound tightened by .09", claim_46("a>=-1", shrink=R(9, 100))))
        for name, (ok, w, n) in wrong:
            chk("WRONG INPUT: " + name + " is rejected", not ok)
        log(f"elapsed {time.time() - t0:.1f}s")
        return finish()
    for k, res in claim_42().items():
        report("Lemma 2.2/4.1 (4.2) " + k, res)
    for name in INTERVALS:
        report(f"Lemma 2.2 (2.5)-(2.6): dlogM/dlogq >= m0 on {name}", claim_m0(name))
    for name in INTERVALS:
        if QUICK and name == "T1":
            continue
        report(f"Lemma 4.1 (4.3)+(4.4) on {name} {INTERVALS[name]}", claim_table(name, TAB))
    for key in TAB46:
        if QUICK and key == "a>=-1":
            continue
        report(f"Lemma 4.1 (4.6) {key} on {TAB46[key][0]}", claim_46(key))
    for k, v in exact_checks().items():
        chk("exact: " + k, v)
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
