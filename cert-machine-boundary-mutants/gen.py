"""Minimal boundary mutants for the EinsteinArena circles-rectangle rows. Exact-rational checker (fractions)."""

import hashlib
import itertools
import json
import math
import os
from fractions import Fraction as F

import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# usage: python gen.py <cert-machine>/corpus/easota/circles-rectangle
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "src")


def dec(s):  # a JSON literal read as the decimal it prints (cert-machine's convention)
    return F(s)


def load_text(path):
    raw = open(path, encoding="utf8").read()
    d = json.loads(raw, parse_float=str, parse_int=str)
    return d["circles"]  # list of [str,str,str]


def constraints(C, reader):
    """every constraint's exact slack under the STRICT reading; C is a list of text triples"""
    Q = [[reader(v) for v in c] for c in C]
    out = {}
    xs0 = [x - r for x, y, r in Q]
    xs1 = [x + r for x, y, r in Q]
    ys0 = [y - r for x, y, r in Q]
    ys1 = [y + r for x, y, r in Q]
    out["box"] = 2 - ((max(xs1) - min(xs0)) + (max(ys1) - min(ys0)))
    for i, j in itertools.combinations(range(len(Q)), 2):
        out[f"pair{i}-{j}"] = (
            (Q[i][0] - Q[j][0]) ** 2
            + (Q[i][1] - Q[j][1]) ** 2
            - (Q[i][2] + Q[j][2]) ** 2
        )
    for i, (x, y, r) in enumerate(Q):
        out[f"r{i}>0"] = r
    return out


def failing(C, reader=dec):
    return sorted(
        k
        for k, v in constraints(C, reader).items()
        if (v <= 0 if k.endswith(">0") else v < 0)
    )


def platform(C):
    """the platform's 1e-9 rule (2026-03-26..2026-08-24), in binary64 exactly as the notebook / ts verifier"""
    c = [tuple(map(float, t)) for t in C]
    w = max(x + r for x, y, r in c) - min(x - r for x, y, r in c)
    h = max(y + r for x, y, r in c) - min(y - r for x, y, r in c)
    if w + h > 2 + 1e-9:
        return "FAIL"
    for a, b in itertools.combinations(c, 2):
        if math.hypot(a[0] - b[0], a[1] - b[1]) < a[2] + b[2] - 1e-9:
            return "FAIL"
    return "PASS"


def _ord(x):
    i = int(np.array([x], dtype=np.float64).view(np.int64)[0])
    return i if i >= 0 else -(i & 0x7FFFFFFFFFFFFFFF)

def _unord(i):
    b = i if i >= 0 else (-i) | (1 << 63)
    return float(np.array([b], dtype=np.uint64).view(np.float64)[0])

def step(s, k, toward):
    """k ulps toward +/-inf: k applications of numpy.nextafter, done on the ordered bit pattern
    (cross-checked against literal nextafter iteration for k <= 64)"""
    x = float(s)
    y = _unord(_ord(x) + (k if toward > 0 else -k))
    if k <= 64:
        z = np.float64(x)
        for _ in range(k): z = np.nextafter(z, toward)
        assert float(z) == y, (s, k)
    return repr(y)

def mutate(C, ci, comp, toward, k):
    D = [list(t) for t in C]
    D[ci][comp] = step(C[ci][comp], k, toward)
    return D


def smallest(C, ci, comp, toward, target):
    """smallest k >= 1 of nextafter steps on C[ci][comp] making `target` fail"""
    fails = lambda k: target in failing(mutate(C, ci, comp, toward, k))
    hi = 1
    while not fails(hi):
        hi *= 2
        if hi > 1 << 40:
            return None
    lo = hi // 2  # fails(lo) false (or lo=0)
    while hi - lo > 1:
        m = (lo + hi) // 2
        if fails(m):
            hi = m
        else:
            lo = m
    return hi


def repaired_base(C):
    """the cert-machine REPAIR idea in binary64: every radius scaled down, rounded down, until the strict reading holds"""
    Q = [[dec(v) for v in c] for c in C]
    lam = min(
        math.sqrt(
            ((Q[i][0] - Q[j][0]) ** 2 + (Q[i][1] - Q[j][1]) ** 2)
            / (Q[i][2] + Q[j][2]) ** 2
        )
        for i, j in itertools.combinations(range(len(Q)), 2)
    )
    lam = min(lam, 1.0)
    while True:
        D = [
            [
                c[0],
                c[1],
                repr(float(np.nextafter(np.float64(float(F(c[2]) * F(lam))), 0))),
            ]
            for c in C
        ]
        if not failing(D):
            return D, lam
        lam = float(np.nextafter(lam, 0)) - 1e-16


def candidates(C):
    """per base: the two tightest pairs, and the box pushed out in x and in y"""
    Q = [[dec(v) for v in c] for c in C]
    s = constraints(C, dec)
    pairs = sorted((v, k) for k, v in s.items() if k.startswith("pair"))[:2]
    out = []
    for v, k in pairs:
        i, j = map(int, k[4:].split("-"))
        ax = 0 if abs(Q[i][0] - Q[j][0]) >= abs(Q[i][1] - Q[j][1]) else 1
        toward = math.inf if Q[i][ax] > Q[j][ax] else -math.inf
        out.append(
            (
                k,
                "pair overlap",
                [(j, ax, toward), (i, ax, -toward), (j, 2, math.inf), (i, 2, math.inf)],
            )
        )
    for ax, name in ((0, "x"), (1, "y")):
        hi = max(range(len(Q)), key=lambda t: Q[t][ax] + Q[t][2])
        lo = min(range(len(Q)), key=lambda t: Q[t][ax] - Q[t][2])
        out.append(
            (
                "box",
                f"box overflow in {name}",
                [(hi, ax, math.inf), (lo, ax, -math.inf), (hi, 2, math.inf)],
            )
        )
    return out


def main():
    rows = {}
    a = load_text(os.path.join(SRC, "alphaevolve_2025.json"))
    o = load_text(os.path.join(SRC, "ours_2026.json"))
    assert failing(a) == [], failing(a)
    of = failing(o)
    print("ours_2026 strict failures as published:", len(of))
    ob, lam = repaired_base(o)
    rows["circles__alphaevolve_2025"] = (a, "pinned bytes (strictly feasible)")
    rows["circles__ours_2026"] = (
        ob,
        f"pinned bytes fail strictly on {len(of)} constraints, so the base is their binary64 repair: every radius times lambda={lam!r}, rounded down; strictly feasible",
    )
    man = []
    for row, (C, basenote) in rows.items():
        base_path = os.path.join(HERE, row + "__BASE.json")
        write(base_path, C)
        n = 0
        for target, kind, moves in candidates(C):
            for ci, comp, toward in moves:
                k = smallest(C, ci, comp, toward, target)
                if k is None:
                    continue
                M = mutate(C, ci, comp, toward, k)
                f = failing(M)
                fb = failing(
                    M, reader=lambda s: F(float(s))
                )  # second reading: the binary64 the text rounds to
                if f != [target] or fb != [target]:
                    continue
                if k > 1:  # minimality re-check: k-1 must not fail the target
                    assert target not in failing(mutate(C, ci, comp, toward, k - 1))
                n += 1
                mid = f"m{n:02d}"
                p = os.path.join(HERE, row + "__" + mid + ".json")
                write(p, M)
                man.append(
                    {
                        "id": f"{row}/{mid}",
                        "row": row.replace("__", "/"),
                        "base": basenote,
                        "constraint": target,
                        "kind": kind,
                        "change": {
                            "circle": ci,
                            "field": "xyr"[comp],
                            "from": C[ci][comp],
                            "to": M[ci][comp],
                            "nextafterSteps": k,
                            "direction": "+" if toward > 0 else "-",
                        },
                        "strictSlack": str(constraints(M, dec)[target]),
                        "strictVerdict": "FAIL (only " + target + ")",
                        "platformVerdict_1e-9_rule_2026-03-26_to_2026-08-24": platform(
                            M
                        ),
                        "platformVerdict_exact_rule_since_2026-08-24": "FAIL",
                        "sha256": hashlib.sha256(open(p, "rb").read()).hexdigest(),
                    }
                )
                break
    json.dump(man, open(os.path.join(HERE, "MANIFEST.json"), "w"), indent=1)
    print(len(man), "mutants")


def write(p, C):
    # his input format: {"circles": [[x, y, r], ...]} with JSON numbers, literal text preserved
    body = ",\n".join("    [" + ", ".join(c) + "]" for c in C)
    open(p, "w", encoding="utf8", newline="\n").write(
        '{\n  "circles": [\n' + body + "\n  ]\n}\n"
    )


if __name__ == "__main__":
    main()
