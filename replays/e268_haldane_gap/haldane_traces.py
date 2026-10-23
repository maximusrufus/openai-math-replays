"""Thermal tables of Lemma 5.1 of 'The periodic spin-one Haldane gap' (module of replay_haldane.py):
Z_n(b,g) = Tr exp[-b (H_n(g) + a n)] for the periodic spin-one Heisenberg ring with twisted seam (g = I, P (theta = pi), C (theta = 2 pi/3)), a = 700741/500000.
Written from the manuscript text only: weight basis, diagonal sum_j x_{j-1} x_j, steps changing two neighbouring letters by opposite units
with coefficient one, seam steps with coefficient exp(i theta (y0 - x0)) (eq. (2.5)).
Method: exact diagonalisation in float64 (NOT a rigorous enclosure).  The seam phase is gauged into a uniform hopping phase exp(i theta/n); S^z
and momentum blocks are then diagonalised (dense, or lowest eigenvalues by Lanczos below a cut E + a n < 2.5; the discarded eigenvalues
number at most 3^n and each contributes less than e^(-12.25*2.5) at b = 49/4)."""

import itertools
import sys
import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

A_SHIFT = 700741 / 500000


# ---------------------------------------------------------------------------------------- direct (ungauged) reference, small n
def direct_spectrum(n, theta):
    """eigenvalues of H_n(g) built literally from (2.5): returns array of all 3^n eigenvalues"""
    words = list(itertools.product((-1, 0, 1), repeat=n))
    idx = {w: i for i, w in enumerate(words)}
    N = len(words)
    H = np.zeros((N, N), dtype=complex)
    for w in words:
        i = idx[w]
        H[i, i] += sum(w[(j - 1) % n] * w[j] for j in range(n))
        for j in range(n):  # bond (j, j+1 mod n); seam bond is j = n-1
            j2 = (j + 1) % n
            for ds in (1, -1):
                # column y = w, row x = w + ds at j, -ds at j2
                x = list(w)
                x[j] += ds
                x[j2] -= ds
                if not (-1 <= x[j] <= 1 and -1 <= x[j2] <= 1):
                    continue
                coef = 1.0
                if (
                    j == n - 1
                ):  # seam step: coefficient exp(i theta (y0 - x0)), site 0 is j2
                    coef = np.exp(1j * theta * (w[0] - x[0]))
                H[idx[tuple(x)], i] += coef
    assert np.allclose(H, H.conj().T)
    return np.linalg.eigvalsh(H)


# ---------------------------------------------------------------------------------------- symmetry-reduced spectrum
def all_codes(n):
    N = 3**n
    codes = np.arange(N, dtype=np.int64)
    digs = np.zeros((N, n), dtype=np.int8)
    c = codes.copy()
    for j in range(n):
        digs[:, j] = (c % 3) - 1
        c //= 3
    return codes, digs


def reduced_spectrum(n, theta, cutoff=None, verbose=True, lanczos=False):
    """returns dict: m -> list of eigenvalue arrays (all momentum blocks), for m = 0..n.  Eigenvalues above `cutoff` are dropped (if given)."""
    phi = theta / n
    N = 3**n
    pw = 3 ** np.arange(n, dtype=np.int64)
    codes, digs = all_codes(n)
    sums = digs.sum(axis=1).astype(np.int32)
    # rotation: (T^t x)_j = x_{(j-t) mod n}; code of rotated word
    rot = np.zeros((n, N), dtype=np.int64)
    for t in range(n):
        d = np.roll(digs, t, axis=1)  # d[:, j] = digs[:, (j - t) % n]
        rot[t] = ((d.astype(np.int64) + 1) * pw[None, :]).sum(axis=1)
    rep = rot.min(axis=0)
    t0 = rot.argmin(axis=0)  # first t with rot[t] == rep :  T^{t0} y = rep
    # orbit size of each representative: smallest t > 0 with rot[t][code] == code
    period = np.full(N, n, dtype=np.int64)
    for t in range(n, 0, -1):
        period = (
            np.where(rot[t % n] == codes, t if t < n else period, period)
            if t < n
            else period
        )
    # (the loop above leaves period = n unless a smaller t reproduces the word; iterate t descending so the smallest wins)
    period = np.full(N, n, dtype=np.int64)
    for t in range(n - 1, 0, -1):
        period = np.where(rot[t] == codes, t, period)
    out = {}
    for m in range(n + 1):
        sel = np.where((sums == m) & (rep == codes))[0]
        reps = codes[sel]
        R = period[sel]
        nrep = len(reps)
        pos = -np.ones(N, dtype=np.int64)
        pos[reps] = np.arange(nrep)
        # build the action of H' on every representative: lists (row_rep_index, col_rep_index, coefficient, shift t0 of the target word)
        d = digs[sel].astype(np.int64)  # (nrep, n)
        rows, cols, cf, shifts = [], [], [], []
        diag = np.zeros(nrep)
        for j in range(n):
            j2 = (j + 1) % n
            diag += d[:, j] * d[:, j2]
        for j in range(n):
            j2 = (j + 1) % n
            for ds, ph in ((1, np.exp(1j * phi)), (-1, np.exp(-1j * phi))):
                ok = (
                    (d[:, j] + ds >= -1)
                    & (d[:, j] + ds <= 1)
                    & (d[:, j2] - ds >= -1)
                    & (d[:, j2] - ds <= 1)
                )
                a = np.where(ok)[0]
                newcode = reps[a] + ds * pw[j] - ds * pw[j2]
                b = pos[rep[newcode]]
                assert (b >= 0).all()
                rows.append(b)
                cols.append(a)
                cf.append(np.full(len(a), ph))
                shifts.append(t0[newcode])
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        cf = np.concatenate(cf)
        shifts = np.concatenate(shifts)
        eigs = []
        for jm in range(n):
            allowed = (jm * R) % n == 0
            if not allowed.any():
                continue
            amap = -np.ones(nrep, dtype=np.int64)
            amap[allowed] = np.arange(allowed.sum())
            keep = allowed[rows] & allowed[cols]
            r_, c_ = rows[keep], cols[keep]
            # matrix element <b,k|H|a,k> = coef * exp(+i k t0) * sqrt(R_a / R_b)    (standard momentum-block formula)
            val = (
                cf[keep]
                * np.exp(2j * np.pi * jm * shifts[keep] / n)
                * np.sqrt(R[c_] / R[r_])
            )
            M = sp.coo_matrix(
                (val, (amap[r_], amap[c_])), shape=(allowed.sum(), allowed.sum())
            ).tocsr()
            M = M + sp.diags(diag[allowed].astype(complex))
            dim = M.shape[0]
            if lanczos and dim > 400:
                import scipy.sparse.linalg as spl

                kk = 40
                while True:
                    kk = min(kk, dim - 2)
                    w = np.sort(spl.eigsh(M, k=kk, which="SA", tol=1e-13, ncv=min(dim - 1, max(2 * kk + 1, 80)), return_eigenvectors=False))
                    if w[-1] >= cutoff or kk >= dim - 2:
                        break
                    kk *= 2
                w = w[w < cutoff]
            else:
                Md = M.toarray()
                herm = np.abs(Md - Md.conj().T).max()
                assert herm < 1e-9, herm
                if cutoff is None:
                    w = sla.eigvalsh(Md)
                else:
                    w = sla.eigvalsh(Md, subset_by_value=(-np.inf, cutoff))
            eigs.append(w)
        out[m] = np.concatenate(eigs)
        if verbose:
            print(
                f"   n={n} theta={theta:.4f} m={m}: sector dim {int((sums == m).sum())}, reps {nrep}, eigenvalues kept {len(out[m])}",
                flush=True,
            )
    return out


def partition(out, n, b, weights=True):
    tot = 0.0
    for m, ev in out.items():
        mult = 1 if m == 0 else 2
        tot += mult * np.exp(-b * (ev + A_SHIFT * n)).sum()
    return tot


TABLE_21_2 = {  # b = 21/2: n -> (Z, ZP, ZC) in 1e-6
    6: (8945303, 346734, 939987),
    8: (3741281, 567167, 984108),
    10: (2327547, 721090, 1002613),
}
TABLE_49_4 = {  # b = 49/4: n -> (Z, ZP) in 1e-6
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
