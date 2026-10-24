# Independent replays of computational claims in the openai/math October 2026 release

Purpose: independently re-run the finite, checkable computations that are stated in the text of manuscripts in
[openai/math](https://github.com/openai/math), using exact arithmetic and code written from the paper text only.
This is not a review of the full proofs. A replay shows that the stated computations are consistent with the
manuscript; it does not show that a theorem is true.

Credit: the manuscripts, results and all displayed formulas belong to their original authors (OpenAI, as listed in each
manuscript). This repository contains only independent re-computation code and its output.

## Replays

| Entry (manuscript, as listed in openai/math) | Replay | Result | Coverage notes |
|---|---|---|---|
| Wall's D(2) problem counterexample (2026-10-06) | [`replays/wall_d2/`](replays/wall_d2/) | 21 of 21 replayed claims pass | |
| 059 Ambiently homeomorphic isolated hypersurface germs in C^4 with multiplicities four and five (2026-09-27) | [`replays/e059_zariski_multiplicity/`](replays/e059_zariski_multiplicity/) | 17 of 17 | homological and topological arguments not replayed |
| 067 The Campana-Peternell conjecture in dimension six (2026-09-25) | [`replays/e067_campana_peternell/`](replays/e067_campana_peternell/) | 33 of 33 | geometric derivation of the constraints not replayed; slow (3-6 min) |
| 101 A sharp entropy bound and the simplex inequality for isotropic constants (2026-10-05) | [`replays/e101_simplex_isotropic/`](replays/e101_simplex_isotropic/) | 29 of 29 (full), 27 of 27 (quick subset) | range a in [-6, 10] only; the row a >= -1 of (4.6) only on [-1, 6]; slow (about 10-15 min) |
| 138 Subset Sum in Time O(2^{0.49n}) (2026-10-04) and A Low-Space Algorithm for Worst-Case Subset Sum (2026-09-26) | [`replays/e138_subset_sum/`](replays/e138_subset_sum/) | 26 of 26 | numerical integrals, not a rigorous quadrature enclosure |
| 153 Arithmetic classification and non-Pisot singularity for Bernoulli convolutions (2026-10-03) | [`replays/e153_bernoulli_non_pisot/`](replays/e153_bernoulli_non_pisot/) | 48 of 48 | the large integers of (A.17)/(A.23) are matched at value level only |
| 268 The periodic spin-one Haldane gap (2026-09-24) | [`replays/e268_haldane_gap/`](replays/e268_haldane_gap/) | 68 of 68 | thermal tables replayed by float64 diagonalisation (not a rigorous enclosure) |
| 318 A rational obstruction to strong chromatic splitting at height three (2026-09-25) | [`replays/e318_chromatic_splitting/`](replays/e318_chromatic_splitting/) | 15 of 15 | only the Chevalley-Eilenberg certificate; no other finite computation was identified in the entry's manuscripts |
| 352 A finite-time singularity of Calabi flow on projective space (2026-09-24) | [`replays/e352_calabi_flow/`](replays/e352_calabi_flow/) | 28 of 28 | the shooting estimate is checked at 8 sample points of the parameter square, not over the whole square |
| Released Ramsey lower-bound graphs of arXiv 2603.09172 v2 (Google, 18 files) | [`ramsey_alphaevolve/`](ramsey_alphaevolve/) | 13 of 18 graphs verified by two methods (A, B); 2 more verified on the independent-set side by one solver family; 3 not verified (too large) | the graphs are Google's and are copied under Apache-2.0; only the stated clique / independent-set property is re-checked, not how the graphs were found; R(4,15): paper text says 159, released file is 163 (see README there) |
| MUB N(6)=3 | separate, already public: [maximusrufus/mub6-certified-rerun](https://github.com/maximusrufus/mub6-certified-rerun) (not copied here) | see that repository | |

## Pinned source (Wall D(2))

- Repository: `https://github.com/openai/math`
- Commit: `1581d4cbdd69a0faf34b2089614fe2a814316a6a` (2026-10-12, "Merge pull request #2 from openai/codex/lock-prs")
- Manuscript path: `preprints/A-Counterexample-to-Walls-D2-Problem-October-6-2026/wall-d2-counterexample.pdf`
- Manuscript sha256: `e34de3716e64bba49880fff61808d5065477785596703bc0debf2379ad246c0c` (10 pages)

Check: `git show 1581d4cbdd69a0faf34b2089614fe2a814316a6a:preprints/A-Counterexample-to-Walls-D2-Problem-October-6-2026/wall-d2-counterexample.pdf | sha256sum`

The manuscript is not redistributed here. The presentation (10), character (14) and matrix (15) were transcribed from
pages 7-9 of the PDF (read as rendered images).

## Wall D(2): claims replayed

Notation follows the manuscript: P is the presentation (10) on x1, a1, x2, a2, s; rho is the character (14) into Q(w),
w^2+w+1=0; t1 = x1^2, t2 = x2^3.

| Claim | Method | Result | What is NOT checked |
|---|---|---|---|
| C1, C1b, C2, C3: rho sends the four relators to 1, kills t1 and t2, and rho(x1) = -1 | prefix-product evaluation in Q(w) (Fractions) | PASS | that x1^2 lies in K holds by definition only; nothing about G beyond the character |
| C4: s -> 1, others -> 0 gives a surjection onto Z | exponent sums | PASS | nothing beyond the stated map |
| C5, C5b, C5c: the twisted Fox matrix at rho equals matrix (15) entry by entry, and the closed-form entries quoted in the text | own free-group word engine and Fox calculus, exact comparison | PASS | the transcription of (15) is from the page image; a transcription error on our side would show as a FAIL, not a silent pass |
| C6, C6b: rank 4 over Q(w); the minor on columns (x1,a1,x2,s) equals -4(1-w^2) = -8-4w | Gaussian elimination in Q(w) | PASS | rank over other fields or rings |
| C7a, C7b, C7c: d1 d2 = 0; twisted Betti numbers h0 = h1 = h2 = 0; 14 cells and Euler characteristic 0 | rank computation; arithmetic | PASS | that the complex has the claimed homotopy type |
| C8a-C8d: the four word identities of (11) in P | explicit certificates as products of relator conjugates, compared as freely reduced words | PASS | the word problem of P in general |
| C8e: 4^2-1 = 15 = 3*5, 3^3-1 = 26 = 2*13, gcd(15,26) = 1 | arithmetic | PASS | the module argument of Lemma 4.1 that uses them |
| C9: the abelianised relator matrix of P is derived from the relator words (exponent sums); its invariant factors are exactly [1, 1, 1, 6] (rank 1, torsion Z/6) and adding t1 = x1^2 and t2 = x2^3 does not change them, so G^ab = P^ab, a necessary consequence of K perfect | Smith normal form (sympy) | PASS | perfectness of K (a necessary condition only) |
| C10: with monodromy rho(g) z^phi(g) the Fox matrix has rank 4 and the minor equals -4(1-w^2) at z=1 | sympy determinants modulo w^2+w+1 | PASS | Lemma 2.1 in general |
| C11: in every homomorphism P -> AGL(1,p), p in {2,3,5,7} (1568 of them), t1 and t2 map to 1 | exhaustive enumeration | PASS | perfectness of K (the search can refute, not establish it) |

Not checked by computation: Theorem 2.2 and Corollary 2.3 (a topological argument), Lemma 2.1 in general,
Proposition 3.1 (the finite plus construction), the splitting (8), the module argument of Lemma 4.1, and the background
citations. These are not computations and were not replayed.

## One-command reproduction

```
python -m venv .venv && . .venv/bin/activate && pip install -r replays/wall_d2/requirements.txt && python replays/wall_d2/replay_wall.py | diff -u replays/wall_d2/expected_output.txt - && echo REPRODUCED
```

On Windows use `.venv\Scripts\activate`. The run takes roughly 10 to 50 seconds on one core, depending on the machine.
Exit code 0 only if all 21 claims pass. Standard output is deterministic (timing goes to standard error), so it can be
compared byte for byte with `replays/wall_d2/expected_output.txt`; hashes are in `replays/wall_d2/SHA256SUMS`
(`expected_output.txt` sha256 `de0f02bd69524713b81221d2adfa7ea031d014fa0f2a517bdfd4c1c9e585a721`).

Tests (perturbing an input must make the replay fail; also checks the baseline against the expected output):

```
python -m unittest tests.test_adversarial -v
```

The tests perturb one entry of matrix (15), one relator exponent, and one value of the character, and assert that the
replay exits 1 with the corresponding FAIL line. CI (`.github/workflows/ci.yml`) runs the replay, the diff against the
expected output, the hash check and the tests on every push.

## Independence

`replay_wall.py` was written from the manuscript text only. It does not import or reuse code from the manuscript's
authors. Self-contained pieces: Q(w) over Python `Fraction`, a free-group word engine, a free-differential-calculus
routine. sympy is used only for the Smith normal form (C9) and the Laurent-polynomial cross-check (C10).

## Limitations

- C11 is a bounded search over affine groups of small prime fields; it can only refute, not establish, perfectness of K.
- C8 verifies identities of P through relator certificates; it does not decide the word problem of P.
- Matrix (15) was read from the page image; the computed matrix agrees with it exactly.
- Passing replays say the stated computations are consistent, nothing more.
- The "stated constants" checks in e101 (S23–S27) and e352 (CF04) compare transcribed constants with a second copy in the same file; they guard against edits, not against a transcription error.
- Dependencies are pinned (`sympy==1.14.0`, `mpmath==1.3.0`). The committed output was produced with CPython 3.14.5;
  CI also runs 3.12.

## Additional entries: pinned sources

All eight additional replays use the same pinned source commit as the Wall replay:
`1581d4cbdd69a0faf34b2089614fe2a814316a6a` of `https://github.com/openai/math` (2026-10-12). Check any file with
`git show 1581d4cbdd69a0faf34b2089614fe2a814316a6a:<path> | sha256sum`. The manuscripts are not redistributed here.

| Entry | Manuscript path (under `preprints/`) | sha256 | Pages |
|---|---|---|---|
| 059 | `Ambiently-homeomorphic-isolated-hypersurface-germs-in-C4-with-multiplicities-four-and-five-September-27-2026/paper.pdf` | `7b9c31dcfb1ba9b7c68b7ec31829963eebffcd6171a8139c26c20ee6c60e3b16` | 32 |
| 067 | `The-Campana-Peternell-conjecture-in-dimension-six-September-25-2026/main.pdf` | `c6c7151f9581ddd5b4d4bea3e6c1e96bf1fa776b97f6e0f2e2de1b952d5a28a1` | 27 |
| 101 | `A-sharp-entropy-bound-and-the-simplex-inequality-for-isotropic-constants-October-5-2026/isotropic-simplex.pdf` | `b6cb0115cd04d4eb0c24db7e600ba8dd9519848a92e87fa8a661f435e41e95a0` | 42 |
| 138 | `Subset-Sum-in-Time-2-power-0-49n-October-4-2026/subset-sum.pdf` | `69f9398e325a2628d32edc73109e2d008466170e35d8f7015636dc4f3f457566` | 33 |
| 138 (companion) | `A-Low-Space-Algorithm-for-Worst-Case-Subset-Sum-September-26-2026/paper.pdf` | `88f951cb233cf765b963cc3881f761846175bf5932649404d6ef3f69756d78ce` | 36 |
| 153 | `Arithmetic-classification-and-non-Pisot-singularity-for-Bernoulli-convolutions-October-3-2026/paper.pdf` | `d7c67b840e4422308a315ddebfa30f7eaddf8dde07affa3954afa84ac63034a2` | 36 |
| 268 | `The-periodic-spin-one-Haldane-gap-September-24-2026/paper.pdf` | `989f14096403897c693467f4007eb8997899648fd50029d5bee50ca9f8471ba6` | 30 |
| 318 | `A-rational-obstruction-to-strong-chromatic-splitting-at-height-three-September-25-2026/paper.pdf` | `6d0bd23f35dc42a1f645f99f7054f2e2f8d0033b78f3e3803b86774fc0558184` | 49 |
| 352 | `A-finite-time-singularity-of-Calabi-flow-on-projective-space-September-24-2026/paper.pdf` | `486e9e5ff47e92a4db90c7904eb5842085de78f02f825f5b68bcc5bd2719f25d` | 54 |

Each script was written from the manuscript text only; displayed matrices, tables and exponents were read from the rendered pages where a table or
formula was transcribed (for example Table 1 of entry 153, the matrix of (3.6) of entry 352). None of the authors' code is imported or reused.

## Additional entries: claims replayed and what is NOT checked

Standard output of every replay is deterministic (timing goes to standard error), so it can be compared byte for byte with `expected_output.txt`;
`SHA256SUMS` in each directory pins the script, the pinned requirements and the expected output. `--perturb` changes exactly one INPUT (a constant or a transcribed table entry) and then runs the unchanged main
checks, which must exit 1 with a FAIL line; `--selftest` runs the built-in wrong-input controls (each must be rejected). Floating-point values in the output are printed
to at most six significant digits.

### 059 Ambiently homeomorphic germs (`replays/e059_zariski_multiplicity/`, 17 claims)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| Lemma 3.1 and 9.1 arithmetic for 24 concrete (r, n0) with r = 4..11: the shift n1 = n0 + 5 * 2^r, congruences, gcd, n0 < n1 < 2 n0, alpha + beta = 4, a alpha = b beta, a, b = 2 mod 5, g >= 2^(r-1) > 5, existence of positive i, j with i = j mod 4 and (b i + a j)/m = q/5 for q = 1, 4 (the found pair is re-verified against these conditions). Every admissible n0 must satisfy every one of these properties; only the two strict lower bounds (a - 2) alpha/4 > 4 and (b - 2) beta/4 > 4, which the lemma states for n0 sufficiently large, are used to skip n0 below that threshold | exact integers and fractions, brute-force search for i, j | the lemmas for all large r (only instances) |
| d0 = h(pi_0) and d1 = h(pi_1) = psi (x - psi^2)(x - 2 psi^2) checked against h(X,T) = T (X - T^2)(X - 2 T^2); expansions phi(0,w), psi(0,w), q2(u,0), q2(v w^2, w); disc(v^2 - v + 1) = -3; delta(d_i) = 28, mu = 51 as arithmetic | polynomial arithmetic (python-flint) | the individual branch intersection multiplicities (their consequence mu = 51 is checked by another route) |
| mu(d0) = mu(d1) = 51 and mu(h) = 6 | exact rank of the truncated Jacobian algebra O/(J + m^N), N = 55 and 60 | |
| gradient estimate (3.1); rank/trace equations a + 4b = 6, a - b = 1; initial forms of degree 5 and 4 (Proposition 3.2) | symbolic | Milnor-lattice homology, Seifert-form congruences, strong approximation, Kato's classification: not finite computations |

### 067 Campana-Peternell in dimension six (`replays/e067_campana_peternell/`, 33 claims; several minutes)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| the ring grading (581 monomials of weight 9), the row families (4.11)-(4.16) (1785 nonzero rows), rank 516 and kernel dimension 65 of A over Q, AK = 0 | exact rational arithmetic (python-flint) | the geometric derivation of the rows, in particular that D_i vanishes for i >= 7 |
| selected index set (129 pairs), P of size 129 x 65 and rank 60, N with 69 rows, N P = 0, denominators prime to 10007 | exact | |
| |C2| = 1140, |C3| = 6820, dim K2 = 163, the 8280 x 6820 matrix of rank 6814 and dim K3 = 6 over F_10007 (and over F_1000003) | modular rank; a rank mod p bounds the rational rank from below when no denominator is divisible by p (checked) | a rational rank computation of the 8280 x 6820 matrix |
| the six witness columns lie in K3 (N kills all 120 double slices), the 6 x 6 minor -16/5, the functional (5.2) vanishes on them, (u - 4)(3u - 8) = 3u^2 - 20u + 32 | exact | |
| Section 3 moment identities and lambda = 9/4 | sympy | the geometric argument around them |

### 101 Simplex inequality for isotropic constants (`replays/e101_simplex_isotropic/`, 29 claims; about 10-15 min, `--quick` 27 claims in about a minute)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| the set (4.2) (0 < b < 1/2, l > 0, p > 0, |h0| <= .67 sqrt p, p >= 1/(7+a^2), p >= .25, p >= .46), the step bound table (2.5)-(2.6), the bounds (4.3) with the 2x2 PSD test (4.4) on the six intervals, (4.6) for a <= -1, the exact identities, the stated constants (alpha, beta, c = 3/2, eps, kappa, Q = [[.43, .13], [.13, 2.72]], the bounds .67, .25, .46 and the step table) and the tiling of [-6, 10] by the six intervals; (4.7) with max_b (b - beta)^2 = .28^2 and the value .11003 < .111 | ball arithmetic (python-flint Arb, 150 bits) with Taylor-model enclosures inside a bisection prover: a claim passes only if decided on every ball of the interval | **range**: a in [-6, 10] only; the tails a < -6 and a > 10 (the manuscript's analytic arguments) |
| (4.6) for a >= -1 | same, on **[-1, 6] only** | (6, 10]: the method does not finish there in reasonable time (node counts grow about tenfold per unit step of 2) |

`--quick` skips the two longest claims (the T1 row of the (4.3)/(4.4) table and (4.6) for a >= -1).

### 138 Subset Sum in time O(2^{0.49n}) (`replays/e138_subset_sum/`, 26 claims)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| parameters Delta, Gamma; |r| <= 2.433, |s| <= 4.010; centres of Appendix A.2; Table 1 (nodes 0..32) within the printed errors; node symmetry; (A.3) 64-node means < 1.4295 < 2^0.5168 | ball arithmetic (python-flint Arb) | |
| Proposition 3.1: both integrals < 2^0.5168, computed directly from definitions (3.1)-(3.2); D, D-A, D-E >= 0; the identities (A.2) D = B - L, u = L + (kappa-1) B - H, v = (2-kappa) B - 2 L + H and L = log Q, each against the continuous quantities at the 64 nodes to better than 1e-9 | float64 trapezoid sums on N = 256..2048 nodes (geometric convergence) | **a rigorous enclosure of the quadrature remainder** ((A.6) is an analytic estimate that is not re-derived); the finite-P statement of Lemma 3.2 |
| companion paper Section 2.2: 3^12/2^19 > 1 + 1/74, 3^53 > 2^84, log 3 > 84/53, h* < .8114, h(1/3) < 14/15, h(1/6) > 5/8, the (39) arithmetic | exact rationals | the algorithms, their correctness and running-time analyses |

### 153 Non-Pisot singular Bernoulli convolution (`replays/e153_bernoulli_non_pisot/`, 48 claims)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| identity (6.2); the 31 root disks of Appendix A.1 (16 transcribed centres plus conjugates): |z| < 2, ||p'(z)||_inf > 1, 2e35 ||p(z)||_1 < ||p'(z)||_inf, separation, Rouche remainder | exact Gaussian-integer Horner evaluation | |
| Proposition 6.2 numbers (beta, |alpha|, moduli of the other roots, simplicity, not-Pisot arithmetic, (6.4)); Corollary 5.4 quartics | independent 80-digit roots (mpmath); sympy | |
| C0 as Laurent coefficients of 1/p, (A.10) maxima, the weights (A.7), |R| = 753, weight sums 198442946/W and 11816095/W (exact integer equality), the eight small weights (exact equality), (A.22) values compared with the printed .010862, .014732, .003156 (ours must lie below them and within 1%), the stated t = .06, (A.14)-(A.16), (A.19), (A.20), (A.22), the final bound .0000056968 < .000007 | float64 cosine products (m in [-8000, 150000]) | **the 36-digit integers of (A.17) and (A.23) bit for bit** (checked at value level: all quantities lie inside the printed brackets; the weight sums match exactly); the error analysis (A.8)-(A.13) of the manuscript's own rounding scheme |

### 268 Periodic spin-one Haldane gap (`replays/e268_haldane_gap/`, 68 claims)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| 32 exact scalar certificates of Sections 4-5: (4.9), the (4.18) table of parameters and the (4.10) hypotheses for both rows, the gap bounds log(1/r)/b0 against the printed closed forms, mI, mO, mN margins, filters, (5.5), R(t) (exact), the two printed comparisons .02093518 and .002777621984 (exact equality at the printed precision), B_J, B_N, FJ/FN margins, m_J, m_N, p0, q0, and the six rounded updates reproducing the printed table of p_j, q_j | exact rationals (fractions, sympy); the Lemma 5.1 tables are inputs here | |
| Appendix C: matrices from Tables 1-2, block dimensions (C.6), Lemma C.1 on a 4-site ring, and V_N > 0, -c_N V_N <= 10^12 U_N < (-c_N + 1) V_N, U_N/V_N < -a at N = 72 and 120 | exact big integers (python-flint) | |
| Lemma 5.1 thermal tables: b = 21/2, n = 6, 8, 10, g in {I, P, C}; b = 49/4, n = 4..12, g in {I, P}; reduced spectra equal literal diagonalisation of (2.5) for n <= 6 | **float64 exact diagonalisation** by S^z and momentum blocks (dense for n <= 9, lowest eigenvalues by Lanczos with the cut E + a n < 2.5 for n >= 10; the claim that the discarded eigenvalues change Z by less than 1e-6 is checked), compared with the printed centres within the printed errors | **a rigorous enclosure** (the manuscript's packed-integer matrix-exponential algorithm of Appendices A and B is not reimplemented); Propositions 3.1, 4.3, 5.4, 5.5 and the purity recursion |

The 68 claims of `replay_haldane.py` are 32 scalar certificates, 7 Appendix C claims, and 29 in the thermal part (the reduction check, the discarded-eigenvalue bound and the 27 table entries).

### 318 Chromatic splitting at height three (`replays/e318_chromatic_splitting/`, 15 claims)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| sl2: d1 equal to the printed matrix entry by entry (signs and order) with determinant 4, d2 = 0; sl3: cochain dimensions 1, 8, 28, 56, 70; d_{k+1} d_k = 0; the minors det d1[R1,C1] = -1, det d2[R2,C2] = 2916, det d3[R3,C3] = 1024; ranks 8, 20, 35; H^1 = H^2 = 0, H^3 = Q; the trace cocycle is closed, is not a coboundary, omega(E12, E21, E11 - E22) = 2 | exact integer matrices (python-flint) | the homotopy-theoretic content; the entry's other four manuscripts were searched for finite computations by keyword and context and none was identified (not a line-by-line read) |

### 352 Finite-time singularity of Calabi flow (`replays/e352_calabi_flow/`, 28 claims)

| Claims replayed | Method | What is NOT checked |
|---|---|---|
| (3.5) from Q(D); the transcription of (3.6) (k0, d0, the matrix, R, S0); (3.7) with the bounds computed from the matrix; Lemma 3.3 exact evaluations, root bounds, the projection matrix Gamma, M_0..M_2, E = 46876434629/15680000000; Lemma 3.2 constants; Lemma 3.4 residual (sampled) and (3.19); the origin series solves (3.5) (ODE residual at s = 1/2 with numerically differentiated series, independent of the coefficient recurrence) | exact fractions, sympy, mpmath (40 digits) | the rigorous supremum in Lemma 3.4 |
| Proposition A.1 / (3.18): |q| <= .25 R, ||(Re xi, Im xi) - z*||_inf <= .13 R and 1 + p > 0 on [1, 56] | high-precision Taylor ODE integration at **8 sample points** of the parameter square Z (centre, four corners, three interior points) | 1 + p > 0 is checked on the integer values s = 1..56 for the shot and on a 50-point grid of (0, 1] for the origin series, not on a continuum; **the enclosure for every z* in Z** (the manuscript's integer-interval certificate over 288 blocks is not reimplemented); the half-line contraction and the matching argument |

## Additional entries: reproduction and tests

Python 3.12 or 3.14; each directory has its own pinned `requirements.txt` (`python-flint==0.9.0`, `sympy==1.14.0`, `mpmath==1.3.0`,
`numpy==2.5.4`, `scipy==1.18.1`):

```
python -m venv .venv && . .venv/bin/activate && pip install -r replays/e059_zariski_multiplicity/requirements.txt
cd replays/e059_zariski_multiplicity && python replay_zariski.py | diff -u expected_output.txt - && sha256sum -c SHA256SUMS
```

Scripts: `replay_zariski.py` (059), `replay_campana_peternell.py` (067), `replay_simplex.py` (101; `--quick` against `expected_output_quick.txt`),
`replay_subset_sum.py` (138), `replay_bernoulli.py` (153), `replay_haldane.py` (268), `replay_chromatic.py` (318), `replay_calabi.py` (352).
Approximate run times on one core: 059 0.5 min, 138 and 318 a few seconds, 153 1 min, 352 2 min, 268 3-4 min, 067 3-6 min, 101 quick 1 min and full about 10-15 min.

Tests (one perturbation per entry: an input is changed, the unchanged main checks run, and the replay must exit 1 with the corresponding FAIL line; the tests run the full main path, 101 in its `--quick` subset, and take several minutes):

```
python -m unittest tests.test_adversarial -v              # Wall D(2) (needs replays/wall_d2/requirements.txt)
python -m unittest tests.test_replays_adversarial -v      # the eight additional entries (needs their pinned requirements)
```

CI (`.github/workflows/ci.yml`) runs on every push, on Python 3.12 and 3.14: the Wall replay and tests; the seven fast replays (the simplex entry in its `--quick`
subset) with the diff against the expected output, the hash check and the built-in wrong-input controls (`--selftest`); a separate slow job with the Campana-Peternell replay and the full simplex run (same steps, including `--selftest`); and the
perturbation tests for the eight additional entries.

## Checker hash and one-command check

Every replay script prints `checker-sha256: <sha256 of its own source file(s)>` as the opening line of its output, and
`CHECKERS.sha256` lists the same hashes (several files of one replay are concatenated in the order listed; CRLF is read as LF).
`python check_all.py` re-runs every replay, compares stdout with the committed expected output byte for byte, verifies each
`SHA256SUMS` and `CHECKERS.sha256`, and exits nonzero on any difference (`--quick` uses the documented quick subset where one exists,
for the simplex entry; a full run takes about 10 minutes with the replays in parallel). It runs as one job in CI.

## Self-audit

Question: does each checker reject a wrong input? Rules were written down before any mutant: [`selfaudit/PREREG.md`](selfaudit/PREREG.md).
Each mutant changes one stated value or one data element (or one edge of a graph), or moves a stated value by one unit past an
inequality, in a temporary copy of the replay; it counts as killed if the replay exits nonzero or its output differs from the expected output.

Result ([`selfaudit/RESULT.md`](selfaudit/RESULT.md)): 77 mutants over the 10 replays, 74 killed (73 by a nonzero exit, 1 by an output
difference only), 3 survivors: 2 equivalent (e067-02, e318-03) and 1 mis-specified mutant (e059-05, whose mutated hypothesis is also true), argued in writing in RESULT.md. The initial run had 5 survivors; 2 were findings about our
checkers and are fixed (commit ee3a574): the claim that the Wall group surjects onto Z did not assert surjectivity (`wall_d2`), and the stated
constant 4.44605e-13 in the Calabi-flow replay was accepted with a tolerance of 1e-17, so a value off by one unit in its last printed digit
passed (`e352`). Both mutants are killed after the fixes. Amendment history: [`selfaudit/PREREG_AMENDMENTS.md`](selfaudit/PREREG_AMENDMENTS.md).

| replay | killed / total | survivors |
|---|---|---|
| wall_d2 | 7 / 7 | none |
| e059 | 8 / 9 | 1 (equivalent) |
| e067 | 7 / 8 | 1 (equivalent) |
| e101 | 8 / 8 | none |
| e138 | 8 / 8 | none |
| e153 | 6 / 6 | none |
| e268 | 7 / 7 | none |
| e318 | 7 / 8 | 1 (equivalent) |
| e352 | 9 / 9 | none |
| ramsey_alphaevolve | 7 / 7 | none |

Re-run (about 30 minutes on 4 cores; the e101 mutants use the `--quick` subset, ramsey_alphaevolve uses `--cheap`):

```
python selfaudit/run_selfaudit.py -j 4
```

## License

Code: MIT. Text: CC-BY-4.0. See [LICENSE](LICENSE).
