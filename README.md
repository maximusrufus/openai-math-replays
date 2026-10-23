# Independent replays of computational claims in the openai/math October 2026 release

Purpose: independently re-run the finite, checkable computations that are stated in the text of manuscripts in
[openai/math](https://github.com/openai/math), using exact arithmetic and code written from the paper text only.
This is not a review of the full proofs. A replay shows that the stated computations are consistent with the
manuscript; it does not show that a theorem is true.

Credit: the manuscripts, results and all displayed formulas belong to their original authors (OpenAI, as listed in each
manuscript). This repository contains only independent re-computation code and its output.

## Replays

| Entry | Replay | Result |
|---|---|---|
| Wall's D(2) problem counterexample (2026-10-06) | [`replays/wall_d2/`](replays/wall_d2/) | 21 of 21 replayed claims pass |
| MUB N(6)=3 | separate, already public: [maximusrufus/mub6-certified-rerun](https://github.com/maximusrufus/mub6-certified-rerun) (not copied here) | see that repository |

## Pinned source

- Repository: `https://github.com/openai/math`
- Commit: `1581d4cbdd69a0faf34b2089614fe2a814316a6a` (2026-10-12, "Merge pull request #2 from openai/codex/lock-prs")
- Manuscript path: `preprints/A-Counterexample-to-Walls-D2-Problem-October-6-2026/wall-d2-counterexample.pdf`
- Manuscript sha256: `e34de3716e64bba49880fff61808d5065477785596703bc0debf2379ad246c0c` (10 pages)

Check: `git show 1581d4cbdd69a0faf34b2089614fe2a814316a6a:preprints/A-Counterexample-to-Walls-D2-Problem-October-6-2026/wall-d2-counterexample.pdf | sha256sum`

The manuscript is not redistributed here. The presentation (10), character (14) and matrix (15) were transcribed from
pages 7-9 of the PDF (read as rendered images; plain text extraction drops the omega glyphs).

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
| C9: G^ab = P^ab (rank 1, torsion Z/6), a necessary consequence of K perfect | Smith normal form (sympy) | PASS | perfectness of K (a necessary condition only) |
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
(`expected_output.txt` sha256 `23d6cef85f5cef533f585d4f4048001988e49a203287711e35895b12135c852c`).

Tests (perturbing an input must make the replay fail; also checks the baseline against the expected output):

```
python -m unittest discover -s tests -v
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
- Dependencies are pinned (`sympy==1.14.0`, `mpmath==1.3.0`). The recorded output was produced with CPython 3.14.5;
  CI also runs 3.12.

## License

Code: MIT. Text: CC-BY-4.0. See [LICENSE](LICENSE).
