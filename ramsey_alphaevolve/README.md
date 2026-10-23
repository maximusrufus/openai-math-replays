# Re-check of the released Ramsey lower-bound graphs of arXiv 2603.09172 (v2)

The paper (arXiv 2603.09172, v2) links a set of adjacency-matrix files that are meant to certify lower bounds
R(s,t) >= m. A file for R(s,t) >= m holds a graph on m-1 vertices; it certifies the bound if the graph has no clique
of size s and no independent set of size t. This folder re-checks exactly that property for each released file and
nothing else. It does not re-check how the graphs were found, and it does not re-check any other claim of the paper.

Credit: the graphs, and the results they encode, are Google's. They are copied unmodified into `graphs/` (file names
only changed to remove spaces and comparison signs) under the Apache License 2.0; see `NOTICE` and
`LICENSE-APACHE-2.0`. This folder adds a checker, its output and this table.

## Pinned source

- Repository: `https://github.com/google-research/google-research`
- Commit: `0d94f8c4d92678ce5f7a0eacae725f79174cfb03` (2026-10-22), directory `ramsey_number_bounds/improved_bounds`
  (18 adjacency-matrix `.txt` files, a README and a notebook; only the 18 `.txt` files are copied).
- `SHA256SUMS` lists the sha256 of every copied file; the files are byte-identical to the pinned commit.
- The arXiv listing for v2 has no ancillary files.

## One-command reproduction

```
python -m pip install -r requirements.txt
python check_ramsey.py --selftest            # planted K_s and I_t must be detected; exit 0
sha256sum -c SHA256SUMS                      # the graph files are the ones tabulated below
python check_ramsey.py --cheap > actual_output.txt && diff expected_output.txt actual_output.txt   # nine graphs, under a minute on a 24-thread machine
```

Any single file: `python check_ramsey.py graphs/R_4_15_ge_163.txt`. Exit status is nonzero for a malformed file, a found
K_s or independent t-set (the witness is printed), or a failed self-test. The four largest verified graphs
(R(4,18) to R(4,21)) take about 3 minutes in total on 24 threads (A and B both) and are not in `--cheap`.

## What is checked

- Format: square, symmetric, entries 0/1, zero diagonal, and exactly m-1 vertices (from the file name).
- Method A: bitset depth-first clique search (uint64 words, popcount pruning, numba, parallel over the smallest vertex).
- Method B: a different algorithm and data structure: Python big-integer bitsets, branch and bound with a greedy-colouring bound
  (Tomita style), multiprocessing over root vertices.
- An independent set of size t is a clique of size t in the complement; both methods are run on the graph (K_s) and its complement (I_t).
- Self-test (`--selftest`): on the R(3,12) graph, a planted K_3, K_4 and a planted independent 12-set must be found by both
  methods with a valid witness, a file with a planted K_3 must be rejected, and malformed matrices must be rejected.

## Results

Verdicts below were produced by the original run of these two methods (and, for the second table, by the solvers named there).
Runtimes are wall-clock seconds from that run on a 24-core/32-thread machine (A parallel, B 12 processes), excluding the self-test.
`check_ramsey.py` reproduces the A+B verdicts for the first 13 rows (the nine `--cheap` rows were re-run for `expected_output.txt`,
R(4,18) to R(4,21) were re-run once more with the same verdicts).

### Fully checked by two methods (A and B)

| bound | n | sha256 | method(s) | verdict | runtime s (A+B) | reproduce here |
|---|---|---|---|---|---|---|
| R(3,12)>=53 | 52 | `f255d720051434018a9095438816c06f2151196e7b9dbf82d55528d691bd84bf` | A, B | VERIFIED | 4.1 | `python check_ramsey.py graphs/R_3_12_ge_53.txt` |
| R(3,13)>=61 | 60 | `d619f0bfed32af5e2fbc8c724d9c65dbf388a19d5647d20e7835b501bc03227a` | A, B | VERIFIED | 5.6 | `python check_ramsey.py graphs/R_3_13_ge_61.txt` |
| R(3,15)>=74 | 73 | `ae94e8e6b3c1a26b90cf55d754c742074efd62449c8227a285bd288eb093d26e` | A, B | VERIFIED | 5.6 | `python check_ramsey.py graphs/R_3_15_ge_74.txt` |
| R(3,17)>=93 | 92 | `71810667cd353e766375041e7968596e6ba7f67cd7d722158b09ab9a30e2336b` | A, B | VERIFIED | 12.1 | `python check_ramsey.py graphs/R_3_17_ge_93.txt` |
| R(3,18)>=100 | 99 | `3e9d16e29111f3a2f25ae3b992235804c2612de2f6ada5570901d17ceedfca45` | A, B | VERIFIED | 10.5 | `python check_ramsey.py graphs/R_3_18_ge_100.txt` |
| R(4,13)>=139 | 138 | `c50a55389d2c69b4cc7c1304639b2af67f16fd1a4a866870b87f907d41b69c54` | A, B | VERIFIED | 6.4 | `python check_ramsey.py graphs/R_4_13_ge_139.txt` |
| R(4,14)>=148 | 147 | `00563a5e5e82289ffd1bf62ed2422a8c2865c11da9eaab1c101c4569a5ff2c8a` | A, B | VERIFIED | 7.5 | `python check_ramsey.py graphs/R_4_14_ge_148.txt` |
| R(4,15)>=163 | 162 | `56f3cf76b87efe7d845b2bd59505d44638a9b49e04dfd33dcb4e86426cb4ec23` | A, B | VERIFIED | 10.9 | `python check_ramsey.py graphs/R_4_15_ge_163.txt` |
| R(4,16)>=174 | 173 | `def3d53ed5879f8984434b93861b3a6a66049f2491fa04da212d704504b82f6f` | A, B | VERIFIED | 17.1 | `python check_ramsey.py graphs/R_4_16_ge_174.txt` |
| R(4,18)>=209 | 208 | `5317e46e94f27bca829f0aa7b7f1bb2dcdf56e040ea9b98c199b00b7685f9945` | A, B | VERIFIED | 31.8 | `python check_ramsey.py graphs/R_4_18_ge_209.txt` |
| R(4,19)>=219 | 218 | `8d2378420e281725f4fa6d510b2be0ec458dcb2a155ef7a61ec5c4b85bfd0f2c` | A, B | VERIFIED | 69.9 | `python check_ramsey.py graphs/R_4_19_ge_219.txt` |
| R(4,20)>=237 | 236 | `d33838caf70a91612889c26d8ac63a9bd96513a9d5e95803042361ed7a999e3e` | A, B | VERIFIED | 55.7 | `python check_ramsey.py graphs/R_4_20_ge_237.txt` |
| R(4,21)>=244 | 243 | `3cdf3083ccf1dde21e41b512d93443cfe55e28e9071976336fff3f9d27d7599e` | A, B | VERIFIED | 108.4 | `python check_ramsey.py graphs/R_4_21_ge_244.txt` |

### Five large graphs

For all five, the K_s side was found free by two methods (the plain bitset search A and a colouring-bound search). The
rows below concern the independent-set side. TOO-LARGE means: not finished within the time spent (progress shown); such a
row is NOT verified. Solver code for this table (a SAT cube-and-conquer run and a dedicated max-clique solver on the complement)
is not part of this folder, and `check_ramsey.py` is not expected to finish on these graphs.

| bound | n | sha256 | method(s): independent-set side (SAT / C2 / C2-deg) | verdict | runtime | reproduce here |
|---|---|---|---|---|---|---|
| R(7,13)>=513 | 512 | `e8160ebe4e877ac22bb2689e2fbb0eca8f517c65d5516e11745d4ef93829d70c` | SAT: TOO-LARGE: 127/512 in 1805 s, hardest cube done 1117 s; C2: stopped at 4100 s: 119/499 roots, 5.7e10 nodes, recent rate about 1 root per 2 min; C2-deg: not run | TOO-LARGE | see method column | not re-run by `check_ramsey.py` (solvers not included) |
| R(8,12)>=511 | 510 | `5db94f26576014f66792c665096e7dd84dd0835a3c048d0b82af5bf0d735a2df` | SAT: TOO-LARGE: 125/510 in 918 s, hardest cube 522 s; C2: NONE (no independent 12-set), 498/498 roots, 596 s, 1.38e10 nodes; C2-deg: NONE, 488/488 roots, 841 s, 2.32e10 nodes | VERIFIED (two clique-solver runs agree; SAT leg TOO-LARGE) | see method column | not re-run by `check_ramsey.py` (solvers not included) |
| R(9,11)>=617 | 616 | `ce12b6cb8becabd658780db0ab8d3248b89f25e3841ed749c55c80902be5afa1` | SAT: TOO-LARGE: 272/616 in 918 s, hardest cube 172 s; C2: NONE (no independent 11-set), 605/605 roots, 582 s, 1.01e10 nodes; C2-deg: NONE, 581/581 roots, 415 s, 1.08e10 nodes | VERIFIED (two clique-solver runs agree; SAT leg TOO-LARGE) | see method column | not re-run by `check_ramsey.py` (solvers not included) |
| R(10,12)>=856 | 855 | `4d270464208e92dcda82232a2618879810c3e11b29af9d001a7484df665f79e2` | SAT: TOO-LARGE: 241/855 in 902 s, hardest cube 224 s; C2: stopped at 1930 s: 51/843 roots, 2.18e10 nodes; C2-deg: not run | TOO-LARGE | see method column | not re-run by `check_ramsey.py` (solvers not included) |
| R(8,18)>=1052 | 1051 | `c5a4706f05436f1d74b56fa401b7ed0e0348422c92f14c03aa522752a72bf1ad` | SAT: TOO-LARGE: 85/1051 in 909 s, hardest cube 633 s; C2: stopped at about 1900 s: 0/1033 roots completed; C2-deg: not run | TOO-LARGE | see method column | not re-run by `check_ramsey.py` (solvers not included) |

Stated plainly:
- R(7,13), R(10,12) and R(8,18) are not verified: the independent-set side did not finish.
- For R(8,12) and R(9,11), the verdict VERIFIED rests on two exhaustive runs of the same colouring-bound solver under different
  vertex orders (different search trees, same bounding family) plus the K_s side by two methods; the SAT leg did not finish, so these two
  are verified by one algorithm family on the independent-set side.
- Time caps were shorter than the 2 hours requested; whether more time would finish the TOO-LARGE rows is untested.
- In both tables, a self-test with a planted independent set of size t (and, for the first table, a planted K_s) was detected on every graph.

## Note on R(4,15)

The paper's summary text states R(4,15) >= 159; the released file is named `R(4,15) >= 163` and holds 162 vertices. Both facts are
recorded here as found. The file as released was checked (162 vertices, no K_4, no independent 15-set, by A and B). This folder does
not decide which statement is the intended one.

## Files

- `check_ramsey.py`, `requirements.txt` (numpy 2.5.0, numba 0.68.0)
- `graphs/` the 18 released files, `SHA256SUMS`, `NOTICE`, `LICENSE-APACHE-2.0`
- `expected_output.txt` output of `python check_ramsey.py --cheap` (stdout; timings go to stderr)

A passing re-check shows that the released files have the stated graph property (as computed by these two programs); it is not a
review of the paper.
