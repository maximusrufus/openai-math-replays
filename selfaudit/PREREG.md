# Pre-registration: mutation self-audit of the replay checkers

Status: written and committed before any mutant is written or run. The mutant list (`mutants.json`) is committed
separately, after this file. Thresholds and rules below are not changed after results are seen; any change is appended
as a dated amendment at the end of this file.

## Question

Do the replay checkers reject a wrong input? Each replay encodes the data of a manuscript (constants, table entries,
witness entries, polynomial coefficients, graph edges) and asserts claims about it. If one stated value or one data element is
changed by an amount large enough to break the claim, an adequate checker must notice. A mutant that is not noticed
is a finding about the checker (it asserts less than its output suggests).

## Replays covered

wall_d2, e059, e067, e101, e138, e153, e268, e318, e352, ramsey_alphaevolve.

## Mutant classes (per replay)

Class D (input-data mutant): change one stated value or one data element to a different value, with a size chosen so that the
claim that depends on it is false for the mutated value (checked by hand arithmetic in the mutant's `why` field).
Class B (boundary mutant): where the claim is a float or exact inequality, move one stated value by one unit in the last
place (floating point) or by one unit in the last printed digit / one integer (exact) past the constraint, so that the
inequality fails by the smallest possible amount. Boundary mutants are written only where such an inequality exists.

Every replay gets at least 5 mutants, and at least 2 of them boundary mutants wherever the replay has an inequality claim
on stated values. Planned data per replay:

| replay | class D data elements | class B (boundary) candidates |
|---|---|---|
| wall_d2 | one exponent of a relator of presentation (10); one character value; one entry of matrix (15); one generator index in a certificate word | none expected (claims are exact identities, no inequalities) |
| e059 | one residue or exponent in the Lemma 3.1 instances; one exponent of a plane-curve polynomial; one stated Milnor number | the strict inequality of the gradient estimate, if stated numerically |
| e067 | one witness entry; one weight; one index of the ideal/base split | none expected (exact rank / identity claims) |
| e101 | one stated constant; one interval endpoint; one table row entry; one step bound | a constant moved one unit in its last stated digit across the bound it must satisfy |
| e138 | one entry of R6 or S6; one table entry of Table 1; the constant KAPPA | KAPPA and one table entry moved one unit past the stated tolerance or bound |
| e153 | one centre of the root-disk table; one polynomial coefficient; one Salem quartic coefficient | a table centre moved to just outside the stated disk radius |
| e268 | one MPS parameter; one entry of each printed thermal table | a table entry moved to exactly the smallest value outside the printed error (10e-6 at b = 21/2, 30e-6 at b = 49/4) |
| e318 | one row or column index of the minor index sets R1..R3, C1..C3; one structure constant | none expected |
| e352 | one of the stated constants (B, C, R, S0, K0, D0); one entry of the printed matrix | K0 or D0 moved by one unit in the last printed digit across a stated bound |
| ramsey_alphaevolve | one edge added (creating a clique of the forbidden size) or removed (creating an independent set of the forbidden size) in a graph file | none (the claim is a combinatorial property, no float inequality) |

"None expected" means no boundary mutants are required; if the audit author finds an inequality while writing the
mutants, boundary mutants are added for it and listed in `mutants.json`.

## Execution rule

A mutant is a list of exact string replacements (file, find, replace; `find` must occur exactly `count` times, default 1) applied to a
temporary copy of the replay directory. The working tree is never modified. The mutated replay is run with its normal
arguments (replay e101 is run with `--quick`; every other replay with no argument, ramsey_alphaevolve with `--cheap`).
A mutant whose `find` string does not match is a runner error, reported separately, and is repaired before the real run;
it is never counted as killed or as a survivor. A run that times out is reported as an error and re-run with a longer limit.

## Pass rule

Output comparison ignores only the opening line (`checker-sha256: ...`), because any source mutation changes it by
construction and counting that as a kill would make every source mutant killed for free.

- Killed (hard): the mutated run exits with a nonzero status.
- Killed (soft): the mutated run exits 0 but its output, excluding the opening line, differs from the expected output.
  This still counts as killed, but is reported separately, because it means the checker's own claims all passed on
  the mutated data and only the textual diff caught it.
- Survivor: exit status 0 and output (excluding the opening line) identical to the expected output.

Pass rule of the audit: every non-equivalent mutant is killed. An audit that has a non-equivalent survivor does not pass,
and the survivor is reported as a finding about our checker.

## Equivalent mutants

A mutant may be declared equivalent (it does not change any stated claim) only by a written argument in `RESULT.md`
stating why the mutated value cannot affect any checked claim. A survivor without such an argument is non-equivalent.
Equivalent mutants are listed with their argument and excluded from the pass rule; they are still counted in the
totals.

## Non-goals and constraints

- Checkers are not fixed in this task; survivors are reported as findings.
- Mutants are not edited or weakened after results are seen. If a mutant is found to be mis-specified (it does not
  actually break the claim), that is noted in `RESULT.md` and the mutant is kept in the totals.
- Additional mutants added after the initial run are labelled `post-hoc` and reported separately.
