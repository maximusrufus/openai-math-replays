# Amendments and history of the pre-registration

Times are local (UTC-7), 2026-10-23, taken from the local reflog. Rules in `PREREG.md` are unchanged by every entry below
except where stated; no threshold or pass rule was changed after any result existed.

## History of PREREG.md

| time | event |
|---|---|
| 20:02:21 | `PREREG.md` committed as 3279db5. No mutant, `mutants.json` or runner existed yet. |
| 20:02:34 | The same commit was amended (now 71be1da) to remove wording only; this is the only rewrite of the file. No mutant, `mutants.json` or runner existed yet. |
| after | `PREREG.md` has not been modified since 20:02:34 (`git log -- selfaudit/PREREG.md` lists the one commit). |

Exact change between 3279db5 and 71be1da (all wording, no rule or number):
- "exactly the first value outside the printed error" became "exactly the smallest value outside the printed error" (e268 row);
- "ignores only the first line" became "ignores only the opening line"; "excluding the first line" became "excluding the opening line" (two places);
- "that is recorded in `RESULT.md`" became "that is noted in `RESULT.md`";
- "added after the first run" became "added after the initial run".

Reason: editorial wording only. An earlier note said the PREREG commit was rewritten twice; it was amended once. A later history rewrite (message wording only, 20:55:29) covered the commits
after 71be1da and did not touch 71be1da.

## Other history rewrites (not PREREG.md)

| time | event |
|---|---|
| 20:24:43 | The commit adding `mutants.json` and the runner was amended (b80ba43 to 30765a3, now a791cee). Changes: one mutant `desc` wording ("in the first graph" became "in the graph with coordinates (u, w)") and one runner message wording. This happened after the full mutant run had started (the run had already loaded the list); no mutation, locator or class was changed. |
| 20:55:29 | Commit messages of the four commits after 71be1da were reworded (one phrase), content unchanged. |

## Amendments

- 2026-10-23 after 20:59 (commit following ee3a574): mutants e352-07 and e352-08 keep the same mutation (stated value 4.44605e-13 changed to 4.44606e-13 and
  4.44615e-13) but their `find` text now targets the fixed comparison `abs(kmax - Fr("4.44605e-13"))`, because the checker line they originally
  matched was replaced by ee3a574. Reason: the fix for the finding changed the source line; the mutation, class and expected effect are unchanged.
  The mutant list is otherwise as committed in a791cee; wall_d2-07 is unchanged and applies to the fixed source as before.

The ordering above rests on local git history; PREREG.md was not externally timestamped.
