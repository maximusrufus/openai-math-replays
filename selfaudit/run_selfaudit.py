#!/usr/bin/env python3
"""Mutation self-audit of the replay checkers (rules: selfaudit/PREREG.md).

  python selfaudit/run_selfaudit.py              # run every mutant, write selfaudit/RESULT.md and selfaudit/results.json
  python selfaudit/run_selfaudit.py --only e153  # restrict to replays whose name contains the text
  python selfaudit/run_selfaudit.py -j 6         # mutants run in parallel (default 8)
  python selfaudit/run_selfaudit.py --strict     # exit 1 if any non-equivalent mutant survives

Each mutant is applied to a temporary copy of its replay directory; the working tree is never modified.
The mutated replay is run with its normal arguments (e101: --quick, ramsey_alphaevolve: --cheap) and its stdout is compared
with the committed expected output, ignoring only the opening `checker-sha256:` line.

  killed (hard)  the run exits nonzero
  killed (soft)  the run exits 0 but the output differs from the expected output
  survived       exit 0 and identical output
  error          the mutant could not be applied or the run timed out (not counted; fix the mutant and re-run)
"""

import argparse
import concurrent.futures
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# name -> (directory, script, extra arguments, expected output file)
REPLAYS = {
    "wall_d2": ("replays/wall_d2", "replay_wall.py", [], "expected_output.txt"),
    "e059": (
        "replays/e059_zariski_multiplicity",
        "replay_zariski.py",
        [],
        "expected_output.txt",
    ),
    "e067": (
        "replays/e067_campana_peternell",
        "replay_campana_peternell.py",
        [],
        "expected_output.txt",
    ),
    "e101": (
        "replays/e101_simplex_isotropic",
        "replay_simplex.py",
        ["--quick"],
        "expected_output_quick.txt",
    ),
    "e138": (
        "replays/e138_subset_sum",
        "replay_subset_sum.py",
        [],
        "expected_output.txt",
    ),
    "e153": (
        "replays/e153_bernoulli_non_pisot",
        "replay_bernoulli.py",
        [],
        "expected_output.txt",
    ),
    "e268": (
        "replays/e268_haldane_gap",
        "replay_haldane.py",
        [],
        "expected_output.txt",
    ),
    "e318": (
        "replays/e318_chromatic_splitting",
        "replay_chromatic.py",
        [],
        "expected_output.txt",
    ),
    "e352": ("replays/e352_calabi_flow", "replay_calabi.py", [], "expected_output.txt"),
    "ramsey_alphaevolve": (
        "ramsey_alphaevolve",
        "check_ramsey.py",
        ["--cheap"],
        "expected_output.txt",
    ),
}


def read_text(path):
    with open(path, "rb") as f:
        return f.read().decode("utf-8")


def write_text(path, text):
    with open(path, "wb") as f:
        f.write(text.encode("utf-8"))


def apply_edit(work, edit):
    path = os.path.join(work, edit["file"])
    text = read_text(path)
    want = edit.get("count", 1)
    got = text.count(edit["find"])
    if got != want:
        raise ValueError(
            "%s: find string occurs %d times, expected %d: %r"
            % (edit["file"], got, want, edit["find"][:80])
        )
    write_text(path, text.replace(edit["find"], edit["replace"]))


def apply_graph_edit(work, ge):
    """set symmetric entries (i, j) of a 0/1 adjacency matrix file; the layout of the file is rewritten in the same bracket style."""
    path = os.path.join(work, ge["file"])
    tok = read_text(path).replace("[", " ").replace("]", " ").split()
    n = int(round(len(tok) ** 0.5))
    assert n * n == len(tok)
    for i, j, v in ge["set"]:
        if tok[i * n + j] == str(v):
            raise ValueError("entry (%d,%d) already equals %d" % (i, j, v))
        tok[i * n + j] = tok[j * n + i] = str(v)
    rows = ["[" + " ".join(tok[r * n : (r + 1) * n]) + "]" for r in range(n)]
    write_text(path, "[" + "\n ".join(rows) + "]\n")


def run_mutant(m, timeout):
    d, script, args, expected = REPLAYS[m["replay"]]
    src = os.path.join(ROOT, d)
    tmp = tempfile.mkdtemp(prefix="mut_")
    t0 = time.time()
    try:
        work = os.path.join(tmp, "w")
        shutil.copytree(
            src, work, ignore=shutil.ignore_patterns("__pycache__", "actual_output.txt")
        )
        for e in m.get("edits", []):
            apply_edit(work, e)
        if "graph_edit" in m:
            apply_graph_edit(work, m["graph_edit"])
        try:
            p = subprocess.run(
                [sys.executable, script] + args,
                cwd=work,
                capture_output=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return dict(
                id=m["id"],
                replay=m["replay"],
                status="error",
                detail="timeout %ds" % timeout,
                secs=time.time() - t0,
            )
        want = read_text(os.path.join(src, expected)).split("\n", 1)[1]
        got_all = p.stdout.decode("utf-8", "replace")
        got = got_all.split("\n", 1)[1] if "\n" in got_all else ""
        fails = [ln for ln in got.splitlines() if ln.startswith("[FAIL]")]
        if p.returncode != 0:
            status = "killed-hard"
        elif got != want:
            status = "killed-soft"
        else:
            status = "survived"
        detail = fails[0][:200] if fails else ""
        if not fails and p.returncode != 0:
            detail = (
                got.strip().splitlines() or [p.stderr.decode("utf-8", "replace")[-200:]]
            )[-1][:200]
        if status == "killed-soft":
            a = got.splitlines()
            b = want.splitlines()
            diff = [(x, y) for x, y in zip(a, b) if x != y]
            detail = (
                ("differing line: %r" % (diff[0][0][:160],))
                if diff
                else "output length differs"
            )
        return dict(
            id=m["id"],
            replay=m["replay"],
            status=status,
            exit=p.returncode,
            detail=detail,
            secs=time.time() - t0,
        )
    except Exception as ex:  # mutant could not be applied
        return dict(
            id=m["id"],
            replay=m["replay"],
            status="error",
            detail=str(ex)[:300],
            secs=time.time() - t0,
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def md_cell(s):
    return re.sub(r"\s+", " ", s).replace("|", "\\|")


def write_report(mutants, results, equiv, path, notes=None):
    by = {m["id"]: m for m in mutants}
    lines = [
        "# Self-audit result",
        "",
        "Rules: `PREREG.md`. Generated by `run_selfaudit.py`; per-mutant data in `results.json`.",
        "",
    ]
    lines += [
        "| replay | mutants | killed (hard) | killed (soft) | killed / total | survivors | errors |",
        "|---|---|---|---|---|---|---|",
    ]
    tot = dict(n=0, h=0, s=0, sv=0, er=0)
    surv = []
    for name in REPLAYS:
        rs = [r for r in results if r["replay"] == name]
        if not rs:
            continue
        h = sum(r["status"] == "killed-hard" for r in rs)
        s = sum(r["status"] == "killed-soft" for r in rs)
        sv = [r for r in rs if r["status"] == "survived"]
        er = sum(r["status"] == "error" for r in rs)
        lines.append(
            "| %s | %d | %d | %d | %d / %d | %s | %d |"
            % (
                name,
                len(rs),
                h,
                s,
                h + s,
                len(rs),
                ", ".join(r["id"] for r in sv) or "none",
                er,
            )
        )
        tot["n"] += len(rs)
        tot["h"] += h
        tot["s"] += s
        tot["sv"] += len(sv)
        tot["er"] += er
        surv += sv
    lines.append(
        "| **all** | %d | %d | %d | %d / %d | %d | %d |"
        % (
            tot["n"],
            tot["h"],
            tot["s"],
            tot["h"] + tot["s"],
            tot["n"],
            tot["sv"],
            tot["er"],
        )
    )
    lines += [
        "",
        "Soft kills: the replay exited 0 and reported all claims as passing on the mutated data; only the text comparison with the expected output caught it.",
        "",
    ]
    lines += ["## Survivors", ""]
    if not surv:
        lines.append("None.")
    for r in surv:
        m = by[r["id"]]
        lines.append(
            "### %s (%s, class %s)" % (r["id"], r["replay"], m.get("cls", "?"))
        )
        lines.append("")
        lines.append("- mutation: %s" % m["desc"])
        lines.append("- why it should break the claim: %s" % m["why"])
        if r["id"] in equiv:
            lines.append("- classification: equivalent. Argument: %s" % equiv[r["id"]])
        else:
            lines.append(
                "- classification: non-equivalent survivor (a finding about the checker; no equivalence argument)"
            )
        lines.append("")
    if notes:
        lines += ["## Notes", ""]
        for k, v in notes.items():
            lines.append("- %s: %s" % (k, v))
        lines.append("")
    lines += [
        "## All mutants",
        "",
        "| id | replay | class | status | evidence |",
        "|---|---|---|---|---|",
    ]
    for r in results:
        m = by[r["id"]]
        lines.append(
            "| %s | %s | %s | %s | %s |"
            % (
                r["id"],
                r["replay"],
                m.get("cls", "?"),
                r["status"],
                md_cell(r.get("detail", ""))[:170],
            )
        )
    write_text(path, "\n".join(lines) + "\n")
    return tot, surv


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--only", default="")
    ap.add_argument("-j", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--report-only", action="store_true", help="regenerate RESULT.md from results.json and equivalence.json without running mutants")
    a = ap.parse_args()
    mutants = json.loads(read_text(os.path.join(HERE, "mutants.json")))["mutants"]
    ids = [m["id"] for m in mutants]
    assert len(ids) == len(set(ids)), "duplicate mutant id"
    sel = [m for m in mutants if a.only in m["replay"]]
    equiv_path = os.path.join(HERE, "equivalence.json")
    eq = json.loads(read_text(equiv_path)) if os.path.isfile(equiv_path) else {}
    equiv, notes = eq.get("equivalent", {}), eq.get("notes", {})
    results = []
    if a.report_only:
        results = json.loads(read_text(os.path.join(HERE, "results.json")))
        sel = [m for m in mutants if m["id"] in {r["id"] for r in results}]
        tot, surv = write_report(sel, results, equiv, os.path.join(HERE, "RESULT.md"), notes)
        print("report regenerated: %d mutants, %d survivors" % (tot["n"], tot["sv"]))
        return 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.j)) as ex:
        futs = {ex.submit(run_mutant, m, a.timeout): m for m in sel}
        for f in concurrent.futures.as_completed(futs):
            r = f.result()
            results.append(r)
            print(
                "%-14s %-10s %-12s %5.0fs"
                % (r["id"], r["replay"], r["status"], r["secs"]),
                flush=True,
            )
    order = {m["id"]: k for k, m in enumerate(mutants)}
    results.sort(key=lambda r: order[r["id"]])
    if a.only:
        out_md, out_js = (
            os.path.join(HERE, "RESULT.partial.md"),
            os.path.join(HERE, "results.partial.json"),
        )
    else:
        out_md, out_js = (
            os.path.join(HERE, "RESULT.md"),
            os.path.join(HERE, "results.json"),
        )
    write_text(out_js, json.dumps(results, indent=1) + "\n")
    tot, surv = write_report(sel, results, equiv, out_md, notes)
    print(
        "\n%d mutants: %d killed (hard %d, soft %d), %d survived, %d errors -> %s"
        % (
            tot["n"],
            tot["h"] + tot["s"],
            tot["h"],
            tot["s"],
            tot["sv"],
            tot["er"],
            os.path.basename(out_md),
        )
    )
    bad = [r for r in surv if r["id"] not in equiv]
    if tot["er"]:
        return 2
    return 1 if (a.strict and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
