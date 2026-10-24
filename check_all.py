#!/usr/bin/env python3
"""Re-run every replay and compare with the committed expected output, byte for byte.

  python check_all.py             # full runs (the two longest replays take about 10-15 minutes on one core)
  python check_all.py --quick     # use the documented quick subset where a replay has one (e101); others run in full
  python check_all.py --only e153 # restrict to replays whose directory name contains the given text
  python check_all.py -j 4        # number of replays run in parallel (default: up to 8)

Checks, in order:
  1. every per-directory SHA256SUMS verifies (file bytes unchanged);
  2. CHECKERS.sha256 matches the source file(s) of every replay script (CRLF normalised to LF);
  3. each replay's stdout equals its expected output file byte for byte (this includes the opening line,
     `checker-sha256: <hash>`, which must also equal the CHECKERS.sha256 entry).
Exit status 0 only if every check agrees; any difference gives a nonzero exit status.
"""

import argparse
import concurrent.futures
import difflib
import hashlib
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))

# (directory, script, expected file, quick arguments, quick expected file)
REPLAYS = [
    ("replays/wall_d2", "replay_wall.py", "expected_output.txt", None, None),
    (
        "replays/e059_zariski_multiplicity",
        "replay_zariski.py",
        "expected_output.txt",
        None,
        None,
    ),
    (
        "replays/e067_campana_peternell",
        "replay_campana_peternell.py",
        "expected_output.txt",
        None,
        None,
    ),
    (
        "replays/e101_simplex_isotropic",
        "replay_simplex.py",
        "expected_output.txt",
        "--quick",
        "expected_output_quick.txt",
    ),
    (
        "replays/e138_subset_sum",
        "replay_subset_sum.py",
        "expected_output.txt",
        None,
        None,
    ),
    (
        "replays/e153_bernoulli_non_pisot",
        "replay_bernoulli.py",
        "expected_output.txt",
        None,
        None,
    ),
    (
        "replays/e268_haldane_gap",
        "replay_haldane.py",
        "expected_output.txt",
        None,
        None,
    ),
    (
        "replays/e318_chromatic_splitting",
        "replay_chromatic.py",
        "expected_output.txt",
        None,
        None,
    ),
    ("replays/e352_calabi_flow", "replay_calabi.py", "expected_output.txt", None, None),
    ("ramsey_alphaevolve", "check_ramsey.py", "expected_output.txt", None, None),
]
# check_ramsey.py needs an explicit mode; the other scripts run with no arguments.
FIXED_ARGS = {"check_ramsey.py": ["--cheap"]}


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()


def verify_sums(directory):
    """sha256sum -c equivalent for <directory>/SHA256SUMS; returns a list of problems."""
    problems = []
    sums = os.path.join(directory, "SHA256SUMS")
    if not os.path.isfile(sums):
        return ["%s: SHA256SUMS missing" % directory]
    for line in read_bytes(sums).decode("utf-8").splitlines():
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        name = name.lstrip("*").strip()
        path = os.path.join(directory, name)
        if not os.path.isfile(path):
            problems.append("%s: listed file %s missing" % (directory, name))
        elif sha256_bytes(read_bytes(path)) != digest:
            problems.append("%s: %s differs from SHA256SUMS" % (directory, name))
    return problems


def checker_hash(files):
    h = hashlib.sha256()
    for p in files:
        h.update(read_bytes(p).replace(b"\r\n", b"\n"))
    return h.hexdigest()


def load_checkers():
    """CHECKERS.sha256 lines: '<hash>  <path>[+<path>...]'; several paths are concatenated in the order listed."""
    out = {}
    path = os.path.join(ROOT, "CHECKERS.sha256")
    if not os.path.isfile(path):
        return out, ["CHECKERS.sha256 missing"]
    for line in read_bytes(path).decode("utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        digest, spec = line.split(None, 1)
        out[spec.strip().split("+")[0]] = (digest, spec.strip().split("+"))
    return out, []


def run_one(entry, quick):
    d, script, exp, qargs, qexp = entry
    args = list(FIXED_ARGS.get(script, []))
    expected = exp
    if quick and qargs:
        args.append(qargs)
        expected = qexp
    cwd = os.path.join(ROOT, d)
    t0 = time.time()
    p = subprocess.run([sys.executable, script] + args, cwd=cwd, capture_output=True)
    got = p.stdout
    want = read_bytes(os.path.join(cwd, expected))
    return d, script, args, expected, p.returncode, got, want, time.time() - t0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("-j", type=int, default=min(8, os.cpu_count() or 1))
    a = ap.parse_args()
    entries = [e for e in REPLAYS if a.only in e[0]]
    if not entries:
        print("no replay matches --only %r" % a.only)
        return 2
    problems = []

    for e in entries:
        problems += verify_sums(os.path.join(ROOT, e[0]))
    print(
        "SHA256SUMS: %s" % ("OK" if not problems else "%d problem(s)" % len(problems)),
        flush=True,
    )

    checkers, errs = load_checkers()
    problems += errs
    hashes = {}
    nbad = len(problems)
    for d, script, _, _, _ in entries:
        key = d + "/" + script
        if key not in checkers:
            problems.append("CHECKERS.sha256: no entry for %s" % key)
            continue
        digest, paths = checkers[key]
        actual = checker_hash([os.path.join(ROOT, p) for p in paths])
        hashes[d] = digest
        if actual != digest:
            problems.append(
                "CHECKERS.sha256: %s lists %s, source hashes to %s"
                % (key, digest, actual)
            )
    print(
        "CHECKERS.sha256: %s" % ("OK" if len(problems) == nbad else "problem(s)"),
        flush=True,
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.j)) as ex:
        futs = [ex.submit(run_one, e, a.quick) for e in entries]
        for f in concurrent.futures.as_completed(futs):
            d, script, args, expected, rc, got, want, dt = f.result()
            ok = rc == 0 and got == want
            head = got.split(b"\n", 1)[0].decode("utf-8", "replace")
            if d in hashes and head != "checker-sha256: " + hashes[d]:
                ok = False
                problems.append(
                    "%s: opening output line %r does not match CHECKERS.sha256"
                    % (d, head)
                )
            print(
                "%-40s %-4s exit=%d %6.0fs  vs %s"
                % (d, "OK" if ok else "FAIL", rc, dt, expected),
                flush=True,
            )
            if rc != 0:
                problems.append("%s: exit status %d" % (d, rc))
            if got != want:
                diff = difflib.unified_diff(
                    want.decode("utf-8", "replace").splitlines(),
                    got.decode("utf-8", "replace").splitlines(),
                    "expected",
                    "actual",
                    lineterm="",
                )
                problems.append(
                    "%s: output differs from %s\n    %s"
                    % (d, expected, "\n    ".join(list(diff)[:20]))
                )

    if problems:
        print("\nFAILURES:")
        for p in problems:
            print(" -", p)
        return 1
    print("\nALL CHECKS PASSED (%d replays)" % len(entries))
    return 0


if __name__ == "__main__":
    sys.exit(main())
