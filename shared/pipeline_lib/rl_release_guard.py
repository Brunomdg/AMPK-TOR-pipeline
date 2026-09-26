#!/usr/bin/env python3
# >>> deposit header
# Purpose: Before a job writes anything, check that every output path it will write lies inside the working
#          copy of the project (or the second allowed folder) and outside the protected folders, and that
#          PROJECT_ROOT points at the working copy; exit with status 2 otherwise. It writes nothing itself.
# Methods: several steps (shared helper)
# Inputs:  output paths on the command line or in a list file; the environment variables `PROJECT_ROOT`,
#          `RELEASE_ROOT`, `BASELINE_ROOT`, `LIVE_ROOT`, `R1_ROOT` and `STAGING_ROOT`
# Outputs: one line per path (OK, FORBIDDEN or OUTSIDE), printed; exit status 0 or 2
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
# Written for the working copy of the project, `TOR-AMPK_orthologs_10target_R1`, set up on 16 September
# 2026 as a full copy of the project in which the ten-target tables were rebuilt. Adapted from the same
# check in an earlier working copy (`TOR-AMPK_remediation_R1`, 13 September 2026): the allowed and
# forbidden folders changed, the logic did not. The script prints its docstring as its usage message,
# so the docstring is kept exactly as it ran.
"""OUTPUT-PATH GUARD for the sister working project TOR-AMPK_orthologs_10target_R1.

Every remediation job lists the paths it intends to write and calls this guard first. The guard
refuses (exit 2) if any path resolves outside the working release or the staging area, or inside
the canonical project (left untouched by Bruno's ruling), the frozen baseline, the remediation_R1
release, the export mirrors, or a protected folder. It also refuses if PROJECT_ROOT is not the
sister, so a script that reads PROJECT_ROOT cannot fall back to the canonical project. Nothing here writes anything.

    rl_release_guard.py --check PATH [PATH ...]        check the listed output paths
    rl_release_guard.py --check-file LIST.txt          one path per line, '#' comments allowed
    rl_release_guard.py --env                          check PROJECT_ROOT only
    rl_release_guard.py --selftest                     the guard's own tests (no filesystem writes)

Env: RELEASE_ROOT (default /gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1),
     BASELINE_ROOT, LIVE_ROOT (the canonical project), R1_ROOT, STAGING_ROOT
     (default ~/q0_stage/remediation/sister_10target_R1). Relative paths resolve against RELEASE_ROOT.
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import os
import sys

RELEASE = os.path.realpath(os.environ.get("RELEASE_ROOT") or
                           (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1"))
BASELINE = os.path.realpath(os.environ.get("BASELINE_ROOT") or
                            (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_baseline_20260913"))
LIVE = os.path.realpath(os.environ.get("LIVE_ROOT") or
                        (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_across_210_eukaryotes"))
R1 = os.path.realpath(os.environ.get("R1_ROOT") or
                      (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_remediation_R1"))
STAGING = os.path.realpath(os.environ.get("STAGING_ROOT") or
                           os.path.expanduser("~/q0_stage/remediation/sister_10target_R1"))
HOME = os.path.realpath(os.path.expanduser("~"))
FORBIDDEN = [BASELINE, LIVE, R1,
             os.path.join(HOME, "mirror_export"), os.path.join(HOME, "full_export"),
             os.path.join(HOME, "q0_stage", "plm", "withdrawn_5163525")]
ALLOWED = [RELEASE, STAGING]


def _under(p, root):
    return p == root or p.startswith(root.rstrip("/") + "/")


def classify(path):
    p = os.path.realpath(path if os.path.isabs(path) else os.path.join(RELEASE, path))
    for f in FORBIDDEN:
        if _under(p, f):
            return "FORBIDDEN", p, f
    for a in ALLOWED:
        if _under(p, a):
            return "OK", p, a
    return "OUTSIDE", p, "-"


def check_env():
    pr = os.environ.get("PROJECT_ROOT")
    if not pr:
        return False, "PROJECT_ROOT is unset (scripts would fall back to the live project)"
    if os.path.realpath(pr) != RELEASE:
        return False, "PROJECT_ROOT=%s is not the working release %s" % (pr, RELEASE)
    return True, "PROJECT_ROOT is the working release"


def selftest():
    cases = [(os.path.join(RELEASE, "04_final_tables", "x.tsv"), "OK"),
             (os.path.join(STAGING, "step1", "log.txt"), "OK"),
             (os.path.join(BASELINE, "04_final_tables", "x.tsv"), "FORBIDDEN"),
             (os.path.join(LIVE, "PURGE_LOG.tsv"), "FORBIDDEN"),
             (os.path.join(LIVE, "04_final_tables", "FINAL_ORTHOLOG_TABLE__replacement_layers.tsv"), "FORBIDDEN"),
             (os.path.join(R1, "07_findings", "x.tsv"), "FORBIDDEN"),
             (os.path.join(HOME, "mirror_export", "a"), "FORBIDDEN"),
             (os.path.join(RELEASE, "..", os.path.basename(BASELINE), "y"), "FORBIDDEN"),
             ((_TA210_SCRATCH + "/ms210_phase2/other"), "OUTSIDE"),
             ("07_findings/publication_tables/T1_ortholog_table.tsv", "OK")]
    bad = 0
    for path, want in cases:
        got = classify(path)[0]
        print("%-9s %-9s %s" % (want, got, path))
        bad += got != want
    print("SELFTEST", "PASS" if bad == 0 else "FAIL (%d)" % bad)
    return bad == 0


def main():
    a = sys.argv[1:]
    if not a or a[0] not in ("--check", "--check-file", "--env", "--selftest"):
        sys.exit(__doc__)
    if a[0] == "--selftest":
        sys.exit(0 if selftest() else 2)
    ok, msg = check_env()
    print("ENV:", msg)
    if a[0] == "--env":
        sys.exit(0 if ok else 2)
    paths = a[1:] if a[0] == "--check" else [l.strip() for l in open(a[1]) if l.strip() and not l.startswith("#")]
    bad = 0
    for p in paths:
        kind, real, root = classify(p)
        print("%-9s %s" % (kind, real))
        bad += kind != "OK"
    if not ok or bad:
        print("GUARD REFUSED: %d path(s) outside the sister or in a forbidden root; env ok=%s" % (bad, ok))
        sys.exit(2)
    print("GUARD OK: %d path(s) resolve inside the working release or staging" % len(paths))


if __name__ == "__main__":
    main()
