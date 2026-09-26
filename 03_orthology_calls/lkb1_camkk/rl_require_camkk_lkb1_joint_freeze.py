#!/usr/bin/env python3
# >>> deposit header
# Purpose: Checks, before any shared downstream step, that the final CAMKK and LKB1 species states are fixed and
#          unchanged (marker file, SHA-256 and recounted state counts of each state file, rule-registry checksum,
#          FINAL status of each target), and refuses otherwise; imported as a module or run as a command.
# Methods: several steps (shared helper)
# Inputs:  16_control_repair/`CAMKK_LKB1_POST_REPAIR_FREEZE.ready.json` (the marker file, or the path in
#          CAMKK_LKB1_FREEZE_SENTINEL); the CAMKK and LKB1 state files it names; LKB1_PROPAGATION/
#          `LKB1_AUTHORITY_STATUS.tsv` and CAMKK_AUTHORITY/`CAMKK_AUTHORITY_STATUS.tsv` (status tables);
#          05_method_and_rules/`CANONICAL_RULE_REGISTRY.tsv` (checksum only)
# Outputs: no file; as a command, exit status 0 or 1 with a short report; when imported, the checked contents of
#          the marker file, or a JointFreezeNotComplete exception
# Tools:   Python 3 standard library only
# <<< deposit header
"""The one joint check of the final CAMKK and LKB1 states.

Shared downstream analysis that uses CAMKK or LKB1 is blocked until the species states of both targets are final
and fixed ("frozen"), and this module is the only implementation of that rule.

It fails closed by design. The existence of the marker file (the "sentinel") means that both targets are final
and that regenerating the downstream analyses was approved. A missing marker file is therefore a refusal, not
permission, and a marker file that fails any check is a refusal too; no caller ever proceeds because a file was
absent. Every failure here, including input and output errors and malformed inputs, is raised as
JointFreezeNotComplete: callers catch that one class, so any other exception type would be a crash where a clean
refusal belongs.

Why one file: the downstream steps are written in Python, R and shell, and three implementations of "are both
targets final" would drift apart, the very kind of defect this repair exists to fix. Python callers import this
module; R and shell run it as a command and stop on a non-zero exit status.

    from rl_require_camkk_lkb1_joint_freeze import require_camkk_lkb1_joint_freeze
    require_camkk_lkb1_joint_freeze()            # raises JointFreezeNotComplete unless every check passes

    Rscript:  if (system2("python3", c(GUARD)) != 0) stop("joint freeze not complete")
    shell:    python3 rl_require_camkk_lkb1_joint_freeze.py || exit 1

What is checked, all of it, on every call:
  1 the marker file exists and parses as a JSON object;
  2 its STATUS reads JOINT_FREEZE_COMPLETE;
  3 for each of CAMKK and LKB1, the state file it names is a regular file,
  3a that lies inside PROJECT_ROOT, both paths resolved through symbolic links, so that a byte-identical copy in
     another project is refused (chosen on 16 September 2026);
  4 the file's SHA-256 equals the one in the marker file;
  5 its state counts, recounted from the file, are not empty and equal the counts declared in the marker file,
    which must be a non-empty object of non-negative JSON integers;
  6 the rule-registry checksum stored in the marker file matches the current CANONICAL_RULE_REGISTRY.tsv;
  7 for each target, the last row of its append-only status table reads FINAL and names, by SHA-256, the same
    state file as the marker file.
    LKB1: 16_control_repair/LKB1_PROPAGATION/LKB1_AUTHORITY_STATUS.tsv (written by lk_status.py).
    CAMKK: 16_control_repair/CAMKK_AUTHORITY/CAMKK_AUTHORITY_STATUS.tsv (written by ck_status.py). This was
    added when LKB1 became final and the CAMKK half was the only part of the marker file that nothing else
    confirmed; requiring it was approved on 16 September 2026, when CAMKK was declared final.
Step 5 is the point of the exercise: a checksum proves the bytes, but recounting proves that the file is still
the state table it claims to be, so an accidental substitution cannot pass. A file with no rows for the target
counts as a substitution: an empty declaration used to match an empty count, which let the CAMKK entry name the
LKB1 file (found on 16 September 2026 in a deliberate attempt to break the check). Step 7 exists because a
target can be reopened after it was declared final (LKB1 was, on 16 September 2026, for the second version of
the contaminant screen); a marker file written against the earlier file would otherwise still pass, because its
checksum and counts are unchanged and only the decision behind them is not.

A status table or state table is malformed, and so refused, if its header repeats a column name or any data row
has a different number of fields from its header.

Before the joint check passes, target-specific repair and quality checks may run: they read the provisional LKB1
states (`LKB1_POST_REPAIR_STAGED_STATES.tsv`) through `load_staged_lkb1_states("target_specific")`. While the
LKB1 status is not FINAL, those states are the result of the first version of the screen. A shared analysis that
asks for LKB1 states is refused unless the joint check passes, and even then it receives exactly the bytes that
were checked and cannot name another file.
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import hashlib
import json
import os
import sys

ROOT = os.environ.get("PROJECT_ROOT") or \
    (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1")
CR = os.path.join(ROOT, "16_control_repair")
SENTINEL = os.environ.get("CAMKK_LKB1_FREEZE_SENTINEL") or \
    os.path.join(CR, "CAMKK_LKB1_POST_REPAIR_FREEZE.ready.json")
STAGED_LKB1 = os.path.join(CR, "LKB1_PROPAGATION", "LKB1_POST_REPAIR_STAGED_STATES.tsv")
LKB1_STATUS = os.environ.get("LKB1_AUTHORITY_STATUS") or \
    os.path.join(CR, "LKB1_PROPAGATION", "LKB1_AUTHORITY_STATUS.tsv")
CAMKK_STATUS = os.environ.get("CAMKK_AUTHORITY_STATUS") or \
    os.path.join(CR, "CAMKK_AUTHORITY", "CAMKK_AUTHORITY_STATUS.tsv")
STATUS_PATHS = {"LKB1": LKB1_STATUS, "CAMKK": CAMKK_STATUS}
RULE_REGISTRY = os.path.join(ROOT, "05_method_and_rules", "CANONICAL_RULE_REGISTRY.tsv")
TARGETS = ("CAMKK", "LKB1")
STATUS_REQUIRED_COLS = ("SEQ", "STATUS", "STATE_SOURCE_SHA256")


class JointFreezeNotComplete(RuntimeError):
    """Raised whenever the joint check cannot be shown to pass, for any reason. Inside this module it
    is caught only to convert another exception type into this one."""


def _parse_tsv(lines, what):
    """(header, rows) of a tab-separated table with '##' comment lines, refusing the malformed shapes
    described in the module docstring."""
    hdr, rows = None, []
    for n, line in enumerate(lines, 1):
        s = line.rstrip("\r\n")
        if s.startswith("##") or not s.strip():
            continue
        parts = s.split("\t")
        if hdr is None:
            if len(set(parts)) != len(parts):
                raise JointFreezeNotComplete("%s is malformed: its header repeats a column name"
                                             % what)
            hdr = parts
            continue
        if len(parts) != len(hdr):
            raise JointFreezeNotComplete(
                "%s is malformed: line %d has %d fields, the header has %d"
                % (what, n, len(parts), len(hdr)))
        rows.append(dict(zip(hdr, parts)))
    return hdr, rows


def _read_tsv(path, what):
    if not isinstance(path, str) or not path or not os.path.isfile(path):
        raise JointFreezeNotComplete("%s missing or not a regular file: %r" % (what, path))
    with open(path, encoding="utf-8") as fh:
        return _parse_tsv(fh, "%s %s" % (what, path))


def authority_status(target, path=None):
    """The current (last) row of a target's append-only status table, or raise.

    A table that is missing, empty, malformed or out of sequence is a refusal, not a default:
    a status nobody can read is not a FINAL status.
    """
    what = "%s authority status record" % target
    try:
        path = path or STATUS_PATHS[target]
        if not os.path.exists(path):
            raise JointFreezeNotComplete("%s missing: %s" % (what, path))
        hdr, recs = _read_tsv(path, what)
        if hdr is None or not recs or any(c not in hdr for c in STATUS_REQUIRED_COLS):
            raise JointFreezeNotComplete("%s is malformed or empty: %s" % (what, path))
        if any(not r[c] for r in recs for c in STATUS_REQUIRED_COLS):
            raise JointFreezeNotComplete("%s is malformed: a row leaves %s empty: %s"
                                         % (what, "/".join(STATUS_REQUIRED_COLS), path))
        if [r["SEQ"] for r in recs] != [str(i) for i in range(1, len(recs) + 1)]:
            raise JointFreezeNotComplete("%s is malformed: SEQ is not 1..%d in order: %s"
                                         % (what, len(recs), path))
        return recs[-1]
    except JointFreezeNotComplete:
        raise
    except Exception as exc:
        raise JointFreezeNotComplete("%s could not be read (%s: %s)"
                                     % (what, type(exc).__name__, exc))


def lkb1_authority_status(path=None):
    """LKB1's current status row (kept for the callers written before CAMKK had a status table)."""
    return authority_status("LKB1", path)


def _status_for_message(path=None, target="LKB1"):
    try:
        r = authority_status(target, path)
        return "%s (row %s)" % (r["STATUS"], r["SEQ"])
    except Exception as exc:
        return "unreadable: %s" % exc


def _sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _count_states(path, target):
    """State counts for one target, counted from the file's own columns, never taken from a stated number."""
    hdr, rows = _read_tsv(path, "%s state file" % target)
    if hdr is None or "TARGET" not in hdr:
        raise JointFreezeNotComplete("%s has no TARGET column; it is not a state table" % path)
    state_col = next((c for c in ("STATE", "NEW_STATE") if c in hdr), None)
    if state_col is None:
        raise JointFreezeNotComplete("%s has no STATE column" % path)
    counts = {}
    for r in rows:
        if r["TARGET"] == target:
            counts[r[state_col]] = counts.get(r[state_col], 0) + 1
    return counts


def governance_version():
    """Checksum of the rule registry in force, so that a marker file cannot outlive the rules it was written under."""
    return "rule-registry:%s" % _sha256(RULE_REGISTRY)[:16]


def _declared_counts(block, target):
    sc = block.get("state_counts")
    if (not isinstance(sc, dict) or not sc
            or any(not isinstance(k, str) or type(v) is not int or v < 0 for k, v in sc.items())):
        raise JointFreezeNotComplete(
            "%s state_counts must be a non-empty object of non-negative integers, got %r"
            % (target, sc))
    return dict(sc)


def _verify(sentinel, status, camkk_status=None):
    overrides = {"LKB1": status, "CAMKK": camkk_status}
    if not os.path.exists(sentinel):
        raise JointFreezeNotComplete(
            "CAMKK/LKB1 joint freeze is not complete: no sentinel at %s. Shared downstream "
            "analysis is blocked until both targets are final and the freeze is approved; a "
            "shared model may not mix a finished target with an unfinished one. LKB1 authority "
            "status: %s. CAMKK authority status: %s."
            % (sentinel, _status_for_message(status, "LKB1"),
               _status_for_message(camkk_status, "CAMKK")))
    try:
        with open(sentinel, encoding="utf-8") as fh:
            meta = json.load(fh)
    except (OSError, ValueError) as exc:
        raise JointFreezeNotComplete("sentinel %s does not parse: %s" % (sentinel, exc))
    if not isinstance(meta, dict):
        raise JointFreezeNotComplete("sentinel %s is not a JSON object" % sentinel)

    if meta.get("STATUS") != "JOINT_FREEZE_COMPLETE":
        raise JointFreezeNotComplete(
            "sentinel STATUS is %r, not JOINT_FREEZE_COMPLETE" % meta.get("STATUS"))

    want_gov = governance_version()
    if meta.get("GOVERNANCE_VERSION") != want_gov:
        raise JointFreezeNotComplete(
            "sentinel governance version %r does not match the live rule registry %r"
            % (meta.get("GOVERNANCE_VERSION"), want_gov))

    verified = {}
    for target in TARGETS:
        block = meta.get(target)
        if not isinstance(block, dict):
            raise JointFreezeNotComplete("sentinel has no block for %s" % target)
        path = block.get("path")
        if not isinstance(path, str) or not path or not os.path.isfile(path):
            raise JointFreezeNotComplete("%s frozen state file missing: %r" % (target, path))
        # A state file must lie inside PROJECT_ROOT: both paths are resolved through symbolic links and
        # compared at a directory boundary, so the marker file cannot name a byte-identical copy in
        # another project (chosen on 16 September 2026).
        root_real = os.path.realpath(ROOT)
        path_real = os.path.realpath(path)
        if not path_real.startswith(os.path.join(root_real, "")):
            raise JointFreezeNotComplete(
                "%s frozen state file is outside the project root %s: %s resolves to %s"
                % (target, root_real, path, path_real))
        observed = _sha256(path)
        if observed != block.get("sha256"):
            raise JointFreezeNotComplete(
                "%s state file hash mismatch: expected %s, observed %s"
                % (target, block.get("sha256"), observed))
        declared = _declared_counts(block, target)
        counted = _count_states(path, target)
        if not counted:
            raise JointFreezeNotComplete(
                "%s state file %s contains no %s rows, so it is not the %s state vector"
                % (target, path, target, target))
        if declared != counted:
            raise JointFreezeNotComplete(
                "%s state counts do not match the file: sentinel says %s, the file contains %s"
                % (target, declared, counted))
        # The marker-file field `gate_a_pass` is deliberately not checked: it was added to this check on
        # 16 September 2026 without being asked for, and it was removed again the same day.
        verified[target] = {"path": path, "sha256": observed, "counts": counted}

    for target in TARGETS:
        current = authority_status(target, overrides[target])
        if current["STATUS"] != "FINAL":
            raise JointFreezeNotComplete(
                "%s authority status is %s (status record row %s), not FINAL, so %s is not "
                "frozen and may not enter shared analysis"
                % (target, current["STATUS"], current["SEQ"], target))
        if current["STATE_SOURCE_SHA256"] != verified[target]["sha256"]:
            raise JointFreezeNotComplete(
                "the sentinel names %s state file sha256 %s, but the FINAL %s status (row %s) "
                "pins %s" % (target, verified[target]["sha256"], target, current["SEQ"],
                             current["STATE_SOURCE_SHA256"]))
    meta["_VERIFIED"] = verified
    return meta


def require_camkk_lkb1_joint_freeze(sentinel=SENTINEL, status=None, camkk_status=None):
    """Return the checked marker file, or raise JointFreezeNotComplete; callers need not check a
    return value. The returned dict carries _VERIFIED: for each target, the path, the SHA-256
    observed and the counts recounted during this call."""
    try:
        return _verify(sentinel, status, camkk_status)
    except JointFreezeNotComplete:
        raise
    except Exception as exc:
        raise JointFreezeNotComplete("the joint freeze could not be verified (%s: %s)"
                                     % (type(exc).__name__, exc))


def load_staged_lkb1_states(caller_class, path=None):
    """Read LKB1 states.

    caller_class must be 'target_specific' or 'shared'.
    target_specific: reads `path`, by default the provisional LKB1 states file, without the joint check.
    shared: refused unless the joint check passes; then returns the rows of exactly the bytes checked
    for LKB1, checksummed again at read time. A shared caller may not pass `path`, because choosing the
    file would let it read something that nobody checked.
    """
    if caller_class not in ("target_specific", "shared"):
        raise ValueError("caller_class must be 'target_specific' or 'shared', got %r"
                         % caller_class)
    if caller_class == "target_specific":
        return _read_tsv(path or STAGED_LKB1, "LKB1 state file")[1]
    if path is not None:
        raise JointFreezeNotComplete(
            "a shared caller may not choose the LKB1 state file (%r); it reads the file the "
            "verified joint freeze pins" % (path,))
    meta = require_camkk_lkb1_joint_freeze()
    v = meta["_VERIFIED"]["LKB1"]
    try:
        with open(v["path"], "rb") as fh:
            data = fh.read()
        if _sha256_bytes(data) != v["sha256"]:
            raise JointFreezeNotComplete(
                "LKB1 state file %s changed between verification and read" % v["path"])
        return _parse_tsv(data.decode("utf-8").splitlines(True),
                          "LKB1 state file %s" % v["path"])[1]
    except JointFreezeNotComplete:
        raise
    except Exception as exc:
        raise JointFreezeNotComplete("verified LKB1 state file could not be read (%s: %s)"
                                     % (type(exc).__name__, exc))


def main(argv):
    quiet = "--quiet" in argv
    try:
        meta = require_camkk_lkb1_joint_freeze()
    except JointFreezeNotComplete as exc:
        if not quiet:
            sys.stderr.write("REFUSED: %s\n" % exc)
        return 1
    if not quiet:
        print("JOINT FREEZE VERIFIED  created=%s" % meta.get("CREATED_AT", "?"))
        for t in TARGETS:
            v = meta["_VERIFIED"][t]
            print("   %-6s %s  %s" % (t, v["counts"], v["path"]))
        print("   LKB1 authority status %s" % _status_for_message(None, "LKB1"))
        print("   CAMKK authority status %s" % _status_for_message(None, "CAMKK"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
