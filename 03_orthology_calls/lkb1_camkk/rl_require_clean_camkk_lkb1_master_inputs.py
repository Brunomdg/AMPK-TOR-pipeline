#!/usr/bin/env python3
# >>> deposit header
# Purpose: Refuse to let later steps read the combined ortholog table while an outdated CAMKK or LKB1
#          input could still reach it, by checking the table of its inputs against the files on disk.
# Methods: Orthology calls
# Inputs:  `CAMKK_LKB1_MASTER_ORTHOLOG_INPUT_AUDIT.tsv` (the inputs of the combined ortholog table, with
#          checksums), `CAMKK_LKB1_POST_REPAIR_FREEZE.ready.json`, `CANONICAL_RULE_REGISTRY.tsv`,
#          `ACTIVE_PUBLICATION_TARGETS.tsv`
# Outputs: none; exit status 0 and a short summary when every condition holds, otherwise 1 and the reason
# Tools:   Python 3 standard library; runs the check in rl_require_camkk_lkb1_joint_freeze.py;
#          versions: see envs/
# <<< deposit header
"""Check that the combined ortholog table no longer depends on outdated CAMKK or LKB1 inputs.

The combined ortholog table (the master ortholog layer: the ortholog calls of the ten target components
from which the publication tables are built) must be rebuilt from the final CAMKK and LKB1 results before
publication tables T1 and T2 (the ortholog table and the species presence table) are written. This order
was fixed on 16 September 2026 and is enforced here with checksums. No step that reads the combined table
may run while an outdated CAMKK or LKB1 input could reach it, and this module is the only implementation
of that check.

The check fails closed: a missing table of inputs is a refusal, and so is one that does not meet every
condition below. The refusal is raised before any output file is opened.

    from rl_require_clean_camkk_lkb1_master_inputs import require_clean_camkk_lkb1_master_inputs
    require_clean_camkk_lkb1_master_inputs()        # raises unless every condition below holds

    Rscript:  if (system2("python3", c(GUARD)) != 0) stop("master inputs not clean")
    shell:    python3 rl_require_clean_camkk_lkb1_master_inputs.py || exit 1

This is not the joint-freeze check, and the two mean different things. The joint-freeze marker file
(`CAMKK_LKB1_POST_REPAIR_FREEZE.ready.json`) states that the CAMKK and LKB1 analyses are final; this check
states that the combined ortholog table has been rebuilt from them. Both are required, so this check also
runs the joint-freeze check rather than assuming it ran. The joint-freeze marker on its own never permits
writing T1 or T2.

Conditions checked on every call:
  1 the table of inputs (`CAMKK_LKB1_MASTER_ORTHOLOG_INPUT_AUDIT.tsv`) exists, parses, and has the
    required columns;
  2 no outdated input of the combined table lacks an existing replacement, counted from the rows. The
    summary lines ACTIVE_STALE_CAMKK_LKB1_SOURCES and MASTER_ORTHOLOG_INPUT_STATUS must agree with that
    count, so editing the summary to read CLEAN cannot pass the check;
  3 every replacement named in a row the check relies on exists on disk;
  4 every file the check relies on still has the checksum stored for it, so a table of inputs written
    before an input changed cannot let the table-writing steps proceed;
  5 no row the check relies on is marked REVIEW_REQUIRED, AMBIGUOUS or MISSING, and every file those rows
    name exists on disk;
  6 the CAMKK and LKB1 checksums stored in the table of inputs equal those in the joint-freeze marker,
    the marker file itself still has its stored checksum, and the joint-freeze check passes;
  7 the stored version of the project's rule table (GOVERNANCE_VERSION) equals the checksum prefix of the
    current `CANONICAL_RULE_REGISTRY.tsv`;
  8 `ACTIVE_PUBLICATION_TARGETS.tsv` lists exactly the approved ten target components and still has its
    stored checksum, so the set of targets cannot be widened after the table of inputs was written.
Condition 8 stops the ten deferred target components from being brought back into the publication tables
from old sources after the CAMKK and LKB1 inputs were replaced.
"""
import hashlib
import json
import os
import sys

def _project_root_from_file(start):
    """Return the project root: the nearest folder above this file that contains 05_method_and_rules.

    Locating the root from the file, rather than from the working directory, stops the check from
    inspecting whatever folder it was started in while still printing MASTER INPUTS CLEAN. The path is
    made absolute but symbolic links are not resolved, so in a project copy built from links the root
    found is the tree that contains the link, not the tree the link points to.
    """
    d = os.path.dirname(os.path.abspath(start))
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "05_method_and_rules")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("REFUSED: PROJECT_ROOT unset and no project root found above %s" % start)


ROOT = os.environ.get("PROJECT_ROOT") or _project_root_from_file(__file__)
CR = os.path.join(ROOT, "16_control_repair")
AUDIT = os.environ.get("CAMKK_LKB1_MASTER_AUDIT") or os.path.join(
    CR, "phase2_discriminant", "master_rebuild", "CAMKK_LKB1_MASTER_ORTHOLOG_INPUT_AUDIT.tsv")
SENTINEL = os.environ.get("CAMKK_LKB1_FREEZE_SENTINEL") or os.path.join(
    CR, "CAMKK_LKB1_POST_REPAIR_FREEZE.ready.json")
REGISTRY = os.path.join(ROOT, "05_method_and_rules", "CANONICAL_RULE_REGISTRY.tsv")
ACTIVE = os.path.join(ROOT, "05_method_and_rules", "ACTIVE_PUBLICATION_TARGETS.tsv")
APPROVED_TEN = ("TOR", "RAPTOR", "RICTOR", "LST8", "AMPKa", "AMPKb", "AMPKg", "SIN1",
                "CAMKK", "LKB1")
BAD_STATES = ("REVIEW_REQUIRED", "AMBIGUOUS", "MISSING")


class MasterInputsNotClean(RuntimeError):
    """Raised whenever the inputs cannot be shown to be current. Within this module only main() catches
    it, to print the reason and return exit status 1."""


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _read(path):
    meta, hdr, rows = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##"):
                t = line[2:].strip()
                if "=" in t:
                    k, v = t.split("=", 1)
                    meta[k.strip()] = v.strip()
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            if len(f) != len(hdr):
                raise MasterInputsNotClean(
                    "audit row has %d fields for a %d-column header: %s" % (len(f), len(hdr), f[:2]))
            rows.append(dict(zip(hdr, f)))
    if hdr is None:
        raise MasterInputsNotClean("audit %s has no header" % path)
    if len(set(hdr)) != len(hdr):
        raise MasterInputsNotClean("audit header repeats a column name")
    return meta, hdr, rows


def require_clean_camkk_lkb1_master_inputs(audit=None, root=None):
    """Return the summary fields of the table of inputs once every condition holds, or raise
    MasterInputsNotClean. Callers need not check the return value."""
    audit = audit or AUDIT
    root = root or ROOT
    if not os.path.exists(audit):
        raise MasterInputsNotClean(
            "no master-input audit at %s. The master ortholog layer is not proven free of stale "
            "CAMKK/LKB1 sources, so publishing from it is blocked. The joint freeze alone does not "
            "authorise it: the sentinel says the target analyses are final, not that the master "
            "layer has incorporated them." % audit)
    meta, hdr, rows = _read(audit)
    for col in ("FILE", "ROLE", "IN_ACTIVE_DEPENDENCY_GRAPH", "STALE", "REPLACEMENT_SOURCE",
                "CURRENT_HASH", "ACTION"):
        if col not in hdr:
            raise MasterInputsNotClean("audit has no %s column; it is not a dependency manifest" % col)

    # Condition 2: count the outdated inputs from the rows; the summary line is never trusted.
    def is_master_input(r):
        return (r["ROLE"] == "master_ortholog_source" and r["STALE"] == "yes"
                and r["IN_ACTIVE_DEPENDENCY_GRAPH"] == "yes")

    relied, blocking = [], []
    for r in rows:
        if r["IN_ACTIVE_DEPENDENCY_GRAPH"] == "yes" or is_master_input(r):
            relied.append(r)
        if not is_master_input(r):
            continue
        rep = r["REPLACEMENT_SOURCE"]
        if rep in ("", "-") or not os.path.exists(os.path.join(root, rep)):
            blocking.append(r["FILE"])
    if blocking:
        raise MasterInputsNotClean(
            "%d stale master ortholog source(s) have no verified replacement: %s. Regenerate the "
            "master layer from the frozen packages before publishing from it."
            % (len(blocking), ", ".join(sorted(blocking)[:5])))

    declared = meta.get("ACTIVE_STALE_CAMKK_LKB1_SOURCES")
    if declared is None or not declared.strip().isdigit() or int(declared) != len(blocking):
        raise MasterInputsNotClean(
            "the audit declares ACTIVE_STALE_CAMKK_LKB1_SOURCES=%r but the rows recount to %d; the "
            "declaration is never trusted over the rows" % (declared, len(blocking)))
    status = meta.get("MASTER_ORTHOLOG_INPUT_STATUS")
    if status != "CLEAN":
        raise MasterInputsNotClean("audit MASTER_ORTHOLOG_INPUT_STATUS is %r, not CLEAN" % status)

    # Conditions 3 to 5: the rows the check relies on must still match the files on disk.
    for r in relied:
        if r["ACTION"] in BAD_STATES or r["STALE"] in BAD_STATES:
            raise MasterInputsNotClean("%s is %s in the audit" % (r["FILE"], r["ACTION"]))
        p = os.path.join(root, r["FILE"])
        if not os.path.exists(p):
            raise MasterInputsNotClean("audited file is missing from disk: %s" % r["FILE"])
        got = _sha(p)
        if got != r["CURRENT_HASH"]:
            raise MasterInputsNotClean(
                "%s changed since the audit was written (audit %s, disk %s). A stale audit must "
                "not keep the writers open." % (r["FILE"], r["CURRENT_HASH"][:16], got[:16]))
        rep = r.get("REPLACEMENT_SOURCE", "-")
        if rep not in ("", "-") and not os.path.exists(os.path.join(root, rep)):
            raise MasterInputsNotClean("replacement named for %s is missing: %s" % (r["FILE"], rep))

    # Condition 6: the joint-freeze marker is checked, not assumed.
    if not os.path.exists(SENTINEL):
        raise MasterInputsNotClean("no joint-freeze sentinel at %s" % SENTINEL)
    sent = json.load(open(SENTINEL, encoding="utf-8"))
    for t in ("CAMKK", "LKB1"):
        want = meta.get("%s_FREEZE_SHA256" % t)
        if want != sent.get(t, {}).get("sha256"):
            raise MasterInputsNotClean(
                "the audit pins %s freeze %s but the sentinel now says %s; the audit predates the "
                "current freeze" % (t, (want or "none")[:16],
                                    (sent.get(t, {}).get("sha256") or "none")[:16]))
    # Compare the checksum of the marker file itself, not only the two checksums it carries. The
    # table of inputs always stored JOINT_FREEZE_SENTINEL_SHA256, but until 16 September 2026 nothing
    # compared it, so a rewritten marker file passed unnoticed and the check still printed CLEAN.
    # A stored checksum that is never compared is not a check.
    want_sent = meta.get("JOINT_FREEZE_SENTINEL_SHA256")
    got_sent = _sha(SENTINEL)
    if want_sent != got_sent:
        raise MasterInputsNotClean(
            "the audit pins the joint-freeze sentinel at %s but it now hashes to %s; the sentinel "
            "was rewritten after the audit" % ((want_sent or "none")[:16], got_sent[:16]))

    sys.path.insert(0, os.path.join(CR, "scripts"))
    try:
        from rl_require_camkk_lkb1_joint_freeze import require_camkk_lkb1_joint_freeze
        require_camkk_lkb1_joint_freeze()
    except ImportError:
        raise MasterInputsNotClean("the joint-freeze check is not importable; cannot prove the freeze")
    except Exception as exc:
        raise MasterInputsNotClean("the joint freeze does not verify: %s" % exc)

    # Condition 7: the stored version of the project's rule table matches the current file.
    want_gov = "rule-registry:%s" % _sha(REGISTRY)[:16]
    if meta.get("GOVERNANCE_VERSION") != want_gov:
        raise MasterInputsNotClean(
            "audit governance version %r does not match the live rule registry %r"
            % (meta.get("GOVERNANCE_VERSION"), want_gov))

    # Condition 8: the list of target components for publication is exactly the approved ten.
    if not os.path.exists(ACTIVE):
        raise MasterInputsNotClean("no ACTIVE_PUBLICATION_TARGETS.tsv at %s" % ACTIVE)
    got_targets = []
    with open(ACTIVE, encoding="utf-8", errors="replace") as fh:
        hdr = None
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            got_targets.append(f[hdr.index("TARGET")])
    if sorted(got_targets) != sorted(APPROVED_TEN):
        extra = sorted(set(got_targets) - set(APPROVED_TEN))
        missing = sorted(set(APPROVED_TEN) - set(got_targets))
        raise MasterInputsNotClean(
            "the active publication scope is not the approved ten: %d rows%s%s"
            % (len(got_targets),
               ", unexpected %s" % extra if extra else "",
               ", missing %s" % missing if missing else ""))
    got_hash = _sha(ACTIVE)
    if meta.get("ACTIVE_PUBLICATION_TARGETS_SHA256") != got_hash:
        raise MasterInputsNotClean(
            "ACTIVE_PUBLICATION_TARGETS.tsv changed since the audit (audit %s, disk %s); the scope "
            "cannot be widened after the audit was written"
            % ((meta.get("ACTIVE_PUBLICATION_TARGETS_SHA256") or "none")[:16], got_hash[:16]))
    return meta


def main(argv):
    quiet = "--quiet" in argv
    try:
        meta = require_clean_camkk_lkb1_master_inputs()
    except MasterInputsNotClean as exc:
        if not quiet:
            sys.stderr.write("REFUSED: %s\n" % exc)
        return 1
    if not quiet:
        print("MASTER INPUTS CLEAN  stale=0  scope=%s" % meta.get("ACTIVE_PUBLICATION_TARGETS"))
        print("   superseded master inputs with replacement: %s"
              % meta.get("SUPERSEDED_MASTER_INPUTS_WITH_REPLACEMENT", "?"))
        print("   downstream outputs awaiting regeneration:  %s"
              % meta.get("DOWNSTREAM_OUTPUTS_AWAITING_REGENERATION", "?"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
