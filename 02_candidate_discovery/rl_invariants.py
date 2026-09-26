#!/usr/bin/env python3
# >>> deposit header
# Purpose: Shared helper that keeps the SHA-256 checksum of each reference input fixed in advance (controls, decoys, pilot cuts, clustering, HMMER panel and others) and stops a calling script if any of them has changed without an explicit acceptance.
# Methods: Candidate discovery (shared helper)
# Inputs:  the reference files listed in `FROZEN` (control declarations, ranked decoys, pilot cuts, SonicParanoid clustering, HMMER panel and its declaration, curated negatives, target registry, replacement panel, table of the supplementary HMMER search)
# Outputs: `05_method_and_rules/INPUT_MANIFEST.tsv` (name, path, SHA-256 checksum, size, date and reason for each input)
# Tools:   Python 3, standard library only, version: see envs/
# <<< deposit header
"""Shared checks on the reference inputs. Import this module and call check_inputs() first.

Why this exists. Two rules of the project were stated but enforced nowhere:

  - every step that uses decoys must use the same corrected decoy pool;
  - parameters are fixed before the calls are made, and each run checks the bytes of its
    inputs.

They are the dangerous kind of unenforced rule, promises about the inputs whose breach is
silent. A script that quietly read an outdated decoy file, or a declaration edited during a
run, would produce a table that looks exactly like a correct one. Nothing detected that; each
script simply opened whatever path its own constant named.

What this does. It writes the SHA-256 checksum of each reference input into a table of inputs
(`INPUT_MANIFEST.tsv`), and on every later run checks that the bytes are unchanged. Every script
that uses controls or decoys calls check_inputs() and either continues on identical bytes or
stops. That makes the single decoy pool a fact about the run rather than a convention: there is
one pool because there is one checksum, and every script that uses it has confirmed it.

Changing an input on purpose is not blocked; it is made explicit. Running with
`--accept-input-change` (or with RL_ACCEPT_INPUT_CHANGE=1) writes the table again, giving the
changed input its new checksum and the date of the change, so a later reader can see that the
inputs moved rather than discovering it from an unexplained shift in the numbers.
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import hashlib
import os
import sys

def _project_root_from_file(start):
    """Find the project root: the nearest folder above this script that holds 05_method_and_rules.

    A script in a copy of the project tree must never fall back to, or hard-code, the path of the
    original tree (chosen on 16 September 2026). Refusing to guess is the point: a hard-coded root
    made a run against the wrong tree look exactly like a correct one, and in its unconditional
    form even setting PROJECT_ROOT did not help.

    os.path.abspath is used rather than os.path.realpath, on purpose. The copied tree holds 2,014
    symbolic links. Walking up from a link finds the tree the script is being run as part of;
    walking up from the link's target would find the tree the file was copied from, which for a
    link into the original tree is exactly the wrong-root error this function prevents.
    """
    d = os.path.dirname(os.path.abspath(start))
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "05_method_and_rules")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("REFUSED: PROJECT_ROOT unset and no project root found above %s" % start)


ROOT = os.environ.get("PROJECT_ROOT") or _project_root_from_file(__file__)
MANIFEST = ROOT + "/05_method_and_rules/INPUT_MANIFEST.tsv"

# The reference inputs, fixed in advance. Every step that uses controls or decoys
# must use exactly these bytes, and each run checks that it does.
FROZEN = {
    "controls_declaration_v3":
        ROOT + "/02_reference_data/RAW_identity-controls_declaration_v3.tsv",
    "decoys_ranked":
        ROOT + "/02_reference_data/decoys_same_species/DECOYS_ranked.tsv",
    "threshold_pilot_pertarget":
        ROOT + "/06_calibration_and_checks/THRESHOLD_PILOT_PERTARGET.tsv",
    # The SonicParanoid clustering is an input in use despite the BACKUP_ in its
    # folder name. rl_arm1_selection.py lists it as the "cov0.75/I1.5 (CURRENT)"
    # candidate beside the five runs under panel/sp2_out/runs, and the jointly
    # chosen settings use it for six of the eight targets. It was left out of this
    # check, the one place designed to catch a silently changed input, only because
    # its folder is named BACKUP. It is not renamed: R59_union_settings__set{A,B}.tsv
    # hold the settings chosen in advance against this exact path, and renaming it
    # for a cosmetic reason would invalidate that choice.
    "sonic_clustering_full02":
        (_TA210_SCRATCH + "/panel/BACKUP_full02_20260819/" "ortholog_groups/flat.ortholog_groups.tsv"),
    # Added on 9 September 2026, when the candidate replacements for lost components were
    # added as targets. Everything below had been outside this check until then, including
    # the panel against which every HMMER score is computed.
    "controls_declaration_current":
        ROOT + "/02_reference_data/RAW_identity-controls_declaration_v7.tsv",
    "panel_hmm":
        ROOT + "/02_reference_data/hmmer_panel/panel_v2.hmm",
    "panel_declaration":
        ROOT + "/02_reference_data/hmmer_panel/PANEL_DECLARATION_v2.tsv",
    "curated_negatives":
        ROOT + "/02_reference_data/CURATED_negatives.tsv",
    "target_registry":
        ROOT + "/05_method_and_rules/TARGET_REGISTRY.tsv",
    "replacement_layer_panel":
        ROOT + "/14_remodeling/panels/REPLACEMENT_LAYER_PANEL.tsv",
    "hmmer_supplement_manifest":
        ROOT + "/03_evidence_by_method/hmmer_genomewide_supplement/SUPPLEMENT_MANIFEST.tsv",
    "threshold_pilot_pertarget__replacement_layers":
        ROOT + "/06_calibration_and_checks/THRESHOLD_PILOT_PERTARGET__replacement_layers.tsv",
    "decoys_ranked__replacement_layers":
        ROOT + "/02_reference_data/decoys_same_species/DECOYS_ranked__replacement_layers.tsv",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_manifest():
    out = {}
    if not os.path.exists(MANIFEST):
        return out
    for line in open(MANIFEST, errors="replace"):
        if line.startswith("##") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[0] != "NAME":
            # The `RECORDED` (date) and `WHY` (reason) columns are read back, not
            # dropped. Reading only path, checksum and size meant that every
            # --accept-input-change wrote a new date and reason on every row,
            # including inputs that had not changed, so adding one input silently
            # dated three others as though they had changed that day. The table
            # is the history of the inputs: a row that did not change must keep
            # the date on which it was first fixed.
            out[f[0]] = dict(path=f[1], sha=f[2], size=int(f[3]),
                             recorded=(f[4] if len(f) > 4 else ""),
                             why=(f[5] if len(f) > 5 else ""))
    return out


def write_manifest(state, why):
    import datetime
    with open(MANIFEST, "w") as fh:
        fh.write("## R26: the frozen reference inputs and their bytes. R20: every "
                 "decoy-consuming layer asserts THESE hashes, so there is one pool "
                 "because there is one hash. Written by rl_invariants.py; a change "
                 "requires --accept-input-change and is dated here.\n")
        fh.write("NAME\tPATH\tSHA256\tSIZE\tRECORDED\tWHY\n")
        for k in sorted(state):
            v = state[k]
            fh.write("%s\t%s\t%s\t%d\t%s\t%s\n"
                     % (k, v["path"], v["sha"], v["size"],
                        v.get("recorded", datetime.date.today().isoformat()),
                        v.get("why", why)))


def check_inputs(quiet=False):
    """Check every fixed input against the table of inputs. Returns the current state of the
    inputs; stops the run (SystemExit) if an input changed without --accept-input-change."""
    accept = ("--accept-input-change" in sys.argv
              or os.environ.get("RL_ACCEPT_INPUT_CHANGE") == "1")
    have = read_manifest()
    cur, drift, missing = {}, [], []
    for name, path in sorted(FROZEN.items()):
        if not os.path.exists(path):
            missing.append("%s (%s)" % (name, path))
            continue
        cur[name] = dict(path=path, sha=sha256(path),
                         size=os.path.getsize(path))
    if missing:
        raise SystemExit("R26: frozen input missing: " + "; ".join(missing))

    if not have:
        write_manifest(cur, "first record")
        if not quiet:
            print("  R20/R26: input manifest created, %d frozen inputs recorded"
                  % len(cur))
        return cur

    # An unchanged input keeps the date and reason it was first given; only a row
    # whose bytes actually changed gets today's date and the new reason.
    for name, v in cur.items():
        prev = have.get(name)
        if prev and prev.get("sha") == v["sha"] and prev.get("path") == v["path"]:
            if prev.get("recorded"):
                v["recorded"] = prev["recorded"]
            if prev.get("why"):
                v["why"] = prev["why"]

    for name, v in cur.items():
        old = have.get(name)
        if old is None:
            drift.append("%s is not in the manifest" % name)
        elif old["sha"] != v["sha"]:
            drift.append("%s CHANGED: %s -> %s (%d -> %d bytes)"
                         % (name, old["sha"][:12], v["sha"][:12],
                            old["size"], v["size"]))
        elif old["path"] != v["path"]:
            drift.append("%s path moved: %s -> %s"
                         % (name, old["path"], v["path"]))
    if drift and not accept:
        raise SystemExit(
            "R26 violation: a frozen reference input changed since the manifest "
            "was written, so calls made before and after are not comparable:\n  "
            + "\n  ".join(drift)
            + "\n  If this change is intended, re-run with --accept-input-change "
              "so the manifest records it; do not edit the manifest by hand.")
    if drift and accept:
        # RL_ACCEPT_WHY carries the reason for the accepted change into the WHY
        # column of the table, ahead of the automatic description of the change,
        # so the table says why the bytes changed and not only that they did.
        why = (os.environ.get("RL_ACCEPT_WHY") or "").strip()
        write_manifest(cur, ((why + " | ") if why else "")
                       + "accepted change: " + "; ".join(drift)[:200])
        print("  R20/R26: input change ACCEPTED and recorded:\n    "
              + "\n    ".join(drift))
        return cur
    if not quiet:
        print("  R20/R26: %d frozen inputs verified unchanged" % len(cur))
    return cur


if __name__ == "__main__":
    check_inputs()
