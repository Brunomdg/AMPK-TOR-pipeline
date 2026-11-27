# >>> deposit header
# Purpose: Give every script the ordered target list of the chosen target set (the eight core targets by
#          default) from the target table, and the set-specific file names that keep the sets apart.
# Methods: several steps (shared helper)
# Inputs:  `05_method_and_rules/TARGET_REGISTRY.tsv`; the environment variables `TARGET_SET`, `PROJECT_ROOT`,
#          `R59_CONTROL_SUBSET`, `PHYLO_SET` and `LOO_BASE`; the companion tables of a non-core set
# Outputs: target lists, suffixes, file names and table rows, returned (printed when run directly)
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
# Ran as: rl_targets.py
"""The single source of the target list. Every script that used to carry the literal
["TOR", "RAPTOR", "RICTOR", "LST8", "AMPKa", "AMPKb", "AMPKg", "SIN1"] reads it from here, and this
module reads it from 05_method_and_rules/TARGET_REGISTRY.tsv.

Why. Adding the twelve candidate replacement targets (layers 3 to 5) to a project whose eight-target
tables were fixed at MD5 fd1803704cf1253d4b126965a043ba9e must not change any existing call. The
safe way is a set: the environment variable TARGET_SET selects which rows of the registry are used,
the default is core8, and the default must reproduce the old list in the old order, which targets()
checks on every call. A run under any other set writes to files and folders whose names carry the
set's suffix (suffix() below), so a run on another set never rewrites the core8 tables.

What is not done here. No target list is written into the code except the eight core targets, and
that list exists only as a check, never as a fallback: a missing or damaged registry stops the run
instead of silently returning the eight.

    targets(set)            ordered target list of the set (default TARGET_SET, otherwise core8)
    suffix(set)             "" for core8, otherwise "__<set>"; appended to every set-specific output
    tagged(path, tag)       base__<tag><suffix>.ext, the one helper for files tagged with a control
                            set (the outputs of the joint choice of settings)
    report_dir(root, tag)   09_per_method_report__<tag><suffix>
    companion(path)         base<suffix>.ext, the set's companion of an input fixed in advance
    pilot_rows(path)        rows of THRESHOLD_PILOT_PERTARGET.tsv plus those of the set's companion
    decoys_ranked_rows(p)   rows of DECOYS_ranked.tsv plus those of the set's companion
    check_track_agreement() PHYLO_SET and LOO_BASE must carry the suffix exactly when the set is
                            not core8
    inverted_from(t), layer(t), controls_declared_in(t)
"""
import csv
import os
import re
import sys

CORE8 = ["TOR", "RAPTOR", "RICTOR", "LST8", "AMPKa", "AMPKb", "AMPKg", "SIN1"]
SETS = ("core8", "replacement_layers")
REGISTRY_REL = "05_method_and_rules/TARGET_REGISTRY.tsv"
# The file-name pattern that the project's layout check (rl_check_layout.py, not part of this
# deposit) refuses. A set name that would produce a forbidden file name is refused here, before
# any file is named.
BAD_NAME = re.compile(r"(^|_)(MASTER|RAW|TABLE1|v[23])(_|\.|$)|20\d{6}")
_ROWS_CACHE = {}


def _project_root_from_file(start):
    """Project root found from the module's own location: the nearest parent folder that holds
    05_method_and_rules.

    A script in the working copy of the project must never fall back to, or hard-code, the path
    of the original project tree (decided on 16 September 2026). A hard-coded root made a run
    against the wrong tree look exactly like a correct one, and where the root was set
    unconditionally even setting PROJECT_ROOT did not help, so the function refuses to guess.

    os.path.abspath is used rather than os.path.realpath, on purpose. The working copy holds 2,014
    symbolic links repointed into it. Walking up from the link finds the tree the script is being
    run in; walking up from the link's target would find the tree the file was copied from, which
    for a link into the original project is exactly the wrong-root error this avoids.
    """
    d = os.path.dirname(os.path.abspath(start))
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "05_method_and_rules")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("REFUSED: PROJECT_ROOT unset and no project root found above %s" % start)


def project_root():
    return os.environ.get("PROJECT_ROOT") or _project_root_from_file(__file__)


def set_name(explicit=None):
    s = (explicit or os.environ.get("TARGET_SET") or "core8").strip()
    if s not in SETS:
        raise SystemExit("REFUSING: TARGET_SET=%r is not one of %s" % (s, ", ".join(SETS)))
    if BAD_NAME.search(s):
        raise SystemExit("REFUSING: TARGET_SET=%r matches the forbidden-name pattern "
                         "rl_check_layout.py enforces" % s)
    if s != "core8" and os.environ.get("R59_CONTROL_SUBSET"):
        raise SystemExit(
            "REFUSING: TARGET_SET=%s together with R59_CONTROL_SUBSET=%r. The subset "
            "files hold core8 accessions only, so a subset solve of a non-core set "
            "would be solved against an empty control set." % (s, os.environ.get("R59_CONTROL_SUBSET")))
    return s


def registry_path(root=None):
    return os.path.join(root or project_root(), REGISTRY_REL)


def rows(root=None):
    p = registry_path(root)
    if p in _ROWS_CACHE:
        return _ROWS_CACHE[p]
    if not os.path.exists(p):
        raise SystemExit("REFUSING: %s is missing. Write it with "
                         "01_scripts/rl_target_registry.py; there is no fallback list." % p)
    out = []
    seen = set()
    with open(p, errors="replace") as fh:
        rd = csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t")
        need = {"TARGET", "TARGET_SET", "LAYER", "ORDER"}
        if not need <= set(rd.fieldnames or []):
            raise SystemExit("REFUSING: %s lacks columns %s"
                             % (p, sorted(need - set(rd.fieldnames or []))))
        for r in rd:
            t = (r.get("TARGET") or "").strip()
            s = (r.get("TARGET_SET") or "").strip()
            if not t:
                continue
            if t in seen:
                raise SystemExit("REFUSING: target %s appears twice in %s" % (t, p))
            seen.add(t)
            if s not in SETS:
                raise SystemExit("REFUSING: %s has TARGET_SET=%r for %s" % (p, s, t))
            try:
                int(r.get("ORDER") or "")
            except ValueError:
                raise SystemExit("REFUSING: non-integer ORDER for %s in %s" % (t, p))
            out.append({k: (v or "").strip() for k, v in r.items()})
    _ROWS_CACHE[p] = out
    return out


def targets(explicit=None, root=None):
    s = set_name(explicit)
    keep = sorted((r for r in rows(root) if r["TARGET_SET"] == s),
                  key=lambda r: int(r["ORDER"]))
    out = [r["TARGET"] for r in keep]
    if not out:
        raise SystemExit("REFUSING: TARGET_SET=%s selects no target in %s"
                         % (s, registry_path(root)))
    if s == "core8" and out != CORE8:
        raise SystemExit("REFUSING: the registry's core8 rows %s differ from the frozen "
                         "eight %s; the eight-target deliverable would no longer be "
                         "reproducible" % (out, CORE8))
    return out


def suffix(explicit=None):
    s = set_name(explicit)
    return "" if s == "core8" else "__" + s


def tagged(path, tag, explicit=None):
    """base__<tag><suffix>.ext. The tag names the control set (R59_TAG) and the suffix the
    target set; the two are independent and are never folded into one token."""
    if not tag:
        raise SystemExit(
            "REFUSING: R59_TAG is not set. R59 outputs are per control set "
            "(setA = the full binding declaration, setB = the clean literature "
            "subset); there is no untagged default because a default IS a fallback.")
    base, dot, ext = path.rpartition(".")
    if dot:
        return "%s__%s%s.%s" % (base, tag, suffix(explicit), ext)
    return "%s__%s%s" % (path, tag, suffix(explicit))


def report_dir(root, tag, explicit=None):
    if not tag:
        raise SystemExit("REFUSING: R59_TAG is not set; the report directory is per control set")
    return os.path.join(root or project_root(), "09_per_method_report__" + tag + suffix(explicit))


def companion(path, explicit=None):
    """The set-specific companion of an input fixed in advance: base<suffix>.ext. For core8
    it is the fixed file itself, so a core8 run reads exactly those bytes."""
    sfx = suffix(explicit)
    if not sfx:
        return path
    base, dot, ext = path.rpartition(".")
    return "%s%s.%s" % (base, sfx, ext) if dot else path + sfx


def _tsv_rows(path):
    with open(path, errors="replace") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"):
            yield r


def pilot_rows(path, explicit=None):
    """Rows of THRESHOLD_PILOT_PERTARGET: the fixed file, then the set's companion. A
    companion row for a target outside the set stops the run: the companion exists only to
    carry the new targets' measurements, never to replace a row of the fixed file."""
    yield from _tsv_rows(path)
    s = set_name(explicit)
    if s == "core8":
        return
    comp = companion(path, s)
    if not os.path.exists(comp):
        raise SystemExit("REFUSING: TARGET_SET=%s needs the pilot companion %s and it "
                         "does not exist" % (s, comp))
    allowed = set(targets(s))
    for r in _tsv_rows(comp):
        t = (r.get("FOR_TARGET") or "").strip()
        if t not in allowed:
            raise SystemExit("REFUSING: %s carries a row for %s, outside TARGET_SET=%s"
                             % (comp, t, s))
        yield r


def decoys_ranked_rows(path, explicit=None):
    """Rows of DECOYS_ranked: the fixed file, then the set's companion (checked on its TARGET column)."""
    yield from _tsv_rows(path)
    s = set_name(explicit)
    if s == "core8":
        return
    comp = companion(path, s)
    if not os.path.exists(comp):
        raise SystemExit("REFUSING: TARGET_SET=%s needs the decoy companion %s and it "
                         "does not exist" % (s, comp))
    allowed = set(targets(s))
    for r in _tsv_rows(comp):
        t = (r.get("TARGET") or "").strip()
        if t not in allowed:
            raise SystemExit("REFUSING: %s carries a row for %s, outside TARGET_SET=%s"
                             % (comp, t, s))
        yield r


def check_track_agreement(explicit=None):
    """The gene-tree output folder (PHYLO_SET) and the leave-one-out folder (LOO_BASE), when
    set, must carry the set suffix exactly when the set is not core8. A core8 run writing into
    a suffixed folder, or a run on another set writing into the core folder, is the collision
    this mechanism exists to prevent, so it stops the run instead of printing a warning."""
    s = set_name(explicit)
    sfx = suffix(s)
    for var in ("PHYLO_SET", "LOO_BASE"):
        v = (os.environ.get(var) or "").strip()
        if not v:
            continue
        body = v[:-4] if v.endswith("_loo") else v
        if sfx:
            if not body.endswith(sfx) and not re.search(re.escape(sfx) + r"(_loo)?(/|$)", v):
                raise SystemExit("REFUSING: TARGET_SET=%s but %s=%s does not carry the "
                                 "suffix %s; a set run must write into its own track"
                                 % (s, var, v, sfx))
        elif "__" in v:
            raise SystemExit("REFUSING: TARGET_SET is core8 (or unset) but %s=%s names a "
                             "set-suffixed track; set TARGET_SET to match" % (var, v))
    return s


def _row(t, root=None):
    for r in rows(root):
        if r["TARGET"] == t:
            return r
    raise SystemExit("REFUSING: target %s is not in %s" % (t, registry_path(root)))


def inverted_from(t, root=None):
    return _row(t, root).get("INVERTED_FROM", "")


def layer(t, root=None):
    return _row(t, root).get("LAYER", "")


def controls_declared_in(t, root=None):
    return _row(t, root).get("CONTROLS_DECLARED_IN", "")


def label(explicit=None):
    return "TARGET_SET=%s" % set_name(explicit)


if __name__ == "__main__":
    s = set_name()
    print("%s -> %s" % (label(s), " ".join(targets(s))))
    print("suffix=%r" % suffix(s))
    check_track_agreement(s)
    sys.exit(0)
