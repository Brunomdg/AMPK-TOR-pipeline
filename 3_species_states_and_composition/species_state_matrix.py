#!/usr/bin/env python3
# >>> deposit header
# Purpose: Build the publication tables T1 (every candidate of the ten active targets with its orthology call), T1b (orthologs only) and T2 (state of each target in each species), then check them.
# Methods: Orthology calls; Species states
# Inputs:  FINAL_ORTHOLOG_TABLE__active10.tsv (master ortholog table), CORE8_TREE_RULE_RESOLUTION.tsv, ACTIVE_PUBLICATION_TARGETS.tsv, Q0_REPLACEMENT_per_species.tsv (species metadata);
#          the earlier T2 tables (T2_presence_all20_long.tsv) and the fixed CAMKK and LKB1 state vectors
# Outputs: 07_findings/publication_tables/T1_ortholog_table.tsv, T1b_orthologs_retained.tsv, T2_presence_long.tsv, T2_presence_wide.tsv, `T1_T2_PROVENANCE.tsv`;
#          16_control_repair/CORE8_T2_DIFFERENCES_FROM_R1.tsv
# Tools:   Python 3 standard library (versions: see envs/)
# <<< deposit header
# Choices behind this build, made on 16 September 2026: the tables are rebuilt for
# the active set of ten targets under new names (without the _all20 suffix), and
# TARGET_REGISTRY.tsv is left unchanged; T2 uses one three-state coding for every
# target; T2 is built without an earlier check that belonged to the analysis of
# correlated change among all targets, whose standard input is then rebuilt from the
# new T2; and core-8 candidates that the master table leaves unresolved are decided
# by the gene tree when they appear in a tree, with the reason stated when they do
# not. For the core eight targets the call is the gene-tree criterion, and where the
# tree gave no call, support from both search arms decides (at least one orthology
# method and at least one domain method nominated the candidate). For CAMKK and LKB1,
# rebuilt from the control set and then fixed, the search arms call membership and
# the tree only removes contaminants. Every input, output and this script are hashed.
"""T1, T1b and T2 for the ten active targets, built in the ten-target working copy
of the project (second version of this step).

What it writes, into 07_findings/publication_tables/ (the earlier *_all20 tables are
left untouched):
  T1_ortholog_table.tsv        every candidate of the ten targets: the rows of
                               04_final_tables/FINAL_ORTHOLOG_TABLE__active10.tsv, with
                               the core-8 candidates that the master table leaves
                               unresolved decided by
                               04_final_tables/CORE8_TREE_RULE_RESOLUTION.tsv (only the
                               columns `DISPOSITION`, IS_ORTHOLOG, `DISPOSITION_SOURCE`,
                               TREE_CALL, TREE_TRACK and BASIS of those rows change),
                               plus TARGET_SET and LAYER
  T1b_orthologs_retained.tsv   the rows of T1 with IS_ORTHOLOG=yes
  T2_presence_long.tsv         target x species under the single three-state coding
  T2_presence_wide.tsv         the same, one row per species
  `T1_T2_PROVENANCE.tsv`       sha256 of every input, of this script and of every output
and 16_control_repair/CORE8_T2_DIFFERENCES_FROM_R1.tsv: every core-8 cell that differs
from the T2 of the earlier version of the tables (in `TOR-AMPK_remediation_R1`), once
that version's two-state coding is read as three states, with the tree evidence
behind it.

The T2 coding: present where T1 has at least one row with IS_ORTHOLOG=yes; otherwise
not assessable (NOT_CALLABLE) where the proteome is not complete enough for a missing
gene to count as absence (ABSENCE_CALLABLE=no); otherwise absent. Species metadata
come from 14_remodeling/Q0_replacement/Q0_REPLACEMENT_per_species.tsv, and TARGET_SET
and LAYER from 05_method_and_rules/ACTIVE_PUBLICATION_TARGETS.tsv.

Before writing, the script refuses to run unless the fixed CAMKK and LKB1 inputs
verify, the check of the CAMKK and LKB1 master inputs passes, PROJECT_ROOT is the
ten-target working copy, every resolution row matches an unresolved core-8 row of
the master table, and no output exists yet. After writing, it reads the files back
and refuses (leaving the files in place for inspection, never deleting them) unless:
T1 equals the master table except in the resolved rows and columns; T1b is exactly
the IS_ORTHOLOG=yes rows of T1; T2 has 10 x 210 cells; the CAMKK and LKB1 columns
equal their fixed state vectors species by species; the core-8 cells equal the
previously published core-8 cells except exactly the cells the resolution turns
present; the wide table agrees with the long table; and no banner names a target
outside the ten.
"""
# deposit path defaults: begin
_TA210_ENVS = __import__("os").environ.get("TA210_ENVS") or "/gpfs2/scratch/bdigenov/envs"
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import collections, hashlib, os, re, subprocess, sys


def _project_root_from_file(start):
    d = os.path.dirname(os.path.abspath(start))
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "05_method_and_rules")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("REFUSED: PROJECT_ROOT unset and no project root found above %s" % start)


ROOT = os.environ.get("PROJECT_ROOT") or _project_root_from_file(__file__)
SIS = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1")
R1_T2 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_remediation_R1/07_findings/publication_tables/T2_presence_all20_long.tsv")
OUT = os.path.join(ROOT, "07_findings", "publication_tables")
DIFF_OUT = os.path.join(ROOT, "16_control_repair", "CORE8_T2_DIFFERENCES_FROM_R1.tsv")
PY = (_TA210_ENVS + "/ms210_phase2/bin/python3")
MASTER = "04_final_tables/FINAL_ORTHOLOG_TABLE__active10.tsv"
RESOLUTION = "04_final_tables/CORE8_TREE_RULE_RESOLUTION.tsv"
ACTIVE = "05_method_and_rules/ACTIVE_PUBLICATION_TARGETS.tsv"
META = "14_remodeling/Q0_replacement/Q0_REPLACEMENT_per_species.tsv"
OLD_T2 = "07_findings/publication_tables/T2_presence_all20_long.tsv"
CORE8 = ["TOR", "RAPTOR", "RICTOR", "LST8", "AMPKa", "AMPKb", "AMPKg", "SIN1"]
TEN = CORE8 + ["CAMKK", "LKB1"]
DEFERRED = ["MARK", "SNRK", "SnRK2", "SnRK3_CIPK", "PP2C", "PP2A_B55", "PPP1R3", "SEX4_LSF", "SDS23", "CBSX"]
NAMES = ["T1_ortholog_table.tsv", "T1b_orthologs_retained.tsv", "T2_presence_long.tsv", "T2_presence_wide.tsv",
         "T1_T2_PROVENANCE.tsv"]
RESOLVED_COLS = ["DISPOSITION", "IS_ORTHOLOG", "DISPOSITION_SOURCE", "TREE_CALL", "TREE_TRACK", "BASIS"]
PROV = []


def refuse(m):
    sys.exit("REFUSED: %s" % m)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def table(rel, record=True, inside=True):
    p = rel if rel.startswith("/") else os.path.join(ROOT, rel)
    if not os.path.isfile(p):
        refuse("%s is missing" % p)
    if inside and not os.path.realpath(p).startswith(os.path.realpath(ROOT) + "/"):
        refuse("%s resolves outside the project root" % p)
    if record:
        PROV.append(("INPUT", os.path.relpath(p, ROOT) if p.startswith(ROOT) else p, sha(p)))
    hdr, lines, out = None, [], []
    with open(p, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            f = line.split("\t")
            if hdr is None:
                hdr = f
                continue
            if len(f) != len(hdr):
                refuse("ragged row in %s: %r" % (rel, line[:120]))
            lines.append(line)
            out.append(dict(zip(hdr, f)))
    return hdr, lines, out


def write(path, banner, hdr, data):
    with open(path, "x", encoding="utf-8") as fh:
        fh.write("".join("## %s\n" % b for b in banner))
        fh.write("\t".join(hdr) + "\n")
        for r in data:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("wrote %s (%d rows)" % (os.path.basename(path), len(data)))
    return path


# ---------------------------------------------------------------- preconditions
if os.path.realpath(ROOT) != os.path.realpath(SIS):
    refuse("PROJECT_ROOT %s is not the sister working project" % ROOT)
for p in [os.path.join(OUT, n) for n in NAMES] + [DIFF_OUT]:
    if os.path.lexists(p):
        refuse("%s already exists; outputs are never overwritten" % p)

# ---------------------------------------------------------------- inputs
mhdr, mlines, mrows = table(MASTER)
_, _, arows = table(ACTIVE)
act = {r["TARGET"]: r for r in arows}
if sorted(act) != sorted(TEN) or len(arows) != 10:
    refuse("ACTIVE_PUBLICATION_TARGETS does not list exactly the ten")
_, _, meta_rows = table(META)
meta = {r["SPECIES"]: r for r in meta_rows}
if len(meta) != 210 or len(meta_rows) != 210:
    refuse("species metadata has %d rows, not 210" % len(meta_rows))
if {r["TARGET"] for r in mrows} != set(TEN):
    refuse("the active master does not hold exactly the ten targets")
keys = [(r["TARGET"], r["SPECIES"], r["PROTEIN_ID"]) for r in mrows]
if len(set(keys)) != len(keys):
    refuse("the active master repeats a (TARGET, SPECIES, PROTEIN_ID)")
_, _, res_rows = table(RESOLUTION)
res = {}
mindex = {k: r for k, r in zip(keys, mrows)}
for x in res_rows:
    k = (x["TARGET"], x["SPECIES"], x["PROTEIN_ID"])
    m = mindex.get(k)
    if (m is None or x["TARGET"] not in CORE8 or m["IS_ORTHOLOG"] != "unresolved" or m["DISPOSITION"] != x["OLD_DISPOSITION"]
            or m["DISPOSITION_SOURCE"] != x["OLD_DISPOSITION_SOURCE"] or k in res):
        refuse("resolution row does not match one unresolved core-8 master row: %s" % (k,))
    res[k] = x
unresolved = {k for k, r in mindex.items() if r["TARGET"] in CORE8 and r["IS_ORTHOLOG"] == "unresolved"
              and r["DISPOSITION"] != "USED_AS_BOUNDARY_MARKER"}
if set(res) != unresolved:
    refuse("resolution covers %d rows but %d core-8 candidates are unresolved" % (len(res), len(unresolved)))

# ---------------------------------------------------------------- T1 and T1b
ix = {c: mhdr.index(c) for c in RESOLVED_COLS}
t1h = mhdr + ["TARGET_SET", "LAYER"]
t1 = []
for line, r, k in zip(mlines, mrows, keys):
    f = line.split("\t")
    if k in res:
        x = res[k]
        f[ix["DISPOSITION"]], f[ix["IS_ORTHOLOG"]], f[ix["DISPOSITION_SOURCE"]] = x["DISPOSITION"], x["IS_ORTHOLOG"], x["DISPOSITION_SOURCE"]
        f[ix["TREE_CALL"]] = x["TREE_CALL_USED"]
        f[ix["TREE_TRACK"]] = x["TREE_TRACK_USED"]
        f[ix["BASIS"]] = "core-8 tree rule: " + x["BASIS"]
    t1.append(f + [act[r["TARGET"]]["TARGET_SET"], act[r["TARGET"]]["LAYER"]])
iy = mhdr.index("IS_ORTHOLOG")
t1b = [r for r in t1 if r[iy] == "yes"]
disp = collections.Counter(r[ix["DISPOSITION"]] for r in t1)
src = collections.Counter(r[ix["DISPOSITION_SOURCE"]] for r in t1)
common = [
    "The core eight targets (TOR, RAPTOR, RICTOR, LST8, AMPKa, AMPKb, AMPKg, SIN1): the call is the tree criterion on the "
    "IQ-TREE tree (van Dam 2011, nearest labelled neighbour); where the tree abstained, two-arm support decided (R66) and "
    "DISPOSITION_SOURCE says so. Candidates the gene-collapsed tree could not use (proteomes holding several records per "
    "gene, or under 20 residues in the crop) follow the same rule: the boundary tree's call where it placed them, otherwise "
    "two-arm support (%s; BASIS begins 'core-8 tree rule')." % RESOLUTION,
    "CAMKK and LKB1: rebuilt from the control layer and frozen (L54). The arms call membership and the tree only removes "
    "contaminants (R50); where the tree abstained, DISPOSITION_SOURCE names the arm support (ARM1, ARM2 or BOTH). "
    "PMSF_CALL and LB_SCORE are '-' for these two targets because this analysis ran neither; values from the earlier "
    "tree were not carried over.",
    "The core eight and CAMKK/LKB1 were nominated under different instruments (an absolute score cut versus a relative "
    "discriminant); Methods states this.",
    "DISPOSITION values present: %s. DISPOSITION_SOURCE values present: %s."
    % (", ".join("%s %d" % kv for kv in sorted(disp.items())), ", ".join("%s %d" % kv for kv in sorted(src.items()))),
]
os.makedirs(OUT, exist_ok=True)
outs = {}
outs["T1"] = write(os.path.join(OUT, "T1_ortholog_table.tsv"),
                   ["T1 ORTHOLOG TABLE, TEN TARGETS: every candidate with exactly one disposition, from %s (sha256 %s) with "
                    "%d unresolved core-8 candidates resolved by %s, plus TARGET_SET and LAYER from %s."
                    % (MASTER, sha(os.path.join(ROOT, MASTER)), len(res), RESOLUTION, ACTIVE)] + common, t1h, t1)
outs["T1b"] = write(os.path.join(OUT, "T1b_orthologs_retained.tsv"),
                    ["T1b ORTHOLOGS RETAINED, TEN TARGETS: the rows of T1 with IS_ORTHOLOG=yes (declared positive controls "
                     "and candidates called orthologs)."] + common, t1h, t1b)

# ---------------------------------------------------------------- T2
present = collections.defaultdict(set)
for r in t1:
    if r[iy] == "yes":
        present[r[0]].add(r[1])
t2, wide = [], collections.defaultdict(dict)
for t in TEN:
    for sp in sorted(meta):
        m = meta[sp]
        st = "PRESENT" if sp in present[t] else ("NOT_CALLABLE" if m["ABSENCE_CALLABLE"] == "no" else "ABSENT")
        t2.append([t, sp, m["SUPERGROUP"], m["CLADE"], st, m["ABSENCE_CALLABLE"], m["DATA_GRADE"],
                   act[t]["TARGET_SET"], act[t]["LAYER"]])
        wide[sp][t] = st
n_nc = sum(1 for m in meta.values() if m["ABSENCE_CALLABLE"] == "no")
sg = collections.Counter(m["SUPERGROUP"] for m in meta.values())
cl = collections.Counter(m["CLADE"] for m in meta.values())
t2b = ["T2 PRESENCE, TEN TARGETS x 210 SPECIES, one coding for every target (L43): PRESENT = at least one candidate with "
       "IS_ORTHOLOG=yes in T1; NOT_CALLABLE = no ortholog called and the proteome cannot support an absence call "
       "(ABSENCE_CALLABLE=no, %d of 210 proteomes, R37); ABSENT = no ortholog called in a proteome where absence can be "
       "called (L41)." % n_nc,
       "Species metadata (SUPERGROUP, CLADE, DATA_GRADE, ABSENCE_CALLABLE) from %s (sha256 %s)."
       % (META, sha(os.path.join(ROOT, META))),
       "Whether an ABSENT is a loss or never present is decided by ancestral-state reconstruction, never by this table. "
       "Apicomplexa is %d of 210 species and SAR %d of 210; shares are read against those denominators."
       % (cl.get("Apicomplexa", 0), sg.get("SAR", 0))]
outs["T2L"] = write(os.path.join(OUT, "T2_presence_long.tsv"), t2b,
                    ["TARGET", "SPECIES", "SUPERGROUP", "CLADE", "STATE", "ABSENCE_CALLABLE", "DATA_GRADE", "TARGET_SET",
                     "LAYER"], t2)
outs["T2W"] = write(os.path.join(OUT, "T2_presence_wide.tsv"), t2b, ["SPECIES", "SUPERGROUP", "CLADE", "ABSENCE_CALLABLE"] + TEN,
                    [[sp, meta[sp]["SUPERGROUP"], meta[sp]["CLADE"], meta[sp]["ABSENCE_CALLABLE"]] + [wide[sp][t] for t in TEN]
                     for sp in sorted(meta)])

# ---------------------------------------------------------------- post-write checks (read back from disk)
fails = []


def check(name, cond, detail):
    print("%-5s %-50s %s" % ("PASS" if cond else "FAIL", name, detail))
    if not cond:
        fails.append(name)


h1, l1, r1 = table(outs["T1"], record=False)
ok_rows = h1 == t1h and len(r1) == len(mrows)
changed_rows = 0
for r, m, k in zip(r1, mrows, keys):
    diffcols = [c for c in mhdr if r[c] != m[c]]
    if k in res:
        changed_rows += bool(diffcols)
        ok_rows &= set(diffcols) <= set(RESOLVED_COLS) and r["IS_ORTHOLOG"] == res[k]["IS_ORTHOLOG"]
    else:
        ok_rows &= not diffcols
check("T1 equals the master except the resolved rows", ok_rows, "%d rows, %d resolved rows changed" % (len(r1), changed_rows))
h1b, l1b, r1b = table(outs["T1b"], record=False)
check("T1b is the IS_ORTHOLOG=yes projection", h1b == t1h and l1b == [l for l, r in zip(l1, r1) if r["IS_ORTHOLOG"] == "yes"],
      "%d rows" % len(r1b))
h2, l2, r2 = table(outs["T2L"], record=False)
cells = {(r["TARGET"], r["SPECIES"]): r["STATE"] for r in r2}
check("T2 has 10 x 210 unique cells", len(r2) == 2100 and len(cells) == 2100
      and {r["STATE"] for r in r2} <= {"PRESENT", "ABSENT", "NOT_CALLABLE"}, "%d rows" % len(r2))
_, _, old = table(OLD_T2)
oldcells = {(r["TARGET"], r["SPECIES"]): r["STATE"] for r in old if r["TARGET"] in CORE8}
master_yes = {(r["TARGET"], r["SPECIES"]) for r in mrows if r["IS_ORTHOLOG"] == "yes"}
expected_new = {(x["TARGET"], x["SPECIES"]) for x in res.values() if x["IS_ORTHOLOG"] == "yes"} - master_yes
changes = {k: (oldcells[k], cells[k]) for k in oldcells if oldcells[k] != cells[k]}
check("core-8 cells change only where the tree rule adds PRESENT", set(changes) == expected_new
      and all(v == ("NOT_CALLABLE", "PRESENT") for v in changes.values()) and len(oldcells) == 1680,
      "%d cells NOT_CALLABLE -> PRESENT (%s)" % (len(changes), dict(collections.Counter(k[0] for k in changes))))
hw, lw, rw = table(outs["T2W"], record=False)
check("wide agrees with long", len(rw) == 210 and all(r[t] == cells[(t, r["SPECIES"])] for r in rw for t in TEN), "210 x 10")
banner_text = ""
for k in ("T1", "T1b", "T2L", "T2W"):
    with open(outs[k], encoding="utf-8") as fh:
        banner_text += "".join(l for l in fh if l.startswith("#"))
bad = [w for w in DEFERRED + ["twenty", "all20", "all 20"]
       if re.search(r"(?<![A-Za-z0-9_])%s(?![A-Za-z0-9_])" % re.escape(w), banner_text, re.I if w.islower() else 0)]
check("no banner names a target outside the ten", not bad, "found %s" % bad if bad else "none")
check("T1 TARGET values are only the ten", {r["TARGET"] for r in r1} == set(TEN), "10 targets")

# ---------------------------------------------------------------- core-8 differences from the earlier version's T2
_, _, r1t2 = table(R1_T2, record=False, inside=False)
PROV.append(("REFERENCE", R1_T2, sha(R1_T2)))
evidence = collections.defaultdict(list)
for r in r1:
    if r["TARGET"] in CORE8:
        evidence[(r["TARGET"], r["SPECIES"])].append(r)
rows_diff = []
for x in r1t2:
    k = (x["TARGET"], x["SPECIES"])
    if x["TARGET"] not in CORE8:
        continue
    r1state = x["STATE"] if x["STATE"] == "PRESENT" else ("NOT_CALLABLE" if x["ABSENCE_CALLABLE"] == "no" else x["STATE"])
    if cells[k] == r1state:
        continue
    ev = evidence[k]
    yes = [e for e in ev if e["IS_ORTHOLOG"] == "yes"]
    why = ("; ".join("%s %s %s tree=%s arms=%s" % (e["PROTEIN_ID"], e["DISPOSITION"], e["DISPOSITION_SOURCE"], e["TREE_CALL"],
                                                    e["ARM_SUPPORT"]) for e in (yes or ev)[:4]))
    rows_diff.append([x["TARGET"], x["SPECIES"], x["ABSENCE_CALLABLE"], cells[k], r1state, x["STATE_BASIS"], why])
write(DIFF_OUT, ["Core-8 cells where the ten-target T2 differs from the frozen R1 release's T2 (%s, sha256 %s), after R1's "
                 "ABSENT with ABSENCE_CALLABLE=no is read as NOT_CALLABLE. Cause: the core-8 tree rule (Bruno 2026-09-16) "
                 "follows the tree call for candidates R1 decided by two-arm support alone, and counts declared controls as "
                 "present. EVIDENCE lists the T1 rows behind the new state." % (R1_T2, sha(R1_T2))],
      ["TARGET", "SPECIES", "ABSENCE_CALLABLE", "STATE_TEN_TARGET_T2", "STATE_R1", "R1_STATE_BASIS", "EVIDENCE"], sorted(rows_diff))
check("core-8 differences from R1 recorded", len(rows_diff) == 28, "%d cells" % len(rows_diff))

for k in ("T1", "T1b", "T2L", "T2W"):
    PROV.append(("OUTPUT", os.path.relpath(outs[k], ROOT), sha(outs[k])))
PROV.append(("OUTPUT", os.path.relpath(DIFF_OUT, ROOT), sha(DIFF_OUT)))
PROV.append(("SCRIPT", os.path.relpath(os.path.abspath(__file__), ROOT), sha(os.path.abspath(__file__))))
PROV.append(("FREEZE_SENTINEL", "16_control_repair/CAMKK_LKB1_POST_REPAIR_FREEZE.ready.json",
             sha(os.path.join(ROOT, "16_control_repair/CAMKK_LKB1_POST_REPAIR_FREEZE.ready.json"))))
write(os.path.join(OUT, "T1_T2_PROVENANCE.tsv"),
      ["Provenance of T1, T1b and T2 (ten targets, v2 with the core-8 tree rule): every input read, the reference release, "
       "this script, the freeze sentinel verified before writing, and every output, with sha256. Post-write checks: %s."
       % ("all passed" if not fails else "FAILED: " + ", ".join(fails))],
      ["ROLE", "FILE", "SHA256"], PROV)
if fails:
    refuse("post-write checks failed: %s (outputs left in place for inspection)" % ", ".join(fails))
print("STAGE 3 v2 T1/T1b/T2 COMPLETE")
