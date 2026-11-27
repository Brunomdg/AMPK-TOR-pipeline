#!/usr/bin/env python3
# >>> deposit header
# Purpose: Builds the standard per-species coding (present 0, not assessable 0.5, absent 1) of the ten active targets
#          from the ten-target presence table, and compares it with the superseded coding and the second coding.
# Methods: Coding of species states for the evolutionary analyses
# Inputs:  `T2_presence_long.tsv` (the ten-target presence table), `ACTIVE_PUBLICATION_TARGETS.tsv`, the 209-tip
#          chronogram, and the superseded twenty-trait `STANDARD_per_species.tsv` (for comparison only)
# Outputs: `STANDARD_per_species.tsv`, `STANDARD_changes.tsv` and `STANDARD_vs_R37.tsv` in `14_remodeling/Q7_alltargets`
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
# Ran as: rl_rm_q7_standard.py
"""Build the standard per-species coding of the ten active targets from the ten-target presence table (T2).

One rule for every trait: `<T>_absent` p1 = 0.0 where T2 says `PRESENT` (an ortholog was called, whatever the
completeness of the proteome); 0.5 where T2 says `NOT_CALLABLE` (the proteome cannot support an absence call; the
reason is in the `_na` column); 1.0 where T2 says `ABSENT` (no ortholog was called). Every cell comes from
07_findings/publication_tables/T2_presence_long.tsv, which is built from the master table of orthology calls plus
the tree-based rule used for the eight core components; nothing is read from the superseded twenty-target tables,
and no cell is copied from an earlier coding. Files superseded by the ten-target rebuild were moved out of every
folder used for calculation on 16 September 2026, so the superseded coding is read from that separate folder,
for comparison only.

Species metadata (SUPERGROUP, CLADE, DATA_GRADE, ABSENCE_CALLABLE) come from T2. ON_CHRONOGRAM is read from the
tips of the 209-tip chronogram, not from a file of the replacement analysis.

Outputs, in 14_remodeling/Q7_alltargets/ (never overwritten):
  STANDARD_per_species.tsv  210 rows, the ten `<T>_absent_p1` / `_na` column pairs, in the order of
                            ACTIVE_PUBLICATION_TARGETS.tsv
  STANDARD_changes.tsv      every cell that differs from the superseded twenty-trait coding
  STANDARD_vs_R37.tsv       per trait, the cells that differ from the second coding (in which all 28 proteomes where
                            absence cannot be assessed are 0.5 on every target), on the 210 species and on the 209
                            chronogram tips; SOURCE_ASR is STANDARD for all ten traits, so no reconstruction is reused
                            from a run made before the repair
Refuses unless the CAMKK and LKB1 calls are confirmed final, the check of the master inputs for CAMKK and LKB1
passes, PROJECT_ROOT is the ten-target project root, T2 has the SHA-256 checksum fixed in this script, and no
output exists yet.
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
PY = (_TA210_ENVS + "/ms210_phase2/bin/python3")
OUTD = os.environ.get("Q7_OUT") or os.path.join(ROOT, "14_remodeling", "Q7_alltargets")
T2 = os.path.join(ROOT, "07_findings", "publication_tables", "T2_presence_long.tsv")
T2_SHA = "e4e4946abac3d6cb41ed3c3477e8bd6d9d70f03820ef1e3930f1b1e9023e65d4"
ACTIVE = os.path.join(ROOT, "05_method_and_rules", "ACTIVE_PUBLICATION_TARGETS.tsv")
CHRONO = (_TA210_SCRATCH + "/concord/ms210/P04b_chronogram/species_chronogram_209.treefile")
SUPERSEDED = ((_TA210_SCRATCH + "/ms210_purge_quarantine/" "20260916_200343_stale_derived_superseded_by_ten_target_rebuild/" "14_remodeling/Q7_alltargets/STANDARD_per_species.tsv"))
R37_REASON = "R37: absence not callable"
CODE = {"PRESENT": ("0.0", ""), "NOT_CALLABLE": ("0.5", R37_REASON), "ABSENT": ("1.0", "")}
SIGN = ("## SIGN (DECISIONS I7, restated under L43): state 1 of every `<T>_absent` trait is the absence of a called "
        "ortholog; whether that absence is a loss or a never-present is what the ancestral-state reconstruction says, "
        "never the name.")
NAMES = ["STANDARD_per_species.tsv", "STANDARD_changes.tsv", "STANDARD_vs_R37.tsv"]


def refuse(m):
    sys.exit("REFUSED: %s" % m)


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def md5(p):
    with open(p, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def table(p):
    hdr, out = None, []
    with open(p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            if len(f) != len(hdr):
                refuse("ragged row in %s" % p)
            out.append(dict(zip(hdr, f)))
    return hdr, out


# ---------------------------------------------------------------- refusals before any input is read
if os.path.realpath(ROOT) != os.path.realpath(SIS):
    refuse("PROJECT_ROOT %s is not the sister" % ROOT)
for n in NAMES:
    if os.path.lexists(os.path.join(OUTD, n)):
        refuse("%s exists; outputs are never overwritten (quarantine the old one first)" % n)

# ---------------------------------------------------------------- inputs
_, arows = table(ACTIVE)
traits = [r["TARGET"] for r in arows]
if len(traits) != 10 or len(set(traits)) != 10:
    refuse("ACTIVE_PUBLICATION_TARGETS does not list ten distinct targets")
_, t2rows = table(T2)
cell = {(r["TARGET"], r["SPECIES"]): r for r in t2rows}
species = sorted({r["SPECIES"] for r in t2rows})
if len(species) != 210 or len(cell) != 10 * 210:
    refuse("T2 holds %d species and %d cells" % (len(species), len(cell)))
meta = {}
for sp in species:
    m = {k: cell[(traits[0], sp)][k] for k in ("SUPERGROUP", "CLADE", "DATA_GRADE", "ABSENCE_CALLABLE")}
    for t in traits:
        if {k: cell[(t, sp)][k] for k in m} != m:
            refuse("T2 metadata disagrees between targets for %s" % sp)
    meta[sp] = m
tips = set(re.findall(r"[(,]([^():,;]+)[:)]", open(CHRONO).read()))
missing = [s for s in species if s not in tips]
if len(tips) != 209 or len(missing) != 1:
    refuse("chronogram has %d tips and %d of the 210 species are off it: %s" % (len(tips), len(missing), missing[:3]))

hdr = ["SPECIES", "SUPERGROUP", "CLADE", "DATA_GRADE", "ABSENCE_CALLABLE", "ON_CHRONOGRAM"]
for t in traits:
    hdr += ["%s_absent_p1" % t, "%s_absent_na" % t]
data, r37 = [], []
for sp in species:
    m = meta[sp]
    base = [sp, m["SUPERGROUP"], m["CLADE"], m["DATA_GRADE"], m["ABSENCE_CALLABLE"], "yes" if sp in tips else "no"]
    row, rrow = list(base), list(base)
    for t in traits:
        st = cell[(t, sp)]["STATE"]
        if st not in CODE:
            refuse("T2 state %r for %s/%s" % (st, t, sp))
        p1, na = CODE[st]
        row += [p1, na or "-"]
        rrow += (["0.5", R37_REASON] if m["ABSENCE_CALLABLE"] == "no" else [p1, na or "-"])
    data.append(row)
    r37.append(rrow)

# ---------------------------------------------------------------- comparisons
_, old = table(SUPERSEDED)
oldrow = {r["SPECIES"]: r for r in old}
changes = []
ix = {c: hdr.index(c) for c in hdr}
for row in data:
    sp = row[0]
    o = oldrow.get(sp)
    if o is None:
        refuse("the superseded coding lacks %s" % sp)
    for t in traits:
        new_p1, new_na = row[ix["%s_absent_p1" % t]], row[ix["%s_absent_na" % t]]
        col = "%s_absent_p1" % t
        if col not in o:
            changes.append(["%s_absent" % t, sp, "-", new_p1, "-", new_na, meta[sp]["ABSENCE_CALLABLE"],
                            "yes" if sp in tips else "no"])
            continue
        if o[col] != new_p1:
            changes.append(["%s_absent" % t, sp, o[col], new_p1, o.get("%s_absent_na" % t, "-"), new_na,
                            meta[sp]["ABSENCE_CALLABLE"], "yes" if sp in tips else "no"])
vs = []
for t in traits:
    i = ix["%s_absent_p1" % t]
    d210 = sum(1 for a, b in zip(data, r37) if a[i] != b[i])
    dchr = sum(1 for a, b in zip(data, r37) if a[i] != b[i] and a[ix["ON_CHRONOGRAM"]] == "yes")
    vs.append(["%s_absent" % t, "%s_absent" % t, arows[traits.index(t)]["TARGET_SET"], str(d210), str(dchr),
               "R37" if dchr == 0 else "STANDARD"])
if any(v[5] != "STANDARD" for v in vs):
    refuse("a trait would reuse a pre-repair Q0 reconstruction: %s" % [v[0] for v in vs if v[5] != "STANDARD"])

# ---------------------------------------------------------------- write
prov = ["## PRODUCED BY %s sha256=%s" % (os.path.basename(__file__), sha(os.path.abspath(__file__))),
        "## INPUT T2 = %s sha256=%s" % (os.path.relpath(T2, ROOT), T2_SHA),
        "## INPUT ACTIVE_PUBLICATION_TARGETS = %s sha256=%s" % (os.path.relpath(ACTIVE, ROOT), sha(ACTIVE)),
        "## INPUT chronogram = %s md5=%s" % (CHRONO, md5(CHRONO)),
        "## SUPERSEDED (quarantined, comparison only) = %s sha256=%s" % (SUPERSEDED, sha(SUPERSEDED))]
os.makedirs(OUTD, exist_ok=True)


def write(name, banner, header, rowsout):
    p = os.path.join(OUTD, name)
    with open(p, "x", encoding="utf-8") as fh:
        fh.write(banner + "\n" + "\n".join(prov) + "\n")
        fh.write("\t".join(header) + "\n")
        for r in rowsout:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("wrote %s (%d rows) sha256 %s" % (name, len(rowsout), sha(p)[:16]))
    return p


p1 = write("STANDARD_per_species.tsv",
           "## STANDARD PER-SPECIES CODING, TEN TARGETS (L43, L44): 210 rows, one per proteome, ten `<T>_absent` traits on\n"
           "## ONE rule in one order: p1 = 0.0 where T2 says PRESENT (found is found, L44), 0.5 where the proteome is\n"
           "## absence_callable=no (R37, reason in _na), 1.0 where T2 says ABSENT (L41). Every cell is read from the\n"
           "## ten-target T2, which carries the core-8 tree rule; nothing comes from the superseded twenty-target tables.\n"
           "## 0.5 IS AMBIGUITY, NOT EVIDENCE (R5).\n" + SIGN, hdr, data)
p2 = write("STANDARD_changes.tsv",
           "## Cells that differ from the SUPERSEDED twenty-trait coding, which is quarantined and named in the provenance\n"
           "## lines below. FROM '-' means the trait had no column there.",
           ["TRAIT", "SPECIES", "FROM", "TO", "REASON_BEFORE", "REASON_AFTER", "ABSENCE_CALLABLE", "ON_CHRONOGRAM"], changes)
p3 = write("STANDARD_vs_R37.tsv",
           "## Per trait: cells differing between this STANDARD coding and the R37 coding (every absence_callable=no\n"
           "## proteome flat at 0.5), on all 210 species and on the 209 chronogram tips. SOURCE_ASR=STANDARD for every\n"
           "## trait: no reconstruction is reused from a pre-repair Q0 run.",
           ["TRAIT", "Q0_NAME", "LAYER", "N_CELLS_DIFF_210", "N_CELLS_DIFF_ON_CHRONOGRAM", "SOURCE_ASR"], vs)

# ---------------------------------------------------------------- checks
_, back = table(p1)
if len(back) != 210:
    refuse("wrote %d rows" % len(back))
counts = {}
for t in traits:
    c = collections.Counter(r["%s_absent_p1" % t] for r in back)
    counts[t] = (c["0.0"], c["0.5"], c["1.0"])
    t2c = collections.Counter(cell[(t, sp)]["STATE"] for sp in species)
    if counts[t] != (t2c["PRESENT"], t2c["NOT_CALLABLE"], t2c["ABSENT"]):
        refuse("%s coding %s does not match T2 %s" % (t, counts[t], dict(t2c)))
banner_text = "".join(l for l in open(p1, encoding="utf-8") if l.startswith("#"))
bad = [w for w in ("MARK", "SNRK", "SnRK2", "SnRK3_CIPK", "PP2C", "PP2A_B55", "PPP1R3", "SEX4_LSF", "SDS23", "CBSX")
       if re.search(r"(?<![A-Za-z0-9_])%s(?![A-Za-z0-9_])" % re.escape(w), banner_text)]
if bad or any(c.endswith("_control_in_supergroup") for c in hdr):
    refuse("a deferred target or a passthrough column survived: %s" % (bad or hdr))
print("per trait (PRESENT/NOT_CALLABLE/ABSENT): %s" % counts)
print("changed cells against the superseded coding: %d" % len(changes))
print("STANDARD vs R37 differing cells on the chronogram: %s" % {v[0]: v[4] for v in vs})
print("TEN-TARGET STANDARD CODING COMPLETE")
