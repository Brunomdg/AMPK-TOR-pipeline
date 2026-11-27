#!/usr/bin/env python3
# >>> deposit header
# Purpose: Shared functions of the chain that measures the AMPKα C-terminal region and the AMPKγ pocket contacts
#          on the retained ortholog copies: input paths and checks, ortholog copies, species states, table writing.
# Methods: several steps (shared helper)
# Inputs:  `T1b_orthologs_retained.tsv` (retained ortholog copies of the ten targets) and `T2_presence_long.tsv`
#          (three-state species table), each checked against its expected SHA-256; `proteome_index.tsv`
# Outputs: the tables of the scripts that import it, under `RBX1_ROOT`, each headed by the SHA-256 of its inputs
#          and of the writing script; `11_logs/outputs_<step>.tsv` (every file written, with its SHA-256)
# Tools:   Python 3 (standard library); the two check scripts run with Python 3.12 of the ms210_phase2
#          environment [env record]
# <<< deposit header
## Shared library of this chain, as rerun on 16 September 2026 on the ten-target ortholog set: the ortholog table
## had changed, so the AMPKα C-terminal and AMPKγ pocket measurements were recomputed and the earlier values set
## aside. The functions are those of the earlier version of this library, so the other scripts of the chain run
## unchanged against it. Three things differ:
##   1. Inputs. The ortholog copies come from `T1b_orthologs_retained.tsv` (the retained orthologs of the ten
##      targets) and the species states from `T2_presence_long.tsv`, not from the earlier twenty-target tables.
##   2. Input check. `verify_freezes()` runs two check scripts, which confirm that the CAMKK and LKB1 calls are
##      final and that the combined ortholog table was built from their current versions, and it requires both
##      input tables to have their expected SHA-256.
##   3. Outputs. Every output is written under the root of this rerun. The earlier result folders (the
##      twenty-target tables, the baseline, the earlier publication tables, the main project folder and the
##      earlier run of this chain) are read-only, and a write that resolves into any of them stops the script.
## Two columns that the earlier ortholog table carried and the current one does not are derived here by the rules
## used for the earlier table:
##   DECISION_SOURCE: `CURATION` (a positive control declared before the analysis) becomes `DECLARED_POSITIVE`;
##                    `TREE` (the call was made from the gene tree) and `R66` (the gene tree was ambiguous and at
##                    least one method of each search arm nominated the candidate) are kept as they are.
##   MODEL_SENSITIVE: "yes" when one of the two gene-tree calls (`TREE_CALL`, and `PMSF_CALL` from the PMSF
##                    cross-check) is `WITH_CONTROLS` and the other `CONTAMINANT`: the rule at line 695 of
##                    `rl_r1_three_state_release.py`.
##   SOURCE_TYPE: transcriptome when absence cannot be assessed in the proteome (ABSENCE_CALLABLE "no"; the 28
##                transcriptome-derived proteomes), otherwise genome. At load the set is compared with the
##                absence_callable column of the shared proteome index, and any difference stops the script.
##                No table of an earlier run is read, here or anywhere else in this rerun.
# deposit path defaults: begin
_TA210_ENVS = __import__("os").environ.get("TA210_ENVS") or "/gpfs2/scratch/bdigenov/envs"
_TA210_HOME = __import__("os").environ.get("TA210_HOME") or "/users/b/d/bdigenov"
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import csv, hashlib, os, subprocess, sys

SIS = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1")
OLD_R1 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_remediation_R1")
BASE = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_baseline_20260913")
RB = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_releaseB_R1")
LIVE = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_across_210_eukaryotes")
RBX1_FROZEN = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_RBX1_R1")
FORBIDDEN = (OLD_R1, BASE, RB, LIVE, RBX1_FROZEN)
R1 = SIS                                    # the read-only input root is now the ten-target project folder
X = os.environ.get("RBX1_ROOT") or os.path.join(SIS, "16_control_repair/refit/rbx1")
ST = os.path.join(X, "scripts")
CONCORD = (_TA210_SCRATCH + "/concord")
P02 = os.path.join(CONCORD, "ms210/P02_gamma_pocket")
P03 = os.path.join(CONCORD, "ms210/P03_actd")
CH = os.path.join(CONCORD, "ms210/P04b_chronogram")
PROT = (_TA210_HOME + "/ref/proteomes")
PUB = os.path.join(SIS, "07_findings/publication_tables")
T1B = os.path.join(PUB, "T1b_orthologs_retained.tsv")
T2F = os.path.join(PUB, "T2_presence_long.tsv")
T1B_SHA = "20abcca4265bfec28e1198eb444859bd37d62acee61a5cd14aa281eaa4a82355"
T2_SHA = "e4e4946abac3d6cb41ed3c3477e8bd6d9d70f03820ef1e3930f1b1e9023e65d4"
GUARDS = [os.path.join(SIS, "16_control_repair/scripts/rl_require_camkk_lkb1_joint_freeze.py"),
          os.path.join(SIS, "16_control_repair/scripts/rl_require_clean_camkk_lkb1_master_inputs.py")]
PY = (_TA210_ENVS + "/ms210_phase2/bin/python3")
Q0SP = os.path.join(SIS, "14_remodeling/Q0_replacement/Q0_REPLACEMENT_per_species.tsv")
ASR = os.path.join(SIS, "14_remodeling/Q0_replacement/asr")
ASRLIB = os.path.join(SIS, "14_remodeling/scripts/rl_rm_asr_lib.R")
MAFFT = (_TA210_ENVS + "/orthofinder3/bin/mafft")
PY_MS210 = (_TA210_ENVS + "/ms210/bin/python3")
RELEASE_ID, RB_ID, X_ID = "SISTER-10TARGET-20260916", "-", "RBX1-REFIT-10T-20260916"
METAPREDICT_VERSION = "3.0.2"
CANON_CTD_BOUNDARY = 290
CORE8 = ["TOR", "RAPTOR", "RICTOR", "LST8", "AMPKa", "AMPKb", "AMPKg", "SIN1"]
COMPONENTS = ["TOR", "RAPTOR", "RICTOR"]
SUBUNITS = ["AMPKa", "AMPKb", "AMPKg"]
SITESETS = ["site1", "site3", "site4"]
PRIMARY_SITESET = "site1"
R1_SHA = RB_SHA = None
WRITTEN = []
_T1B_CACHE = _T2_CACHE = None


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rows(p, required=True):
    if not os.path.exists(p):
        if required:
            sys.exit("REFUSED: missing %s" % p)
        return []
    with open(p, errors="replace") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t", quoting=csv.QUOTE_NONE))


def self_sha(path=None):
    return sha(os.path.abspath(path or sys.argv[0]))


def out_ok(p):
    """A write is allowed only inside the root of this rerun (X), never into a folder listed in FORBIDDEN,
    whatever the path looks like."""
    rp = os.path.realpath(p)
    for bad in FORBIDDEN:
        if rp == os.path.realpath(bad) or rp.startswith(os.path.realpath(bad) + os.sep):
            return False
    return rp.startswith(os.path.realpath(X) + os.sep)


def banner(script, extra=()):
    return ["RBX1 %s | inputs %s: T1b sha256 %s, T2 sha256 %s | generated by %s sha256 %s"
            % (X_ID, RELEASE_ID, R1_SHA, RB_SHA, os.path.basename(script), self_sha(script))] + list(extra)


def write(rel, hdr, data, notes=(), script=None):
    p = rel if os.path.isabs(rel) else os.path.join(X, rel)
    if not out_ok(p):
        sys.exit("FATAL: %s resolves outside the refit root" % p)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as fh:
        for b in banner(script or sys.argv[0], notes):
            fh.write("## " + b + "\n")
        fh.write("\t".join(hdr) + "\n")
        for r in data:
            fh.write("\t".join("" if c is None else str(c) for c in r) + "\n")
    os.replace(p + ".tmp", p)
    WRITTEN.append(p)
    return p


def record_outputs(tag, script):
    p = os.path.join(X, "11_logs", "outputs_%s.tsv" % tag)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        fh.write("## %s\n" % banner(script)[0])
        fh.write("PATH\tSHA256\tBYTES\tWRITTEN_BY\tSCRIPT_SHA256\tSLURM_JOB_ID\n")
        for f in WRITTEN:
            fh.write("%s\t%s\t%d\t%s\t%s\t%s\n" % (os.path.relpath(f, X), sha(f), os.path.getsize(f),
                                                   os.path.basename(script), self_sha(script),
                                                   os.environ.get("SLURM_JOB_ID", "")))
    print("recorded %d outputs in 11_logs/outputs_%s.tsv" % (len(WRITTEN), tag))


def verify_freezes():
    """Input check before any computation: the two CAMKK and LKB1 check scripts (GUARDS) must both exit without
    error, and the two publication tables must have their expected SHA-256."""
    global R1_SHA, RB_SHA
    for g in GUARDS:
        if not os.path.exists(g):
            sys.exit("REFUSED: guard %s is missing" % g)
        r = subprocess.run([PY, g], capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("REFUSED: %s exited %d: %s" % (os.path.basename(g), r.returncode,
                                                    (r.stdout + r.stderr).strip()[-400:]))
    bad = []
    for name, path, pin in (("T1b", T1B, T1B_SHA), ("T2", T2F, T2_SHA)):
        got = sha(path)
        if got != pin:
            bad.append("%s %s is %s, pinned %s" % (name, path, got, pin))
    if bad:
        sys.exit("REFUSED: the ten-target publication tables are not the pinned ones: %s" % "; ".join(bad))
    R1_SHA, RB_SHA = T1B_SHA, T2_SHA
    print("SISTER GATE VERIFIED: joint freeze and master inputs clean; T1b %s; T2 %s" % (T1B_SHA[:12], T2_SHA[:12]))
    return R1_SHA, RB_SHA


def fingerprint(tag, root=SIS):
    out = subprocess.run(["find", root, "-type", "f", "-printf", "%P\t%s\t%T@\n"],
                         capture_output=True, text=True, check=True).stdout
    h = hashlib.sha256(out.encode()).hexdigest()
    p = os.path.join(X, "10_validation", "fingerprint_%s_%s.txt" % (os.path.basename(root), tag))
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        fh.write(out)
    print("fingerprint %s %s: sha256 %s over %d files" % (os.path.basename(root), tag, h, out.count("\n")))
    return h


def _load_t1b():
    global _T1B_CACHE
    if _T1B_CACHE is None:
        out = []
        for r in rows(T1B):
            if r.get("IS_ORTHOLOG") != "yes":
                continue
            src = r.get("DISPOSITION_SOURCE", "-")
            r["DECISION_SOURCE"] = "DECLARED_POSITIVE" if src == "CURATION" else src
            r["MODEL_SENSITIVE"] = "yes" if {r.get("TREE_CALL"), r.get("PMSF_CALL")} == {"WITH_CONTROLS", "CONTAMINANT"} else "no"
            r["FINAL_CLASS"] = "ORTHOLOG"
            out.append(r)
        _T1B_CACHE = out
    return _T1B_CACHE


def t1_orthologs(target):
    """Every retained ortholog copy of one target in the ten-target table. No other row enters the sequence analysis."""
    got = [r for r in _load_t1b() if r["TARGET"] == target]
    if not got:
        sys.exit("REFUSED: %s has no retained ortholog row in %s" % (target, T1B))
    return got


def t2():
    """The ten-target three-state table, with SOURCE_TYPE derived from ABSENCE_CALLABLE and checked against the
    proteomes that the proteome index marks as not assessable (absence_callable "no")."""
    global _T2_CACHE
    if _T2_CACHE is None:
        out = []
        for r in rows(T2F):
            r["SOURCE_TYPE"] = "transcriptome" if r.get("ABSENCE_CALLABLE") == "no" else "genome"
            out.append(r)
        new_t = {r["SPECIES"] for r in out if r["SOURCE_TYPE"] == "transcriptome"}
        idx_t = {r["file"].replace(".faa", "") for r in rows(os.path.join(PROT, "proteome_index.tsv"))
                 if r.get("absence_callable") == "no"}
        if new_t != idx_t:
            sys.exit("REFUSED: the re-derived transcriptome set (%d species) is not the proteome index's "
                     "absence_callable=no set (%d); symmetric difference %s"
                     % (len(new_t), len(idx_t), sorted(new_t ^ idx_t)[:6]))
        _T2_CACHE = out
    return _T2_CACHE


def trait01(cell):
    """Binary coding of one species cell, fixed before the analysis: present 0, absent in a proteome where absence
    can be assessed 1, every other cell NA (None)."""
    if cell["STATE"] == "PRESENT":
        return 0
    if cell["STATE"] == "ABSENT" and cell["ABSENCE_CALLABLE"] == "yes":
        return 1
    return None
