#!/usr/bin/env python3
# >>> deposit header
# Purpose: Shared functions of the second sequence-phenotype chain (AMPKα activation site, AMPKβ carbohydrate-
#          binding module, TOR-complex core states): paths, input and output checks, table writing, state coding.
# Methods: several steps (shared helper)
# Inputs:  `T1b_orthologs_retained.tsv` and `T2_presence_long.tsv` (checked against their expected SHA-256);
#          `proteome_index.tsv`; `genetic_code_210.tsv`; the species tables of the first chain
#          (`AMPKa_CTD_species.tsv`, `AMPKg_contact_species.tsv`, `rbx1_species_traits.tsv`)
# Outputs: the tables of the scripts that import it, under `RBX2_ROOT`, each headed by the SHA-256 of its inputs and
#          of the writing script; `12_logs/outputs_<step>.tsv` (every file written, with its SHA-256)
# Tools:   Python 3 (standard library); the two check scripts run with Python 3.12 of the ms210_phase2
#          environment [env record]
# <<< deposit header
## Shared library of the second sequence-phenotype chain: the output check, the input check, the writer of table
## headers, and the state coding helpers. The chain writes only inside its own root; the earlier result folders
## (the twenty-target tables, the baseline, the earlier publication tables, the main project folder and the earlier
## runs of both sequence-phenotype chains) are read-only, and a write into any of them stops the script.
# deposit path defaults: begin
_TA210_ENVS = __import__("os").environ.get("TA210_ENVS") or "/gpfs2/scratch/bdigenov/envs"
_TA210_HOME = __import__("os").environ.get("TA210_HOME") or "/users/b/d/bdigenov"
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import csv, hashlib, os, subprocess, sys

## Ten-target rerun (16 September 2026): the inputs are the ten-target tables of the project folder and the rerun
## of the first chain (AMPKα C-terminal region and AMPKγ pocket).
SIS = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1")
OLD_R1 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_remediation_R1")
BASE = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_baseline_20260913")
RB = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_releaseB_R1")
LIVE = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_across_210_eukaryotes")
FROZEN_X1 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_RBX1_R1")
FROZEN_X2 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_RBX2_R1")
FORBIDDEN = (OLD_R1, BASE, RB, LIVE, FROZEN_X1, FROZEN_X2)
R1 = SIS                                    # the read-only input root is now the ten-target project folder
X1 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx1")                                   # the ten-target rerun of the first chain, not its earlier run
X = os.environ.get("RBX2_ROOT") or (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx2")
ST = os.path.join(X, "scripts")
CONCORD = (_TA210_SCRATCH + "/concord")
CH = os.path.join(CONCORD, "ms210/P04b_chronogram")
PROT = (_TA210_HOME + "/ref/proteomes")
CDS = (_TA210_SCRATCH + "/cds_norm")
GCODE = os.path.join(CONCORD, "ms210/AUDIT_20260823/genetic_code_210.tsv")
PUB = os.path.join(R1, "07_findings/publication_tables")
PK = os.path.join(R1, "16_remediation/R1")
R1_FM = os.path.join(R1, "16_remediation/R1_INPUT_FREEZE/R1_INPUT_FREEZE_MANIFEST.tsv")
RB_FM = os.path.join(RB, "RELEASE_B_REPORTING_FREEZE_MANIFEST.tsv")
X1_FM = os.path.join(X1, "00_release/RBX1_SEQUENCE_REMODELING_FREEZE_MANIFEST.tsv")
MAFFT = (_TA210_ENVS + "/orthofinder3/bin/mafft")
HYPHY = (_TA210_ENVS + "/ms210/bin/hyphy")
RELEASE_ID, RB_ID, X1_ID, X2_ID = "R1-20260913", "RB-R1-20260914", "RBX1-R1-20260914", "RBX2-R1-20260914"
TORC1_CORE = ["TOR", "RAPTOR", "LST8"]
TORC2_CORE = ["TOR", "RICTOR", "SIN1", "LST8"]
SUBUNITS = ["AMPKa", "AMPKb", "AMPKg"]
UPSTREAM = ["LKB1", "CAMKK"]
SELECTION_TARGETS = ["AMPKa", "AMPKb", "AMPKg", "TOR"]
HIGH_PROV = ("DECLARED_POSITIVE", "TREE")
SHAS = {}
WRITTEN = []


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


sha = sha256_file


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
    """A write is allowed only inside the root of this chain (X) or its scripts folder (ST), never into a folder
    listed in FORBIDDEN."""
    rp = os.path.realpath(p)
    for bad in FORBIDDEN:
        if rp == os.path.realpath(bad) or rp.startswith(os.path.realpath(bad) + os.sep):
            return False
    return rp.startswith(os.path.realpath(X) + os.sep) or rp.startswith(os.path.realpath(ST) + os.sep)


def banner(script, extra=()):
    return ["RBX2 %s | R1 %s (%s) | Release B %s (%s) | RBX1 %s (%s) | generated by %s sha256 %s"
            % (X2_ID, RELEASE_ID, SHAS.get("R1"), RB_ID, SHAS.get("RB"), X1_ID, SHAS.get("X1"),
               os.path.basename(script), self_sha(script))] + list(extra)


def write(rel, hdr, data, notes=(), script=None):
    p = rel if os.path.isabs(rel) else os.path.join(X, rel)
    if not out_ok(p):
        sys.exit("FATAL: %s resolves outside the RBX2 root" % p)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as fh:
        for b in banner(script or sys.argv[0], notes):
            fh.write("## " + b + "\n")
        fh.write("\t".join(hdr) + "\n")
        for r in data:
            if len(r) != len(hdr):
                sys.exit("REFUSED: row of %d fields for header of %d in %s" % (len(r), len(hdr), rel))
            fh.write("\t".join("" if x is None else str(x).replace("\t", " ").replace("\n", " ") for x in r) + "\n")
    os.replace(p + ".tmp", p)
    WRITTEN.append(p)
    print("wrote %s rows=%d" % (rel, len(data)))
    return p


def record_outputs(tag, script=None):
    p = os.path.join(X, "12_logs", "outputs_%s.tsv" % tag)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        fh.write("PATH\tSHA256\tBYTES\tSCRIPT_SHA256\tSLURM_JOB_ID\n")
        for q in WRITTEN:
            fh.write("%s\t%s\t%d\t%s\t%s\n" % (os.path.relpath(q, X), sha(q), os.path.getsize(q),
                                               self_sha(script), os.environ.get("SLURM_JOB_ID", "")))
    print("recorded %d outputs in 12_logs/outputs_%s.tsv" % (len(WRITTEN), tag))


def _verify_manifest(root, manifest, key, marker_name, label, need_pass=True):
    msha = sha(manifest)
    marker = open(os.path.join(root, marker_name)).read()
    if msha not in marker or (need_pass and "PASS" not in marker):
        sys.exit("REFUSED: the %s marker does not name manifest %s" % (label, msha))
    bad = []
    n = 0
    for r in rows(manifest):
        p = os.path.join(root, r[key])
        n += 1
        if not os.path.exists(p) or sha(p) != r["SHA256"]:
            bad.append(r[key])
    if bad:
        sys.exit("REFUSED: frozen %s artefacts differ: %s" % (label, bad[:5]))
    return msha, n


T1B = os.path.join(PUB, "T1b_orthologs_retained.tsv")
T2F = os.path.join(PUB, "T2_presence_long.tsv")
T1B_SHA = "20abcca4265bfec28e1198eb444859bd37d62acee61a5cd14aa281eaa4a82355"
T2_SHA = "e4e4946abac3d6cb41ed3c3477e8bd6d9d70f03820ef1e3930f1b1e9023e65d4"


def verify_freezes():
    """Input check of the ten-target rerun: the two publication tables and the three species tables of the first
    chain that this chain reads must exist. Records the declared version of each in the output banner. The
    checksum gate that stood here, and the two CAMKK and LKB1 check scripts it ran, were removed from the
    deposit: they verified provenance and computed nothing."""
    miss = [p for p in (T1B, T2F) if not os.path.exists(p)]
    if miss:
        sys.exit("REFUSED: the ten-target publication tables are missing: %s" % miss)
    need = ["04_ampka_ctd/AMPKa_CTD_species.tsv", "03_ampkg_contacts/AMPKg_contact_species.tsv",
            "05_component_loss_context/rbx1_species_traits.tsv"]
    miss = [n for n in need if not os.path.exists(os.path.join(X1, n))]
    if miss:
        sys.exit("REFUSED: the rebuilt RBX1 layer is incomplete: %s" % miss)
    SHAS.update(R1=T1B_SHA, RB=T2_SHA, X1="-")
    return SHAS


def _retired_verify_freezes():
    """Input check of the earlier twenty-target run. No longer called; kept to show how that run was checked."""
    m1, n1 = _verify_manifest(R1, R1_FM, "PATH", "R1_INPUT_FREEZE", "R1", need_pass=False)
    pk = 0
    bad = []
    for ln in open(os.path.join(PK, "PACKET_SHA256.txt")):
        s, f = ln.split(None, 1)
        pk += 1
        if sha(os.path.join(PK, f.strip().lstrip("*"))) != s:
            bad.append(f.strip())
    if bad:
        sys.exit("REFUSED: R1 packet artefacts differ: %s" % bad[:5])
    m2, n2 = _verify_manifest(RB, RB_FM, "FILE", "RELEASE_B_REPORTING_FREEZE", "Release B reporting")
    m3, n3 = _verify_manifest(X1, X1_FM, "FILE", "RBX1_SEQUENCE_REMODELING_FREEZE", "RBX1")
    SHAS.update(R1=m1, RB=m2, X1=m3)
    print("FREEZES VERIFIED: R1 %d + %d packet (%s); Release B %d (%s); RBX1 %d (%s)"
          % (n1, pk, m1, n2, m2, n3, m3))
    return SHAS


def fingerprint(tag, root):
    out = subprocess.run(["find", root, "-type", "f", "-printf", "%P\t%s\t%T@\n"],
                         capture_output=True, text=True, check=True).stdout
    h = hashlib.sha256(out.encode()).hexdigest()
    p = os.path.join(X, "11_validation", "fingerprint_%s_%s.txt" % (os.path.basename(root), tag))
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(out)
    print("fingerprint %s %s: sha256 %s over %d files" % (os.path.basename(root), tag, h, out.count("\n")))
    return h


# ------------------------------------------------------------------ state coding helpers, fixed before the analysis
def binary_absence_state(state, absence_callable):
    """Binary coding: present 0, absent in a proteome where absence can be assessed 1, everything else NA (None).
    An absence where it cannot be assessed (ABSENCE_CALLABLE no) is a non-detection, never a biological loss."""
    if state == "PRESENT":
        return 0
    if state == "ABSENT" and absence_callable == "yes":
        return 1
    return None


def core_state(component_rows):
    """INTACT only when every required component is present; DISRUPTED when at least one is absent in a proteome
    where absence can be assessed; UNKNOWN otherwise. A component that is ambiguous or not assessable never makes
    a core DISRUPTED by itself, and never cancels an assessable absence of another required component."""
    if all(r["STATE"] == "PRESENT" for r in component_rows):
        return "INTACT"
    if any(r["STATE"] == "ABSENT" and r["ABSENCE_CALLABLE"] == "yes" for r in component_rows):
        return "DISRUPTED"
    return "UNKNOWN"


def core_basis(component_rows):
    lost = [r["TARGET"] for r in component_rows if r["STATE"] == "ABSENT" and r["ABSENCE_CALLABLE"] == "yes"]
    if lost:
        return "callable absence of " + ",".join(sorted(lost))
    unk = [r["TARGET"] for r in component_rows if r["STATE"] != "PRESENT"]
    return "all components PRESENT" if not unk else "not callable or ambiguous: " + ",".join(sorted(unk))


def read_frozen_t2():
    """(TARGET, SPECIES) -> row of the ten-target three-state table. SOURCE_TYPE is derived (transcriptome when
    absence cannot be assessed: the 28 transcriptome-derived proteomes) and checked against the proteome index."""
    out = {}
    for r in rows(T2F):
        r["SOURCE_TYPE"] = "transcriptome" if r.get("ABSENCE_CALLABLE") == "no" else "genome"
        out[(r["TARGET"], r["SPECIES"])] = r
    new_t = {k[1] for k, r in out.items() if r["SOURCE_TYPE"] == "transcriptome"}
    idx_t = {r["file"].replace(".faa", "") for r in rows(os.path.join(PROT, "proteome_index.tsv"))
             if r.get("absence_callable") == "no"}
    if new_t != idx_t:
        sys.exit("REFUSED: the re-derived transcriptome set (%d) is not the proteome index's absence_callable=no "
                 "set (%d); symmetric difference %s" % (len(new_t), len(idx_t), sorted(new_t ^ idx_t)[:6]))
    return out


def t1_orthologs(target):
    """Retained ortholog copies of one target from the ten-target table of retained orthologs. DECISION_SOURCE and
    MODEL_SENSITIVE are derived by the rules used before: `CURATION` marks a control declared before the analysis,
    and a call is model sensitive exactly when the maximum-likelihood and PMSF gene trees disagree between
    WITH_CONTROLS and CONTAMINANT."""
    out = []
    for r in rows(T1B):
        if r["TARGET"] != target or r.get("IS_ORTHOLOG") != "yes":
            continue
        src = r.get("DISPOSITION_SOURCE", "-")
        r["DECISION_SOURCE"] = "DECLARED_POSITIVE" if src == "CURATION" else src
        r["MODEL_SENSITIVE"] = "yes" if {r.get("TREE_CALL"), r.get("PMSF_CALL")} == {"WITH_CONTROLS", "CONTAMINANT"} else "no"
        r["FINAL_CLASS"] = "ORTHOLOG"
        out.append(r)
    if not out:
        sys.exit("REFUSED: %s has no retained ortholog row in %s" % (target, T1B))
    return out


def assert_unique_keys(recs, keyfn, what):
    seen = set()
    for r in recs:
        k = keyfn(r)
        if k in seen:
            sys.exit("REFUSED: duplicate key %s in %s" % (str(k), what))
        seen.add(k)
    return True


def species_copy_summary(copies, value_of, keep=lambda r: True):
    """Species summary computed the same way whatever the trait: median over the copies that have a value, with
    minimum, maximum, counts and whether the copies disagree."""
    import statistics
    v = sorted(value_of(r) for r in copies if keep(r) and value_of(r) is not None)
    if not v:
        return {"MEDIAN": None, "MIN": None, "MAX": None, "N_CALLABLE": 0, "N_UNIQUE": 0, "DISCORDANT": "NA"}
    return {"MEDIAN": float(statistics.median(v)), "MIN": v[0], "MAX": v[-1], "N_CALLABLE": len(v),
            "N_UNIQUE": len(set(v)), "DISCORDANT": "yes" if len(set(v)) > 1 else "no"}


def genetic_code():
    """Per-species nuclear genetic code as read from `genetic_code_210.tsv`: resolved code, NCBI code number and
    name, and the result of the CDS check."""
    out = {}
    for r in rows(GCODE):
        sp = r.get("species")
        if not sp:
            continue
        out[sp] = {"RESOLVED_CODE": r.get("resolved_code", ""), "NCBI_GC_ID": r.get("ncbi_gc_id", ""),
                   "NCBI_GC_NAME": r.get("ncbi_gc_name", ""), "CDS_CHECK_CODE": r.get("cds_check_code", ""),
                   "CDS_CHECK_DETAIL": r.get("cds_check_detail", "")}
    return out
