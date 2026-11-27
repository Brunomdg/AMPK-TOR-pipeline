#!/usr/bin/env python3
# >>> deposit header
# Purpose: Shared helpers of the upstream-kinase refit: folder and file locations, the input check run by every step,
#          the table writer with its banner line, the CAMKK reference list and the species-state helpers.
# Methods: several steps (shared helper)
# Inputs:  the publication tables `T1_ortholog_table.tsv`, `T1b_orthologs_retained.tsv` and `T2_presence_long.tsv`;
#          the rebuilt
#          `rbx1_species_traits.tsv`, `TORC_core_states.tsv` and `upstream_kinase_states.tsv`
# Outputs: none of its own; write() and record_outputs() write each step's tables and `11_logs/outputs_<step>.tsv`
# Tools:   Python 3 (standard library); paths of MAFFT and IQ-TREE 3 (orthofinder3 environment); versions: see envs/
# <<< deposit header
## Shared helpers for the upstream-kinase refit scripts in this folder: folder locations, the check of the input
## tables, the banner line written at the top of every output table, the table writer and the species-state
## helpers, adapted from rbx2c_common.py with the folders changed. The scripts write only inside this refit's own
## output folder; the earlier analysis folders, the baseline and the live project are read-only, and a write into
## any of them stops the run. A subfamily label or a sequence-support status is an annotation added after the
## orthology calls and never decides one.
# deposit path defaults: begin
_TA210_ENVS = __import__("os").environ.get("TA210_ENVS") or "/gpfs2/scratch/bdigenov/envs"
_TA210_HOME = __import__("os").environ.get("TA210_HOME") or "/users/b/d/bdigenov"
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import csv, hashlib, os, re, subprocess, sys

BASE = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_baseline_20260913")
R1 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_remediation_R1")
RB = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_releaseB_R1")
X1 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_RBX1_R1")
X2 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_RBX2_R1")
XC = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_RBX2_clade_robustness_R1")
LIVE = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_across_210_eukaryotes")
## From 16 September 2026 the inputs are the publication tables of the ten-target survey and the rebuilt earlier
## analyses (the `rbx1` and `rbx2` folders); the older folders stay in `FORBIDDEN`, so no write can reach them.
SIS = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1")
OLD_R1, OLD_X1, OLD_X2, OLD_X3 = R1, X1, X2, (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_RBX3_R1")
FORBIDDEN = (OLD_R1, RB, OLD_X1, OLD_X2, OLD_X3, XC, BASE, LIVE)
R1 = SIS
X1 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx1")
X2 = (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx2")
X3 = os.environ.get("RBX3_ROOT") or (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx3")
ST = os.path.join(X3, "scripts")
CONCORD = (_TA210_SCRATCH + "/concord")
CH = os.path.join(CONCORD, "ms210/P04b_chronogram")
COMBINED = os.path.join(CONCORD, "combined")
PROT = (_TA210_HOME + "/ref/proteomes")
PUB = os.path.join(R1, "07_findings/publication_tables")
DECL = os.path.join(R1, "02_reference_data/RAW_identity-controls_declaration_v7.tsv")
NEGS = os.path.join(R1, "02_reference_data/CURATED_negatives.tsv")
DOMTBL = os.path.join(R1, "03_evidence_by_method/hmmer_genomewide/domtbl")
## Tree inputs, as chosen on 17 September 2026, when the tree analysis was rerun as originally designed. The
## membership tree that decides CAMKK orthology is the second-version screening tree of the CAMKK analysis; its
## cropped alignment, tip roles and per-copy calls are copied, under the names the scripts read, into an input
## folder of this refit, with the checksum of every source in `ADAPTER_SOURCES.tsv`. The CAMKK target definition
## is copied unchanged from the ten-target project into the same input folder. The ancestral-state
## reconstructions use the standard three-state coding (files ending `__TRAITS`), the co-transition source is the
## rebuilt table of the correlated-change analysis of all targets, and the clade inventory is the rebuilt
## clade-robustness analysis; `FORBIDDEN` above still lists the older copy of that analysis.
PK = os.path.join(X3, "02_inputs/declarations")
PHYLO = os.path.join(X3, "02_inputs/membership_tree")
ASR = os.path.join(R1, "14_remodeling/Q0_replacement/asr")
Q7_COTRANSITION = os.path.join(R1, "14_remodeling/Q7_alltargets/Q7_COTRANSITION.tsv")
XC = os.path.join(R1, "16_control_repair/refit/rbx2c")
MAFFT = (_TA210_ENVS + "/orthofinder3/bin/mafft")
IQTREE = (_TA210_ENVS + "/orthofinder3/bin/iqtree3")
CHRONO = os.path.join(CH, "species_chronogram_209.treefile")
RELEASE_ID = "R1-20260913"
X3_ID = "RBX3-R1-20260914"
SHAS = {}
WRITTEN = []

## The seven CAMKK subfamily reference sequences, declared here before any candidate is classified;
## rbx3_preflight.py checks each against the declaration table and the proteome FASTA and stops on any mismatch.
CAMKK_REFS = [
    ("CAMKKREF01", "Homo_sapiens", "NP_757343.2", "CAMKK1", "CAMKK_METAZOAN_LIKE", "-"),
    ("CAMKKREF02", "Homo_sapiens", "XP_047284057.1", "CAMKK2", "CAMKK_METAZOAN_LIKE", "-"),
    ("CAMKKREF03", "Saccharomyces_cerevisiae", "NP_011055.3", "SAK1", "SNF1_SNRK1_ACTIVATOR_LIKE",
     "SAK1_TOS3_ELM1_LIKE"),
    ("CAMKKREF04", "Saccharomyces_cerevisiae", "NP_011336.3", "TOS3", "SNF1_SNRK1_ACTIVATOR_LIKE",
     "SAK1_TOS3_ELM1_LIKE"),
    ("CAMKKREF05", "Saccharomyces_cerevisiae", "NP_012876.1", "ELM1", "SNF1_SNRK1_ACTIVATOR_LIKE",
     "SAK1_TOS3_ELM1_LIKE"),
    ("CAMKKREF06", "Arabidopsis_thaliana", "NP_566876.3", "GRIK1", "SNF1_SNRK1_ACTIVATOR_LIKE", "GRIK_SNAK_LIKE"),
    ("CAMKKREF07", "Arabidopsis_thaliana", "NP_200863.2", "GRIK2", "SNF1_SNRK1_ACTIVATOR_LIKE", "GRIK_SNAK_LIKE"),
]
L1_CLASSES = ["CAMKK_METAZOAN_LIKE", "SNF1_SNRK1_ACTIVATOR_LIKE", "OTHER_CAMKK_CATEGORY", "SUBFAMILY_UNRESOLVED"]
L2_CLASSES = ["SAK1_TOS3_ELM1_LIKE", "GRIK_SNAK_LIKE", "SNF1_ACTIVATOR_SUBTYPE_UNRESOLVED"]
NEGATIVE_TARGETS = ["LKB1", "AMPKa", "MARK", "SNRK", "SnRK2", "SnRK3_CIPK"]
TORC1_CORE = ["TOR", "RAPTOR", "LST8"]
TORC2_CORE = ["TOR", "RICTOR", "SIN1", "LST8"]


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


sha = sha256_file


def _read(p, required=True):
    if not os.path.exists(p):
        if required:
            sys.exit("REFUSED: missing %s" % p)
        return []
    with open(p, errors="replace") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t",
                                   quoting=csv.QUOTE_NONE))


## Several scripts still name the retired twenty-target publication tables; rows() reads those names from the
## ten-target tables instead, never from a copy of the old ones. `T1_ortholog_table_all20.tsv` is read from
## `T1_ortholog_table.tsv`: FINAL_CLASS is `ORTHOLOG` when IS_ORTHOLOG is yes (accepted orthologs and curated
## controls), otherwise the value of the `DISPOSITION` column; DECISION_SOURCE is `DECLARED_POSITIVE` for sequences
## decided by manual curation (source `CURATION`) and otherwise the source column as given; REFERENCE_ROLE and
## TRUTH_CLASS read `CONTROL` and `POSITIVE` for the declared controls; MODEL_SENSITIVE is yes when one of the
## maximum-likelihood and site-heterogeneous (PMSF) tree calls is `WITH_CONTROLS` and the other `CONTAMINANT`.
## `T2_presence_all20_long.tsv` is read from `T2_presence_long.tsv`, with SOURCE_TYPE transcriptome when absence
## cannot be assessed and genome otherwise. Every other path is read as it is.
T1_ALL20 = os.path.join(PUB, "T1_ortholog_table_all20.tsv")
T2_ALL20 = os.path.join(PUB, "T2_presence_all20_long.tsv")


def rows(p, required=True):
    ap = os.path.abspath(p)
    if ap == os.path.abspath(T1_ALL20):
        out = []
        for r in _read(T1F):
            src = r.get("DISPOSITION_SOURCE", "-")
            r["FINAL_CLASS"] = "ORTHOLOG" if r.get("IS_ORTHOLOG") == "yes" else r.get("DISPOSITION", "-")
            r["DECISION_SOURCE"] = "DECLARED_POSITIVE" if src == "CURATION" else src
            r["MODEL_SENSITIVE"] = "yes" if {r.get("TREE_CALL"), r.get("PMSF_CALL")} == {"WITH_CONTROLS", "CONTAMINANT"} else "no"
            r["REFERENCE_ROLE"] = "CONTROL" if src == "CURATION" else "-"
            r["TRUTH_CLASS"] = "POSITIVE" if src == "CURATION" else "-"
            out.append(r)
        return out
    if ap == os.path.abspath(T2_ALL20):
        out = _read(T2F)
        for r in out:
            r["SOURCE_TYPE"] = "transcriptome" if r.get("ABSENCE_CALLABLE") == "no" else "genome"
        return out
    return _read(p, required)


def self_sha(path=None):
    return sha(os.path.abspath(path or sys.argv[0]))


def out_ok(p):
    """A write is allowed only inside this refit's output folder or its scripts folder."""
    rp = os.path.realpath(p)
    for bad in FORBIDDEN:
        if rp == os.path.realpath(bad) or rp.startswith(os.path.realpath(bad) + os.sep):
            return False
    return rp.startswith(os.path.realpath(X3) + os.sep) or rp.startswith(os.path.realpath(ST) + os.sep)


def banner(script, extra=()):
    return ["RBX3 %s | R1 %s (%s) | Release B (%s) | RBX1 (%s) | RBX2 (%s) | RBX2C (%s) | generated by %s "
            "sha256 %s" % (X3_ID, RELEASE_ID, SHAS.get("R1"), SHAS.get("RB"), SHAS.get("X1"), SHAS.get("X2"),
                           SHAS.get("XC"), os.path.basename(script), self_sha(script))] + list(extra)


def write(rel, hdr, data, notes=(), script=None):
    p = rel if os.path.isabs(rel) else os.path.join(X3, rel)
    if not out_ok(p):
        sys.exit("FATAL: %s resolves outside the RBX3 root" % p)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as fh:
        for b in banner(script or sys.argv[0], notes):
            fh.write("## " + b + "\n")
        fh.write("\t".join(hdr) + "\n")
        for r in data:
            if len(r) != len(hdr):
                sys.exit("REFUSED: row of %d fields for header of %d in %s" % (len(r), len(hdr), rel))
            fh.write("\t".join("" if x is None else str(x).replace("\t", " ").replace("\n", " ")
                               for x in r) + "\n")
    os.replace(p + ".tmp", p)
    WRITTEN.append(p)
    print("wrote %s rows=%d" % (rel, len(data)))
    return p


def record_outputs(tag, script=None):
    p = os.path.join(X3, "11_logs", "outputs_%s.tsv" % tag)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        fh.write("PATH\tSHA256\tBYTES\tSCRIPT_SHA256\tSLURM_JOB_ID\n")
        for q in WRITTEN:
            fh.write("%s\t%s\t%d\t%s\t%s\n" % (os.path.relpath(q, X3), sha(q), os.path.getsize(q),
                                               self_sha(script), os.environ.get("SLURM_JOB_ID", "")))
    print("recorded %d outputs in 11_logs/outputs_%s.tsv" % (len(WRITTEN), tag))


T1F = os.path.join(PUB, "T1_ortholog_table.tsv")
T1_SHA = "3a71dbcf05408d7102d375cf49bfe6c5633a88216b338b03a1ea860e41041b23"
T1B = os.path.join(PUB, "T1b_orthologs_retained.tsv")
T2F = os.path.join(PUB, "T2_presence_long.tsv")
T1B_SHA = "20abcca4265bfec28e1198eb444859bd37d62acee61a5cd14aa281eaa4a82355"
T2_SHA = "e4e4946abac3d6cb41ed3c3477e8bd6d9d70f03820ef1e3930f1b1e9023e65d4"


def verify_freezes():
    """Input check run by every step: the publication tables T1, T1b and T2 and the three rebuilt tables this
    refit reads must exist. Fills `SHAS` for the banner line of the output tables. The checksum gate that stood
    here, and the two CAMKK and LKB1 check scripts it ran, were removed from the deposit: they verified
    provenance and computed nothing."""
    miss = [p for p in (T1F, T1B, T2F) if not os.path.exists(p)]
    if miss:
        sys.exit("REFUSED: the publication tables are missing: %s" % miss)
    need = [(X1, "05_component_loss_context/rbx1_species_traits.tsv"), (X2, "03_torc_core/TORC_core_states.tsv"),
            (X2, "04_ampka_activation/upstream_kinase_states.tsv")]
    miss = ["%s/%s" % (os.path.basename(r), n) for r, n in need if not os.path.exists(os.path.join(r, n))]
    if miss:
        sys.exit("REFUSED: the rebuilt layers are incomplete: %s" % miss)
    SHAS.update(R1=T1B_SHA, RB=T2_SHA, X1="-", X2="-", XC="-")
    return SHAS


_NEWICK_SAFE = re.compile(r"[^A-Za-z0-9_-]")


def newick_safe(s):
    return _NEWICK_SAFE.sub("_", s)


def tip_in(d, tip):
    """Look up a tip name built from the domain table in a dict keyed by the tree's own tip names."""
    if tip in d:
        return d[tip]
    k = newick_safe(tip)
    hits = [v for kk, v in d.items() if newick_safe(kk) == k]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        sys.exit("REFUSED: %s matches %d tree tips once Newick label sanitisation is undone; the mapping is "
                 "not one to one and no call may be taken from it" % (tip, len(hits)))
    return None


def binary_absence_state(state, absence_callable):
    """Present 0, absent in a proteome where absence can be assessed 1, anything else None (missing). An absence
    with ABSENCE_CALLABLE no is a non-detection by the searches and never a biological loss."""
    if state == "PRESENT":
        return 0
    if state == "ABSENT" and absence_callable == "yes":
        return 1
    return None


def core_state(component_rows):
    if all(r["STATE"] == "PRESENT" for r in component_rows):
        return "INTACT"
    if any(r["STATE"] == "ABSENT" and r["ABSENCE_CALLABLE"] == "yes" for r in component_rows):
        return "DISRUPTED"
    return "UNKNOWN"


def read_frozen_t2():
    """The ten-target species-state table (T2) keyed by (target, species), with SOURCE_TYPE derived again
    (transcriptome when absence cannot be assessed, genome otherwise)."""
    out = {}
    for r in rows(T2F):
        r["SOURCE_TYPE"] = "transcriptome" if r.get("ABSENCE_CALLABLE") == "no" else "genome"
        out[(r["TARGET"], r["SPECIES"])] = r
    return out


def t1_orthologs(target):
    return [r for r in rows(os.path.join(PUB, "T1_ortholog_table_all20.tsv"))
            if r["TARGET"] == target and r["FINAL_CLASS"] == "ORTHOLOG"]


def clade_of():
    return {r["SPECIES"]: r["CLADE"] for r in rows(os.path.join(PUB, "T2_presence_all20_long.tsv"))}


def supergroup_of():
    return {r["SPECIES"]: r["SUPERGROUP"] for r in rows(os.path.join(PUB, "T2_presence_all20_long.tsv"))}


def read_fasta(p):
    seqs, name, buf = {}, None, []
    op = open
    with op(p, errors="replace") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].split()[0], []
            else:
                buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def tool_version(binary, args=("--version",)):
    try:
        r = subprocess.run([binary] + list(args), capture_output=True, text=True, timeout=120)
        return (r.stdout + r.stderr).strip().splitlines()[0]
    except Exception as e:
        return "UNAVAILABLE: %s" % e
