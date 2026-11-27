#!/usr/bin/env python3
# >>> deposit header
# Purpose: Derive each method's candidate calls for the eight targets across the 210 proteomes from the raw
#          output of five methods at the jointly chosen settings, then build the candidate list, run the
#          winner-take-all contaminant screen and write the per-method and final lists with a run report.
# Methods: Candidate discovery
# Inputs:  OrthoFinder `Orthogroups.txt` and graph-rule member tables; SonicParanoid `flat.ortholog_groups.tsv`
#          and pairwise member tables; the genome-wide HMMER tblout files; the CDD and InterProScan scans in
#          `concord/combined`; the proteomes in `concord/null_split/fwd`; `callers.tsv` (signatures allowed to
#          call); `THRESHOLD_PILOT_PERTARGET.tsv` (zero-decoy cuts); `R59_union_settings` and `R59_arm1_selected`
#          (the jointly chosen settings, named by the control set given in `R59_TAG`); `PANEL_DECLARATION_v2.tsv`
#          (target-side and contaminant-side models); the current declared-control table
# Outputs: in the report folder of the control set: `01_orthofinder_full_calls.tsv` to `05_interpro_calls.tsv`
#          and `04b_cdd_LST8_raw_hits.tsv` (calls per method), `06_candidate_list.tsv`,
#          `06b_candidate_provenance.tsv` (which methods nominated each protein), `07_contaminant_screen.tsv`,
#          `08_final_retained_orthologs.tsv`, `09_final_lists_per_method.tsv` and `RUN_REPORT.md`
# Tools:   Python 3 standard library and the project modules `accessions`, `control_declaration`,
#          `call_vocabulary`, `targets`, `evidence_roots` and `join`; reads the output of OrthoFinder,
#          SonicParanoid, HMMER, CDD and InterProScan without running them; versions: see envs/
# <<< deposit header
"""Per-method ortholog calls derived from raw method output, for 8 targets in 210 species.

Five methods in two independent arms, each method reporting only for itself:

    orthology arm   OrthoFinder (full run) and SonicParanoid
    domain arm      the HMMER profile panel, CDD and InterProScan

The target list and the output folder names come from `lib/rl_targets`, and the
HMMER inputs from the evidence roots in `lib/rl_evidence_roots`; with the default
set of eight targets the output is byte-identical to the original.

The governing principle. A protein nominated by several methods is not stronger
evidence than a protein nominated by one. The number of methods that nominate a
protein tells where the call came from; it is not evidence. Everything structural
here implements that one sentence: the long format, the separate table of which
methods nominated each protein, the forbidden column names and the absence of any
count.

No derived view is read. `EVIDENCE_MAP.tsv` exists because reading a derived view
as though it were the evidence once led to a wrong root cause being reported. The
inputs are raw method output only; the files of the candidate-pool layers are not
used.

What these files are:

    raw method output      ->  membership rule of that method  ->  MEMBER list
    (what the tool wrote)      (shared controls, published cut)     (what is reported)

Every output here is the third object, calls derived from raw method output,
never the raw output itself. Absence from a MEMBER list means "not nominated under
this rule", not "the method did not see it". Questions of coverage go to
`rl_evidence.py` against the raw source.

The candidate list is not evidence. It is a working list naming the proteins that
enter the contaminant screen (since 8 September 2026 `rl_phylo_prepare.py` also
builds the gene-tree tip set from it). It has no biological interpretation and
carries no column saying which methods nominated a protein: that information is
kept in a separate long table, so that a support count cannot be obtained with a
split().
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import collections
import csv
import glob
import gzip
import os
import re
import sys

ROOT = os.environ.get("PROJECT_ROOT") or (
    (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1"))
sys.path.insert(0, os.path.join(ROOT, "01_scripts", "lib"))

import rl_accession                       # noqa: E402
import rl_declaration                     # noqa: E402
import rl_vocab                           # noqa: E402
import rl_targets                         # noqa: E402
import rl_evidence_roots                  # noqa: E402
from rl_join import JoinFailed            # noqa: E402


# --- The joint-selection outputs are named by control set; a tag is required ---
# The untagged files were solved against a different control set and have been
# deleted. A reader that silently fell back to one could not tell a current result
# from an outdated one; this once let a job run for 32 minutes on the wrong input
# while its preliminary check reported success. A missing tagged input stops the
# script and is never replaced by another file, and there is no default tag,
# because a default is a fallback.
_R59_TAG = os.environ.get("R59_TAG") or ""


def _r59(p, must_exist=True):
    if not _R59_TAG:
        raise SystemExit(
            "REFUSING: R59_TAG is not set. R59 outputs are per control set "
            "(setA = declaration v5, setB = the clean literature subset); there "
            "is no untagged default. Re-run with R59_TAG=setA or R59_TAG=setB.")
    # base__<tag><set suffix>.ext: the tag names the control set and the suffix the
    # target set (empty for the default eight targets). One helper names every
    # tagged file of the joint selection.
    q = rl_targets.tagged(p, _R59_TAG)
    if must_exist and not os.path.exists(q):
        raise SystemExit("REFUSING: %s does not exist. Solve that control set "
                         "first; this script will not read an untagged file." % q)
    return q
# ------------------------------------------------------------------------------


CONCORD = (_TA210_SCRATCH + "/concord")
PANELDIR = (_TA210_SCRATCH + "/panel/BACKUP_full02_20260819")

# --- Raw evidence sources; each is listed as KIND=RAW in EVIDENCE_MAP.tsv ------
SRC = {
    "orthofinder_full": ROOT + "/03_evidence_by_method/orthofinder_full/Orthogroups.txt",
    "sonic":            PANELDIR + "/ortholog_groups/flat.ortholog_groups.tsv",
    "hmmer_tblout":     ROOT + "/03_evidence_by_method/hmmer_genomewide/tblout",
    "hmmer_domtbl":     ROOT + "/03_evidence_by_method/hmmer_genomewide/domtbl",
    "combined":         CONCORD + "/combined",
    "proteomes_plain":  CONCORD + "/null_split/fwd",
}
# CDD and InterProScan both come from `combined`, the raw source registered as
# `interproscan_combined`: all 17 signature databases (CDD and 16 InterProScan
# member databases) per protein, with the score kind, score and e-value. The other
# raw source, cddfull/out/fwd, stores signatures as numeric PSSM identifiers
# (gnl|CDD|270713) rather than the cd05169-style accessions that callers.tsv
# uses, so reading it would need an identifier-mapping step that can silently
# drop signatures. The two were checked to agree: on Homo_sapiens, cd05169 rows
# carry identical bit scores and e-values in both.

# --- Declared rules: not evidence, and not derived views ---------------------
CALLERS = CONCORD + "/callers.tsv"                       # which signatures may call
PILOT = ROOT + "/06_calibration_and_checks/THRESHOLD_PILOT_PERTARGET.tsv"
# Jointly chosen settings. These are not evidence: they are solved settings,
# produced by `rl_r59_union_optimum.py` from the same raw sources that this script
# reads. One row per target and method; the union over a target's rows is what
# recovers every declared control. The earlier per-method settings (each method
# required to recover every control alone) are superseded and no longer read.
R59_SETTINGS = _r59(ROOT + "/07_findings/R59_union_settings.tsv")
PANEL = ROOT + "/02_reference_data/hmmer_panel/PANEL_DECLARATION_v2.tsv"
NCBI_CUT = CONCORD + "/bitscore_specific.txt"            # reference only
DECLDIR = ROOT + "/02_reference_data"

# The report folder is tagged too. Two control sets produce two sets of tables,
# and an untagged folder would let them overwrite each other, the same kind of
# collision that made the first parallel solve useless.
OUT = rl_targets.report_dir(ROOT, _R59_TAG)

# The single source of the target list: the environment variable TARGET_SET selects
# the rows of the target registry in use; the default set is checked to equal the
# eight targets fixed in advance, in their fixed order.
TARGETS = rl_targets.targets()
ARM = {"orthofinder_full": "1_orthology", "sonic": "1_orthology",
       "hmmer": "2_domain", "cdd": "2_domain", "interpro": "2_domain"}

NEG_INF = float("-inf")
EXPECT_SPECIES = 210

# Trial runs only. Limits how many of the 210 per-species files are read, so one
# target can be followed end to end before the full batch. The size checks in
# `gate()` still require all 210 files to be present, and every output written
# under this limit is marked DRY-RUN, so a trial can never be taken for a
# complete run.
LIMIT = int(os.environ.get("MS210_LIMIT_SPECIES") or 0)


def corpus(pattern_dir, pattern):
    f = sorted(glob.glob(os.path.join(pattern_dir, pattern)))
    return f[:LIMIT] if LIMIT else f


def hmmer_corpus():
    """The HMMER tblout files this run reads, from the declared evidence roots
    (03_evidence_by_method/HMMER_EVIDENCE_ROOTS.tsv). For the default eight targets
    that is exactly the 210 files of the primary root, in the sorted order that
    `corpus()` returned; a supplementary root (models searched separately with
    identical options) is read only for another target set, so the tables for the
    eight targets are built from identical bytes."""
    f = rl_evidence_roots.tblout_files()
    return f[:LIMIT] if LIMIT else f

# Models built from the target's own proteins may not judge them; copied
# unchanged from rl_phase2_panelv2_screen.py. A "contaminant" model built from the
# target's own in-paralogs is the target, not a look-alike. Both TOR models are
# PIKK-family models that contain TOR; without this exclusion the screen removed
# 25 of 30 TOR controls.
SELF_CONTAM = {"RAPTOR__RAPTOR_s_own_lineage-specific_in-paralogs": "RAPTOR",
               "KOG0891.faa.final_tree.fa": "TOR",
               "TOR__ATR": "TOR"}

# InterProScan member databases that may not call: predictions of disorder, coiled
# coils, signal peptides and transmembrane segments are not protein families.
IPS_BARRED = {"IPS:MobiDBLite", "IPS:Coils", "IPS:SignalP", "IPS:TMHMM",
              "IPS:Phobius"}

# Verification 5: column names that would carry a support statistic are refused.
# A count of how often methods differ is a support statistic reached from the
# other direction, so the terms for disagreement are refused too.
BANNED_COL = re.compile(
    r"SUPPORT|CONSENSUS|N_METHODS|RANK|TIER|VOTE|AGREE|MAJORITY|SCORE_COMBINED"
    r"|CONFIDENCE|ARM_CALL|DISAGREE|OVERLAP|CONCORDAN|N_SHARED|N_BOTH")
FRACTION = re.compile(r"^\d+\s*/\s*\d+$")

DERIVED_MARKERS = ("layer_calls_domain_scored", "layer_calls_orthology",
                   "/inputs/", "04_final_tables")

norm = rl_accession.normalise
LOG = []


def say(msg):
    print(msg, flush=True)
    LOG.append(msg)


def rows(path, delim="\t"):
    """Read a table, skipping `##` banner lines; the first other line is the header."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        yield from csv.DictReader(
            (l for l in fh if not l.startswith("##")), delimiter=delim)


# ------------------------------------------------------- step 0  input checks
def gate():
    """Check, before any row is read, that every source is raw method output and has its expected size."""
    say("STEP 0  input gate")
    for name, p in SRC.items():
        low = p.lower()
        for m in DERIVED_MARKERS:
            if m in low:
                raise SystemExit(
                    "REFUSED: %s resolves inside a DERIVED area (%s). "
                    "EVIDENCE_MAP: only KIND=RAW may be cited for what an "
                    "instrument saw." % (name, m))
        if not os.path.exists(p):
            raise SystemExit("missing RAW source %s: %s" % (name, p))

    checks = [
        ("orthofinder_full", "lines", 643339, SRC["orthofinder_full"]),
        ("sonic", "lines", 134540, SRC["sonic"]),
        ("combined", "files:*.csv.gz", EXPECT_SPECIES, SRC["combined"]),
        ("proteomes_plain", "files:*.faa", EXPECT_SPECIES, SRC["proteomes_plain"]),
        # HMMER is the one method not listed in EVIDENCE_MAP.tsv, so
        # `rl_check_evidence_provenance.py` never checks its size. The check is
        # made here instead, and the omission is reported.
        ("hmmer_tblout", "files:*.gz", EXPECT_SPECIES, SRC["hmmer_tblout"]),
    ]
    for name, kind, expect, path in checks:
        if kind == "lines":
            with open(path, errors="replace") as fh:
                got = sum(1 for _ in fh)
        else:
            got = len(glob.glob(os.path.join(path, kind.split(":")[1])))
        flag = "ok" if got == expect else "MISMATCH"
        say("  %-18s %-16s expected=%-8d actual=%-8d %s"
            % (name, kind, expect, got, flag))
        if got != expect:
            raise SystemExit(
                "REFUSED: %s is %d, expected %d. A source that shrank or moved "
                "must fail rather than be read as complete." % (name, got, expect))
    say("  NOTE  hmmer_genomewide is absent from EVIDENCE_MAP.tsv; asserted here "
        "at %d+%d. Recommend registering it." % (EXPECT_SPECIES, EXPECT_SPECIES))
    # Every declared HMMER evidence root that this target set may read is checked
    # for 210 domtbl and 210 tblout files by rl_evidence_roots.roots() itself, which
    # stops on an incomplete root, so a half-written supplementary root cannot be read.
    for r in rl_evidence_roots.roots():
        say("  %-18s evidence root   %s (scope %s)"
            % ("hmmer_root", r["ROOT"], r.get("SET_SCOPE") or "all"))


# ------------------------------------------------------------- shared loaders
def load_controls():
    """Map (target, normalised accession) to CONTROL or CONTROL_HELDOUT, and return the controls in use.

    A control is in use (`binding` in the code) when its REVIEW field is empty; any row with a
    REVIEW entry is held out, as in rl_phase2_panelv2_screen.py and
    MASTER_freeze_parameters.py.
    """
    path, label = rl_declaration.current(DECLDIR)
    role, binding = {}, collections.defaultdict(set)
    for r in rows(path):
        t = (r.get("TARGET") or "").strip()
        a = (r.get("PROTEIN_ID") or "").strip()
        if t not in TARGETS or not a:
            continue
        held = bool((r.get("REVIEW") or "").strip())
        role[(t, norm(a))] = "CONTROL_HELDOUT" if held else "CONTROL"
        if not held:
            binding[t].add(norm(a))
    say("  controls %s: %d rows, binding per target %s"
        % (label, len(role), {t: len(binding[t]) for t in TARGETS}))
    return role, binding, os.path.basename(path)


def species_index():
    """Map each normalised accession to its species, from the plain proteomes (headers `Species||ACC`)."""
    idx = {}
    for fn in sorted(glob.glob(os.path.join(SRC["proteomes_plain"], "*.faa"))):
        sp = os.path.basename(fn)[:-4]
        with open(fn, errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    idx[norm(line[1:].split()[0])] = sp
    say("  species index: %d accessions over %d proteomes"
        % (len(idx), len(glob.glob(os.path.join(SRC["proteomes_plain"], "*.faa")))))
    return idx


def load_callers():
    """For each target and signature database, the kept signatures that may call, with score kind and threshold."""
    keep = collections.defaultdict(lambda: collections.defaultdict(dict))
    n = 0
    for r in rows(CALLERS):
        if (r.get("status") or "").strip() != "kept":
            continue
        fam = (r.get("family") or "").strip()
        inst = (r.get("instrument") or "").strip()
        if fam not in TARGETS or inst in IPS_BARRED:
            continue
        thr = (r.get("threshold") or "").strip()
        try:
            thr = float(thr)
        except ValueError:
            continue
        keep[fam][inst][(r.get("signature") or "").strip()] = (
            (r.get("kind") or "").strip(), thr)
        n += 1
    say("  callers.tsv kept and callable: %d signatures" % n)
    return keep


def published_cuts():
    """Every model that can be made free of decoys calls, at its zero-decoy cut or a stricter published cut.

    What this replaces, and why. The earlier rule kept, for each target, the
    models with the highest control recall at each model's published cut. It had
    two faults, both measured:

      (a) It ignored contamination. At equal recall it kept the most promiscuous
          model: TOR used FRB_dom (0.00% decoy false-positive rate), FATC (0.31%)
          and PI3_PI4_kinase (1.24%) together at 100% recall, and PI3_PI4_kinase
          is the shared PIKK domain that matches ATM, ATR, DNA-PK, SMG1 and TRRAP.
      (b) Adding a contamination tie-break fixed (a) but created a worse problem:
          it discarded clean models to keep less clean ones. AMPKβ kept K07199
          (94% recall, 1.92% false-positive rate) and dropped AMPKBI (91%, 0.00%)
          for three points of recall. Choosing a single model was the mistake.

    The rule now. Contamination is a property not of a model but of a model at
    a cut, and the pilot measured the cut at which each model calls no decoy. So:

      1. Take every target-side model that was not built from the controls
         themselves (a control may not recover itself).
      2. Give each its own zero-decoy cut (`CUTSEQ_FPR0` / `CUTDOM_FPR0`), raised to
         its published cut where that is stricter (see the loop below).
      3. Keep every model that still recovers at least one control at that cut.
      4. Take the union of their calls. No model is chosen over another.

    Contamination is zero by construction rather than by ranking, and recall is
    that of the union rather than of the best single model. Nothing is discarded
    for being slightly less sensitive, and no threshold is invented: the cuts are
    the pilot's own measurements.

    This also lets in models that had no published cut at all, which the earlier
    rule excluded outright. Measured consequence: AMPKγ was calling on CBS, which
    has 94% recall at a 23.35% decoy false-positive rate and falls to 17% when
    made decoy-free, whereas KOG1764 reaches 81% recall with zero decoys and was
    excluded only for lacking a Pfam gathering threshold. Likewise LST8 PTHR19842
    (96% against 87% for K08266) and AMPKα LECA_OG047_AMPKA2 (86% against 69% for
    K07198).

    Limitation: "zero decoys" means zero among the declared decoys, which come
    from five species (Arabidopsis, human, yeast, T. gondii, T. brucei). It is a
    lower bound on contamination, not a complete census.
    """
    raw = collections.defaultdict(list)
    skipped = collections.defaultdict(list)
    # The pilot fixed in advance, then its companion for the target set
    # (THRESHOLD_PILOT_PERTARGET__<set>.tsv) when another target set is in use; the
    # default eight targets read exactly the original file.
    for r in rl_targets.pilot_rows(PILOT):
        if r.get("POLE") != "target" or r.get("FOR_TARGET") not in TARGETS:
            continue
        t = r["FOR_TARGET"]
        if r.get("SELF_BUILT") != "no":                       # built from the controls
            skipped[t].append((r["MODEL"], "self-built (R28)"))
            continue
        s_, d_ = r.get("CUTSEQ_FPR0", ""), r.get("CUTDOM_FPR0", "")
        if s_ in ("", "none") and d_ in ("", "none"):
            skipped[t].append((r["MODEL"], "no zero-decoy cut exists"))
            continue
        try:
            rec0 = float((r.get("REC_FPR0") or "0%").rstrip("%"))
        except ValueError:
            rec0 = 0.0
        if rec0 <= 0:
            skipped[t].append(
                (r["MODEL"], "recovers 0%% of controls at its zero-decoy cut -- "
                             "it cannot separate this target from its decoys"))
            continue
        cs = NEG_INF if s_ in ("", "none") else float(s_)
        cd = NEG_INF if d_ in ("", "none") else float(d_)
        # A model is not clean or unclean in itself; a model at a cut is. Every
        # model that can be made clean calls, and the calls are combined as a union.
        # Take the stricter of the published cut and the zero-decoy cut, separately
        # for the sequence score and the domain score. A zero-decoy cut alone says
        # little about specificity when the decoy set is small: AMPKβ has 104
        # declared decoys, none scoring above 1.9 bits on AMPKBI, so its zero-decoy
        # cut is 1.9 bits against a published 24.3, 13 times looser. Applied alone,
        # it produced 80,662 AMPKβ calls (384 per species) in one run. The published
        # cut shows that the domain is present; the zero-decoy cut shows
        # specificity; a model must meet both. This keeps every gain (TOR
        # PI3_PI4_kinase stays at a domain cut of 237.8, AMPKγ CBS at 118.4/35.8),
        # because there the zero-decoy cut is the stricter one. Models with no
        # published cut use the zero-decoy cut alone.
        p_ = (r.get("PUBLISHED_CUT") or "").strip()
        if p_ and p_ != "none":
            parts = [x.strip() for x in p_.split("/")]
            try:
                ps = float(parts[0])
                pd = float(parts[1]) if len(parts) > 1 else ps
                cs, cd = max(cs, ps), max(cd, pd)
            except ValueError:
                pass
        try:
            pubfpr = float((r.get("PUB_FPR") or "").rstrip("%") or "nan")
        except ValueError:
            pubfpr = float("nan")
        raw[t].append((r["MODEL"], cs, cd, rec0, pubfpr,
                       (r.get("PUBLISHED_SOURCE") or "zero-decoy cut").strip()))

    cuts = {}
    say("  HMMER: EVERY model that can be made clean calls, at the STRICTER of "
        "its published and zero-decoy cuts; calls are UNIONED:")
    for t in TARGETS:
        cuts[t] = raw[t]
        if not raw[t]:
            say("    %-7s NO model stays clean with any control -- HMMER cannot "
                "call for this target" % t)
            continue
        best = max(x[3] for x in raw[t])
        say("    %-7s %d models call (union). best single-model recall at zero "
            "decoys %.0f%%  [%s]"
            % (t, len(raw[t]), best,
               ", ".join("%s@%.0f%%" % (x[0], x[3]) for x in sorted(
                   raw[t], key=lambda z: -z[3]))))
    return cuts, skipped


def r59_settings():
    """Map (target, method) to the jointly chosen setting, or to None where the method makes no call.

    The settings are chosen jointly: the fewest proteins called across the 210
    proteomes, subject to the union of all methods recovering every declared
    control. A row here may therefore cover only part of its target's control
    set, and that is correct: the union over the target's rows is what must
    recover every control.

    Why the earlier per-method settings are no longer read. They required each
    method to reach every control alone, which misreads the requirement on the
    union as a requirement on each method. That forced every method to loosen
    until it caught the single most divergent control by itself: InterProScan
    for AMPKα ended at IPS:SMART SM00220, the generic serine/threonine kinase
    domain, calling 57,468 proteins, 274 per species for a family expected at
    one to three copies. That rule was replaced by the joint choice, and this
    function is why the code no longer applies it either.
    """
    out = {}
    if not os.path.exists(R59_SETTINGS):
        raise SystemExit(
            "REFUSED: %s is missing. The joint settings are not optional -- "
            "without them this script would fall back to the R52 cuts everywhere "
            "and look like it worked, which is the R58 failure. Run "
            "rl_r59_union_optimum.py first." % R59_SETTINGS)
    for r in rows(R59_SETTINGS):
        t = (r.get("TARGET") or "").strip()
        inst = (r.get("INSTRUMENT") or "").strip()
        if t not in TARGETS or (r.get("SOURCE") or "").strip() != "r59_joint":
            continue                      # orthology-arm rows are carried, not re-applied
        setting = (r.get("SETTING") or "").strip()
        if not setting or setting == "(no setting)":
            out[(t, inst)] = None         # this method makes no call here
            continue
        out[(t, inst)] = (setting, (r.get("CUT") or "").strip(),
                          (r.get("CUT_DOM") or "").strip(),
                          (r.get("SCORE_KIND") or "").strip(),
                          (r.get("CONTROLS_COVERED") or "").strip(),
                          (r.get("COST") or "").strip(),
                          (r.get("PER_SPECIES") or "").strip())
    return out


def apply_r59(cuts, keep):
    """Apply the jointly chosen settings, per target and per method.

    A rule takes effect only where code that writes the tables reads it, and this
    function is where the joint settings are read. The table of settings it adds
    to RUN_REPORT.md lets a reader check that from the output rather than from
    this file.

    A method may make no call. If the joint solution assigns it no controls and
    the union is cheaper without it, its cell is emptied and marked
    `R59_NO_SETTING`. That is a decision, not a gap, and it differs both from a
    method silenced because optimisation could not improve it and from a
    fallback to the earlier rule (`FALLBACK_R52`).
    """
    ledger = []
    sel = r59_settings()
    say("\n  R59 in force: settings chosen JOINTLY so the UNION recovers every "
        "declared control at the lowest corpus-wide cost. A cell covering part of "
        "the control set is correct, not a shortfall.")
    for t in TARGETS:
        got = sel.get((t, "hmmer"), "absent")
        if got == "absent":
            ledger.append(("hmmer", t, "FALLBACK_R52",
                           "%d model(s) at the stricter of published and "
                           "zero-decoy cuts" % len(cuts.get(t, [])),
                           "no R59 row for this cell"))
            say("    hmmer   %-7s FALLBACK to R52 -- no R59 row" % t)
        elif got is None:
            cuts[t] = []
            ledger.append(("hmmer", t, "R59_NO_SETTING", "makes no call",
                           "the joint solution covers this target's controls "
                           "more cheaply without it"))
            say("    hmmer   %-7s NO CALL by the joint solution" % t)
        else:
            m, cs, cd, _k, ncov, cost, persp = got
            cuts[t] = [(m, float(cs), float(cd), 100.0, float("nan"),
                        "R59 joint")]
            ledger.append(("hmmer", t, "R59",
                           "%s seq>=%s dom>=%s" % (m, cs, cd),
                           "covers %s controls, %s non-targets, %s per species"
                           % (ncov, cost, persp)))
            say("    hmmer   %-7s R59  %s@%s/%s  covers %s controls, cost %s "
                "(%s/species)" % (t, m, cs, cd, ncov, cost, persp))

    for inst in ("cdd", "interpro"):
        for t in TARGETS:
            got = sel.get((t, inst), "absent")
            n_prior = sum(len(v) for k, v in keep[t].items()
                          if (k == "CDD") == (inst == "cdd"))
            keep[t] = collections.defaultdict(dict, keep[t])
            if got == "absent":
                ledger.append((inst, t, "FALLBACK_R52",
                               "%d kept signature(s) at declared thresholds"
                               % n_prior, "no R59 row for this cell"))
                say("    %-7s %-7s FALLBACK to R52 -- no R59 row (%d kept "
                    "signatures still call)" % (inst, t, n_prior))
                continue
            for k in [k for k in keep[t] if (k == "CDD") == (inst == "cdd")]:
                keep[t][k] = {}
            if got is None:
                ledger.append((inst, t, "R59_NO_SETTING", "makes no call",
                               "the joint solution covers this target's controls "
                               "more cheaply without it"))
                say("    %-7s %-7s NO CALL by the joint solution" % (inst, t))
                continue
            setting, cut, _cd, kind, ncov, cost, persp = got
            layer, sig = setting.split(":", 1)
            keep[t][layer] = {sig: (kind, float(cut))}
            ledger.append((inst, t, "R59", "%s:%s %s %s" % (layer, sig, kind, cut),
                           "covers %s controls, %s non-targets, %s per species"
                           % (ncov, cost, persp)))
            say("    %-7s %-7s R59  %s:%s %s %s  covers %s controls, cost %s "
                "(%s/species)" % (inst, t, layer, sig, kind, cut, ncov, cost,
                                  persp))
    return cuts, keep, ledger


def load_poles():
    """For each target, the target-side and contaminant-side models of panel v2.

    Models built from the target's own in-paralogs are left off the contaminant side.
    """
    tgt = collections.defaultdict(set)
    con = collections.defaultdict(set)
    for r in rows(PANEL):
        t = (r.get("FOR_TARGET") or "").strip()
        m = (r.get("MODEL") or "").strip()
        pole = (r.get("POLE") or "").strip()
        if t not in TARGETS or not m:
            continue
        if pole == "target":
            tgt[t].add(m)
        elif pole == "contaminant":
            if SELF_CONTAM.get(m) == t:                        # built from the target itself
                continue
            con[t].add(m)
    say("  panel v2 poles (contaminant pole excludes %d self-in-paralog models): %s"
        % (len(SELF_CONTAM),
           {t: "%d/%d" % (len(tgt[t]), len(con[t])) for t in TARGETS}))
    return tgt, con


# ------------------------------------------------------------ output writers
CALL_HDR = ("TARGET\tARM\tMETHOD\tSPECIES\tPROTEIN_ID\tROLE\tCALL\tSCOPE\t"
            "ABSENCE_MEANS\tEVIDENCE_1\tEVIDENCE_2\tEVIDENCE_3\tRAW_SOURCE\n")


def emit(path, banner, header, body):
    """Write a call file, banner first, refusing column names or values that would carry a support count."""
    cols = header.rstrip("\n").split("\t")
    for c in cols:
        if BANNED_COL.search(c.upper()):
            raise SystemExit(
                "REFUSED: column %r in %s matches the anti-consensus pattern. "
                "This design treats multiplicity of nomination as provenance, "
                "not evidence." % (c, os.path.basename(path)))
    with open(path, "w") as fh:
        if LIMIT:
            fh.write("## DRY-RUN: only %d of %d per-species corpora were streamed. "
                     "NOT A COMPLETE RESULT; do not cite.\n" % (LIMIT, EXPECT_SPECIES))
        fh.write(banner if banner.endswith("\n") else banner + "\n")
        fh.write(header)
        n = 0
        for row in body:
            vals = [("" if v is None else str(v)) for v in row]
            if len(vals) != len(cols):
                raise SystemExit("%s: %d fields, header has %d: %r"
                                 % (os.path.basename(path), len(vals), len(cols), vals))
            for v in vals:
                if FRACTION.match(v):
                    raise SystemExit(
                        "REFUSED: %s holds %r, a k/n fraction. That is a support "
                        "statistic in disguise." % (os.path.basename(path), v))
            fh.write("\t".join(vals) + "\n")
            n += 1
    say("  wrote %-34s %8d rows" % (os.path.basename(path), n))
    return n


# ------------------ step 1a  orthology-arm calls (OrthoFinder and SonicParanoid)
# The joint selection. The earlier per-target FLOOR/SELECT file is superseded: it
# chose each orthology method in isolation under a per-method floor on control
# recovery, which no longer applies now that the floor is tested at the union.
ARM1_SEL = os.environ.get("ARM1_SEL") or (
    _r59(ROOT + "/07_findings/R59_arm1_selected.tsv", must_exist=False))
PW_MEMBERS = ROOT + "/07_findings/sonic_pairwise_members.tsv"


def arm1_selection():
    """Map (method, target) to its setting, membership rule and source, from the orthology-arm selection file.

    The clustering setting of the orthology methods is chosen per target, because
    the floor on control recovery is written per target: control recovery for a
    target must be at least the recovery at the current setting. The eight
    targets are independent genes that already carry their own controls, decoys,
    HMM models and bit cuts; a single global clustering setting would make every
    target pay for the worst-behaved one. The criterion itself is fixed and was
    declared before any clustering existed; only its application is per target,
    and `rl_arm1_selection.py` reports the full scored table, including every
    setting that lost and why.

    Returns {} when no selection file exists, in which case both builders fall
    back to the single declared raw source and behave exactly as before. That
    keeps this script runnable on its own rather than silently depending on
    another script's output.
    """
    if not os.path.exists(ARM1_SEL):
        say("  NOTE  no per-target Arm 1 selection on disk; using the single "
            "declared RAW source for every target")
        return {}
    # The two orthology methods are chosen per target and independently of each
    # other, and that is checked here again rather than trusted.
    # `rl_arm1_selection.py` checks it when it writes this file, but a reader that assumes the writer was
    # correct is how a broken rule reaches a published table. The file carries
    # CUR_CONTROLS and CUR_SPECIES so that the floor test can be reproduced from
    # the file alone.
    NAMESPACE = {"orthofinder_full": "of-graph-", "sonic": "pairwise-"}
    sel = {}
    arm1_cover, arm1_species = {}, {}
    for r in rows(ARM1_SEL):
        m = (r.get("METHOD") or "").strip()
        t = (r.get("TARGET") or "").strip()
        if not m or t not in TARGETS:
            continue
        m = "orthofinder_full" if m == "orthofinder" else m
        if (m, t) in sel:
            raise SystemExit("R53: %s/%s appears twice in %s" % (m, t, ARM1_SEL))
        setting = (r.get("SETTING") or "").strip()
        for other, pfx in NAMESPACE.items():
            if other != m and setting.startswith(pfx):
                raise SystemExit(
                    "R53: %s/%s selected %r, which is %s's namespace"
                    % (m, t, setting, other))
        # The floor on control recovery is now tested at the union, not per
        # method. The old check rejected any cell recovering fewer controls than
        # its predecessor. Under the joint choice that is exactly what a correct
        # cell looks like: OrthoFinder for TOR covers 2 of 30 controls at zero
        # non-targets, because SonicParanoid covers 29 and the domain methods the
        # rest, against 30 of 30 at 1,593 non-targets before. Applying the old
        # test here would reject the better setting. What still matters, reported
        # below and tested at the union, is that species are not given up to make
        # the calls cheaper.
        try:
            arm1_cover[(m, t)] = set(
                x.strip() for x in (r.get("CONTROLS_COVERED_LIST") or "").split(";")
                if x.strip())
            arm1_species[(m, t)] = int(r.get("SPECIES_COVERED") or 0)
        except ValueError:
            arm1_species[(m, t)] = 0
        sel[(m, t)] = dict(setting=setting,
                           rule=(r.get("MEMBERSHIP_RULE") or "").strip(),
                           source=(r.get("SOURCE") or "").strip(),
                           not_informative=(r.get("NOT_INFORMATIVE_R55")
                                            or "").strip())
    # A cell where optimisation could not improve the method, and where the
    # unoptimised default calls more declared negatives than positives, makes no
    # call for that target. The cell is reported with its reason (it is not
    # deleted), but its candidate set is not carried into the target's table. The
    # absence of a call from one method is not evidence against the target, and
    # it does not strengthen the other method's call.
    mute = {k: v["not_informative"] for k, v in sel.items() if v["not_informative"]}
    for (m, t), why in sorted(mute.items()):
        say("  R55  %s makes NO CALL for %s -- %s" % (m, t, why))
    per_method = collections.Counter(m for m, _t in sel)
    for m, n in per_method.items():
        if n != len(TARGETS):
            raise SystemExit("R53: %s has %d of %d target decisions in %s"
                             % (m, n, len(TARGETS), ARM1_SEL))
    say("  R53 re-checked at the consumer: %d decisions over %d methods, "
        "namespaces clean" % (len(sel), len(per_method)))
    say("  R54 is a UNION test now: per-instrument control coverage is EXPECTED "
        "to be partial. Per-cell coverage: %s"
        % ", ".join("%s/%s=%d sp" % (m, t, arm1_species.get((m, t), 0))
                    for (m, t) in sorted(sel)))
    return sel


RULE_MEMBER_FILES = [
    _r59(ROOT + "/07_findings/R59_arm1_members.tsv", must_exist=False),
    ROOT + "/07_findings/sonic_pairwise_members.tsv",
    ROOT + "/07_findings/of_graph_members.tsv",
    ROOT + "/07_findings/arm1_threshold_sweep_members.tsv",
]


def pairwise_members():
    """Map (target, rule) to its member accessions, for every membership rule of either orthology method.

    Both methods now have rules read below their own clustering: SonicParanoid's
    pairwise co-orthologs, OrthoFinder's similarity-graph thresholds, and the
    exhaustive sweep over both. Reading only the SonicParanoid table was the
    fault that broke the first full run: three OrthoFinder targets had selected
    graph rules, found no members here, and went on to open the SOURCE string as
    a file path.
    """
    out = collections.defaultdict(set)
    for path in RULE_MEMBER_FILES:
        if not os.path.exists(path):
            continue
        for r in rows(path):
            out[((r.get("TARGET") or "").strip(),
                 (r.get("RULE") or "").strip())].add(
                     norm((r.get("PROTEIN_ID") or "").strip()))
    return out


def build_orthofinder(role, binding, spidx, sel=None):
    """OrthoFinder calls. MEMBER means sharing an orthogroup with at least one declared
    control in use or, for a target on a graph rule, being a member under that rule.

    OG0000000 holds 15,239 proteins, and `EVIDENCE_MAP.tsv` describes it as nearly
    uninformative for kinases. The group size is carried on every row, so a call
    that rests on it is visible. Nothing is filtered on size, which would invent
    a threshold.
    """
    sel = sel or {}
    pw = pairwise_members()
    out, seen = [], collections.defaultdict(set)

    # ---- targets on a graph rule ---------------------------------------------
    # OrthoFinder membership can be read below the clustering, from its own
    # OrthoFinder_graph.txt, at connectivity thresholds taken from the data or at
    # a swept threshold. Those settings name a member table rather than a
    # clustering file.
    rule_targets = []
    for t in TARGETS:
        e = sel.get(("orthofinder_full", t)) or {}
        if e.get("not_informative"):          # the method makes no call for this target
            continue
        lab = e.get("setting", "")
        if lab.startswith("of-graph-"):
            members = pw.get((t, lab), set())
            if not members:
                raise SystemExit(
                    "orthofinder: %s selected rule %s but no members were found "
                    "in %s" % (t, lab, ", ".join(RULE_MEMBER_FILES)))
            # A rule's member table lists the partners of the controls, so the
            # controls themselves are absent by construction: they are the seeds
            # of the query, not candidates. The home group of a clustering, by
            # contrast, contains its controls. Adding the controls back makes both
            # kinds of setting mean the same thing: the target's candidate set
            # includes its own known orthologs. Without this, the control-recovery
            # check before writing correctly reported a broken join (0 controls
            # for RAPTOR and AMPKβ). The selection measures are unaffected:
            # CONTROLS_RECOVERED remains the leave-one-out measure, a control
            # linked to at least one other control.
            members = members | binding[t]
            say("  orthofinder: %s uses %s -- %d members (incl. %d controls)"
                % (t, lab, len(members), len(binding[t])))
            rule_targets.append(t)
            for a in sorted(members):
                seen[t].add(a)
                out.append((t, ARM["orthofinder_full"], "orthofinder_full",
                            spidx.get(a, "UNMAPPED"), a,
                            role.get((t, a), "CANDIDATE"), rl_vocab.MEMBER,
                            "genome_wide_210",
                            "not nominated by orthofinder_full under this rule",
                            lab, len(binding[t]), len(members),
                            "orthofinder_210:OrthoFinder_graph.txt|rule=%s" % lab))

    # ---- targets on a clustering, grouped so that each file is read once ----
    by_src = collections.defaultdict(list)
    for t in TARGETS:
        e = sel.get(("orthofinder_full", t)) or {}
        if e.get("not_informative") or t in rule_targets:
            continue
        src = e.get("source") or SRC["orthofinder_full"]
        if src.startswith("MEMBERS:"):
            raise SystemExit("orthofinder: %s names a member table (%s) but its "
                             "setting %r is not a graph rule"
                             % (t, src, e.get("setting")))
        by_src[src].append(t)
    for src, tgts in sorted(by_src.items()):
        if not os.path.exists(src):
            raise SystemExit("orthofinder source missing for %s: %s"
                             % (",".join(tgts), src))
        groups = []
        for line in open(src, errors="replace"):
            og, _, rest = line.partition(":")
            members = [norm(a) for a in rest.split()]
            if members:
                groups.append((og.strip(), members))
        setting = (sel.get(("orthofinder_full", tgts[0])) or {}).get(
            "setting", "declared RAW")
        say("  orthofinder: %d orthogroups for %s  [%s]"
            % (len(groups), ",".join(tgts), setting))
        for og, members in groups:
            mset = set(members)
            for t in tgts:
                hit = mset & binding[t]
                if not hit:
                    continue
                for a in members:
                    if a in seen[t]:
                        continue
                    seen[t].add(a)
                    out.append((t, ARM["orthofinder_full"], "orthofinder_full",
                                spidx.get(a, "UNMAPPED"), a,
                                role.get((t, a), "CANDIDATE"), rl_vocab.MEMBER,
                                "genome_wide_210",
                                "not nominated by orthofinder_full under this rule",
                                og, len(hit), len(members),
                                "orthofinder_210:Orthogroups.txt|setting=%s"
                                % setting))
    return out


def build_sonic(role, binding, spidx=None, sel=None):
    """SonicParanoid calls: MEMBER under whichever SonicParanoid rule the target was assigned.

    Two rules are possible, and which one a target uses is decided by the
    selection test fixed in advance in `rl_arm1_selection.py`, not here:

    control-sharing-group
        shares an ortholog group with at least one declared control in use,
        judged leave-one-out, so a control alone in its group cannot recover
        itself.

    pairwise
        stands in a direct pairwise ortholog relation with at least one control
        in use in SonicParanoid's own species_to_species_orthologs/ tables, at or
        above the confidence the rule names. MCL discards that confidence when it
        collapses the graph into flat groups; the tables keep it. This is
        SonicParanoid's own raw output read at finer resolution, and no other
        method is consulted.

    Why a target may use the second: the AMPKα group has 15,148 proteins holding
    41 of 42 controls and 341 declared decoys, whereas the pairwise rule at
    confidence 1.0 gives 496 proteins, 42 controls and 11 decoys. Both are
    SonicParanoid's own output; neither is a vote and neither consults the
    domain methods.
    """
    sel = sel or {}
    pw = pairwise_members()
    out, seen = [], collections.defaultdict(set)

    # ---- targets on the pairwise rule ---------------------------------------
    muted = {t for t in TARGETS
             if (sel.get(("sonic", t)) or {}).get("not_informative")}   # no call for this target
    pw_targets = [t for t in TARGETS
                  if t not in muted
                  and (sel.get(("sonic", t)) or {}).get("rule") == "pairwise"]
    for t in pw_targets:
        rule = sel[("sonic", t)]["setting"]
        members = pw.get((t, rule), set())
        if not members:
            raise SystemExit("sonic: %s selected rule %s but no members found in %s"
                             % (t, rule, PW_MEMBERS))
        # A rule's member table lists the partners of the controls, so the
        # controls themselves are absent by construction: they are the seeds
        # of the query, not candidates. The home group of a clustering, by
        # contrast, contains its controls. Adding the controls back makes both
        # kinds of setting mean the same thing: the target's candidate set
        # includes its own known orthologs. Without this, the control-recovery
        # check before writing correctly reported a broken join (0 controls for
        # RAPTOR and AMPKβ). The selection measures are unaffected:
        # CONTROLS_RECOVERED remains the leave-one-out measure, a control linked
        # to at least one other control.
        members = members | binding[t]
        say("  sonic: %s uses %s -- %d members (incl. %d controls)"
            % (t, rule, len(members), len(binding[t])))
        for a in sorted(members):
            seen[t].add(a)
            out.append((t, ARM["sonic"], "sonic",
                        (spidx or {}).get(a, "UNMAPPED"), a,
                        role.get((t, a), "CANDIDATE"), rl_vocab.MEMBER,
                        "genome_wide_210",
                        "not nominated by sonic under this rule",
                        rule, len(binding[t]), len(members),
                        "sonicparanoid_210:species_to_species_orthologs|rule=%s"
                        % rule))

    # ---- targets on the group rule, grouped by their assigned clustering ----
    by_src = collections.defaultdict(list)
    for t in TARGETS:
        if t in pw_targets or t in muted:
            continue
        e = sel.get(("sonic", t))
        by_src[(e or {}).get("source") or SRC["sonic"]].append(t)

    for src, tgts in sorted(by_src.items()):
        if not os.path.exists(src):
            raise SystemExit("sonic source missing for %s: %s"
                             % (",".join(tgts), src))
        setting = (sel.get(("sonic", tgts[0])) or {}).get("setting", "declared RAW")
        ngroups = 0
        with open(src, errors="replace") as fh:
            rdr = csv.reader(fh, delimiter="\t")
            hdr = next(rdr)
            species = [h[:-4] if h.endswith(".faa") else h for h in hdr[1:]]
            for rec in rdr:
                if not rec:
                    continue
                ngroups += 1
                gid = rec[0]
                members = []
                for sp, cell in zip(species, rec[1:]):
                    cell = (cell or "").strip()
                    if not cell or cell == "*":
                        continue
                    for a in cell.split(","):
                        a = a.strip()
                        if a and a != "*":
                            members.append((sp, norm(a)))
                if not members:
                    continue
                accs = {a for _, a in members}
                for t in tgts:
                    hit = accs & binding[t]
                    if not hit:
                        continue
                    for sp, a in members:
                        # Leave-one-out: a control alone in its group cannot
                        # recover itself.
                        if a in hit and len(hit) == 1:
                            continue
                        if a in seen[t]:
                            continue
                        seen[t].add(a)
                        out.append((t, ARM["sonic"], "sonic", sp, a,
                                    role.get((t, a), "CANDIDATE"),
                                    rl_vocab.MEMBER, "genome_wide_210",
                                    "not nominated by sonic under this rule",
                                    gid, len(hit), len(members),
                                    "sonicparanoid_210:flat.ortholog_groups.tsv"
                                    "|setting=%s" % setting))
        say("  sonic: %d groups for %s  [%s]"
            % (ngroups, ",".join(tgts), setting))
    return out


# ------------------------------ step 1b  domain arm: CDD and InterProScan calls
def build_domains_from_combined(role, keep):
    """One pass over the 210 raw combined files, forward leg (leg=fwd) only.

    CDD           MEMBER = a kept CDD signature at or above its declared bit cut.
    InterProScan  MEMBER = a kept IPS:* signature that passes its declared threshold
                  (an e-value, or a score for score-based member databases).

    The `dec` leg is the decoy leg used to estimate false-positive rates: 166,043
    rows for Homo_sapiens alone. Reading it would turn shuffled sequences into calls.
    """
    files = corpus(SRC["combined"], "*.csv.gz")
    cdd_hit = collections.defaultdict(dict)      # target -> accession -> (signature, value, threshold, kind)
    ips_hit = collections.defaultdict(dict)
    spof = {}
    n_dec = 0
    for i, fn in enumerate(files, 1):
        with gzip.open(fn, "rt", errors="replace") as fh:
            for r in csv.DictReader(fh):
                if (r.get("leg") or "") != "fwd":
                    n_dec += 1
                    continue
                inst = r.get("instrument") or ""
                sig = r.get("signature") or ""
                acc = norm(r.get("protein") or "")
                sp = r.get("species") or ""
                for t in TARGETS:
                    spec = keep[t].get(inst)
                    if not spec or sig not in spec:
                        continue
                    kind, thr = spec[sig]
                    # The direction of the test is a property of the score, not
                    # of the database. `bits` and `score` are higher-is-better and
                    # sit in the `score` column; `evalue` is lower-is-better.
                    # Treating a ProSiteProfiles `score` as an e-value would
                    # invert its test.
                    try:
                        val = float(r.get("score") if kind in ("bits", "score")
                                    else r.get("evalue"))
                    except (TypeError, ValueError):
                        continue
                    ok = ((val >= thr) if kind in ("bits", "score")
                          else (val <= thr))
                    if not ok:
                        continue
                    spof[acc] = sp
                    tgt = cdd_hit if inst == "CDD" else ips_hit
                    prev = tgt[t].get(acc)
                    better = prev is None or (
                        val > prev[1] if kind == "bits" else val < prev[1])
                    if better:
                        tgt[t][acc] = (sig, val, thr, kind)
        if i % 40 == 0:
            say("    combined %d/%d" % (i, len(files)))
    say("  combined: skipped %d decoy-leg (leg!=fwd) rows" % n_dec)

    def rowsfor(store, method, src):
        out = []
        for t in TARGETS:
            for acc, (sig, val, thr, kind) in sorted(store[t].items()):
                out.append((t, ARM[method], method, spof.get(acc, "UNMAPPED"),
                            acc, role.get((t, acc), "CANDIDATE"), rl_vocab.MEMBER,
                            "genome_wide_210",
                            "not nominated by %s under this rule" % method,
                            sig, val, "%s>=%s" % (kind, thr)
                            if kind in ("bits", "score")
                            else "%s<=%s" % (kind, thr), src))
        return out

    return (rowsfor(cdd_hit, "cdd", "interproscan_combined:CDD"),
            rowsfor(ips_hit, "interpro", "interproscan_combined:IPS"))


# ----------------------------------------------------------- step 1c  HMMER
def parse_tblout(path, wanted):
    """Yield (protein, model, sequence bits, domain bits) for the models in `wanted`.

    tblout columns are separated by white space: target (0), accession (1),
    query (2), accession (3), then the full-sequence E-value (4), score (5) and
    bias (6), and the best-domain E-value (7), score (8) and bias (9). The bit
    scores are fields 5 and 8, matching the project's convention that both the
    sequence score and the domain score must pass.
    """
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 9 or f[2] not in wanted:
                continue
            try:
                yield f[0], f[2], float(f[5]), float(f[8])
            except ValueError:
                continue


def build_hmmer(role, cuts, skipped):
    """HMMER calls. MEMBER means a calling model clears both its sequence cut and its
    domain cut, as set in `cuts` (the jointly chosen setting where one exists,
    otherwise the stricter of the published and zero-decoy cuts). Requiring both
    follows the structure of the Pfam gathering-threshold pair, so no new rule is
    introduced.
    """
    wanted = {m for t in TARGETS for m, _, _, _, _, _ in cuts[t]}
    if not wanted:
        raise SystemExit("no HMMER model carries a published cut -- refusing")
    best = collections.defaultdict(dict)          # target -> accession -> best hit (model, scores, cuts, source, species, recall, FPR)
    files = hmmer_corpus()
    for i, fn in enumerate(files, 1):
        sp = os.path.basename(fn).split(".")[0]
        for prot, model, sq, dm in parse_tblout(fn, wanted):
            acc = norm(prot)
            for t in TARGETS:
                for m, cs, cd, rec, fpr, src in cuts[t]:
                    if m != model or sq < cs or dm < cd:
                        continue
                    prev = best[t].get(acc)
                    if prev is None or sq > prev[1]:
                        best[t][acc] = (m, sq, dm, cs, cd, src, sp, rec, fpr)
        if i % 40 == 0:
            say("    tblout %d/%d" % (i, len(files)))

    out = []
    for t in TARGETS:
        for acc, (m, sq, dm, cs, cd, src, sp, rec, fpr) in sorted(best[t].items()):
            out.append((t, ARM["hmmer"], "hmmer", sp, acc,
                        role.get((t, acc), "CANDIDATE"), rl_vocab.MEMBER,
                        "genome_wide_210",
                        "not nominated by hmmer under this rule",
                        m, "seq=%.1f/%.1f" % (sq, cs), "dom=%.1f/%.1f" % (dm, cd),
                        "hmmer_genomewide:tblout|cut=%s" % (src or "published")))
    for t in TARGETS:
        for m, why in skipped[t][:200]:
            LOG.append("    hmmer %-7s model not calling: %-46s %s" % (t, m, why))
    return out


# ------------------------------------------------- step 3  candidate list
def candidate_list(percall):
    """The union of every MEMBER call: a working list, not evidence.

    This set has no biological interpretation. It names the proteins that enter
    the contaminant screen (since 8 September 2026 `rl_phylo_prepare.py` also
    builds the gene-tree tip set from it). It carries no column saying which
    methods called a protein: a ';'-joined CALLED_BY field would hold the same
    information but sit one split() away from a support count, and people
    naturally start counting. That information goes to a separate long table.
    """
    cand = set()
    prov = []
    for method, calls in percall.items():
        for r in calls:
            if r[6] != rl_vocab.MEMBER:
                continue
            key = (r[0], r[3], r[4])              # TARGET, SPECIES, PROTEIN_ID
            cand.add(key)
            prov.append((r[0], r[3], r[4], method))
    return sorted(cand), sorted(set(prov))


# ------------------------------------------------ step 4  contaminant screen
def screen(cand, tgt_models, con_models, role=None):
    """Winner-take-all screen, unchanged: CONTAMINANT when the best contaminant-side
    model outscores every target-side model.

    A declared control is never called a contaminant. The declared-control table
    is the authors' statement that the protein is the target, and a later
    bit-score comparison is not entitled to overturn it. Measured case: the AMPKα
    control XP_002676974.1 (Naegleria gruberi, NAEGRDRAFT_2587,
    HIGH_CONFIDENCE_3LAYER) was once removed because LECA_OG048_BRSK1 scored 474.6
    bits against 462.1 for AMPKa_controls, a 12.5-bit margin on a protein that the
    table itself marks TRUNCATED_NOT_TESTABLE, a partial fragment of 253 residues.
    A truncated fragment scoring slightly higher on a model of a sister kinase is
    what a fragment does, not evidence that the table is wrong.

    The verdict is still computed for controls and still written to
    07_contaminant_screen.tsv, as CONTROL_EXEMPT with both scores, so the
    disagreement stays visible. The only change is that the row is not removed.
    The exemption covers the proteins of the declared-control table (in use or
    held out) and nothing else: no candidate gains protection, and nomination
    cannot turn a protein into a control. A test based on presence alone is
    wrong: the corrected earlier screen notes that it returned "TOR 374/374, a
    self-hit, not a screen".

    A verdict is a property of the protein: it depends only on that protein's own
    best target-side and best contaminant-side scores, never on which or how many
    methods nominated it. Screening the union and screening each method's set
    separately therefore give identical verdicts, which Verification 6 checks.

    The universe() of rl_phase2_panelv2_screen.py is deliberately not reproduced:
    it pools MEMBER calls across nine layers spanning both arms plus the gene
    tree, `armB` and KinOrtho, and reads derived layer files to do so.
    """
    want = {(t, sp, acc) for t, sp, acc in cand}
    by_target = collections.defaultdict(set)
    for t, sp, acc in cand:
        by_target[t].add(acc)
    allmodels = set()
    for t in TARGETS:
        allmodels |= tgt_models[t] | con_models[t]

    best = {}                       # (target, accession) -> [target bits, model, contaminant bits, model]
    files = hmmer_corpus()
    for i, fn in enumerate(files, 1):
        for prot, model, sq, dm in parse_tblout(fn, allmodels):
            acc = norm(prot)
            for t in TARGETS:
                if acc not in by_target[t]:
                    continue
                cur = best.setdefault((t, acc), [NEG_INF, "", NEG_INF, ""])
                if model in tgt_models[t] and sq > cur[0]:
                    cur[0], cur[1] = sq, model
                elif model in con_models[t] and sq > cur[2]:
                    cur[2], cur[3] = sq, model
        if i % 40 == 0:
            say("    screen tblout %d/%d" % (i, len(files)))

    out = {}
    n_exempt = []
    for (t, sp, acc) in cand:
        own, ownm, con, conm = best.get((t, acc), [NEG_INF, "", NEG_INF, ""])
        if own == NEG_INF and con == NEG_INF:
            call = "NO_SCORES"
        elif con > NEG_INF and con > own:
            call = "CONTAMINANT"
        else:
            call = "CLEAN"
        if call == "CONTAMINANT" and role is not None \
                and role.get((t, acc)) in ("CONTROL", "CONTROL_HELDOUT"):
            call = "CONTROL_EXEMPT"
            n_exempt.append((t, sp, acc, ownm, own, conm, con))
        out[(t, sp, acc)] = (call, ownm, own, conm, con)
    if n_exempt:
        say("  screen: %d DECLARED CONTROL(s) scored higher on a contaminant-pole "
            "model and are EXEMPT, not removed:" % len(n_exempt))
        for t, sp, acc, om, ov, cm, cv in n_exempt:
            say("    %-7s %-24s %-18s target %s %.1f vs contaminant %s %.1f"
                % (t, sp, acc, om, ov, cm, cv))
            LOG.append("    CONTROL_EXEMPT %s %s %s: target %s %.1f, "
                       "contaminant %s %.1f" % (t, sp, acc, om, ov, cm, cv))
    return out


# ------------------------------------------------------ verification helpers
def verify_screen_equivalence(percall, verdict):
    """Verification 6: every MEMBER row of every method must have a verdict from the
    screen of the union. Because a verdict depends only on the protein, this shows
    that the union is a working list that decides nothing itself."""
    bad = 0
    for method, calls in percall.items():
        for r in calls:
            if r[6] != rl_vocab.MEMBER:
                continue
            key = (r[0], r[3], r[4])
            if key not in verdict:
                bad += 1
    if bad:
        raise SystemExit(
            "REFUSED: %d MEMBER rows have no screen verdict. The union is not a "
            "superset of every method's own set, so it is doing something other "
            "than bookkeeping." % bad)
    say("  Verification 6  screen equivalence: every method's MEMBER set is "
        "covered by the union verdict, identical by construction")


def disagreement_examples(percall):
    """Verification 8: representative examples only. No measure of agreement,
    overlap, concordance, disagreement or support is computed, not even as an
    observation: a count of how often methods differ is a support statistic in
    another form."""
    by = collections.defaultdict(set)
    for method, calls in percall.items():
        for r in calls:
            if r[6] == rl_vocab.MEMBER:
                by[(r[0], r[3], r[4])].add(method)
    ex_arm1, ex_arms = [], []
    for key, ms in by.items():
        if len(ex_arm1) < 5 and ("orthofinder_full" in ms) != ("sonic" in ms):
            ex_arm1.append((key, sorted(ms)))
        a1 = ms & {"orthofinder_full", "sonic"}
        a2 = ms & {"hmmer", "cdd", "interpro"}
        if len(ex_arms) < 5 and bool(a1) != bool(a2):
            ex_arms.append((key, sorted(ms)))
        if len(ex_arm1) >= 5 and len(ex_arms) >= 5:
            break
    say("  Verification 8  preserved disagreement, worked examples "
        "(no metric computed):")
    for label, ex in (("within Arm 1", ex_arm1), ("between arms", ex_arms)):
        for (t, sp, acc), ms in ex:
            say("    %-12s %-7s %-28s %-22s called by %s"
                % (label, t, sp, acc, ",".join(ms)))


def group_size_exposure(percall):
    """Verification 9: calls that rest on each orthology method's largest, nearly
    uninformative group, as a report on that method's own evidence rather than a
    comparison between methods. Rows are shown rather than counted, because a
    count invites tallying."""
    say("  Verification 9  MEMBER calls resting on the near-vacuous groups "
        "(rows shown, not tallied):")
    for method, big in (("orthofinder_full", "OG0000000"), ("sonic", "1")):
        shown = 0
        for r in percall.get(method, []):
            if str(r[9]) == big and shown < 6:
                say("    %-17s %-7s %-26s %-22s group=%s size=%s"
                    % (method, r[0], r[3], r[4], r[9], r[11]))
                shown += 1
        if not shown:
            say("    %-17s no MEMBER call rests on %s" % (method, big))


def control_group_concentration(percall, binding):
    """Where a target's declared controls sit, per orthology method.

    A report on each method's own evidence. It names the group, its size and how
    many controls it holds, so a reader can see why a candidate set has the size
    it has, without any corrective rule being applied.

    No maximum group size, orthogroup specificity score, minimum control density,
    reciprocal-support rule, or agreement with another method or arm is computed
    here or anywhere else. Each would be a new decision layer, which is what this
    design excludes. A method must report faithfully what it found; it is not
    required to be informative.
    """
    say("  Verification 9b  where the declared controls sit, per orthology "
        "method (reporting only, no rule applied):")
    lines = []
    for method in ("orthofinder_full", "sonic"):
        for t in TARGETS:
            groups = {}
            for r in percall.get(method, []):
                if r[0] != t or r[6] != rl_vocab.MEMBER:
                    continue
                gid, nctl, size = r[9], r[10], r[11]
                if gid == "":
                    continue
                groups[gid] = (int(nctl or 0), int(size or 0))
            if not groups:
                continue
            gid, (nctl, size) = max(groups.items(), key=lambda kv: kv[1][0])
            share = 100.0 * nctl / len(binding[t]) if binding[t] else 0.0
            lines.append((method, t, gid, size, nctl, len(binding[t]), share))
            say("    %-17s %-7s group=%-11s size=%-7d holds %d of %d binding "
                "controls (%.0f%%)"
                % (method, t, gid, size, nctl, len(binding[t]), share))
    return lines


def near_vacuous_statement(lines, percall):
    """The plain statement for every method and target, with no threshold anywhere.

    An earlier version limited this to `size >= 5000 and nctl == tot`. Both parts
    were wrong: 5000 was an invented cutoff, and the equality silently left out
    SonicParanoid's AMPKα group (19,867 proteins, the least informative group in
    the table) only because it holds 41 of 42 controls rather than 42. A
    selection rule here would be the same fault that Verification 10 exists to
    avoid, so every row is stated with its own numbers and the reader judges.
    """
    out = []
    for method, t, gid, size, nctl, tot, share in sorted(
            lines, key=lambda r: -r[3]):
        out.append("%-17s %-7s group %-11s holds %d of %d binding controls and "
                   "%d proteins; the control-sharing rule admits all %d."
                   % (method, t, gid, nctl, tot, size, size))

    # Two targets falling into the same group of one method is a factual
    # observation, not a threshold: that method returns identical candidate sets
    # for two different targets and therefore cannot separate them.
    shared = []
    for method in ("orthofinder_full", "sonic"):
        bygroup = collections.defaultdict(list)
        for m, t, gid, size, nctl, tot, share in lines:
            if m == method:
                bygroup[gid].append(t)
        for gid, ts in sorted(bygroup.items()):
            if len(ts) > 1:
                shared.append(
                    "%s resolves %s to the SAME group (%s). It returns identical "
                    "candidate sets for them and cannot separate them under the "
                    "control-sharing rule."
                    % (method, " and ".join(sorted(ts)), gid))

    say("  Where each target's controls sit, every row, no cutoff applied:")
    for s in out:
        say("    " + s)
    for s in shared:
        say("  FINDING  " + s)
    return out, shared


def copy_number(percall):
    """Verification 10: mean copy number per species, per method and target. Reported
    only; no filter, rule or decision uses these values."""
    say("  Verification 10  per-species mean copy number (reporting only):")
    for method in ("orthofinder_full", "sonic", "hmmer", "cdd", "interpro"):
        per = collections.Counter()
        sp = collections.defaultdict(set)
        for r in percall.get(method, []):
            if r[6] == rl_vocab.MEMBER:
                per[r[0]] += 1
                sp[r[0]].add(r[3])
        bits = []
        for t in TARGETS:
            n = per[t]
            s = len(sp[t]) or 1
            bits.append("%s %d (%.1f/sp)" % (t, n, float(n) / s))
        say("    %-17s %s" % (method, "  ".join(bits)))
    say("    Highlighted where values greatly exceed the published expectations "
        "in the plan: TOR (PI3_PI4_kinase matches all PIKKs and PI3K/PI4K) and "
        "AMPKg (CBS matches CLC, CorC, CNNM, IMPDH, cystathionine beta-synthase). "
        "Promiscuity there is for the contaminant screen to remove, not for a "
        "raised cut.")


# -------------------------------------------- step 2  CDD / LST8 abstention
def cdd_lst8_abstention(role, cand, r57_in_force=False):
    """CDD ran on LST8 but has no verdict, which is not the same as a gap.

    LST8 has 0 kept CDD signatures out of 117 candidates: 84 excluded and 33
    singletons. KOG0315 is excluded at 83.3 bits because it also matches APAF1,
    ATG16L1, ATG16L2, BTRC, DAW1 and hundreds more: LST8 is a WD40 β-propeller,
    and every WD40 KOG matches the whole WD40 repertoire.

    So CDD abstains, with its reason stated, rather than being dropped, and the
    raw uncalled hits are written alongside, so the evidence stays visible even
    though no call is made.

    The joint choice of settings can overturn this. The abstention rests on a
    premise, that no CDD signature is specific enough to call, and the joint
    solve tests that premise over CDD's whole parameter space. Under the
    superseded per-method rule, which made CDD reach all 25 LST8 controls alone,
    the answer was KOG0315 at 73.7 bits: looser than the 83.3 bits the curation
    had already rejected, calling 22,791 proteins, 108 per species. The curation
    was right and that rule was wrong. Under the joint choice, CDD is responsible
    only for the controls that the union does not obtain more cheaply elsewhere,
    so it can sit at a far stricter cut and still contribute.

    When a jointly chosen setting is in force for CDD on LST8 (the argument
    `r57_in_force`, despite its name), this function writes nothing: abstention
    rows (NO_LAYER) for the same method and target that also has MEMBER rows would
    put two contradictory verdicts in one file. The reason is still returned and
    printed, so the overturned premise stays visible.
    """
    status = collections.Counter()
    sigs = {}
    for r in rows(CALLERS):
        if (r.get("family") or "").strip() != "LST8":
            continue
        if (r.get("instrument") or "").strip() != "CDD":
            continue
        st = (r.get("status") or "").strip() or "unlabelled"
        status[st] += 1
        sigs[(r.get("signature") or "").strip()] = st
    reason = ("CDD ran; %d CDD signatures examined, %s, 0 kept -- no signature "
              "is specific enough to call, so CDD has no verdict for LST8"
              % (sum(status.values()),
                 ", ".join("%d %s" % (n, s) for s, n in sorted(status.items()))))
    if r57_in_force:
        overturned = (
            "CDD/LST8 abstention NOT APPLIED this run: R59 put a setting in "
            "force for CDD on LST8, so CDD has a verdict. The "
            "curated premise it rested on is recorded and superseded, not "
            "deleted -- premise was: %s" % reason)
        say("  " + overturned)
        LOG.append("  " + overturned)
        return [], [], overturned
    say("  CDD/LST8 abstention: %s" % reason)

    want = {acc for (t, sp, acc) in cand if t == "LST8"}
    say("  CDD/LST8 uncalled-hit scan over %d LST8 candidates" % len(want))

    abst = []
    for (t, sp, acc) in sorted(cand):
        if t != "LST8":
            continue
        abst.append((t, ARM["cdd"], "cdd", sp, acc,
                     role.get((t, acc), "CANDIDATE"), "NO_LAYER",
                     "genome_wide_210", reason, "", "", "",
                     "interproscan_combined:CDD"))

    raw_hits = []
    for i, fn in enumerate(corpus(SRC["combined"], "*.csv.gz"), 1):
        with gzip.open(fn, "rt", errors="replace") as fh:
            for r in csv.DictReader(fh):
                if (r.get("leg") or "") != "fwd" or (r.get("instrument") or "") != "CDD":
                    continue
                acc = norm(r.get("protein") or "")
                if acc not in want:
                    continue
                sig = r.get("signature") or ""
                raw_hits.append((r.get("species") or "", acc, sig,
                                 r.get("name") or "", r.get("score") or "",
                                 r.get("evalue") or "", r.get("start") or "",
                                 r.get("end") or "",
                                 sigs.get(sig, "not-an-LST8-candidate-signature"),
                                 "no"))
    return abst, sorted(set(raw_hits)), reason


# --------------------------------------------------------------------- main
def main():
    # Every method that uses decoys must use the same corrected decoy pool, and
    # the inputs fixed before the analysis are checked against their SHA-256
    # checksums before any call is made, so a table cannot be built on a silently
    # changed control table.
    os.makedirs(OUT, exist_ok=True)
    gate()

    say("\nSTEP 1  declared rules and shared indexes")
    role, binding, declname = load_controls()
    for t in TARGETS:
        if not binding[t]:
            raise SystemExit(
                "REFUSED: %s has no BINDING control. Arm 1's membership rule is "
                "control-sharing, so its output would be empty for a reason that "
                "is not biology." % t)
    spidx = species_index()
    keep = load_callers()
    cuts, skipped = published_cuts()
    # The joint settings are applied after the fallback cuts are built, so the run
    # report can state what each cell would have done and what it does instead.
    cuts, keep, r59ledger = apply_r59(cuts, keep)
    tgt_models, con_models = load_poles()

    say("\nSTEP 2  per-method membership rules applied to RAW")
    percall = {}
    arm1sel = arm1_selection()
    if arm1sel:
        say("  Arm 1 per-target selection loaded: %d (method,target) assignments"
            % len(arm1sel))
        for m in ("orthofinder_full", "sonic"):
            say("    %-16s %s" % (m, ", ".join(
                "%s=%s" % (t, arm1sel[(m, t)]["setting"])
                for t in TARGETS if (m, t) in arm1sel)))
    percall["orthofinder_full"] = build_orthofinder(role, binding, spidx, arm1sel)
    percall["sonic"] = build_sonic(role, binding, spidx, arm1sel)
    cdd_rows, ips_rows = build_domains_from_combined(role, keep)
    percall["cdd"] = cdd_rows
    percall["interpro"] = ips_rows
    percall["hmmer"] = build_hmmer(role, cuts, skipped)

    for method, calls in percall.items():
        rl_vocab.assert_known({r[6] for r in calls}, "%s calls" % method)

    # Verification 3: control recovery, checked before writing; a failure stops the run.
    say("\n  Verification 3  control recovery (blocking):")
    for method, calls in percall.items():
        got = collections.Counter()
        for r in calls:
            if r[5] in ("CONTROL", "CONTROL_HELDOUT") and r[6] == rl_vocab.MEMBER:
                got[r[0]] += 1
        # A cell that makes no call by decision (optimisation could not improve
        # the method there) recovers zero controls on purpose; that is not a
        # broken join.
        muted_here = {t for t in TARGETS
                      if (arm1sel.get((method, t)) or {}).get("not_informative")}
        zero = [t for t in TARGETS
                if binding[t] and got[t] == 0 and t not in muted_here]
        say("    %-17s %s" % (method,
                              " ".join("%s=%d" % (t, got[t]) for t in TARGETS)))
        if method in ("orthofinder_full", "sonic") and zero:
            raise JoinFailed(
                "%s recovered 0 controls for %s. A control-sharing rule that "
                "recovers no control is a broken join, and its zero says nothing "
                "about the data." % (method, ",".join(zero)))

    # The candidate list is built before the files are written, because CDD's
    # LST8 abstention rows are written against the LST8 candidate set. Those rows
    # are NO_LAYER, never MEMBER, so they cannot feed back into the list they
    # come from.
    say("\nSTEP 3  candidate list (bookkeeping) and provenance (separate)")
    cand, prov = candidate_list(percall)
    cdd_lst8_r57 = any(i == "cdd" and t == "LST8" and st == "R59"
                       for i, t, st, _s, _n in r59ledger)
    abst, lst8_raw, lst8_reason = cdd_lst8_abstention(role, cand, cdd_lst8_r57)
    percall["cdd"] = percall["cdd"] + abst
    rl_vocab.assert_known({r[6] for r in percall["cdd"]}, "cdd calls with abstention")

    say("\nSTEP 4  emit RAW-derived method calls")
    files = [("01_orthofinder_full_calls.tsv", "orthofinder_full",
              "MEMBER = shares an OrthoFinder orthogroup with >=1 BINDING declared "
              "control. EVIDENCE_1/2/3 are OG, OG_N_CONTROLS, OG_SIZE. OG0000000 "
              "holds 15239 proteins and is near-vacuous for kinases; size is "
              "carried on every row and is NOT filtered on."),
             ("02_sonic_calls.tsv", "sonic",
              "MEMBER = shares a SonicParanoid ortholog group with >=1 BINDING "
              "declared control, judged LEAVE-ONE-OUT so a control alone in its "
              "group cannot recover itself. EVIDENCE_1/2/3 are GROUP, "
              "GROUP_N_CONTROLS, GROUP_SIZE. The largest group holds 19867 "
              "proteins; size is carried and NOT filtered on."),
             ("03_hmmer_calls.tsv", "hmmer",
              "MEMBER = a non-self-built target-pole model clears BOTH a "
              "sequence cut and a domain cut. WHICH cuts is decided per target by "
              "R59 and is stated per target in the operating-point ledger in "
              "RUN_REPORT.md; the CUT columns on every row carry the cut actually "
              "applied. R59 chooses settings JOINTLY across the instruments so "
              "their UNION recovers every declared control at the lowest "
              "corpus-wide cost, so this cell may cover only PART of the target's "
              "controls -- that is correct under R59, and the ledger names how "
              "many it covers. Where no R59 row exists the cell FALLS BACK to the "
              "rule below and does not go silent. "
              "R52 fallback rule: each cut is the STRICTER of "
              "the model's PUBLISHED cut (Pfam GA pair, KOfam bit cut, member-DB "
              "curated e-value) and its ZERO-DECOY cut measured on that target's "
              "own declared decoys. Every such model calls and the calls are "
              "UNIONED: no model is chosen over another, no vote is taken. "
              "R19: a model calls only if it still recovers >=1 control at its own "
              "cut. R28: models built from the controls do not call. A model with "
              "NO published cut calls at its zero-decoy cut alone -- see the CUT "
              "column, where cut=NONE marks exactly those. LIMITATION, and it is "
              "load-bearing: a zero-decoy cut is certified against that target's "
              "declared decoys only (5 species) and is applied here to 210. Where "
              "a cut does not generalise, 03b_hmmer_model_yield.tsv reports it "
              "per model; SIN1 KOG3739 is the worst case in this run, calling "
              "7,043 times = 33.6 per species for a single-copy gene."),
             ("04_cdd_calls.tsv", "cdd",
              "MEMBER = a CDD signature clears its threshold in the RAW combined "
              "scan, leg=fwd only. WHICH signature and threshold is decided per "
              "target by R59 -- chosen jointly with the other instruments so "
              "their UNION recovers every declared control at the lowest "
              "corpus-wide cost -- and is stated in the operating-point ledger in "
              "RUN_REPORT.md. Targets with no R59 row FALL BACK to the "
              "callers.tsv status=kept set at declared thresholds. LST8 is the "
              "cell to read carefully: it has 0 kept CDD signatures of 117 "
              "candidates (84 excluded, 33 singleton) because every WD40 KOG "
              "admits the whole WD40 repertoire, so under the fallback rule CDD "
              "ABSTAINS on LST8 and 04b carries its uncalled hits. Whether that "
              "abstention applied THIS run is stated in the ledger and in 04b's "
              "own banner -- R59 overturns it if a setting is in force, and "
              "then CDD calls instead of abstaining."),
             ("05_interpro_calls.tsv", "interpro",
              "MEMBER = an IPS:* signature clears its threshold in the RAW "
              "combined scan, leg=fwd only. WHICH member database, signature and "
              "threshold is decided per target by R59, searching InterPro's own "
              "parameter space jointly with the other instruments; it is ONE "
              "parameterisation, never a union across "
              "member databases, which would be a composite instrument InterPro "
              "is not. Ledger in RUN_REPORT.md. Targets with no covering signature "
              "FALL BACK to the callers.tsv status=kept set. "
              "MobiDBLite, Coils, SignalP, TMHMM and Phobius are barred from "
              "calling: disorder and coiled-coil prediction are not families.")]

    banner_head = ("## RAW-DERIVED METHOD CALLS, not raw instrument output. "
                   "RAW output -> method-specific membership rule -> MEMBER list; "
                   "this file is the third object. Absence means NOT NOMINATED "
                   "UNDER THIS RULE, never 'the instrument did not see it' -- ask "
                   "rl_evidence.py against the RAW source for coverage. "
                   "Controls from %s. This method speaks only for itself: no "
                   "vote, no rank, no arm-level call, and multiplicity of "
                   "nomination is provenance, not evidence.\n## %s")
    for fname, method, rule in files:
        emit(os.path.join(OUT, fname), banner_head % (declname, rule),
             CALL_HDR, percall[method])

    emit(os.path.join(OUT, "04b_cdd_LST8_raw_hits.tsv"),
         ("## UNCALLED. Raw CDD hits on LST8 candidates, from the RAW combined "
          "scan, leg=fwd. " + lst8_reason + ". CALLED=no on every row: these are "
          "shown so the evidence stays visible, and NONE of them is a call. "
          "SIGNATURE_STATUS is that signature's status in callers.tsv."
          if lst8_raw else
          "## EMPTY BY DESIGN this run, and the emptiness is the message. " +
          lst8_reason + ". This file exists to show CDD's LST8 evidence when CDD "
          "declines to call; CDD did call, so its LST8 hits are in "
          "04_cdd_calls.tsv as MEMBER rows and duplicating them here would "
          "present the same hits twice under two incompatible verdicts."),
         "SPECIES\tPROTEIN_ID\tSIGNATURE\tSIGNATURE_NAME\tBITS\tEVALUE\t"
         "START\tEND\tSIGNATURE_STATUS\tCALLED\n", lst8_raw)

    say("\nSTEP 5  emit candidate list and provenance")
    emit(os.path.join(OUT, "06_candidate_list.tsv"),
         "## BOOKKEEPING ONLY, NOT EVIDENCE. This set has no biological "
         "interpretation. It names which proteins enter contaminant screening and "
         "is, since 2026-09-08, READ BY rl_phylo_prepare.py, which builds the "
         "phylogenetic tip set from it -- the earlier claim that it is read by "
         "nothing else was false. It carries NO provenance column on purpose: a "
         "joined CALLED_BY field sits one split() away from a support count. "
         "Provenance is in 06b, long form. Never cite this file as a result.",
         "TARGET\tSPECIES\tPROTEIN_ID\n", cand)
    emit(os.path.join(OUT, "06b_candidate_provenance.tsv"),
         "## Provenance, long form, one row per contributing method. "
         "Multiplicity of nomination is PROVENANCE, NOT EVIDENCE: a protein "
         "named by three methods is not stronger evidence than one named by a "
         "single method. Do not GROUP BY (TARGET, SPECIES, PROTEIN_ID) to "
         "produce a size -- that reconstructs the support count this design "
         "exists to prevent.",
         "TARGET\tSPECIES\tPROTEIN_ID\tMETHOD\n", prov)

    say("\nSTEP 6  contaminant screen, winner-take-all")
    verdict = screen(cand, tgt_models, con_models, role)
    emit(os.path.join(OUT, "07_contaminant_screen.tsv"),
         "## Contaminant screen on panel v2's contaminant pole, genome-wide. "
         "WINNER-TAKE-ALL: CONTAMINANT when the best contaminant-pole model "
         "outscores every target-pole model. Presence alone is wrong -- it "
         "returned TOR 374/374, a self-hit, not a screen. R28 excludes 3 "
         "misdeclared self-in-paralog models. The verdict is a property of the "
         "PROTEIN, never of who nominated it, so this screen decides nothing "
         "about which method or arm was correct.",
         "TARGET\tSPECIES\tPROTEIN_ID\tCALL\tBEST_TARGET_MODEL\tBEST_TARGET_BITS\t"
         "BEST_CONTAM_MODEL\tBEST_CONTAM_BITS\n",
         [(t, sp, a, v[0], v[1], "" if v[2] == NEG_INF else "%.1f" % v[2],
           v[3], "" if v[4] == NEG_INF else "%.1f" % v[4])
          for (t, sp, a), v in sorted(verdict.items())])

    say("\nSTEP 7  final lists")
    retained = [(t, sp, a) for (t, sp, a), v in sorted(verdict.items())
                if v[0] != "CONTAMINANT"]
    emit(os.path.join(OUT, "08_final_retained_orthologs.tsv"),
         "## Candidate list minus CONTAMINANT. NO_SCORES rows are RETAINED and "
         "flagged in 07, never silently dropped. No provenance column: join 06b "
         "deliberately if you need to know which methods called a protein.",
         "TARGET\tSPECIES\tPROTEIN_ID\n", retained)

    per_method_final = []
    for method, calls in percall.items():
        for r in calls:
            if r[6] != rl_vocab.MEMBER:
                continue
            v = verdict.get((r[0], r[3], r[4]))
            if v and v[0] == "CONTAMINANT":
                continue
            per_method_final.append((r[0], r[1], method, r[3], r[4], r[5],
                                     v[0] if v else "NO_SCORES"))
    emit(os.path.join(OUT, "09_final_lists_per_method.tsv"),
         "## Five lists stacked, never reconciled. Each method's own MEMBER set "
         "minus CONTAMINANT. Rows from different methods are NOT compared, "
         "counted across, or merged; disagreement between them is a result, not "
         "a problem to be fixed.",
         "TARGET\tARM\tMETHOD\tSPECIES\tPROTEIN_ID\tROLE\tSCREEN\n",
         sorted(per_method_final))

    say("\nSTEP 8  verification")
    verify_screen_equivalence(percall, verdict)
    disagreement_examples(percall)
    group_size_exposure(percall)
    conc = control_group_concentration(percall, binding)
    findings, shared = near_vacuous_statement(conc, percall)
    copy_number(percall)

    with open(os.path.join(OUT, "RUN_REPORT.md"), "w") as fh:
        fh.write("# RAW-derived per-method ortholog calls -- run report\n\n")
        fh.write("Multiplicity of nomination is provenance, not evidence. "
                 "No agreement, overlap, concordance or support metric appears "
                 "anywhere in this report or in any emitted file.\n\n")
        fh.write(
            "## Why CDD and InterPro are read from `concord/combined`\n\n"
            "Both `cdd_dense_scan` (cddfull/out/fwd) and `interproscan_combined`\n"
            "(concord/combined) are KIND=RAW in EVIDENCE_MAP.tsv, so the RAW-only\n"
            "rule does not choose between them. `combined` was chosen because it\n"
            "preserves the declared signature namespace:\n\n"
            "    callers.tsv identifiers -> direct match          (combined)\n"
            "    callers.tsv identifiers -> ID translation -> PSSM (cddfull)\n\n"
            "cddfull stores signatures as PSSM numeric ids (`gnl|CDD|270713`)\n"
            "while callers.tsv declares `cd05169`. Every translation step is an\n"
            "opportunity for silent loss, and a signature that fails to translate\n"
            "disappears without an error. Verified equivalent before choosing: on\n"
            "Homo_sapiens, cd05169 rows carry identical bits and e-values in both\n"
            "sources. This is recorded here because the reason is not recoverable\n"
            "from the output.\n\n"
            "## Why only leg=fwd is read\n\n"
            "`combined` carries a `dec` decoy leg alongside the forward biological\n"
            "proteomes. This is not a biological threshold; it restricts the scan\n"
            "to the declared biological leg. Reading `dec` would turn shuffled\n"
            "control sequences into candidates. The count skipped is logged below.\n\n"
            "## Verification 10 is observational\n\n"
            "Copy numbers are reported, never acted on. No rule, filter or cut is\n"
            "adjusted because a copy number looks wrong -- that would be\n"
            "calibration, and it is exactly what this design excludes.\n\n")
        fh.write("## Operating point of every Arm 2 cell (R59)\n\n"
                 "R59: settings are chosen JOINTLY so the UNION of all\n"
                 "instruments recovers every declared control at the lowest cost\n"
                 "measured across the 210 proteomes. A cell here may cover only\n"
                 "PART of its target's controls; the note column says how many.\n"
                 "That is the point of running five instruments, and it is what\n"
                 "the SUPERSEDED R57 threw away by requiring each instrument to\n"
                 "reach every control alone -- which drove AMPKa InterPro to the\n"
                 "generic kinase domain at 274 calls per species.\n\n"
                 "STATUS=R59_NO_SETTING means the joint solution covers this\n"
                 "target more cheaply without this instrument, so it makes no\n"
                 "call. STATUS=FALLBACK_R52 means no R59 row exists for the cell\n"
                 "and its prior declared rule still applies. Neither is silence\n"
                 "imposed on an instrument that had something to say.\n\n")
        fh.write("| instrument | target | status | setting | note |\n")
        fh.write("|---|---|---|---|---|\n")
        for inst, t, st, setting, note in r59ledger:
            fh.write("| %s | %s | %s | `%s` | %s |\n"
                     % (inst, t, st, setting, note))
        fh.write("\nArm 1 is NOT re-searched by R59. OrthoFinder and SonicParanoid have no\n"
                 "continuous score left in their output to solve on -- the merge\n"
                 "threshold acts upstream of MCL -- so Arm 1 is optimised by its own\n"
                 "exhaustive sweep under the pre-declared FLOOR/SELECT criterion in\n"
                 "07_findings/arm1_selected_per_target.tsv. FLOOR carries a species\n"
                 "term that R59 does not, and for sonic RICTOR and sonic AMPKb the\n"
                 "two criteria disagree: a swept setting recovers 19/19 and 33/33\n"
                 "respectively at zero decoys but covers fewer species (75 vs 95,\n"
                 "158 vs 178), so FLOOR rejects it. FLOOR is kept, because species\n"
                 "are not dropped to buy control recovery that the domain arm\n"
                 "already supplies. Both facts are recorded rather than reconciled.\n\n")
        if findings:
            fh.write("## Where each target's declared controls sit, Arm 1\n\n")
            fh.write("Every method/target row, largest group first. No cutoff is\n"
                     "applied and no row is selected: the numbers speak.\n\n```\n")
            for s in findings:
                fh.write(s + "\n")
            fh.write("```\n\n")
            if shared:
                fh.write("### Targets a method cannot separate\n\n")
                for s in shared:
                    fh.write("- " + s + "\n")
                fh.write("\n")
            fh.write(
                "\nThe output is CORRECT under the declared rule and the method is\n"
                "behaving to specification -- control recovery is exact. What is\n"
                "weak is the discriminatory power of this orthology formulation\n"
                "for these targets in these data, which is a scientific result,\n"
                "not a pipeline malfunction. The design requires each method to\n"
                "report honestly what it found; it does not require every method\n"
                "to be informative.\n\n"
                "No rescue is applied. A maximum group size, an orthogroup\n"
                "specificity score, a minimum control density, a reciprocal-support\n"
                "rule, or agreement with HMMER or the domain arm would each be a\n"
                "new adjudication layer. OG_SIZE and GROUP_SIZE on every row are\n"
                "what let this limitation be communicated without inventing one.\n\n"
                "The contaminant screen does not address this and is not expected\n"
                "to: it removes proteins that score higher against a contaminant\n"
                "pole than a target pole -- wrong families recruited by promiscuous\n"
                "domains -- and is not a general specificity filter. A low removal\n"
                "rate across a large orthology set is consistent behaviour, not a\n"
                "defect.\n\n")
        fh.write("```\n")
        fh.write("\n".join(LOG))
        fh.write("\n```\n")
    say("\nwrote %s/RUN_REPORT.md" % OUT)


if __name__ == "__main__":
    main()
