#!/usr/bin/env python3
# >>> deposit header
# Purpose: Choose jointly, for each target, one setting for each of the five candidate-discovery methods so that
#          together they recover every declared positive control, maximising the controls recovered minus the
#          declared negatives called, under a ceiling on calls per species.
# Methods: Candidate discovery
# Inputs:  the genome-wide HMMER tblout files (through `rl_evidence_roots`); the CDD and InterProScan scans in
#          `concord/combined`; `THRESHOLD_PILOT_PERTARGET.tsv` (the target-side models);
#          `07_findings/R59_arm1_settings__<tag>.tsv` and `R59_arm1_cover__<tag>.tsv` (the OrthoFinder and
#          SonicParanoid settings enumerated in the preceding step); the controls in use, declared decoys
#          and negative controls, read through `rl_arm1_selection`
# Outputs: `07_findings/R59_union_settings__<tag>.tsv` (the chosen setting per target and method),
#          `R59_frontier__<tag>.tsv` and `R59_frontier_cover__<tag>.tsv` (every allowed setting after Pareto
#          pruning, with the controls it covers) and `R59_union_optimum__<tag>.md` (summary table); <tag> is the
#          control set in `R59_TAG` plus the target-set suffix
# Tools:   Python 3 standard library and the project modules `rl_arm1_selection`, `rl_invariants`, `rl_targets`
#          and `rl_evidence_roots`; reads HMMER, CDD, InterProScan, OrthoFinder and SonicParanoid output without
#          running them; versions: see envs/
# <<< deposit header
"""The joint optimum of the search settings, chosen at the union of all methods.

    choose one setting s_i per method
    minimise   non-targets called across the 210 proteomes
    subject to union of cover(s_i) = every declared control

(That is the rule as first stated; the objective the code maximises is given
below.) The earlier rule required each method to reach every control alone. That
misreads the requirement on the union as a requirement on each method, and it
made every method loosen until it caught the single most divergent control by
itself: InterProScan for AMPKα ended at IPS:SMART SM00220, the generic kinase
domain, calling 57,468 proteins, 274 per species. This script uses the
flexibility instead: a method that is sharp for most of the control set runs at
a strict setting, and the hard controls are picked up by whichever method
reaches them most cheaply.

How settings are enumerated. The cost is monotone in the cut, so the only cuts
worth considering for a model or signature are those at a control's own score.
Each gives one (cover, cost) pair, where the cover is the set of controls that
also clear the cut. No threshold is invented anywhere: every candidate cut is a
score that some declared control actually has.

The cost, and the correction that matters most. An intermediate version defined
the cost as every protein called minus the declared controls, counted over all
210 proteomes. That inverts the objective: an undeclared true ortholog, the
thing this project exists to find, counts as cost, so the cheapest setting is
the one that finds nothing. It produced exactly that: seven cells of the
orthology methods came back at cost zero with zero members, covering two
controls each because those two controls link to each other and nothing else
clears the bar. The cost reached its minimum precisely when the method became
useless. The score is now

    score(setting) = |controls recovered| - |declared decoys and curated
                     negative controls called|

and the search maximises it.

Calling nothing scores zero, so silence is no longer optimal and does not have
to be forbidden by a separate rule. Coverage is part of the objective rather
than an added tie-break, and the two things a method can do, find a true
ortholog or call a look-alike, are weighed against each other directly.

The exchange rate is 1:1, and that is a declared position, not a default. One
decoy called cancels one control recovered. A setting recovering 30 controls
while calling 100 decoys scores -70 and loses to one recovering 2 controls with
no decoys, at +2. Calling a look-alike is held to be exactly as bad as missing a
true ortholog. NEG_WEIGHT below is that coefficient, written out so that it can
be changed deliberately rather than discovered later inside a formula.

A copy-number ceiling is applied here. COPY_CEILING = 5 calls per species limits
which settings are allowed: a setting above it is rejected however well it
scores, and the refusal at the CALL_CAP check below enforces it. The full
reasoning is in the comment above the COPY_CEILING assignment.

An earlier decision was that no ceiling applied, on the grounds that a cap
"would exclude every OrthoFinder clustering for TOR (5.2-7.6 calls/species) and
force back the near-empty settings this model exists to eliminate". That reason
was measured and does not hold: 732 of 743 enumerated TOR OrthoFinder settings
are under the cap, and across all 40 cells (one per target and method) every
cell keeps at least one allowed setting, so the cap empties nothing.

N_CALLED and CALLED_PER_SPECIES are still carried and reported.

    maximise   sum over methods of ( controls_i - NEG_WEIGHT x negatives_i )
    subject to union of cover = every declared control

The sum is over methods, not over the union: a method gets credit for reaching
on its own a control that another method already has. That credits each method
with what it finds, and it is what makes running five methods worth anything. It
also means that no method is pushed into silence to avoid duplicating work.

Two passes over the proteomes for each kind of domain search (HMMER, and CDD
with InterProScan): the first fixes the cuts from the controls and scores the
declared negatives, the second counts the total calls for the ceiling test.

The search is over methods, not arms. All five are optimised on the same
footing: OrthoFinder and SonicParanoid alongside HMMER, CDD and InterProScan. An
earlier version kept the orthology methods at their FLOOR/SELECT setting instead
of searching them, and the consequence was immediate: every target's remaining
cost came from the orthology methods, TOR's OrthoFinder cell alone contributing
1,593 of its 1,941 non-targets. Searched properly, OrthoFinder reaches full
control cover for TOR at a cost of 184.

Every setting on disk is a candidate. For OrthoFinder: the clusterings at
inflation 1.2, 1.5, 2.0 and 3.0, the original orthofinder_full clustering, and
every edge-weight threshold in its own weighted graph. For SonicParanoid: the
cov0.75/I1.5 and cov0.90 runs and the inflation 2.0 and 3.0 runs, and every
confidence threshold in its own pairwise tables. `rl_arm1_all_settings.py`
enumerates them with their cover and their cost over the proteomes; 10,082 were
found.

This does not break the requirement that OrthoFinder and SonicParanoid are
optimised per target and independently of each other. The search selects each
method's own setting from that method's own settings; it never merges, weights,
scores or votes across methods, and no method's calls are altered by another's.
What is shared is only the constraint that between them the controls are
covered, which was always a statement about the union. The two arms remain
independent in the report; what is rejected here is letting the boundary between
the arms decide which methods get optimised.
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import bisect
import collections
import csv
import glob
import gzip
import math
import os
import sys

ROOT = os.environ.get("PROJECT_ROOT") or (
    (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1"))
sys.path.insert(0, os.path.join(ROOT, "01_scripts"))
sys.path.insert(0, os.path.join(ROOT, "01_scripts", "lib"))
import rl_arm1_selection as A                                    # noqa: E402
import rl_invariants                                             # noqa: E402
import rl_targets                                                # noqa: E402
import rl_evidence_roots                                         # noqa: E402

TBL = ROOT + "/03_evidence_by_method/hmmer_genomewide/tblout"
COMBINED = (_TA210_SCRATCH + "/concord/combined")
PILOT = ROOT + "/06_calibration_and_checks/THRESHOLD_PILOT_PERTARGET.tsv"
CALLERS = (_TA210_SCRATCH + "/concord/callers.tsv")

# --- Output names carry the control-set tag, so two control sets can be solved at once ---
# This script reads the output of `rl_arm1_all_settings.py`, so untagged parallel
# runs would not merely overwrite each other's files, they would feed each other:
# the union solve of one control set would read the orthology-method enumeration
# of the other.
_TAG = os.environ.get("R59_TAG") or ""


def _tagged(p):
    if not _TAG:
        raise SystemExit(
            "REFUSING: R59_TAG is not set. Writing untagged R59 outputs is how "
            "stale files accumulate: a later reader cannot tell which control "
            "set produced them. Set R59_TAG=setA or R59_TAG=setB.")
    return rl_targets.tagged(p, _TAG)           # base__<tag><set suffix>.ext
# ------------------------------------------------------------------------------
ARM1_SETTINGS = _tagged(ROOT + "/07_findings/R59_arm1_settings.tsv")
ARM1_COVER = _tagged(ROOT + "/07_findings/R59_arm1_cover.tsv")
F = ROOT + "/07_findings/"
OUT_SET = _tagged(F + "R59_union_settings.tsv")
OUT_FRONT = _tagged(F + "R59_frontier.tsv")
OUT_FCOV = _tagged(F + "R59_frontier_cover.tsv")
OUT_MD = _tagged(F + "R59_union_optimum.md")
TARGETS = A.TARGETS
N_SPECIES = 210
# All five methods, optimised equally. The order is only the order in which the
# dynamic programme processes them and carries no precedence: the objective is
# symmetric in the methods.
INSTRUMENTS = ("orthofinder_full", "sonic", "hmmer", "cdd", "interpro")
# The exchange rate between one control recovered and one declared negative
# called; 1.0 means they cancel one for one. It is the model's only free
# parameter, named here rather than left as a literal inside the formula.
NEG_WEIGHT = float(os.environ.get("NEG_WEIGHT") or 1.0)
# The copy-number ceiling. A setting that calls more than COPY_CEILING proteins per
# species on average (more than CALL_CAP = COPY_CEILING x 210 in total) is not
# allowed, however well it scores. An earlier decision applied no ceiling here, on
# the grounds that a cap "would exclude every OrthoFinder clustering for TOR
# (5.2-7.6 calls/species) and force back the near-empty settings this model exists
# to eliminate". That reason was measured and does not hold: 732 of 743 enumerated
# TOR OrthoFinder settings are under the cap, because the graph-edge enumeration
# supplies allowed alternatives that the discrete clusterings did not. Across all
# 40 cells (one per target and method), every cell keeps at least one allowed
# setting, so the cap empties nothing.
#
# The ceiling also brings the code into line with the rule as written in
# METHOD_RULES.md, "minimise |union of the called non-targets| subject to union of
# cover(s_i) = C", which the code had never implemented: it minimises declared
# negatives, so a protein that is neither a control nor a declared negative costs
# nothing. OrthoFinder (orthofinder_full) for AMPKγ at inflation 1.2 called 1,389
# proteins (6.61 per species) and was charged for 25.
#
# A hard cap rather than a weighted term for calls, deliberately: a one-sided
# upper bound cannot recreate the degenerate optimum already met once, when "a
# cost function counting every call made the cheapest setting the one that finds
# NOTHING". The objective still pays for coverage, so silence stays unattractive.
# The ceiling carries no free parameter to tune: it is the authors' biological
# prior, defined identically in `rl_r59_plots.py`, whose docstring calls settings
# above it "INADMISSIBLE however cheap they look".
COPY_CEILING = int(os.environ.get("COPY_CEILING") or 5)   # calls per species
CALL_CAP = COPY_CEILING * N_SPECIES
MAX_STATES = int(os.environ.get("R59_MAX_STATES") or 2000000)
# Cells emptied by the ceiling (decided 9 September 2026). A cell (one target and
# one method) in which every enumerated setting calls more than CALL_CAP has no
# allowed setting. The default, "refuse", keeps the original behaviour: the whole
# run stops, which is correct for the eight core targets (every cell has an
# allowed setting) and leaves their output byte-identical. "report" marks the cell
# as having no allowed setting, writes a row with `SETTING="(no admissible
# setting)"` and `SOURCE=r64_cap_emptied`, and lets the remaining methods carry the
# union. That is the measured CAMKK case: a functional category with no separable
# orthogroup at any clustering setting on disk. A cell with no candidate setting
# at all still stops the run in both modes.
CAP_EMPTIED = (os.environ.get("R59_CAP_EMPTIED") or "refuse").strip()
if CAP_EMPTIED not in ("refuse", "report"):
    raise SystemExit("REFUSING: R59_CAP_EMPTIED=%r is not refuse or report" % CAP_EMPTIED)
CAP_EMPTIED_SETTING = "(no admissible setting)"
EMPTIED = {}                       # (target, method) -> (smallest N_CALLED, number of settings)
NEG_INF = float("-inf")
PRESENCE = {"", "none", "presence"}
LOWER_BETTER = {"evalue"}
IPS_BARRED = {"IPS:MobiDBLite", "IPS:Coils", "IPS:SignalP", "IPS:TMHMM",
              "IPS:Phobius", "IPS:ProSitePatterns"}


def rows(p):
    return csv.DictReader((l for l in open(p, errors="replace")
                           if not l.startswith("##")), delimiter="\t")


# ---------------------------------------------------------------- enumeration
def settings_from_scores(scores, controls, lower=False):
    """Map each cover (a frozenset of controls) to its cut, for every cut at a control's own score.

    scores: {control -> value} for one model or signature, over the controls it
    reaches. A cut at control c passes exactly the controls whose score is at
    least as good as c's. Ties collapse: two controls with the same score give the
    same cover and the same cut, which is returned once.
    """
    out = {}
    for c, v in scores.items():
        if lower:
            cover = frozenset(x for x, w in scores.items() if w <= v)
        else:
            cover = frozenset(x for x, w in scores.items() if w >= v)
        cover &= controls
        if not cover:
            continue
        prev = out.get(cover)
        # Only tied scores give the same cover, and they give the same cut, so
        # one entry per cover is kept.
        if prev is None or (v > prev if lower else v < prev):
            out[cover] = v
    return out


def pareto(items):
    """Keep the settings that are not dominated: drop a setting when another covers a
    superset of its controls at no greater cost. Settings tied on both are
    collapsed to one."""
    items = sorted(items, key=lambda x: (x[0], -len(x[1])))    # ascending cost
    keep = []
    for cost, cover, meta in items:
        if any(cover <= k[1] and k[0] <= cost for k in keep):
            continue
        keep = [k for k in keep if not (k[1] <= cover and cost <= k[0])]
        keep.append((cost, cover, meta))
    return keep


# --------------------------------------------------------------------- HMMER
def hmmer_frontier(allc, neg):
    models = collections.defaultdict(set)
    # The pilot fixed in advance, then its companion for the target set when
    # another target set is in use; the default eight targets read exactly the
    # original file.
    for r in rl_targets.pilot_rows(PILOT):
        if (r.get("POLE") or "").strip() != "target":
            continue
        t = (r.get("FOR_TARGET") or "").strip()
        if t in TARGETS and (r.get("SELF_BUILT") or "").strip() == "no":
            models[t].add((r.get("MODEL") or "").strip())
    allmodels = set().union(*models.values())
    every_control = set().union(*allc.values())

    interest = set(every_control)
    for t in TARGETS:
        interest |= neg[t]
    print("HMMER pass 1: control AND declared-negative scores over %d models, "
          "%d proteins of interest" % (len(allmodels), len(interest)))
    sc = collections.defaultdict(dict)              # model -> accession -> (sequence bits, domain bits)
    # The declared evidence roots: the primary tblout set, plus a supplementary
    # root (models searched separately with identical options) only for another
    # target set.
    files = rl_evidence_roots.tblout_files()
    for n, f in enumerate(files):
        if n % 60 == 0:
            print("  %d/%d" % (n, len(files)), flush=True)
        with gzip.open(f, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 9 or p[2] not in allmodels:
                    continue
                acc = A.norm(p[0])
                if acc not in interest:
                    continue
                try:
                    b, dm = float(p[5]), float(p[8])
                except ValueError:
                    continue
                d = sc[p[2]]
                if b > d.get(acc, (NEG_INF, NEG_INF))[0]:
                    d[acc] = (b, dm)

    # Candidate cuts: (model, sequence cut, domain cut) taken at each control's own scores.
    cand = collections.defaultdict(list)            # (target, model) -> [(cover, sequence cut, domain cut)]
    cuts_by_model = collections.defaultdict(set)
    for t in TARGETS:
        for m in sorted(models[t]):
            d = sc.get(m, {})
            hit = {c: d[c] for c in allc[t] if c in d}
            if not hit:
                continue
            seen = {}
            for c, (vs, vd) in hit.items():
                cover = frozenset(x for x, (ws, wd) in hit.items()
                                  if ws >= vs and wd >= vd)
                cur = seen.get(cover)
                if cur is None or (vs, vd) > cur:      # strictest cut, same cover
                    seen[cover] = (vs, vd)
            for cover, (cs, cd) in seen.items():
                cand[(t, m)].append((cover, cs, cd))
                cuts_by_model[m].add((cs, cd))
    print("  %d (target, model) cells, %d distinct cuts to cost"
          % (len(cand), sum(len(v) for v in cuts_by_model.values())))

    print("HMMER pass 2: corpus cost of each cut")
    # Cuts sorted by sequence score, so a protein is tested only against those it can clear.
    sorted_cuts = {m: sorted(v) for m, v in cuts_by_model.items()}
    seqkeys = {m: [x[0] for x in v] for m, v in sorted_cuts.items()}
    count = collections.Counter()
    for n, f in enumerate(files):
        if n % 60 == 0:
            print("  %d/%d" % (n, len(files)), flush=True)
        seen = collections.defaultdict(set)
        with gzip.open(f, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 9 or p[2] not in sorted_cuts:
                    continue
                try:
                    b, dm = float(p[5]), float(p[8])
                except ValueError:
                    continue
                m = p[2]
                j = bisect.bisect_right(seqkeys[m], b)
                if not j:
                    continue
                acc = A.norm(p[0])
                for cs, cd in sorted_cuts[m][:j]:
                    if dm >= cd:
                        seen[(m, cs, cd)].add(acc)
        for k, accs in seen.items():
            count[k] += len(accs)

    front = collections.defaultdict(list)
    for (t, m), lst in cand.items():
        d = sc.get(m, {})
        for cover, cs, cd in lst:
            # Cost = the declared negatives called. They were scored in pass 1,
            # so this is exact and needs no second lookup.
            cost = sum(1 for a in neg[t]
                       if a in d and d[a][0] >= cs and d[a][1] >= cd)
            ncall = count[(m, cs, cd)]
            front[(t, "hmmer")].append(
                (cost, cover, ("hmmer", m, "%.1f" % cs, "%.1f" % cd, ncall)))
    return front


# ------------------------------------------------------------------- domains
def domain_frontier(allc, neg):
    every_control = set().union(*allc.values())
    interest = set(every_control)
    for t in TARGETS:
        interest |= neg[t]
    print("domains pass 1: control AND declared-negative scores, %d proteins"
          % len(interest))
    sc = collections.defaultdict(dict)
    kinds = {}
    files = sorted(glob.glob(COMBINED + "/*.csv.gz"))
    for n, f in enumerate(files):
        if n % 40 == 0:
            print("  %d/%d" % (n, len(files)), flush=True)
        with gzip.open(f, "rt", errors="replace") as fh:
            for row in csv.DictReader(fh):
                if (row.get("leg") or "") != "fwd":
                    continue
                acc = A.norm((row.get("protein") or "").strip())
                if acc not in interest:
                    continue
                inst = (row.get("instrument") or "").strip()
                if inst in IPS_BARRED:
                    continue
                kind = (row.get("score_kind") or "").strip()
                kinds[inst] = kind
                if kind in PRESENCE:
                    continue
                try:
                    v = float(row.get("score") or "nan")
                except ValueError:
                    continue
                if math.isnan(v):
                    continue
                k = (inst, (row.get("signature") or "").strip())
                cur = sc[k].get(acc)
                if cur is None or (v < cur if kind in LOWER_BETTER else v > cur):
                    sc[k][acc] = v

    cand = collections.defaultdict(list)
    cuts_by_sig = collections.defaultdict(set)
    for t in TARGETS:
        for (inst, sig), d in sc.items():
            lower = kinds.get(inst) in LOWER_BETTER
            hit = {c: d[c] for c in allc[t] if c in d}
            if not hit:
                continue
            for cover, cut in settings_from_scores(hit, allc[t], lower).items():
                cand[(t, inst, sig)].append((cover, cut))
                cuts_by_sig[(inst, sig)].add(cut)
    print("  %d (target, layer, signature) cells, %d distinct cuts"
          % (len(cand), sum(len(v) for v in cuts_by_sig.values())))

    print("domains pass 2: corpus cost of each cut")
    sorted_cuts = {k: sorted(v) for k, v in cuts_by_sig.items()}
    count = collections.Counter()
    for n, f in enumerate(files):
        if n % 40 == 0:
            print("  %d/%d" % (n, len(files)), flush=True)
        seen = collections.defaultdict(set)
        with gzip.open(f, "rt", errors="replace") as fh:
            for row in csv.DictReader(fh):
                if (row.get("leg") or "") != "fwd":
                    continue
                k = ((row.get("instrument") or "").strip(),
                     (row.get("signature") or "").strip())
                if k not in sorted_cuts:
                    continue
                try:
                    v = float(row.get("score") or "nan")
                except ValueError:
                    continue
                if math.isnan(v):
                    continue
                lower = kinds.get(k[0]) in LOWER_BETTER
                acc = A.norm((row.get("protein") or "").strip())
                for cut in sorted_cuts[k]:
                    if (v <= cut) if lower else (v >= cut):
                        seen[(k[0], k[1], cut)].add(acc)
        for kk, accs in seen.items():
            count[kk] += len(accs)

    front = collections.defaultdict(list)
    for (t, inst, sig), lst in cand.items():
        arm = "cdd" if inst == "CDD" else "interpro"
        d = sc[(inst, sig)]
        lower = kinds.get(inst) in LOWER_BETTER
        for cover, cut in lst:
            cost = sum(1 for a in neg[t] if a in d
                       and ((d[a] <= cut) if lower else (d[a] >= cut)))
            ncall = count[(inst, sig, cut)]
            front[(t, arm)].append(
                (cost, cover, (arm, "%s:%s" % (inst, sig), "%.6g" % cut,
                               kinds.get(inst, ""), ncall)))
    return front


# ------------------------------------------- OrthoFinder and SonicParanoid
def arm1_frontier(allc):
    """Every OrthoFinder and SonicParanoid setting, as (cost, cover, meta).

    Read from the output of `rl_arm1_all_settings.py` rather than derived again:
    the cover there is leave-one-out (a control may not recover itself), which is
    the only correct reading and cannot be obtained from a member list, since the
    members of a control-sharing rule are the partners of the controls and
    exclude the controls themselves.
    """
    if not (os.path.exists(ARM1_SETTINGS) and os.path.exists(ARM1_COVER)):
        raise SystemExit(
            "REFUSED: %s / %s missing. Without them OrthoFinder and SonicParanoid "
            "would be absent from the search and the other three instruments would "
            "be made to carry every control -- a milder form of the mistake R59 "
            "exists to correct. Run rl_arm1_all_settings.py first."
            % (ARM1_SETTINGS, ARM1_COVER))
    cover = collections.defaultdict(set)
    for r in rows(ARM1_COVER):
        cover[((r.get("TARGET") or "").strip(),
               (r.get("INSTRUMENT") or "").strip(),
               (r.get("SETTING") or "").strip())].add(
                   A.norm((r.get("PROTEIN_ID") or "").strip()))
    front = collections.defaultdict(list)
    n = 0
    for r in rows(ARM1_SETTINGS):
        t = (r.get("TARGET") or "").strip()
        inst = (r.get("INSTRUMENT") or "").strip()
        setting = (r.get("SETTING") or "").strip()
        if t not in TARGETS:
            continue
        try:
            cost = int(r.get("COST") or 0)
        except ValueError:
            continue
        try:
            ncall = int(r.get("N_CALLED") or 0)
        except ValueError:
            ncall = 0
        cov = frozenset(cover.get((t, inst, setting), set()) & allc[t])
        key = "orthofinder_full" if inst == "orthofinder" else inst
        front[(t, key)].append((cost, cov, (key, setting,
                                            (r.get("SOURCE") or "").strip(), "",
                                            ncall)))
        n += 1
    print("Arm 1: %d settings loaded across %d cells" % (n, len(front)))
    return front


# ----------------------------------------------------------------- joint solve
def solve(controls, frontiers, base_cover, base_cost, skip=frozenset()):
    """Maximise the summed score subject to the union covering `controls`.

    `skip` names the methods whose cell the copy-number ceiling emptied (reported,
    never chosen); the dynamic programme passes over them and the remaining
    methods carry the union.

        score = sum over methods of ( controls_i - NEG_WEIGHT x negatives_i )

    Dominance pruning on (covered set, score): a state is dropped when another
    covers a superset at no worse a score. Values are carried negated, so that
    "smaller is better" holds throughout and the pruning reads the same way.

    Negatives are summed across methods rather than combined as a union, so a
    decoy called by two methods is counted twice. That is deliberate: each method
    is scored on its own behaviour, which is what "optimised equally" means, and
    it keeps the score separable, so the search stays exact.
    """
    # states: covered set -> (negated running score, total calls, chain of
    # settings). The score is maximised and stored negated, so the same "smaller
    # is better" comparisons and dominance pruning work unchanged. Calls are a
    # tie-break, never a term in the score, but they must be carried through the
    # dynamic programme rather than applied only at the end: if a state were kept
    # by score alone, then between two chains reaching the same covered set at
    # equal score whichever came first would be kept, and the chain with the
    # fewest calls could be discarded before the final comparison. The dominance
    # pruning would have the same fault.
    states = {base_cover: (base_cost, 0, ())}
    for inst in INSTRUMENTS:
        # No empty option. Under the earlier cost-only objective a method could
        # score best by calling nothing, and seven cells did exactly that. The
        # score removes that incentive at the root (a setting that finds nothing
        # scores zero, and any setting recovering a control with no decoys beats
        # it), so silence no longer has to be forbidden. The empty option is left
        # out all the same: a method with no setting for a target is a fact to
        # report, not a free choice.
        opts = frontiers.get(inst)
        if not opts:
            if inst in skip:
                print("    %-17s cap-emptied under R64: no admissible setting, "
                      "reported and skipped" % inst, flush=True)
                continue
            raise SystemExit(
                "REFUSED: %s has no setting at all for this target, so it cannot "
                "be made to call. That is a fact about the instrument and must be "
                "reported, not papered over with an empty option." % inst)
        nxt = {}
        for cov, (cost, calls, chain) in states.items():
            for c2, cov2, meta in opts:
                k = cov | cov2
                # c2 is that setting's count of declared negatives; the score it
                # adds is (controls it covers) - NEG_WEIGHT x c2, and the running
                # value is the negated total, so that lower stays better.
                v = cost + (NEG_WEIGHT * c2 - len(cov2))
                n = calls + meta[4]
                cur = nxt.get(k)
                if cur is None or (v, n) < (cur[0], cur[1]):
                    nxt[k] = (v, n, chain + (meta + (c2, len(cov2)),))
        # Prune dominated states: a state whose covered set is the same or a
        # subset, at no better score, cannot lead anywhere better.
        items = sorted(nxt.items(), key=lambda kv: (kv[1][0], kv[1][1], -len(kv[0])))
        keep = {}
        for k, v in items:
            # Dominance compares (score, calls): a state is dropped only when
            # another is at least as good on both; otherwise the tie-break on
            # calls would depend on which state happened to be pruned first.
            if any(k <= k2 and (v2[0], v2[1]) <= (v[0], v[1])
                   for k2, v2 in keep.items()):
                continue
            keep[k] = v
        states = keep
        print("    %-17s %5d settings -> %6d states (from %d)"
              % (inst, len(opts), len(states), len(nxt)), flush=True)
        if len(states) > MAX_STATES:
            # Measured, not assumed. Stopping is better than exhausting memory:
            # the planned alternative is branch-and-bound with the admissible
            # bound V + sum over the remaining methods of max_s f(s).
            raise SystemExit(
                "REFUSED: %d states after %s exceeds MAX_STATES=%d. The exact DP "
                "is no longer tractable for this target; switch to "
                "branch-and-bound rather than letting it run."
                % (len(states), inst, MAX_STATES))
    full = [(v[0], v[1], k, v[2]) for k, v in states.items() if controls <= k]
    if full:
        # Determinism. Values carry the negated score, so min() finds the
        # maximum. Ties are broken first by the fewest calls and then
        # lexicographically on the chosen setting names: without a final key
        # the winner would depend on dictionary order, and two runs on identical
        # inputs could disagree.
        b = min(full, key=lambda x: (x[0], x[1],
                                     tuple(str(m[1]) for m in x[3])))
        return (b[0], b[2], b[3]), True
    best = max(states.items(), key=lambda kv: (len(kv[0]), -kv[1][0]))
    return (best[1][0], best[0], best[1][2]), False


def load_frontier_cache(allc):
    """Rebuild the frontier from `R59_frontier.tsv` and `R59_frontier_cover.tsv`.

    The four passes over the proteomes take about forty minutes and do not depend
    on the selection rule, so repeating them to change a tie-break is wasteful.
    Set `R59_REUSE_FRONTIER=1` to solve from the cached frontier instead. It must be
    asked for and is never the default: an outdated cache would silently describe
    settings measured against different inputs.
    """
    cov = collections.defaultdict(set)
    for r in rows(OUT_FCOV):
        cov[((r.get("TARGET") or "").strip(), (r.get("INSTRUMENT") or "").strip(),
             (r.get("SETTING") or "").strip(),
             (r.get("COST") or "").strip())].add(
                 A.norm((r.get("PROTEIN_ID") or "").strip()))
    front = collections.defaultdict(list)
    n = 0
    for r in rows(OUT_FRONT):
        t = (r.get("TARGET") or "").strip()
        inst = (r.get("INSTRUMENT") or "").strip()
        setting = (r.get("SETTING") or "").strip()
        # The frontier file names this column COST_DECLARED_NEG; the cover file
        # names it COST. Reading only "COST" here got an empty string and failed
        # on int(""), so `R59_REUSE_FRONTIER=1` was broken from the renaming of the
        # column until this fix; it went unnoticed because the option must be
        # asked for and rl_resolve_v5.sbatch unsets it. Either name is accepted.
        cost = ((r.get("COST_DECLARED_NEG") or r.get("COST") or "").strip())
        if not cost:
            raise SystemExit(
                "REFUSING: frontier cache row has no cost column "
                "(looked for COST_DECLARED_NEG and COST): %r" % (r,))
        if t not in TARGETS:
            continue
        c = frozenset(cov.get((t, inst, setting, cost), set()) & allc[t])
        front[(t, inst)].append((int(cost), c,
                                 (inst, setting, (r.get("CUT") or "").strip(),
                                  (r.get("CUT2") or "").strip(),
                                  int(r.get("N_CALLED") or 0))))
        n += 1
    # The ceiling must be applied here too. write_frontier() is skipped on the
    # cached path, so without this a solve with `R59_REUSE_FRONTIER=1` would have no
    # ceiling while a fresh run had one: the same answer reached two ways,
    # differing silently. Filtering a cache that already respects the ceiling
    # changes nothing, so this is safe whether or not the cache predates it.
    kept = {}
    n_drop = 0
    for k, opts in front.items():
        # The copy-number ceiling filters the allowed settings here, and the
        # refusal below enforces it.
        adm = [o for o in opts if o[2][4] <= CALL_CAP]
        n_drop += len(opts) - len(adm)
        if opts and not adm:
            raise SystemExit(
                "REFUSED: the cached frontier for %s has no setting under "
                "CALL_CAP=%d. Re-solve without R59_REUSE_FRONTIER rather than "
                "relaxing the cap against a stale cache." % (k, CALL_CAP))
        kept[k] = adm
    print("frontier cache: %d settings over %d cells (%d dropped by "
          "CALL_CAP=%d)" % (n, len(front), n_drop, CALL_CAP))
    return kept


def main():
    rl_invariants.check_inputs()
    # The control set is the set of controls in use, not every declared row.
    #
    # This once read the declared-control table directly, with no REVIEW filter,
    # so it solved against all 247 declared rows whatever had been held out. Two
    # consequences, both seen in one job before the fault was caught: rows that
    # the curation had contested were still required as controls, and
    # `R59_CONTROL_SUBSET` had no effect at all (the setA and setB tasks produced
    # byte-identical output, because neither read the set it had been given).
    #
    # `binding_controls()` is the single definition, already used by every
    # earlier script, and it is where the subset is applied. A second reading of
    # the table must not be reintroduced here.
    allc_d, label = A.binding_controls()
    allc = {t: frozenset(allc_d.get(t, ())) for t in TARGETS}
    empty = [t for t in TARGETS if not allc[t]]
    if empty:
        raise SystemExit("REFUSING: no binding control for %s" % ", ".join(empty))
    print("control set: %s -- %d controls %s"
          % (label, sum(len(v) for v in allc.values()),
             {t: len(allc[t]) for t in TARGETS}))

    reuse = os.environ.get("R59_REUSE_FRONTIER") == "1"
    if reuse and os.path.exists(OUT_FRONT) and os.path.exists(OUT_FCOV):
        print("R59_REUSE_FRONTIER=1: solving from the cached frontier, corpus "
              "passes SKIPPED")
        front = load_frontier_cache(allc)
        pruned = {(t, i): front.get((t, i), []) for t in TARGETS
                  for i in INSTRUMENTS}
        return report(allc, front, pruned, cached=True)
    front = collections.defaultdict(list)
    dec, negc = A.declared_decoys(), A.negative_controls()
    neg = {t: (dec[t] | negc[t]) for t in TARGETS}
    print("declared negatives: %s" % {t: len(neg[t]) for t in TARGETS})
    for k, v in list(hmmer_frontier(allc, neg).items()) \
            + list(domain_frontier(allc, neg).items()) \
            + list(arm1_frontier(allc).items()):
        front[k].extend(v)
    return report(allc, front, None, cached=False)


def write_frontier(allc, front):
    """Write the frontier and its cover sets, and return the pruned frontier.

    This is a separate function so that a cached solve can skip it: in that mode
    these two files are the input, and rewriting them from themselves could
    truncate the cache if anything failed during the write."""
    n_excl = {}
    pruned = {}
    with open(OUT_FRONT, "w") as fh:
        fh.write("## Every enumerated setting for every Arm 2 instrument, after "
                 "Pareto pruning. COST is proteins called across the 210 "
                 "proteomes minus the controls recovered. A setting covering only "
                 "part of the control set is not a failure under R59; it is the "
                 "point.\n")
        fh.write("TARGET\tINSTRUMENT\tSETTING\tCUT\tCUT2\tCONTROLS_COVERED\t"
                 "N_CONTROLS\tCOST_DECLARED_NEG\tN_CALLED\t"
                 "CALLED_PER_SPECIES\tINSTRUMENT_SCORE\n")
        pruned = {}
        # The cover set itself, not only its size, so that a new solve can reuse
        # this instead of repeating the four passes over the proteomes; the
        # passes are the expensive part and do not change when only the
        # selection rule changes.
        with open(OUT_FCOV, "w") as ch:
            ch.write("## Which controls each frontier setting recovers, long "
                     "form. Companion to R59_frontier.tsv; together they let the "
                     "joint solve be re-run without re-scanning the corpus.\n")
            ch.write("TARGET\tINSTRUMENT\tSETTING\tCUT\tCUT2\tCOST\t"
                     "N_CALLED\tPROTEIN_ID\n")
            for t in TARGETS:
                for inst in INSTRUMENTS:
                    # The order matters: keep the allowed settings first, then
                    # prune. pareto() drops a setting when another covers a
                    # superset at no greater cost. If the ceiling were applied
                    # afterwards, an allowed setting could already have been
                    # dropped because a setting over the ceiling dominated it,
                    # and the solver would never see it. Filtering first makes
                    # the frontier the frontier of allowed settings, which is
                    # also what is cached.
                    raw = front.get((t, inst), [])
                    # The copy-number ceiling, enforced by the refusal below.
                    adm = [o for o in raw if o[2][4] <= CALL_CAP]
                    if raw and not adm:
                        if CAP_EMPTIED != "report":
                            raise SystemExit(
                                "REFUSED: every %s setting for %s calls more than "
                                "CALL_CAP=%d (%d per species). The cap has emptied "
                                "this cell; it is not that the instrument is silent. "
                                "Raise COPY_CEILING or accept that this instrument "
                                "cannot be run within the ceiling."
                                % (inst, t, CALL_CAP, COPY_CEILING))
                        # With `R59_CAP_EMPTIED=report` the cell is reported as a
                        # fact and the remaining methods carry the union (see
                        # CAP_EMPTIED).
                        EMPTIED[(t, inst)] = (min(o[2][4] for o in raw), len(raw))
                        n_excl[(t, inst)] = len(raw)
                        pruned[(t, inst)] = []
                        print("  %-8s %-17s INADMISSIBLE under R64: every one of %d "
                              "settings calls more than CALL_CAP=%d (min N_CALLED %d)"
                              % (t, inst, len(raw), CALL_CAP, EMPTIED[(t, inst)][0]))
                        continue
                    n_excl[(t, inst)] = len(raw) - len(adm)
                    p = pareto(adm)
                    pruned[(t, inst)] = p
                    for cost, cover, meta in sorted(p, key=lambda x: (x[0], -len(x[1]))):
                        fh.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t"
                                 "%.2f\t%.1f\n"
                                 % (t, inst, meta[1], meta[2], meta[3],
                                    len(cover), len(allc[t]), cost, meta[4],
                                    meta[4] / float(N_SPECIES),
                                    len(cover) - NEG_WEIGHT * cost))
                        for a in sorted(cover):
                            ch.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\n"
                                     % (t, inst, meta[1], meta[2], meta[3],
                                        cost, meta[4], a))

    return pruned


def report(allc, front, pruned, cached=False):
    if not cached:
        pruned = write_frontier(allc, front)
    md = ["# R59: the joint optimum, chosen at the union", "",
          "Settings are selected together to minimise proteins called across the "
          "210 proteomes, subject to the UNION recovering every declared control. "
          "An instrument covering part of the control set is correct, not a "
          "shortfall.", "",
          "| target | controls | covered | instrument | setting | cut | covers | "
          "cost | per species |", "|---|---|---|---|---|---|---|---|---|"]
    out = []
    print("\n%-8s %-9s %-15s %s"
          % ("target", "covered", "Control_recovery", "solution"))
    for t in TARGETS:
        fr = {inst: pruned[(t, inst)] for inst in INSTRUMENTS}
        skip = frozenset(i for i in INSTRUMENTS if (t, i) in EMPTIED)
        if cached and CAP_EMPTIED == "report":
            # A cached frontier holds no rows for a cell emptied by the ceiling (it
            # never had an allowed setting), so in report mode an empty cell is
            # read as emptied.
            skip = skip | frozenset(i for i in INSTRUMENTS if not fr[i])
        (cost, cov, chain), ok = solve(allc[t], fr, frozenset(), 0, skip)
        print("%-8s %2d/%-6d score=%-9.1f %s%s"
              % (t, len(cov & allc[t]), len(allc[t]), -cost,
                 "COMPLETE" if ok else "INCOMPLETE",
                 ("  [R64 cap-emptied: %s]" % ", ".join(sorted(skip))) if skip else ""))
        for inst in sorted(skip):
            ncall_min = EMPTIED.get((t, inst), ("", 0))[0]
            md.append("| %s | %d | 0 | %s | `%s` | | 0 | 0 | %s |"
                      % (t, len(allc[t]), inst, CAP_EMPTIED_SETTING,
                         ("INADMISSIBLE under R64: every enumerated setting calls "
                          "more than CALL_CAP=%d (min N_CALLED %s)"
                          % (CALL_CAP, ncall_min))))
            out.append((t, inst, CAP_EMPTIED_SETTING, "", "", "",
                        0, len(allc[t]), 0, ncall_min if ncall_min != "" else 0,
                        ("%.2f" % (ncall_min / float(N_SPECIES))) if ncall_min != "" else "",
                        "0.0", "%.1f" % (-cost),
                        "COMPLETE" if ok else "INCOMPLETE", "r64_cap_emptied"))
        for meta in chain:
            inst, setting, cut, cut2, ncall, c, ncov = meta
            # HMMER carries two cuts, because the report requires both a sequence
            # score and a domain score, so cut2 is CUT_DOM there; for CDD and
            # InterProScan it is the score kind. An earlier version dropped this
            # field, which would have given the report a sequence cut with no
            # domain cut and silently weakened the rule.
            md.append("| %s | %d | %d | %s | `%s` | %s | %d | %d | %.2f |"
                      % (t, len(allc[t]), len(cov & allc[t]), inst, setting,
                         ("%s / %s" % (cut, cut2)) if inst == "hmmer" else cut,
                         ncov, c, c / float(N_SPECIES)))
            out.append((t, inst, setting, cut,
                        cut2 if inst == "hmmer" else "",
                        "" if inst == "hmmer" else cut2,
                        ncov, len(allc[t]), c, ncall,
                        "%.2f" % (ncall / float(N_SPECIES)),
                        "%.1f" % (ncov - NEG_WEIGHT * c),
                        "%.1f" % (-cost),          # CONTROL_RECOVERY is the negated solver value
                        "COMPLETE" if ok else "INCOMPLETE", "r59_joint"))
        if not ok:
            miss = sorted(allc[t] - cov)
            md.append("| %s | | | **UNCOVERED** | %s | | | | |"
                      % (t, ", ".join("`%s`" % x for x in miss[:10])))

    with open(OUT_SET, "w") as fh:
        fh.write("## R59 joint settings: one row per (target, instrument). "
                 "CONTROLS_COVERED is what THAT instrument recovers; the union "
                 "over the rows of a target is what R56 requires. COST is "
                 "corpus-wide proteins called minus controls recovered. "
                 "SOURCE=arm1_fixed means the setting came from FLOOR/SELECT and "
                 "was carried, not re-searched.\n")
        fh.write("## INSTRUMENT_SCORE = CONTROLS_COVERED - NEG_WEIGHT x "
                 "COST_DECLARED_NEG, where COST_DECLARED_NEG is the declared "
                 "decoys and adjudicated negative controls the setting admits. "
                 "CONTROL_RECOVERY is the sum of INSTRUMENT_SCORE over that "
                 "target's five rows and is the quantity maximised. N_CALLED is "
                 "bounded: a setting calling more than CALL_CAP is INADMISSIBLE "
                 "by ruling -- broken species are cut downstream. "
                 "CONTROLS_COVERED is this instrument's share and is expected to "
                 "be partial; the UNION over a target's rows is what R56 "
                 "requires.\n")
        if EMPTIED:
            fh.write("## R64 CAP-EMPTIED CELLS (R59_CAP_EMPTIED=report): %s. In each, "
                     "EVERY enumerated setting called more than CALL_CAP=%d, so the "
                     "instrument has no admissible operating point for that target; "
                     "the row carries SETTING=\"%s\", SOURCE=r64_cap_emptied and "
                     "INSTRUMENT_SCORE 0, is a reported fact and not a chosen setting, "
                     "and the remaining instruments carry the union.\n"
                     % (", ".join("%s/%s (min N_CALLED %d of %d settings)"
                                  % (t, i, v[0], v[1]) for (t, i), v in sorted(EMPTIED.items())),
                        CALL_CAP, CAP_EMPTIED_SETTING))
        fh.write("TARGET\tINSTRUMENT\tSETTING\tCUT\tCUT_DOM\tSCORE_KIND\t"
                 "CONTROLS_COVERED\tN_CONTROLS\tCOST_DECLARED_NEG\tN_CALLED\t"
                 "CALLED_PER_SPECIES\tINSTRUMENT_SCORE\tCONTROL_RECOVERY\t"
                 "UNION_STATUS\tSOURCE\n")
        for r in out:
            fh.write("\t".join(str(x) for x in r) + "\n")
    open(OUT_MD, "w").write("\n".join(md) + "\n")
    print("\nwrote %s\n      %s\n      %s" % (OUT_SET, OUT_FRONT, OUT_MD))


if __name__ == "__main__":
    main()
