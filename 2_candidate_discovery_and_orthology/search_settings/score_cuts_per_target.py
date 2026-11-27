#!/usr/bin/env python3
# >>> deposit header
# Purpose: For each target and each of its panel models, find the pair of sequence and domain score cuts that
#          keeps the false-positive rate on the target's same-species decoys at or below 5, 1, 0.5, 0.1 and 0%,
#          and the recall of the target's own controls at each.
# Methods: Candidate discovery
# Inputs:  the current declared-control table (through `control_declaration`); the decoy sets
#          `02_reference_data/decoys_same_species/<target>_same-species-decoys.faa`; the 210 proteomes
#          (`concord/null_split/fwd/*.faa`); `panel_v2.hmm`; `PANEL_DECLARATION_v2.tsv`; `THRESHOLD_SOURCES.tsv`
# Outputs: `06_calibration_and_checks/THRESHOLD_PILOT_PERTARGET.tsv`, or `THRESHOLD_PILOT_PERTARGET__<set>.tsv`
#          for a non-core target set; the benchmark FASTA and tblout in `threshold_pilot_pertarget<suffix>/`
# Tools:   HMMER hmmsearch (busco6 environment; --max -E 1e5 --domE 1e5); Python 3 standard library and the
#          project modules `control_declaration` and `targets`; versions: see envs/
# <<< deposit header
"""Threshold pilot by target: every target measured against its own positives and negatives.

Why this exists. rl_threshold_pilot.py sets POS = {"POSITIVE_AMPK"}: it scored all 229 models
against the 24 verified AMPKα proteins, whatever target a model belonged to. Every row of
THRESHOLD_PILOT.tsv for a target other than AMPKα therefore reports N_POS=24 and asks, in effect,
whether a SIN1 model recovers AMPKα proteins, which is why twelve models report exactly 0% recall.
Those zeros come from the wrong positive set, and seven of the eight targets had never had a
threshold measured at all.

What is measured, for each target:
   positives   the target's own declared identity controls (222 across the eight targets, never
               pooled)
   negatives   the target's own same-species decoy set: real proteins from the same proteomes
               (for example >DECOY__Homo_sapiens__NP_000446.1), so the false-positive rate is a
               rate on hard negatives, not on shuffled sequence. No decoy is a control of its own
               target (checked).
   models      only the models that PANEL_DECLARATION_v2 declares for that target, on either pole.

The script chooses nothing. For each model and target it reports the score cuts needed to reach
each of several false-positive rates, and the recall obtained there, so that the rate to adopt can
be chosen on a measured error rate rather than arbitrarily. The cut proposed before, the lowest bit
score of any control, was fitted to the positives and had no measured error rate at all.

Two safeguards:

No self-recovery. Four target-pole models were built from the controls they would be scored on:
AMPKa_controls, AMPKb_controls, SIN1_controls and AMPKalpha_ours (SOURCES project_controls or
project_ampka, as the declaration's own NOTE says). Their recall on those controls is not a
measurement. They are marked SELF_BUILT=yes, and their recall must not be compared with any other
model's. Their false-positive rate on decoys is valid, because decoys were not used to build them.
This matters: AMPKalpha_ours is the AMPKα model with the lowest false-positive rate, and its recall
is circular.

Resolution of the rate. A rate cannot be measured more finely than 1/N_NEG. AMPKβ has 104 decoys
(0.96%) and AMPKγ 167 (0.60%), so a 0.1% cut cannot be measured for either. FINEST_FPR gives this
limit per target instead of printing a misleading 0.00%.

Under a target set other than the core eight (environment variable TARGET_SET), the script works in
its own directory and writes the companion table THRESHOLD_PILOT_PERTARGET__<set>.tsv; the core
table, fixed in advance, is never rewritten.
"""
# deposit path defaults: begin
_TA210_ENVS = __import__("os").environ.get("TA210_ENVS") or "/gpfs2/scratch/bdigenov/envs"
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import csv, collections, glob, os, subprocess, sys

R     = os.environ.get("PROJECT_ROOT") or (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1")
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "lib"))
import control_declaration                                            # noqa: E402
import targets                                                # noqa: E402
PANEL = R + "/02_reference_data/hmmer_panel/panel_v2.hmm"
DECL  = R + "/02_reference_data/hmmer_panel/PANEL_DECLARATION_v2.tsv"
CTRLS = control_declaration.current(R + "/02_reference_data")[0]     # the current (highest-numbered) control table
DECOY = R + "/02_reference_data/decoys_same_species"
THR   = R + "/05_method_and_rules/THRESHOLD_SOURCES.tsv"
PROT  = (_TA210_SCRATCH + "/concord/null_split/fwd")
SFX   = targets.suffix()
WORK  = (_TA210_SCRATCH + "/ms210_phase2/threshold_pilot_pertarget") + SFX
BIN   = (_TA210_ENVS + "/busco6/bin")
OUT   = R + "/06_calibration_and_checks"

TARGETS   = targets.targets()          # from TARGET_SET; default the core eight, in their fixed order
SELF_SRC  = {"project_controls", "project_ampka"}   # declared, not inferred; see "No self-recovery" above
# A "contaminant" model built from the target's own in-paralogs is the target, not a confusable
# family. It scores 1,580 bits on a real Arabidopsis RAPTOR against 1,543 for the target pole and
# would reject 14 of the 23 RAPTOR controls if counted. The first candidate screen excluded it;
# this pilot marks it (SELF_BUILT=yes-own-inparalogs) so that no cut is ever read from it.
SELF_CONTAM = {"RAPTOR__RAPTOR_s_own_lineage-specific_in-paralogs": "RAPTOR"}
FPR_LEVELS = [0.05, 0.01, 0.005, 0.001, 0.0]

def norm(a): return a.split("||")[-1].split(":")[0].split()[0]

def load_controls():
    """accession -> set of targets it is a declared control for; rows held out in REVIEW are excluded."""
    pos = collections.defaultdict(set); sp = {}
    with open(CTRLS, errors="replace") as fh:
        rd = csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t")
        for r in rd:
            t = (r.get("TARGET") or "").strip()
            a = norm((r.get("PROTEIN_ID") or "").strip())
            if t in TARGETS and a and not (r.get("REVIEW") or "").strip():
                pos[a].add(t); sp[a] = (r.get("SPECIES") or "").strip()
    return pos, sp

def load_decoys():
    """accession -> (set of targets it is a decoy for, species, sequence)."""
    neg = collections.defaultdict(set); sp = {}; seq = {}
    for t in TARGETS:
        f = "%s/%s_same-species-decoys.faa" % (DECOY, t)
        if not os.path.exists(f): sys.exit("missing decoy set: " + f)
        acc = None
        take = False
        for line in open(f, errors="replace"):
            if line[0] == ">":
                parts = line[1:].strip().split("__")
                acc = norm(parts[-1]); neg[acc].add(t)
                sp.setdefault(acc, parts[1] if len(parts) > 2 else "")
                # One copy per accession. An earlier version guarded the dictionary key but not
                # the appending of sequence lines, so a protein that is a decoy for k targets had
                # its residues concatenated k times: human titin reached 215,946 residues across
                # six pools and hmmsearch stopped at its 100,000-residue limit. Where the search
                # does not stop, a concatenated sequence also scores too high, because the
                # full-sequence bit score sums every domain hit.
                take = acc not in seq
                if take:
                    seq[acc] = []
            elif acc is not None and take:
                seq[acc].append(line.strip())
    return neg, sp, {a: "".join(v) for a, v in seq.items()}

def build_fasta(pos, dseq):
    """One FASTA of unique accessions: the controls, taken from the 210 proteomes, and the decoys."""
    os.makedirs(WORK, exist_ok=True)
    fa = WORK + "/pertarget_benchmark.faa"
    if os.path.exists(fa) and os.path.getsize(fa) > 0:
        return fa, sum(1 for l in open(fa) if l[0] == ">")
    found = {}
    want = set(pos)
    for p in sorted(glob.glob(PROT + "/*.faa")):
        keep, acc = False, None
        for line in open(p, errors="replace"):
            if line[0] == ">":
                acc = norm(line[1:]); keep = acc in want and acc not in found
                if keep: found[acc] = []
            elif keep:
                found[acc].append(line.strip())
    missing = want - set(found)
    if missing:
        sys.stderr.write("controls NOT resolved in the 210 proteomes (%d): %s\n"
                         % (len(missing), ", ".join(sorted(missing)[:10])))
    with open(fa, "w") as out:
        for a, s in sorted(found.items()):
            out.write(">%s\n%s\n" % (a, "".join(s)))
        for a, s in sorted(dseq.items()):
            if a not in found:                       # a control is never re-emitted as a decoy
                out.write(">%s\n%s\n" % (a, s))
    return fa, sum(1 for l in open(fa) if l[0] == ">")

def run_panel(fa):
    """Run hmmsearch --max and return two scores per protein per model, read from the tblout:

       column index 5 (0-based)  full-sequence bit score (the first candidate screen used it
                                 for HMMER_BEST_TARGET_BITS)
       column index 8 (0-based)  best single-domain score

    The pilot sets a cut on each (see pair_at_fpr). The full-sequence score sums every domain hit,
    so it measures how much of a protein matches, not how well. Measured here: LST8 is a WD40
    protein of about 330 residues whose declared decoy set is the WD40 clan, and multi-repeat
    proteins of 1,892 to 5,898 residues outscore every real LST8 on a WD40 model purely by carrying
    more repeats; AMPKα against the protein kinase clan behaves the same way. The effect is large:
    at a false-positive rate of 1% or less on those decoys, recall rises from 0% to 100% for LST8's
    KOG0315 model and for AMPKa_controls simply by using the best single-domain score, while the
    targets that already worked (TOR, RAPTOR, RICTOR, SIN1) are unchanged. Both scores are returned
    so that the choice can be checked."""
    tbl = WORK + "/pertarget_benchmark.tblout"
    def complete(f):
        if not (os.path.exists(f) and os.path.getsize(f) > 0): return False
        with open(f, "rb") as fh:
            fh.seek(max(0, os.path.getsize(f) - 4096))
            return b"# [ok]" in fh.read()
    if not complete(tbl):
        if os.path.exists(tbl):
            sys.stderr.write("existing tblout INCOMPLETE (no '# [ok]') -- rerunning\n")
            os.remove(tbl)
        subprocess.run([BIN + "/hmmsearch", "--cpu", os.environ.get("SLURM_CPUS_PER_TASK", "8"),
                        "-E", "1e5", "--domE", "1e5", "--max",
                        "--tblout", tbl, PANEL, fa],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    if not complete(tbl):
        sys.exit("hmmsearch did not finish: %s has no '# [ok]' footer" % tbl)
    dom  = collections.defaultdict(dict)
    full = collections.defaultdict(dict)
    with open(tbl, errors="replace") as fh:
        for line in fh:
            if line[0] == "#": continue
            p = line.split()
            if len(p) < 9: continue
            try: fb, db = float(p[5]), float(p[8])
            except ValueError: continue
            a, m = p[0], p[2]
            if m not in dom[a]  or db > dom[a][m]:  dom[a][m]  = db
            if m not in full[a] or fb > full[a][m]: full[a][m] = fb
    return dom, full

def cut_at_fpr(pv, nv, level):
    """Cut on a single score, kept only to report what a rule on one score alone would have given."""
    n = len(nv)
    if not n or not pv: return None, None
    allowed = int(level * n)
    ns = sorted(nv, reverse=True)
    cut = ns[min(allowed, n - 1)] + 0.1
    if cut == float("-inf"):
        fin = [x for x in pv if x > float("-inf")]
        if not fin: return None, None
        cut = min(fin)
    return cut, sum(1 for x in pv if x >= cut) / len(pv)

MAX_CAND = 400          # candidate domain cuts swept per model per level

def _kth_cut(vals, allowed):
    """Lowest cut letting at most `allowed` of vals through. -inf means no cut is needed."""
    v = sorted(vals, reverse=True)
    return (v[allowed] + 0.1) if len(v) > allowed else float("-inf")

def pair_at_fpr(pos, neg, level):
    """Pfam's convention, made measurable: a hit must clear both its sequence cut and its domain
    cut, and the false-positive rate is the rate of negatives that clear both, not either alone.

    Chosen on 5 September 2026. Pfam gives every family a sequence gathering threshold and a
    domain gathering threshold and requires both.

    pos and neg are lists of (sequence score, domain score). Many pairs of cuts reach a given rate,
    so this sweeps candidate domain cuts and, for each, takes the lowest sequence cut that keeps the
    joint rate at or below `level`, keeping the pair with the best recall.

    Each candidate domain cut is a negative's domain score plus a step of 0.1, not the score
    itself: a cut equal to a negative's score still accepts that negative, which once left the
    domain-only solution outside the search space and made the two-score rule do worse than either
    score alone (LST8 KOG0315 at 0% recall for a rate of 0.1% or less, where the domain score alone
    gives 96%). The two single-score optima are also added explicitly, so the result is never worse
    than the better of them."""
    n = len(neg)
    if not n or not pos: return None, None, None
    allowed = int(level * n)
    doms = [b for _, b in neg]
    cands = sorted({b + 0.1 for b in doms})
    if len(cands) > MAX_CAND:
        step = len(cands) / float(MAX_CAND)
        cands = [cands[min(len(cands) - 1, int(k * step))] for k in range(MAX_CAND)]
    cands.append(float("-inf"))                 # anchor: no domain requirement (sequence-only)
    dom_only = _kth_cut(doms, allowed)          # anchor: the domain-only optimum
    if dom_only != float("-inf"): cands.append(dom_only)

    best = (-1.0, None, None)
    for d in cands:
        passing = [a for a, b in neg if b >= d]
        sc = _kth_cut(passing, allowed)
        rec = sum(1 for a, b in pos if b >= d and a >= sc) / float(len(pos))
        if rec > best[0]: best = (rec, sc, d)
    rec, sc, d = best
    # A score with no effective cut is reported as None ("none" in the table), never as an
    # invented number. Substituting the lowest score of any accepted positive was tried and is
    # wrong: for TOR FRB_dom it gave the pair 9.1 / 51.4, a stricter rule than the one whose
    # false-positive rate was measured, which genome-wide would reject proteins the measured rule
    # accepts. A reported threshold has to be the rule that was measured.
    #
    # A domain cut above the sequence cut is possible. The best single-domain score exceeds the
    # full-sequence score in 204,384 of 793,811 tblout rows (25.7%, by up to 1,754 bits), because
    # HMMER applies its null2 bias correction differently to the whole sequence and to each domain.
    # The two scores are separate quantities and neither bounds the other, which is part of why
    # requiring both is informative.
    if sc == float("-inf"): sc = None
    if d  == float("-inf"): d  = None
    return sc, d, rec

def main():
    pos, psp = load_controls()
    neg, nsp, dseq = load_decoys()
    print("declared controls: %d proteins across %d targets" % (len(pos), len(TARGETS)))
    for t in TARGETS:
        print("   %-7s positives %-3d   decoys %-5d   finest measurable FPR %.3f%%"
              % (t, sum(1 for a in pos if t in pos[a]), sum(1 for a in neg if t in neg[a]),
                 100.0 / max(1, sum(1 for a in neg if t in neg[a]))))
    fa, n = build_fasta(pos, dseq)
    print("\nbenchmark FASTA: %d unique proteins" % n)
    best, fullsc = run_panel(fa)
    print("scored: %d proteins (statistic: BEST SINGLE DOMAIN, tblout col 8)\n" % len(best))

    pole, forT = collections.defaultdict(dict), collections.defaultdict(set)
    src = {}
    for r in csv.DictReader(open(DECL), delimiter="\t"):
        forT[r["FOR_TARGET"]].add(r["MODEL"])
        pole[r["FOR_TARGET"]][r["MODEL"]] = r["POLE"]
        src[r["MODEL"]] = r["SOURCES"]
    pub = {}
    for r in csv.DictReader((l for l in open(THR) if not l.lstrip('"').startswith("##")),
                            delimiter="\t"):
        if r["THRESHOLD_SOURCE"] != "NONE":
            pub.setdefault(r["MODEL"], (r["THRESHOLD_SOURCE"], r["THRESHOLD_VALUE"]))

    rows = []
    for t in TARGETS:
        pacc = sorted(a for a in pos if t in pos[a] and a in best)
        nacc = sorted(a for a in neg if t in neg[a] and a in best)
        pmiss = sum(1 for a in pos if t in pos[a] and a not in best)
        nspec = len({nsp.get(a, "") for a in nacc}) ; pspec = len({psp.get(a, "") for a in pacc})
        finest = 100.0 / max(1, len(nacc))
        for m in sorted(forT[t]):
            pv = [best[a][m] for a in pacc if m in best[a]]
            nv = [best[a][m] for a in nacc if m in best[a]]
            pv += [float("-inf")] * (len(pacc) - len(pv))
            nv += [float("-inf")] * (len(nacc) - len(nv))
            self_built = "yes" if src.get(m, "") in SELF_SRC else "no"
            if SELF_CONTAM.get(m) == t: self_built = "yes-own-inparalogs"
            row = [m, pole[t][m], t, self_built, len(pacc), pspec, pmiss, len(nacc), nspec,
                   "%.3f%%" % finest]
            ps, pvv = pub.get(m, ("NONE", ""))
            pr = pf = ""
            try:
                pc = float(str(pvv).split("/")[0].strip())
                pr = "%.0f%%" % (100 * sum(1 for x in pv if x >= pc) / len(pv))
                pf = "%.2f%%" % (100 * sum(1 for x in nv if x >= pc) / len(nv))
            except (ValueError, ZeroDivisionError, IndexError):
                pass
            row += [ps, pvv, pr, pf]
            ppair = [(fullsc[a].get(m, float("-inf")), best[a].get(m, float("-inf")))
                     for a in pacc]
            npair = [(fullsc[a].get(m, float("-inf")), best[a].get(m, float("-inf")))
                     for a in nacc]
            for lv in FPR_LEVELS:
                sc, dc, rc = pair_at_fpr(ppair, npair, lv)
                row += ["none" if sc is None else "%.1f" % sc,
                        "none" if dc is None else "%.1f" % dc,
                        "" if rc is None else "%.0f%%" % (100 * rc)]
            # what a rule on the domain score alone and on the sequence score alone gives at 1%, for comparison
            _, rdom = cut_at_fpr([b for _, b in ppair], [b for _, b in npair], 0.01)
            _, rseq = cut_at_fpr([a for a, _ in ppair], [a for a, _ in npair], 0.01)
            row += ["" if rdom is None else "%.0f%%" % (100 * rdom),
                    "" if rseq is None else "%.0f%%" % (100 * rseq)]
            rows.append(row)

    os.makedirs(OUT, exist_ok=True)
    p = OUT + "/THRESHOLD_PILOT_PERTARGET%s.tsv" % SFX
    if not SFX and os.path.exists(p):
        sys.exit("REFUSED: %s is a FROZEN input (INPUT_MANIFEST.tsv) and is never rewritten; "
                 "run under a non-core TARGET_SET to write its companion" % p)
    hdr = ["MODEL","POLE","FOR_TARGET","SELF_BUILT","N_POS","N_POS_SPECIES","N_POS_UNSCORED",
           "N_NEG","N_NEG_SPECIES","FINEST_FPR","PUBLISHED_SOURCE","PUBLISHED_CUT",
           "PUB_RECALL","PUB_FPR"]
    for lv in FPR_LEVELS:
        hdr += ["CUTSEQ_FPR%g" % (lv * 100), "CUTDOM_FPR%g" % (lv * 100),
                "REC_FPR%g" % (lv * 100)]
    hdr += ["REC_FPR1_DOMONLY", "REC_FPR1_SEQONLY"]
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if SFX:
            w.writerow(["## COMPANION of the frozen THRESHOLD_PILOT_PERTARGET.tsv for TARGET_SET=%s "
                        "(L13): rows for these targets only, read after the frozen file by "
                        "lib/targets.pilot_rows(). Declaration %s; panel %s."
                        % (targets.set_name(), os.path.basename(CTRLS), os.path.basename(PANEL))])
        w.writerow(["## Each target measured against ITS OWN controls (positives) and ITS OWN "
                    "same-species decoys (negatives). CUT_FPRx is the lowest bit cut whose "
                    "false-positive rate on that target's decoys is <= x%, and REC_FPRx the "
                    "recall there. MEASUREMENTS, NOT RULINGS. SELF_BUILT=yes means the model was "
                    "built from the controls it is scored on (R28): its RECALL is circular and "
                    "must not be compared; its FPR is valid; yes-own-inparalogs marks a "
                    "'contaminant' built from the target's own in-paralogs, which is the target. "
                    "A rate cannot be finer than FINEST_FPR = 1/N_NEG. Scores are the BEST "
                    "SINGLE DOMAIN (tblout col 8) and the FULL SEQUENCE (col 5) TOGETHER: "
                    "Bruno ruled 2026-09-05 for Pfam's own convention, where a call must clear a "
                    "sequence cut AND a domain cut. CUTSEQ_FPRx / CUTDOM_FPRx is the pair whose "
                    "JOINT false-positive rate on that target's decoys is <= x%, chosen to "
                    "maximise recall; REC_FPR1_DOMONLY and REC_FPR1_SEQONLY show what each score "
                    "alone would have given at 1%."])
        w.writerow(hdr); w.writerows(rows)
    print("wrote %s  (%d model x target rows)\n" % (p, len(rows)))

    # Output columns: 14 fixed columns, then CUTSEQ, CUTDOM and REC for each rate 5, 1, 0.5, 0.1, 0%.
    S1, D1, R1 = 17, 18, 19                  # false-positive rate <= 1%
    S01, D01, REC01 = 23, 24, 25             # false-positive rate <= 0.1%
    for t in TARGETS:
        tr = [r for r in rows if r[2] == t and r[1] == "target"]
        if not tr:
            print("=== %s: NO target-pole model produced a row ===\n" % t); continue
        tr.sort(key=lambda r: -float((r[REC01] or "0%").rstrip("%")))
        print("=== %s   pos=%s (%s sp, %s unscored)  neg=%s (%s sp)  finest FPR %s ==="
              % (t, tr[0][4], tr[0][5], tr[0][6], tr[0][7], tr[0][8], tr[0][9]))
        print("  %-40s %-5s %17s %6s %17s %6s"
              % ("target-pole model", "self", "seq/dom cut @1%", "rec",
                 "seq/dom cut @0.1%", "rec"))
        for r in tr[:8]:
            print("  %-40s %-5s %17s %6s %17s %6s"
                  % (r[0][:40], r[3],
                     "%s / %s" % (r[S1] or "-", r[D1] or "-"), r[R1] or "-",
                     "%s / %s" % (r[S01] or "-", r[D01] or "-"), r[REC01] or "-"))
        print()

if __name__ == "__main__":
    main()
