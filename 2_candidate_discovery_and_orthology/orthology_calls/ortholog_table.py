#!/usr/bin/env python3
# >>> deposit header
# Purpose: Builds the final ortholog table, one orthology call per candidate from the gene-tree clade
#          criterion (two-arm support where the tree abstained), and the list of retained orthologs.
# Methods: Orthology calls
# Inputs:  06_candidate_list.tsv, 07_contaminant_screen.tsv, 14_arm_rule_calls.tsv (per-method report
#          folder); LB_SCORES.tsv and, per target, _contaminants_iqtree.tsv, _contaminants_pmsf.tsv,
#          _subtree_calls_iqtree.tsv, _labels.tsv, _excluded_multirecord.tsv, _cropped_out.tsv
#          of the gene-tree track
# Outputs: 04_final_tables/FINAL_ORTHOLOG_TABLE__<track>.tsv and FINAL_ORTHOLOGS__<track>.tsv; for the
#          primary (gene-collapsed) track also FINAL_ORTHOLOG_TABLE.tsv and FINAL_ORTHOLOGS.tsv
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
"""The final ortholog table: every candidate hit, with the non-orthologs removed.

The aim is a trustworthy ortholog table per target across 210 eukaryotic proteomes.
Every other step (the leave-one-out test of the clade criterion, the HMM
cross-check, PMSF, codon composition, OMArk) exists to make this table defensible,
not to replace it. This script builds it.

Two files:
  FINAL_ORTHOLOG_TABLE.tsv  every one of the 2,726 candidates with its orthology call
  FINAL_ORTHOLOGS.tsv       the hits that remain once the non-orthologs are removed

Every candidate receives exactly one call, and the partition is checked:
  CURATED_CONTROL  declared an ortholog before any analysis; not decided here
  ORTHOLOG         the tree places it with the controls (WITH_CONTROLS)
  NON_ORTHOLOG     the tree places it outside the control clade (CONTAMINANT)
  AMBIGUOUS        the criterion did not resolve it
  NOT_PLACEABLE    fewer than 20 residues survived the domain crop, so it was never
                   in the tree. This is not a call: absence from the tree says nothing
  PRESENCE_ONLY    from a proteome with more than one sequence per gene, so its tips
                   are not independent. The hit is real and kept; it cannot be placed
  USED_AS_BOUNDARY_MARKER
                   a candidate promoted into the declared set of boundary markers,
                   the reference the criterion measures against. A tip cannot be
                   judged by an edge it helps to define, so the tree returns no call
                   for it. There are 19 such candidates, of mixed screen status (14
                   removed, 5 retained), so this is not a hidden non-ortholog call
  NOT_TESTED       the target's tree has no declared boundary marker (a hits-only
                   tree), so the criterion was not tested
  PENDING          the target's tree has not finished

The tree criterion decides, and that is not an arbitration between methods. It is the
clade criterion of van Dam et al. (2011), which this analysis reproduces, and it is the
only method with a measured error rate at the boundary on every target (leave-one-out:
sensitivity 88 to 97%, specificity 83 to 100%). The other methods are declared
cross-checks and are carried as columns beside the call, never merged into it. A reader
who wants to apply a different rule has every number needed to do so.

Where the tree abstained, the two-arm rule decides (see rl_ambiguous_arm_rule.py),
except for a target whose orthology arm produced no calls, where the rows stay
AMBIGUOUS. The target list and file names come from the shared helper rl_targets; only
a run on the eight core targets writes the unsuffixed final table, and a run on another
target set never does.

Why the HMM screen may not override the tree. Measured on 8 September 2026: its
specificity at the boundary ranges from 100% (RAPTOR) to 30% (SIN1), it has no measured
error rate at all for the candidates, and for AMPKalpha, AMPKbeta and SIN1 its winning
target model is built from the same curated controls on which the tree criterion is
anchored, so for those three it is not an independent method. Its column is reported;
it does not vote.

The script stops unless every candidate is accounted for exactly once.
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import csv
import os
import sys

ROOT = os.environ.get("PROJECT_ROOT") or (
    (_TA210_SCRATCH + "/ms210_phase2/" "TOR-AMPK_orthologs_10target_R1"))
sys.path.insert(0, os.path.join(ROOT, "01_scripts", "lib"))
import rl_targets                                                # noqa: E402

# Target set. The target registry names the targets in use; the file-name suffix is
# empty for the eight core targets and "__<set>" otherwise. A run on another target
# set must read a suffixed tree track and report folder, and only the eight core
# targets may write the unsuffixed final table.
SFX = rl_targets.suffix()
TRACK = os.environ.get("LOO_BASE") or ("11_phylo_vandam_genecollapsed" + SFX)
rl_targets.check_track_agreement()
BASE = os.path.join(ROOT, TRACK)
TAG = os.environ.get("R59_TAG") or "setA"
RPT = rl_targets.report_dir(ROOT, TAG)
CAND = os.path.join(RPT, "06_candidate_list.tsv")
SCREEN = os.path.join(RPT, "07_contaminant_screen.tsv")
ARM = os.path.join(RPT, "14_arm_rule_calls.tsv")
# The two-arm rule can be switched off, and the setting in use is written into the
# output. ARM_RULE=none reproduces the table as it was before the rule, except for
# the arm-support columns, which are always added because the source of a
# nomination is not a call.
ARM_RULE = (os.environ.get("ARM_RULE") or "R66").strip()
if ARM_RULE not in ("R66", "none"):
    sys.exit("REFUSED: ARM_RULE=%s is not R66 or none" % ARM_RULE)
# The output name carries the tree track, because both tracks used to write this
# one path and whichever ran last determined the final table (on 9 September 2026
# the boundary track wrote it at 01:21, the gene-collapsed track at 07:29 and the
# boundary track again at 08:22). A checksum taken of that file would only have
# fixed one arbitrary moment of a race.
#
# The unsuffixed name is the gene-collapsed copy, and that follows from the methods
# rather than being a choice: a proteome with more than one sequence per gene has
# non-independent tips and is excluded from trees, clustering, copy number and
# paralogy. ORTHOLOG and NON_ORTHOLOG here come from the tree criterion, so
# including those tips would break that rule. The boundary track is kept as the
# secondary track; it holds the boundary markers and the leave-one-out test.
PRIMARY_TRACK = "11_phylo_vandam_genecollapsed" + SFX
SHORT = TRACK.replace("11_phylo_vandam_", "").replace("11_phylo_vandam", "vandam")
OUT = os.path.join(ROOT, "04_final_tables",
                   "FINAL_ORTHOLOG_TABLE__%s.tsv" % SHORT)
HITS = os.path.join(ROOT, "04_final_tables", "FINAL_ORTHOLOGS__%s.tsv" % SHORT)
# For another target set the primary copy carries the suffix; only the eight core
# targets use the bare name.
OUT_PRIMARY = os.path.join(ROOT, "04_final_tables", "FINAL_ORTHOLOG_TABLE%s.tsv" % SFX)
HITS_PRIMARY = os.path.join(ROOT, "04_final_tables", "FINAL_ORTHOLOGS%s.tsv" % SFX)
TARGETS = rl_targets.targets()

TREE_TO_DISP = {"WITH_CONTROLS": "ORTHOLOG",
                "CONTAMINANT": "NON_ORTHOLOG",
                "AMBIGUOUS": "AMBIGUOUS"}
# NON_ORTHOLOG_ARM_RULE is a separate value from NON_ORTHOLOG and is never summed
# with it. One comes from the tree criterion, which has a measured error rate at the
# boundary; the other comes from the two-arm rule, which does not perform well on the
# held-out controls. Merging them would make a call that failed its test on the
# controls look like one that passed it.
IS_ORTHOLOG = {"CURATED_CONTROL": "yes", "ORTHOLOG": "yes",
               "NON_ORTHOLOG": "no", "NON_ORTHOLOG_ARM_RULE": "no",
               "AMBIGUOUS": "unresolved",
               "NOT_PLACEABLE": "unresolved", "PRESENCE_ONLY": "unresolved",
               "USED_AS_BOUNDARY_MARKER": "unresolved", "PENDING": "unresolved",
               "NOT_TESTED": "unresolved"}


def rows(p):
    if not os.path.exists(p):
        return []
    return list(csv.DictReader((l for l in open(p, errors="replace")
                                if not l.startswith("##")), delimiter="\t"))


def main():
    cands = rows(CAND)
    if not cands:
        sys.exit("REFUSED: %s is missing or empty" % CAND)
    hmm = {(r["TARGET"], r["PROTEIN_ID"]): r for r in rows(SCREEN)}
    arm = {(r["TARGET"], r["PROTEIN_ID"]): r for r in rows(ARM)}
    if not arm:
        sys.exit("REFUSED: %s is missing. Run rl_ambiguous_arm_rule.py first -- the\n"
                 "arm provenance columns are not optional, they are how a reader\n"
                 "sees why a protein was called." % ARM)
    lb = {}
    for r in rows(os.path.join(BASE, "LB_SCORES.tsv")):
        if r.get("LB_SCORE"):
            lb[(r["TARGET"], r["PROTEIN_ID"])] = r["LB_SCORE"]

    out, per = [], {}
    for t in TARGETS:
        d = os.path.join(BASE, t)
        finished = os.path.exists(os.path.join(d, "iqtree", "%s.iqtree" % t))
        tree = {r["PROTEIN_ID"]: r for r in rows(
            os.path.join(d, "%s_contaminants_iqtree.tsv" % t))}
        # CLASS=BOUNDARY is in _subtree_calls_iqtree.tsv (all tips), not in the
        # contaminants file (per-tip calls only); that is why these 19 candidates had
        # no call and were first counted as PENDING.
        bnd = {r["PROTEIN_ID"] for r in rows(
            os.path.join(d, "%s_subtree_calls_iqtree.tsv" % t))
            if (r.get("CLASS") or "").strip() == "BOUNDARY"}
        pmsf = {r["PROTEIN_ID"]: r["TREE_CALL"] for r in rows(
            os.path.join(d, "%s_contaminants_pmsf.tsv" % t))}
        labels = {r["PROTEIN_ID"]: r["ROLE"] for r in rows(
            os.path.join(d, "%s_labels.tsv" % t))}
        excl = {r["PROTEIN_ID"] for r in rows(
            os.path.join(d, "%s_excluded_multirecord.tsv" % t))}
        crop = {r["PROTEIN_ID"] for r in rows(
            os.path.join(d, "%s_cropped_out.tsv" % t))}
        # Hits-only track (LKB1 and SNRK). With no declared boundary marker,
        # rl_phylo_subtree.py writes no per-tip call at all (its call block sits under
        # `if bnd:` and the else branch prints NOT TESTED), so a finished tree here
        # has no _contaminants file and no BOUNDARY class. That is not "UNEXPLAINED":
        # the criterion never examined the candidate, and it must not fall through to
        # PENDING, which the presence matrix reads as a finished target and then codes
        # as absent. Absence is evidence only where the sensitivity of the search is
        # known, and a method with no opinion is never a vote against, so an untested
        # absence is not allowed.
        hits_only = (finished and not bnd and not tree and
                     not any(v.startswith("BOUNDARY") for v in labels.values()))

        for c in [x for x in cands if x["TARGET"] == t]:
            pid, sp = c["PROTEIN_ID"], c["SPECIES"]
            role = labels.get(pid, "")
            tr = tree.get(pid, {})
            call = (tr.get("TREE_CALL") or "").strip()
            if role in ("CONTROL", "CONTROL_HELDOUT"):
                disp, basis = "CURATED_CONTROL", "declared before analysis"
            elif pid in excl:
                disp, basis = "PRESENCE_ONLY", ">1 record per gene; tips not independent"
            elif pid in crop:
                disp, basis = "NOT_PLACEABLE", "<20 residues survived the domain crop"
            elif not finished:
                disp, basis = "PENDING", "tree not finished"
            elif hits_only:
                disp = "NOT_TESTED"
                basis = ("hits-only tree (L35): no declared boundary marker, so the "
                         "criterion was NOT TESTED; neither an ortholog nor an absence")
            elif pid in bnd:
                disp = "USED_AS_BOUNDARY_MARKER"
                basis = ("promoted into the declared boundary set; the criterion "
                         "measures against it")
            elif call in TREE_TO_DISP:
                disp, basis = TREE_TO_DISP[call], "tree criterion (van Dam 2011)"
            else:
                disp = "PENDING"
                basis = "UNEXPLAINED: finished tree, no call, not a boundary marker"
            a = arm.get((t, pid), {})
            src = "TREE" if disp in ("ORTHOLOG", "NON_ORTHOLOG") else \
                  "CURATION" if disp == "CURATED_CONTROL" else "BOOKKEEPING"
            # ---- The two-arm rule, applied only where the tree abstained. The tree
            # is never overruled: TREE_CALL keeps the tree's own call beside this one.
            if disp == "AMBIGUOUS" and ARM_RULE == "R66":
                if (a.get("ARM1_INCOMPLETE") or "").strip() == "yes":
                    # Decided on 9 September 2026: the two-arm rule does not apply
                    # where the orthology arm has no acceptable setting for the
                    # target (every setting called more proteins than the copy-number
                    # ceiling allows). The rule is consistent with "a method with no
                    # opinion is never a vote against" only because absence from a
                    # call file means "not nominated under this method's rule";
                    # where no setting was chosen no rule is in use, that reading
                    # fails, and the tree's abstention stands. These rows stay
                    # AMBIGUOUS.
                    basis = ("R66 not applicable: arm 1 has no admissible operating "
                             "point (R64); a muted arm is no opinion (R40)")
                    src = "R66_NOT_APPLICABLE"
                else:
                    disp = ("ORTHOLOG" if a.get("ARM_SUPPORT") == "BOTH"
                            else "NON_ORTHOLOG_ARM_RULE")
                    basis = ("R66 arm support (%s); the tree returned no call"
                             % (a.get("ARM_SUPPORT") or "?"))
                    src = "R66"
            out.append([t, sp, pid, disp, IS_ORTHOLOG[disp], basis, call or "-",
                        TRACK,
                        (hmm.get((t, pid)) or {}).get("CALL", "-"),
                        pmsf.get(pid, "-"),
                        lb.get((t, pid), "-"),
                        a.get("ARM1_INSTRUMENTS", "-"),
                        a.get("ARM2_INSTRUMENTS", "-"),
                        a.get("ARM_SUPPORT", "-"),
                        a.get("ARM2_EVIDENCE", "-"),
                        src,
                        a.get("ARM2_INCOMPLETE", "-")])
            per.setdefault(t, {}).setdefault(disp, 0)
            per[t][disp] += 1

    # ---- Stop unless every candidate appears exactly once.
    if len(out) != len(cands):
        sys.exit("REFUSED: %d rows emitted for %d candidates -- every candidate must "
                 "appear exactly once." % (len(out), len(cands)))
    seen = set()
    for r in out:
        k = (r[0], r[2])
        if k in seen:
            sys.exit("REFUSED: %s/%s appears twice" % k)
        seen.add(k)

    # ---- The two-arm rule may change nothing but the tree's abstentions. Applied
    # more widely it would demote every declared boundary marker and almost every
    # NOT_PLACEABLE candidate, so this is checked rather than assumed.
    stray = [(r[0], r[2], r[3]) for r in out
             if r[15] == "R66" and r[3] not in ("ORTHOLOG", "NON_ORTHOLOG_ARM_RULE")]
    if stray:
        sys.exit("REFUSED: R66 wrote a disposition outside its scope: %s" % stray[:5])
    n66 = sum(1 for r in out if r[15] == "R66")
    # A row left AMBIGUOUS because the two-arm rule does not apply (no orthology-arm
    # calls for the target) is the intended outcome, not a leftover; any other
    # AMBIGUOUS row still stops the run.
    n_na = sum(1 for r in out if r[15] == "R66_NOT_APPLICABLE")
    namb = sum(1 for r in out if r[3] == "AMBIGUOUS" and r[15] != "R66_NOT_APPLICABLE")
    if ARM_RULE == "R66" and namb:
        sys.exit("REFUSED: ARM_RULE=R66 but %d rows are still AMBIGUOUS" % namb)
    if ARM_RULE == "none" and n66:
        sys.exit("REFUSED: ARM_RULE=none but %d rows claim R66 as their source" % n66)

    with open(OUT, "w") as fh:
        fh.write("## THE FINAL ORTHOLOG TABLE. Every candidate hit across 210 "
                 "eukaryotic proteomes, with its disposition. IS_ORTHOLOG=no is a "
                 "hit REMOVED as a non-ortholog; yes is a retained ortholog; "
                 "unresolved is neither, and the DISPOSITION says why.\n")
        fh.write("## The call is the TREE CRITERION (van Dam 2011), the method this "
                 "project reproduces and the only instrument with a measured error "
                 "rate at the boundary on every target. HMM_CALL, PMSF_CALL and "
                 "LB_SCORE are DECLARED CROSS-CHECKS carried beside the call and "
                 "never merged into it -- no vote is taken.\n")
        fh.write("## NOT_PLACEABLE and PRESENCE_ONLY ARE NOT VERDICTS. The first "
                 "carried too little of the domain to be placed; the second comes "
                 "from a proteome with more than one record per gene, so its tips "
                 "are not independent observations. Both hits are real and retained "
                 "on the record.\n")
        fh.write("## Written by 01_scripts/rl_final_ortholog_table.py, which refuses "
                 "unless every candidate is accounted for exactly once.\n")
        fh.write("## ARM1_INSTRUMENTS / ARM2_INSTRUMENTS / ARM_SUPPORT say HOW THE "
                 "PROTEIN WAS NOMINATED, on every row. ARM2_EVIDENCE carries each "
                 "domain instrument's score and the cut it already cleared, "
                 "verbatim; nothing compares it. DISPOSITION_SOURCE names the layer "
                 "that wrote DISPOSITION, and TREE_CALL always keeps the tree's own "
                 "word whatever that layer did.\n")
        fh.write("## ARM_RULE IN FORCE: %s. Under R66 a candidate the tree left "
                 "AMBIGUOUS is ORTHOLOG when at least one instrument on EACH arm "
                 "nominated it, and NON_ORTHOLOG_ARM_RULE otherwise. THAT IS A "
                 "DECLARED, SCOPED EXCEPTION TO R50, author-instructed 2026-09-08 "
                 "for the tree's abstention residue only. It recovers 12 of 39 "
                 "HELD-OUT curated controls and promotes 0 of 323 declared boundary "
                 "markers (09_per_method_report/14_arm_rule_benchmark.tsv); it does "
                 "NOT validate and may not be quoted without that number.\n"
                 % ARM_RULE)
        fh.write("## NON_ORTHOLOG_ARM_RULE IS NEVER SUMMED WITH NON_ORTHOLOG. The "
                 "first is R66, unvalidated; the second is the tree criterion, which "
                 "has a measured error rate at the boundary. Two criteria, two "
                 "values, never one count.\n")
        if n_na:
            fh.write("## %d ROW(S) STAY AMBIGUOUS WITH DISPOSITION_SOURCE="
                     "R66_NOT_APPLICABLE (DECISIONS D9): their target's arm 1 has no "
                     "admissible operating point under R64 (ARM1_INCOMPLETE=yes in "
                     "14_arm_rule_calls.tsv), so R66's two-arm support cannot be "
                     "BOTH by construction. A muted arm is no opinion (R40), not a "
                     "vote against; the tree's abstention stands.\n" % n_na)
        if SFX:
            fh.write("## TARGET SET %s (DECISIONS L2, L9): the replacement-layer "
                     "targets of layers 3-5, scored with the identical method and "
                     "written into this set-suffixed table; the unsuffixed "
                     "FINAL_ORTHOLOG_TABLE.tsv is the eight core targets and is NOT "
                     "touched by a set run. LAYER-5 SIGN INVERSION (DECISIONS I7): "
                     "SDS23, CBSX, SNRK, SnRK2, SnRK3_CIPK and MARK are detected by "
                     "HMMs built as AMPKa/AMPKg CONTAMINANT poles; the measurement is "
                     "unchanged and is read here as evidence FOR the family. Their "
                     "tree calls are anchored on their own declared controls. Arm 2 "
                     "runs at two-thirds strength for this set exactly as for the "
                     "eight (G-h1, D4).\n" % rl_targets.set_name())
        fh.write("TARGET\tSPECIES\tPROTEIN_ID\tDISPOSITION\tIS_ORTHOLOG\tBASIS\t"
                 "TREE_CALL\tTREE_TRACK\tHMM_CALL\tPMSF_CALL\tLB_SCORE\t"
                 "ARM1_INSTRUMENTS\tARM2_INSTRUMENTS\tARM_SUPPORT\tARM2_EVIDENCE\t"
                 "DISPOSITION_SOURCE\tARM2_INCOMPLETE\n")
        for r in out:
            fh.write("\t".join(str(x) for x in r) + "\n")

    keep = [r for r in out if r[4] == "yes"]
    with open(HITS, "w") as fh:
        fh.write("## THE FINAL HITS TABLE -- non-orthologs removed. Every retained "
                 "ortholog: curated controls plus the candidates the tree criterion "
                 "places with them. Hits called NON_ORTHOLOG, and those left "
                 "unresolved, are in FINAL_ORTHOLOG_TABLE.tsv with the reason.\n")
        fh.write("## TREE_CALL, DISPOSITION_SOURCE and ARM_SUPPORT are carried here "
                 "because without them an R66 promotion is indistinguishable from a "
                 "tree-decided ortholog in the list people actually use.\n")
        fh.write("TARGET\tSPECIES\tPROTEIN_ID\tDISPOSITION\tBASIS\tHMM_CALL\t"
                 "PMSF_CALL\tLB_SCORE\tTREE_CALL\tDISPOSITION_SOURCE\t"
                 "ARM_SUPPORT\n")
        for r in keep:
            fh.write("\t".join([r[0], r[1], r[2], r[3], r[5], r[8], r[9], r[10],
                                r[6], r[15], r[13]]) + "\n")

    # Only the primary track writes the unsuffixed name. A run on another track leaves
    # it alone, which stops the two tracks from overwriting each other.
    if TRACK == PRIMARY_TRACK:
        import shutil
        shutil.copyfile(OUT, OUT_PRIMARY)
        shutil.copyfile(HITS, HITS_PRIMARY)
        print("primary track (%s): also wrote %s"
              % (rl_targets.label(), os.path.basename(OUT_PRIMARY)))
    else:
        print("NON-PRIMARY track %s: the un-suffixed FINAL_ORTHOLOG_TABLE.tsv was "
              "NOT touched. The primary is %s (METHODS section 9)."
              % (SHORT, PRIMARY_TRACK))

    order = ["CURATED_CONTROL", "ORTHOLOG", "NON_ORTHOLOG", "NON_ORTHOLOG_ARM_RULE",
             "AMBIGUOUS",
             "NOT_PLACEABLE", "PRESENCE_ONLY", "USED_AS_BOUNDARY_MARKER", "PENDING"]
    print("%-8s %6s | %s" % ("TARGET", "cand", "  ".join("%-14s" % o for o in order)))
    print("-" * 128)
    tot = {}
    for t in TARGETS:
        p = per.get(t, {})
        n = sum(p.values())
        print("%-8s %6d | %s" % (t, n, "  ".join("%-14d" % p.get(o, 0) for o in order)))
        for o in order:
            tot[o] = tot.get(o, 0) + p.get(o, 0)
    print("-" * 128)
    print("%-8s %6d | %s" % ("ALL", len(out),
                             "  ".join("%-14d" % tot.get(o, 0) for o in order)))
    print("\nRETAINED ORTHOLOGS: %d of %d candidates (%d curated controls + %d "
          "placed with them)" % (len(keep), len(out), tot.get("CURATED_CONTROL", 0),
                                 tot.get("ORTHOLOG", 0)))
    print("REMOVED AS NON-ORTHOLOGS: %d" % tot.get("NON_ORTHOLOG", 0))
    unres = sum(tot.get(o, 0) for o in ("AMBIGUOUS", "NOT_PLACEABLE",
                                        "PRESENCE_ONLY", "USED_AS_BOUNDARY_MARKER",
                                        "PENDING"))
    print("UNRESOLVED: %d (ambiguous %d, not placeable %d, presence-only %d, "
          "boundary-marker %d, pending %d)"
          % (unres, tot.get("AMBIGUOUS", 0), tot.get("NOT_PLACEABLE", 0),
             tot.get("PRESENCE_ONLY", 0), tot.get("USED_AS_BOUNDARY_MARKER", 0),
             tot.get("PENDING", 0)))
    print("\nwrote %s\n      %s" % (OUT, HITS))
    if tot.get("PENDING", 0):
        print("\nTABLE IS NOT FINAL: %d candidate(s) PENDING on unfinished trees."
              % tot["PENDING"])


if __name__ == "__main__":
    main()
