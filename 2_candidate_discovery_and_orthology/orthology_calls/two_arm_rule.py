#!/usr/bin/env python3
# >>> deposit header
# Purpose: For every candidate, lists which orthology-arm and domain-arm methods nominated it, writes the
#          call the two-arm rule would make from that support, and scores the rule on held-out controls.
# Methods: Orthology calls
# Inputs:  06_candidate_list.tsv; the per-method call files 0[1-9]_*_calls.tsv of the per-method report
#          folder; the boundary markers in <target>_subtree_calls_iqtree.tsv of the boundary tree track
# Outputs: 14_arm_rule_calls.tsv (arm support and two-arm call per candidate); 14_arm_rule_benchmark.tsv
#          (held-out controls recovered and boundary markers wrongly promoted, per target)
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
"""Support from the two search arms for every candidate, and the call it implies.

The gene tree returns AMBIGUOUS for a candidate when the smallest clade that encloses
it and contains any declared tip holds both a positive control and a boundary marker,
that is, when the topology does not separate the two reference sets at that point. This script writes, for every candidate, which methods of each arm nominated it:
the orthology arm (OrthoFinder, SonicParanoid) and the domain arm (CDD, the HMMER
profile panel, InterProScan). It also writes the call the two-arm rule would make from
that support.

The script itself decides nothing. It only writes the column `R66_CALL`. Whether that
call is applied is decided in rl_final_ortholog_table.py, only for candidates the tree
left AMBIGUOUS, and it can be switched off with ARM_RULE=none. The list of targets and
the names of the output folders come from the shared helper targets.

The two-arm rule (chosen on 8 September 2026):
    nominated by at least one orthology-arm method
    and by at least one domain-arm method              -> ORTHOLOG
    nominated by one arm only                          -> NON_ORTHOLOG_ARM_RULE

The rule is a limited exception to the general design, in which the two arms are
independent, every method is reported, and phylogeny is a third line of evidence that
never decides membership on its own. Here the arms decide, but only where the tree
abstained, so the rule sets an order of priority; it is a k-of-n threshold across the
arms. It applies to candidates the tree left unresolved and never to the assembly of
the candidate list. Two further tensions are stated rather than argued away:

  - A line of evidence may decide only if it performs well on the declared controls,
    and lines of evidence are never tested against each other. The two-arm rule does
    not pass that test: of the 39 held-out curated controls it recovers 12 (30.8%).
    The script prints that number on every run so that it travels with the calls.
  - A method with no opinion is never a vote against. The two-arm rule turns absence
    from one arm into a demotion. The two principles agree only under the definition
    the per-method call files state for themselves: absence from a call file means
    "not nominated under this method's membership rule", a defined negative opinion,
    not "the method did not see the sequence". That definition is the one used here.

No score is compared here. The bit-score condition of the original instruction was
dropped after it was measured: 274 of the 275 AMBIGUOUS candidates already had the
target model scoring above the contaminant model, because the upstream winner-take-all
screen had filtered on exactly that, so as a yes/no condition it changed no call. A cut
on the size of the margin would separate the classes, but it would be a new free
parameter fitted to the same tree that abstained. The column `ARM2_EVIDENCE` carries
each domain-arm method's score and the cutoff it already passed, copied unchanged, and
nothing in this script compares them.

There are deliberately no integer columns counting the methods per arm: a count column
invites a k-of-n rule across methods, which the design does not use. The method names
are joined by semicolons; count them if a number is needed.
"""
import collections
import csv
import glob
import os
import sys

def _project_root_from_file(start):
    """Find the project root: the nearest folder above this file that holds 05_method_and_rules.

    A script in a copy of the project tree must never fall back to, or hard-code, the original
    tree (chosen on 16 September 2026). The function refuses rather than guesses: with a
    hard-coded root, a run against the wrong tree could not be told apart from a correct one, and
    in the unconditional form even setting PROJECT_ROOT did not help.

    The path is resolved with abspath, not realpath, on purpose. The copy holds 2,014 symbolic
    links that point into the original tree. Walking up from the link finds the tree the script is
    being run in; walking up from the link's target would find the tree the file was copied from,
    which for a link into the original tree is exactly the wrong-root error this function prevents.
    """
    d = os.path.dirname(os.path.abspath(start))
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "05_method_and_rules")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("REFUSED: PROJECT_ROOT unset and no project root found above %s" % start)


ROOT = os.environ.get("PROJECT_ROOT") or _project_root_from_file(__file__)
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "lib"))
import targets                                                # noqa: E402

TAG = os.environ.get("R59_TAG") or "setA"
SFX = targets.suffix()
targets.check_track_agreement()
RPT = targets.report_dir(ROOT, TAG)
CAND = os.path.join(RPT, "06_candidate_list.tsv")
OUT = os.path.join(RPT, "14_arm_rule_calls.tsv")
BENCH = os.path.join(RPT, "14_arm_rule_benchmark.tsv")

# The methods of each arm are listed explicitly: the orthology arm (OrthoFinder,
# SonicParanoid) and the domain arm (CDD, the HMMER profile panel, InterProScan). A
# method that appears in the call files but is not named here stops the run instead of
# joining an arm by accident; a new method enters the two-arm rule only on purpose,
# after its effect on the calls has been measured.
ARM1 = ["orthofinder_full", "sonic"]
ARM2 = ["cdd", "hmmer", "interpro"]
TARGETS = targets.targets()


def rows(p):
    if not os.path.exists(p):
        return []
    return list(csv.DictReader((l for l in open(p, errors="replace")
                                if not l.startswith("##")), delimiter="\t"))


def main():
    cands = rows(CAND)
    if not cands:
        sys.exit("REFUSED: %s is missing or empty" % CAND)

    # ---- Read the call file of each method. Support means CALL == MEMBER, tested by
    # name, and never ROLE: ROLE says what the protein is (CANDIDATE, CONTROL or
    # CONTROL_HELDOUT), not what the method said. Testing the CALL value by name means
    # that a future value such as NOT_MEMBER or UNDETERMINED cannot silently count as
    # support.
    inst = collections.defaultdict(set)
    ev = collections.defaultdict(dict)
    role = {}
    seen_methods, delivered = set(), collections.Counter()
    delivered_per_target = collections.Counter()      # MEMBER rows per (target, method)
    for p in sorted(glob.glob(os.path.join(RPT, "0[1-9]_*_calls.tsv"))):
        for r in rows(p):
            m = (r.get("METHOD") or "").strip()
            seen_methods.add(m)
            delivered[m] += 1
            k = (r["TARGET"], r["PROTEIN_ID"])
            role.setdefault(k, (r.get("ROLE") or "").strip())
            if (r.get("CALL") or "").strip() != "MEMBER":
                continue
            inst[k].add(m)
            delivered_per_target[(r["TARGET"], m)] += 1
            if m in ARM2:
                bits = ";".join(x for x in (r.get("EVIDENCE_2"), r.get("EVIDENCE_3"))
                                if x and x != "-")
                if bits:
                    ev[k][m] = bits

    # ---- Every method found in the call files must belong to one of the two arms.
    unrostered = sorted(seen_methods - set(ARM1) - set(ARM2))
    if unrostered:
        sys.exit("REFUSED: method(s) %s appear in the call files but are on neither "
                 "arm roster. A new instrument enters R66 deliberately, never by "
                 "being absorbed. Add it to ARM1 or ARM2 and re-measure the "
                 "consequence." % ", ".join(unrostered))

    # ---- A listed domain-arm method that produced no rows at all has failed; it has
    # not merely found nothing. This is marked on every row (ARM2_INCOMPLETE=yes) but
    # does not stop the run: stopping here would block the chosen rule on a condition
    # that held for the whole analysis.
    dead2 = [m for m in ARM2 if delivered.get(m, 0) == 0]
    arm2_incomplete = "yes" if dead2 else "no"

    # ---- `ARM1_INCOMPLETE`, per target (decided on 9 September 2026). A target for
    # which every orthology-arm method produced zero MEMBER rows has no usable
    # orthology-arm setting at all: every clustering setting tried called more
    # proteins than the copy-number ceiling allows (a setting calling more than the
    # ceiling times 210 proteins is not acceptable, however few non-targets it calls),
    # so no setting was chosen for it. CAMKK is the measured case: no setting on disk
    # gives a separable orthogroup. For such a target `ARM_SUPPORT` can never be
    # `BOTH`, so the two-arm rule would demote every candidate the tree left
    # unresolved, by construction. The target is marked here, and
    # rl_final_ortholog_table.py leaves those rows AMBIGUOUS, because an arm with no
    # calls has no opinion and is not a vote against. Each of the eight core targets
    # has OrthoFinder and SonicParanoid rows and reads "no".
    arm1_incomplete = {t: ("yes" if all(delivered_per_target.get((t, m), 0) == 0
                                        for m in ARM1) else "no")
                       for t in TARGETS}
    muted1 = sorted(t for t in TARGETS if arm1_incomplete[t] == "yes")

    def support(k):
        ms = inst.get(k, set())
        a1 = sorted(ms & set(ARM1))
        a2 = sorted(ms & set(ARM2))
        s = ("BOTH" if a1 and a2 else "ARM1_ONLY" if a1 else
             "ARM2_ONLY" if a2 else "NEITHER")
        return a1, a2, s

    out = []
    for c in cands:
        k = (c["TARGET"], c["PROTEIN_ID"])
        a1, a2, s = support(k)
        # ---- Support from neither arm cannot occur, because the candidate list is
        # the union of the call files. If it does occur, the list and the call files
        # have diverged, and a plausible-looking NON_ORTHOLOG call here would hide
        # that error.
        if s == "NEITHER":
            sys.exit("REFUSED: %s/%s is a candidate with no MEMBER row in any call "
                     "file. The candidate list is built as the union of those files, "
                     "so this means they have diverged. Fix the build; do not emit a "
                     "call for it." % k)
        r66 = "ORTHOLOG" if s == "BOTH" else "NON_ORTHOLOG_ARM_RULE"
        if arm1_incomplete[k[0]] == "yes":
            r66 = "NOT_APPLICABLE_ARM1_MUTED"
        out.append([k[0], k[1], c.get("SPECIES", "-"), role.get(k, "-"),
                    ";".join(a1) or "-", ";".join(a2) or "-", s,
                    "|".join("%s:%s" % (m, ev[k][m]) for m in sorted(ev[k])) or "-",
                    r66, arm2_incomplete, arm1_incomplete[k[0]]])

    with open(OUT, "w") as fh:
        fh.write("## ARM SUPPORT PER CANDIDATE, and the call R66 would make from it. "
                 "R66_CALL is written for EVERY candidate so the rule's opinion is "
                 "auditable everywhere; it is APPLIED only where the tree said "
                 "AMBIGUOUS, by rl_final_ortholog_table.py, and only while "
                 "ARM_RULE=R66.\n")
        fh.write("## SUPPORT means CALL=MEMBER in that method's own call file. Per "
                 "those files' banner, absence means NOT NOMINATED UNDER THIS RULE -- "
                 "a defined negative opinion of the method's membership rule, not "
                 "'the instrument did not see it'. That reading is what reconciles "
                 "R66 with R40, and it is the only one that does.\n")
        fh.write("## ARM2_EVIDENCE carries each arm-2 instrument's score AND THE CUT "
                 "IT ALREADY CLEARED, verbatim. NOTHING HERE COMPARES IT. Arm-2 "
                 "support already IS a bit-score test at the R59 operating point; a "
                 "second cut above it would be a new free parameter.\n")
        fh.write("## R66 IS A DECLARED, SCOPED EXCEPTION TO R50, author-instructed "
                 "2026-09-08 for the tree's abstention residue only. It does not "
                 "validate on the controls (12 of 39 held out; see "
                 "14_arm_rule_benchmark.tsv) and may not be quoted without that "
                 "number.\n")
        if dead2:
            fh.write("## ARM2_INCOMPLETE=yes ON EVERY ROW: rostered arm-2 "
                     "instrument(s) %s delivered zero rows. Arm 2 is running below "
                     "strength and every single-arm call rests on that.\n"
                     % ", ".join(dead2))
        if muted1:
            fh.write("## ARM1_INCOMPLETE=yes FOR %s: no arm-1 instrument delivered a "
                     "MEMBER row for that target, because every enumerated clustering "
                     "setting exceeded the R64 ceiling and the cell was reported as "
                     "inadmissible (r64_cap_emptied). R66 IS NOT APPLICABLE there "
                     "(DECISIONS D9): R66_CALL reads NOT_APPLICABLE_ARM1_MUTED and "
                     "rl_final_ortholog_table.py leaves those rows AMBIGUOUS. A muted "
                     "arm is no opinion (R40), not a vote against.\n"
                     % ", ".join(muted1))
        fh.write("TARGET\tPROTEIN_ID\tSPECIES\tROLE\tARM1_INSTRUMENTS\t"
                 "ARM2_INSTRUMENTS\tARM_SUPPORT\tARM2_EVIDENCE\tR66_CALL\t"
                 "ARM2_INCOMPLETE\tARM1_INCOMPLETE\n")
        for r in out:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # ================================================================ benchmark of the rule
    # The error rate is written next to the calls. A line of evidence is tested only
    # on the declared controls and never against another line of evidence, so the one
    # test used is this: the held-out controls (declared orthologs, kept out of the
    # settings of both arms) and the declared boundary markers (declared non-orthologs).
    ho = [k for k, v in role.items() if v == "CONTROL_HELDOUT"]
    bndp = {}
    # For the eight core targets the default is the boundary tree track, which holds
    # the declared boundary markers; a run on another target set has no boundary track
    # and must name its own track with LOO_BASE.
    track = os.environ.get("LOO_BASE") or "11_phylo_vandam_boundary"
    if SFX and not os.environ.get("LOO_BASE"):
        sys.exit("REFUSED: %s needs LOO_BASE set to the set's track (e.g. "
                 "11_phylo_vandam_genecollapsed%s)" % (targets.label(), SFX))
    for t in TARGETS:
        for r in rows(os.path.join(ROOT, track, t,
                                   "%s_subtree_calls_iqtree.tsv" % t)):
            if (r.get("CLASS") or "").strip() == "BOUNDARY":
                bndp[(t, r["PROTEIN_ID"])] = t

    lines = []
    tho = tre = tb = tfp = 0
    for t in TARGETS:
        hk = [k for k in ho if k[0] == t]
        bk = [k for k in bndp if k[0] == t]
        rec = sum(1 for k in hk if support(k)[2] == "BOTH")
        fp = sum(1 for k in bk if support(k)[2] == "BOTH")
        tho += len(hk); tre += rec; tb += len(bk); tfp += fp
        lines.append([t, len(hk), rec,
                      "%.1f" % (100.0 * rec / len(hk)) if hk else "-",
                      len(bk), fp,
                      "%.1f" % (100.0 * fp / len(bk)) if bk else "-"])

    with open(BENCH, "w") as fh:
        fh.write("## THE ONLY TEST R19 PERMITS. R19: a layer decides only if it "
                 "validates on the CONTROLS, and layers are never validated against "
                 "each other -- so R66 is NOT scored against the tree here. It is "
                 "scored on the HELD-OUT curated controls (declared orthologs, out of "
                 "sample for both arms) and the DECLARED BOUNDARY MARKERS (declared "
                 "non-orthologs).\n")
        fh.write("## R66 DOES NOT PASS THIS TEST AND IS NOT CLAIMED TO. It is "
                 "recorded so the rule can never be quoted without its error rate.\n")
        fh.write("## Concordance with the tree is NOT reported here and is not a "
                 "test: every candidate is arm-nominated by construction, so a "
                 "sensitivity against the tree restates the arm-support cross-tab it "
                 "was computed from.\n")
        fh.write("TARGET\tN_HELDOUT_CONTROL\tRECOVERED\tRECALL_PCT\t"
                 "N_BOUNDARY_MARKER\tWRONGLY_PROMOTED\tFP_PCT\n")
        for r in lines:
            fh.write("\t".join(str(x) for x in r) + "\n")
        fh.write("\t".join(str(x) for x in
                 ["ALL", tho, tre, "%.1f" % (100.0 * tre / tho) if tho else "-",
                  tb, tfp, "%.1f" % (100.0 * tfp / tb) if tb else "-"]) + "\n")

    n = collections.Counter(r[6] for r in out)
    print("arm support over %d candidates: %s" % (len(out), dict(n)))
    if muted1:
        print("  ARM1_INCOMPLETE=yes for %s: R66 not applicable there (D9); those rows "
              "stay AMBIGUOUS in the final table." % ", ".join(muted1))
    if dead2:
        print("  WARNING: rostered arm-2 instrument(s) %s delivered ZERO rows. "
              "Arm 2 is below strength and every ARM1_ONLY call rests on that."
              % ", ".join(dead2))
    print("\nR19-legal benchmark -- held-out controls and declared boundary markers:")
    print("  %-8s %9s %10s %9s %11s %9s" % ("TARGET", "N_HELDOUT", "RECOVERED",
                                            "RECALL", "N_BOUNDARY", "PROMOTED"))
    for r in lines:
        print("  %-8s %9d %10d %8s%% %11d %9d" % (r[0], r[1], r[2], r[3], r[4], r[5]))
    print("  %-8s %9d %10d %8s%% %11d %9d"
          % ("ALL", tho, tre, "%.1f" % (100.0 * tre / tho) if tho else "-", tb, tfp))
    print("\n  R66 recovers %d of %d held-out curated controls. It does not validate; "
          "this is recorded, not passed." % (tre, tho))
    print("wrote %s\n      %s" % (OUT, BENCH))


if __name__ == "__main__":
    main()
