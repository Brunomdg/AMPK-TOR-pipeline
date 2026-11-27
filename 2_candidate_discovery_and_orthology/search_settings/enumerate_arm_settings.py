#!/usr/bin/env python3
# >>> deposit header
# Purpose: Enumerate every OrthoFinder and SonicParanoid setting held on disk, with the declared controls each
#          recovers (leave-one-out) and the declared negatives it admits, so that the joint choice of one setting
#          per method and target has a complete, non-sampled frontier to solve over.
# Methods: Candidate discovery
# Inputs:  the OrthoFinder run folder and its saved similarity graph; the SonicParanoid run folders and their
#          pairwise ortholog tables; the controls, declared decoys and negative controls read through
#          `rl_arm1_selection`
# Outputs: `07_findings/R59_arm1_settings__<tag>.tsv` (every setting, with its cover and cost) and
#          `R59_arm1_cover__<tag>.tsv` (which controls each setting recovers, long form); <tag> is the control set
#          named in `R59_TAG` plus the target-set suffix
# Tools:   Python 3 standard library and the project modules `rl_arm1_selection` and `rl_targets`; reads
#          OrthoFinder and SonicParanoid results without running them; versions: see envs/
# <<< deposit header
# Step 2 of method 1, and the step `rl_r59_union_optimum.py` reads: that script solves the joint choice over the
# frontier enumerated here and stops with "Run rl_arm1_all_settings.py first" if this output is absent.
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
## GOVERNED BY: R59, L15 (output names carry the TARGET_SET suffix through lib/rl_targets.tagged)
"""Every OrthoFinder and SonicParanoid setting we hold, with its cover and cost.

WHY THIS EXISTS. R59 chooses settings jointly so the UNION recovers every declared
control at the lowest corpus-wide cost. That search is over INSTRUMENTS, not arms:
OrthoFinder and SonicParanoid are instruments like HMMER, CDD and InterPro and are
optimised on the same footing. The first R59 run carried them at their FLOOR/SELECT
setting instead of searching them, and the consequence was visible immediately --
every target's remaining cost was Arm 1, TOR's OrthoFinder cell alone contributing
1,593 of its 1,941 non-targets, 7.6 per species.

WHAT COUNTS AS A SETTING. Everything on disk:

  OrthoFinder, group rules   the base I=1.2 clustering, the original
                             orthofinder_full clustering, and the I=1.5, I=2.0 and
                             I=3.0 re-clusterings from the saved DIAMOND results
  OrthoFinder, graph rule    every edge-weight threshold in OrthoFinder's own
                             weighted graph
  SonicParanoid, group rules the cov0.75/I1.5 run, cov0.90 (full03), and the
                             inflation 2.0 and 3.0 runs
  SonicParanoid, pairwise    every confidence threshold in its own
                             species_to_species_orthologs tables, per run

THE THRESHOLD ENUMERATION IS COMPLETE, NOT SAMPLED. Cover changes only at a
CONTROL-TO-CONTROL edge weight: between two such weights the set of controls
linked to another control cannot change. And for a given cover the cheapest
setting is the HIGHEST threshold achieving it, which is that control-to-control
weight itself. So enumerating at exactly those weights enumerates every
non-dominated setting the instrument offers, and nothing in between is missed.

COVER IS LEAVE-ONE-OUT (R28). A control counts as recovered when it links to at
least one OTHER control of the same target. This is why cover cannot be read off
the member files: a control-sharing rule's member list holds the PARTNERS of the
controls and excludes the controls themselves, so intersecting it with the control
set answers ~0 by construction -- measured on the sweep members: 3 controls listed
where 42 are recovered. That mistake has been made in this project before and is
not repeated here.

COST IS DECLARED NEGATIVES ADMITTED, NOT EVERY CALL. An earlier version defined
cost as members minus the target's own controls -- every protein called that was
not already a declared control. That makes an UNDECLARED TRUE ORTHOLOG, the thing
this project exists to find, count as cost, so the cheapest setting is the one
that finds nothing: seven Arm 1 cells came back emitting zero members at cost
zero. Cost is now the declared decoys and adjudicated negative controls the
setting admits -- the false positives we can actually name. Finding an ortholog
costs nothing.

N_CALLED is carried separately, because it is not a cost but a PLAUSIBILITY
CHECK: a family expected at one copy per genome, with a ceiling of five, cannot
legitimately call 274 per species however few declared decoys it admits.
"""
import bisect
import collections
import csv
import glob
import os
import sys

ROOT = os.environ.get("PROJECT_ROOT") or (
    (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1"))
sys.path.insert(0, os.path.join(ROOT, "01_scripts"))
sys.path.insert(0, os.path.join(ROOT, "01_scripts", "lib"))
import rl_arm1_selection as A                                    # noqa: E402
import rl_invariants                                             # noqa: E402
import rl_targets                                                # noqa: E402

OF2 = (_TA210_SCRATCH + "/ofrun2/I1.2/Results_Sep06")
WD = OF2 + "/WorkingDirectory"
SP = (_TA210_SCRATCH + "/panel/sp2_out/runs")

# --- output tagging, so two control sets can be solved concurrently ------------
# rl_r59_union_optimum.py READS rl_arm1_all_settings.py's output, so untagged
# parallel runs would not merely overwrite each other's files, they would
# cross-feed: one control set's union solve reading the other's Arm 1 enumeration.
_TAG = os.environ.get("R59_TAG") or ""


def _tagged(p):
    if not _TAG:
        raise SystemExit(
            "REFUSING: R59_TAG is not set. Writing untagged R59 outputs is how "
            "stale files accumulate: a later reader cannot tell which control "
            "set produced them. Set R59_TAG=setA or R59_TAG=setB.")
    return rl_targets.tagged(p, _TAG)           # base__<tag><set suffix>.ext (L15)
# ------------------------------------------------------------------------------
OUT_SET = _tagged(ROOT + "/07_findings/R59_arm1_settings.tsv")
OUT_COV = _tagged(ROOT + "/07_findings/R59_arm1_cover.tsv")
TARGETS = A.TARGETS

GROUPS_OF = [
    ("orthofinder", "I1.2 (orthofinder_full)",
     ROOT + "/03_evidence_by_method/orthofinder_full/Orthogroups.txt"),
    ("orthofinder", "I1.2 (ofrun2)", OF2 + "/Orthogroups/Orthogroups.txt"),
    ("orthofinder", "I1.5",
     WD + "/OrthoFinder/Results_I1.5/Orthogroups/Orthogroups.txt"),
    ("orthofinder", "I2.0",
     WD + "/OrthoFinder/Results_I2.0/Orthogroups/Orthogroups.txt"),
    ("orthofinder", "I3.0",
     WD + "/OrthoFinder/Results_I3.0/Orthogroups/Orthogroups.txt"),
]
GROUPS_SP = [
    ("sonic", "cov0.75/I1.5",
     _TA210_SCRATCH + "/panel/BACKUP_full02_20260819/"
                      "ortholog_groups/flat.ortholog_groups.tsv"),
    ("sonic", "cov0.90 (full03)",
     SP + "/full03/ortholog_groups/flat.ortholog_groups.tsv"),
    ("sonic", "infl2.0", SP + "/infl20/ortholog_groups/flat.ortholog_groups.tsv"),
    ("sonic", "infl3.0", SP + "/infl30/ortholog_groups/flat.ortholog_groups.tsv"),
]


def rows(p):
    return csv.DictReader((l for l in open(p, errors="replace")
                           if not l.startswith("##")), delimiter="\t")


def load_controls():
    """The binding set, twice. `allc` was formerly EVERY declared row.

    Leave-one-out cover computed against every declared row asks each setting to
    recover controls the declaration has held out -- including, after the
    2026-09-07 curation, rows adjudicated as the wrong protein. It also ignored
    R59_CONTROL_SUBSET entirely. Both callers now get the binding set, which is
    the only set the project asserts ARE the target.
    """
    binding, label = A.binding_controls()
    binding = {t: set(binding.get(t, ())) for t in TARGETS}
    empty = [t for t in TARGETS if not binding[t]]
    if empty:
        raise SystemExit("REFUSING: no binding control for %s" % ", ".join(empty))
    print("control set: %s -- %d controls %s"
          % (label, sum(len(v) for v in binding.values()),
             {t: len(binding[t]) for t in TARGETS}))
    return {t: set(v) for t, v in binding.items()}, binding


def emit(fh, ch, t, inst, setting, source, members, cover, allc, spidx,
         neg=None):
    cand = members - allc[t]
    cost = len(cand & (neg or set()))          # DECLARED negatives admitted
    spp = len({spidx.get(a) for a in members} - {None})
    fh.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%d\t%.3f\n"
             % (t, inst, setting, source, len(members), len(cover), len(allc[t]),
                cost, spp, cost / 210.0, len(cand), len(cand) / 210.0))
    for a in sorted(cover):
        ch.write("%s\t%s\t%s\t%s\n" % (t, inst, setting, a))


def do_groups(fh, ch, spec, allc, binding, spidx, neg):
    inst, setting, path = spec
    if not os.path.exists(path):
        print("  MISSING, skipped and recorded: %s" % path)
        return
    gid = {}
    members_of = collections.defaultdict(set)
    if path.endswith("Orthogroups.txt"):
        for line in open(path, errors="replace"):
            g, _, rest = line.partition(":")
            g = g.strip()
            for tok in rest.split():
                a = A.norm(tok.split("|")[-1])
                gid[a] = g
                members_of[g].add(a)
    else:
        first = True
        for line in open(path, errors="replace"):
            f = line.rstrip("\n").split("\t")
            if first:
                first = False
                continue
            if len(f) < 5:
                continue
            g = f[0].strip()
            for col in f[4:]:
                for tok in col.replace(",", " ").split():
                    if tok in ("*", "-"):
                        continue
                    a = A.norm(tok.split("|")[-1])
                    gid[a] = g
                    members_of[g].add(a)
    for t in TARGETS:
        gs = {gid[c] for c in binding[t] if c in gid}
        members = set()
        for g in gs:
            members |= members_of[g]
        # leave-one-out: a control counts only if its group holds ANOTHER control
        cover = set()
        for c in allc[t]:
            g = gid.get(c)
            if g and len((members_of[g] & allc[t]) - {c}) >= 1:
                cover.add(c)
        emit(fh, ch, t, inst, setting, os.path.basename(path), members, cover,
             allc, spidx, neg[t])
    print("  %s %s: done" % (inst, setting))


def do_graph(fh, ch, allc, binding, spidx, neg):
    """OrthoFinder's own weighted graph, every edge-weight threshold."""
    seq = ROOT and (WD + "/SequenceIDs.txt")
    if not (os.path.exists(seq) and os.path.exists(WD + "/OrthoFinder_graph.txt")):
        print("  graph or SequenceIDs missing -- skipped")
        return
    idx2acc, acc2idx = {}, {}
    with open(seq, errors="replace") as fhh:
        for i, line in enumerate(fhh):
            if ":" not in line:
                continue
            a = A.norm(line.split(":", 1)[1].split()[0])
            idx2acc[i] = a
            acc2idx[a] = i
    cidx = {t: {acc2idx[c] for c in allc[t] if c in acc2idx} for t in TARGETS}
    bidx = {t: {acc2idx[c] for c in binding[t] if c in acc2idx} for t in TARGETS}
    owner = collections.defaultdict(set)
    for t in TARGETS:
        for j in cidx[t]:
            owner[j].add(t)
    print("  scanning OrthoFinder_graph.txt (8.2 GB, control rows only)")
    cc = collections.defaultdict(list)          # t -> [(w, i, j)] control-control
    nbr = collections.defaultdict(dict)         # t -> partner idx -> best w to a
    n = 0                                       #      BINDING control
    with open(WD + "/OrthoFinder_graph.txt", errors="replace") as fhh:
        for line in fhh:
            if not line or line[0] not in "0123456789":
                continue
            n += 1
            if n % 2000000 == 0:
                print("    %d rows" % n, flush=True)
            head, _, rest = line.partition(" ")
            try:
                row = int(head)
            except ValueError:
                continue
            ts = owner.get(row)
            if not ts:
                continue
            for tok in rest.split():
                if tok == "$" or ":" not in tok:
                    continue
                a_, _, b_ = tok.partition(":")
                try:
                    j, w = int(a_), float(b_)
                except ValueError:
                    continue
                for t in ts:
                    if j in cidx[t]:
                        cc[t].append((w, row, j))
                    if row in bidx[t] and j not in cidx[t]:
                        if w > nbr[t].get(j, -1.0):
                            nbr[t][j] = w
    print("    scanned %d rows" % n)
    for t in TARGETS:
        if not cc[t]:
            print("  graph %s: no control-to-control edge -- no setting" % t)
            continue
        ws = sorted({w for w, _i, _j in cc[t]}, reverse=True)
        # partners sorted by best weight to a BINDING control, descending, so the
        # member set at threshold w is a PREFIX of this list
        parts = sorted(nbr[t].items(), key=lambda kv: -kv[1])
        pw = [w for _j, w in parts]
        for w in ws:
            cover = set()
            for ww, i, j in cc[t]:
                if ww >= w:
                    cover.add(idx2acc.get(i))
                    cover.add(idx2acc.get(j))
            cover = {c for c in cover if c in allc[t]}
            k = bisect.bisect_left([-x for x in pw], -w)
            members = {idx2acc[j] for j, _w in parts[:k] if j in idx2acc}
            emit(fh, ch, t, "orthofinder", "of-graph-edge>=%.4f" % w,
                 "OrthoFinder_graph.txt", members, cover, allc, spidx, neg[t])
        print("  graph %s: %d thresholds" % (t, len(ws)))


def parse_side(toks):
    out = []
    for i in range(0, len(toks) - 1, 2):
        try:
            out.append((A.norm(toks[i].split("|")[-1]), float(toks[i + 1])))
        except ValueError:
            continue
    return out


def do_pairwise(fh, ch, run, allc, binding, spidx, neg):
    base = SP + "/" + run + "/species_to_species_orthologs"
    files = sorted(glob.glob(base + "/*/*"))
    if not files:
        print("  pairwise %s: no tables -- skipped" % run)
        return
    print("  scanning %d pairwise tables for %s" % (len(files), run))
    cc = collections.defaultdict(list)          # t -> [(conf, c1, c2)]
    nbr = collections.defaultdict(dict)         # t -> acc -> best conf to binding
    every = set().union(*allc.values())
    for n, path in enumerate(files):
        if n % 5000 == 0:
            print("    %d/%d" % (n, len(files)), flush=True)
        try:
            fhh = open(path, errors="replace")
        except OSError:
            continue
        with fhh:
            fhh.readline()
            for line in fhh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                sides = parse_side(f[2].split()) + parse_side(f[3].split())
                if not any(a in every for a, _ in sides):
                    continue
                for t in TARGETS:
                    hits = [(a, s) for a, s in sides if a in allc[t]]
                    if not hits:
                        continue
                    if len(hits) > 1:
                        m = min(s for _a, s in hits)
                        for i in range(len(hits)):
                            for j in range(i + 1, len(hits)):
                                cc[t].append((m, hits[i][0], hits[j][0]))
                    if any(a in binding[t] for a, _ in hits):
                        for a, s in sides:
                            if a in allc[t]:
                                continue
                            if s > nbr[t].get(a, -1.0):
                                nbr[t][a] = s
    for t in TARGETS:
        if not cc[t]:
            print("  pairwise %s %s: no control-to-control relation" % (run, t))
            continue
        ws = sorted({w for w, _a, _b in cc[t]}, reverse=True)
        for w in ws:
            cover = set()
            for ww, a, b in cc[t]:
                if ww >= w:
                    cover.add(a)
                    cover.add(b)
            members = {a for a, s in nbr[t].items() if s >= w}
            emit(fh, ch, t, "sonic", "pairwise-conf>=%.4f@%s" % (w, run),
                 "species_to_species_orthologs", members, cover, allc, spidx,
                 neg[t])
        print("  pairwise %s %s: %d thresholds" % (run, t, len(ws)))


def main():
    rl_invariants.check_inputs()
    allc, binding = load_controls()
    spidx = A.species_index()
    # The declared negatives: the ranked decoy pools unioned with the corpus-wide
    # pools, orthogroup-flagged rows excluded because a decoy must be a
    # non-member, plus the adjudicated negative controls.
    dec = A.declared_decoys()
    negc = A.negative_controls()
    neg = {t: (dec[t] | negc[t]) for t in TARGETS}
    print("declared negatives: %s" % {t: len(neg[t]) for t in TARGETS})
    print("controls: %s" % {t: len(allc[t]) for t in TARGETS})
    with open(OUT_SET, "w") as fh, open(OUT_COV, "w") as ch:
        fh.write("## Every OrthoFinder and SonicParanoid setting on disk. COST is "
                 "the DECLARED decoys and negative controls the setting admits -- "
                 "not every call, because an undeclared true ortholog must not "
                 "count as cost. N_CALLED is carried separately as the "
                 "plausibility check against the copy-number ceiling. Controls "
                 "recovered is LEAVE-ONE-OUT (R28).\n")
        fh.write("TARGET\tINSTRUMENT\tSETTING\tSOURCE\tN_MEMBERS\t"
                 "CONTROLS_COVERED\tN_CONTROLS\tCOST\tSPECIES\tPER_SPECIES\t"
                 "N_CALLED\tCALLED_PER_SPECIES\n")
        ch.write("## Which controls each setting recovers, long form.\n")
        ch.write("TARGET\tINSTRUMENT\tSETTING\tPROTEIN_ID\n")
        for spec in GROUPS_OF + GROUPS_SP:
            do_groups(fh, ch, spec, allc, binding, spidx, neg)
        do_graph(fh, ch, allc, binding, spidx, neg)
        for run in ("full03", "infl20", "infl30"):
            do_pairwise(fh, ch, run, allc, binding, spidx, neg)
    print("wrote %s\n      %s" % (OUT_SET, OUT_COV))


if __name__ == "__main__":
    main()
