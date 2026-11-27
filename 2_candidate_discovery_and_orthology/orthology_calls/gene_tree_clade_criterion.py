#!/usr/bin/env python3
# >>> deposit header
# Purpose: Apply the clade criterion of van Dam et al. (2011), anchored on curated controls, to one target's gene tree, and make the per-tip contaminant call on that tree.
# Methods: Orthology calls
# Inputs:  the target's tree (RAxML_bipartitions.<target>, <target>.treefile, or the tree named in PHYLO_TREE), <target>_labels.tsv and <target>_cropped_out.tsv under 11_phylo_vandam<PHYLO_SET>/<target>/;
#          declared negatives and decoys from select_arm_settings
# Outputs: <target>_subtree_calls<suffix>.tsv, <target>_subtree_summary<suffix>.tsv, and <target>_contaminants<suffix>.tsv when the tree carries boundary markers
# Tools:   Python 3 with DendroPy (versions: see envs/)
# <<< deposit header
# Ran as: rl_phylo_subtree.py
"""The gene-tree clade criterion of van Dam et al. (2011), anchored on the curated
positive controls, and the per-tip contaminant call made on the same tree.

The paper's criterion, quoted. It is stated once, for Rheb:

    "From the resulting tree a sub-tree was selected that contained Rheb and
     as many genes from the searched species as possible without including other
     known small GTPases."

It has three parts: (1) the sub-tree contains the target, (2) it maximises the
species represented, and (3) it excludes the known related genes.

The formal name for this optimisation. It is not ad hoc. Yang and Smith (2014;
Mol Biol Evol 31:3081-3092, PMID 25158799) call it maximum inclusion: extract the
subtree containing the largest number of taxa subject to a purity constraint,
with no bootstrap or branch-length threshold anywhere in it. Their purity
constraint is "no repeated species"; here it is "no boundary marker". What this
script computes is therefore

    maximum inclusion (Yang and Smith 2014) anchored on curated positive controls

and only the anchoring is new.

What the script computes, for one target and one tree:
  CONTROL_CLADE      the smallest side of any edge that contains every binding
                     control, and how many withheld controls it holds;
  MAXIMUM_INCLUSION  the side with the most species that holds every binding
                     control and no boundary marker, with a check of whether the
                     boundary excluded anything beyond itself. It has two
                     declared failure states, both REPORTED and never worked
                     around: UNCONSTRAINED, where no declared negative and no
                     declared decoy is present in this tree, so the purity
                     clause excludes nothing and the satisfying sub-tree is the
                     whole tree; and NO_SOLUTION, where every split holding all
                     controls also holds a boundary marker, which means the
                     family and the known related genes are NOT separated by
                     this tree, and that is the result. Neither state changes a
                     per-tip call, because the call does not come from this
                     object;
  THE CALL, per tip  this is the object the ortholog table consumes, and it is
                     neither clade object above. For every tip other than the
                     binding controls and the boundary markers: the nearest
                     labelled neighbour (TREE_CALL)
                     and the nearer, by patristic distance, of the closest control
                     and the closest boundary marker (DIST_CALL); this is the
                     contaminant call, written to <T>_contaminants<suffix>.tsv;
  VANDAM_SUBTREE     van Dam's own rule (at least one control, most species, no
                     boundary marker), kept as a diagnostic of whether the family
                     can be separated as a clade at all;
  a two-by-two table of candidate position (inside or outside the control clade)
                     against the contaminant screen's decision, with the number of
                     outside species that also contribute a binding control.

Three tip classes, which must never be merged:

  CONTROL           binding controls of the named control set (R59_TAG). The
                    anchor: the control clade must contain all of them.

  CONTROL_HELDOUT   controls deliberately withheld from the optimisation. A
                    test, never a boundary. Finding one inside the control clade
                    is a pass: leave-one-out recovery performed on topology
                    instead of on scores, so a control never recovers itself.
                    Recovery is reported; it never constrains the search.

  BOUNDARY          declared negatives and declared decoys present in the tree.
                    These, and only these, play the role of van Dam's "other
                    known small GTPases". Measured disjoint from CONTROL and
                    CONTROL_HELDOUT for all eight targets; the script warns if a
                    withheld control is also a negative or decoy.

Why the classes are kept apart. An earlier measurement pooled "a declared decoy,
a curated negative, or a curation-contested control" into one set of boundary
markers and used it as the purity constraint. That was a design error, because
the third member is a curated true positive: for AMPKg, 21 of those 29 pooled
markers were withheld controls, so the criterion would have had to build a clade
excluding 21 proteins curated as genuine AMPKg, which inverts the test.

Measured on the tip sets of search hits only, before the boundary set was added:

    target  tips  CONTROL  HELD-OUT  neg  decoy  BOUNDARY
    TOR      336    29        1       0     0      0
    RAPTOR   307    23        0       0     0      0
    RICTOR   198    18        1       0     0      0
    LST8     213    21        4       0     4      4
    AMPKa    591    31        8       6    11     15
    AMPKb    293    31        2       0     0      0
    AMPKg    651    33       21       5    46     47
    SIN1     137    14        2       0     0      0

On those tip sets van Dam's purity clause binds for three targets (LST8, AMPKa
and AMPKg) and is unconstrained for the other five. That is a result about what
the screen had already removed, and it is reported per target as such. For the
five, the tree still answers two real questions: do the controls form a coherent
clade at all, and does that clade recover the withheld controls. With no boundary
marker in the tree no contaminant can be called, and an empty result means not
tested, not clean.

No support threshold anywhere. The criterion maximises species under a purity
constraint; it contains no bootstrap term, as neither van Dam's nor Yang and
Smith's does. Support is reported for the clade, never used to define it. The
70% figure of Hillis and Bull (1993) assumes equal rates of change; these tip
sets put Giardia, trypanosomatids and Tetrahymena beside slowly evolving
opisthokonts, the unequal-rate situation in which that paper finds bootstrap
proportions to be overestimates. Importing the threshold would import simulation
conditions these data violate.

Unrooted. The paper designates no outgroup and RAxML's output is unrooted, so
every edge is tested and both of its sides are candidates. The answer does not
depend on where the file happens to be rooted for printing.

Nothing is subtracted. This script writes an assessment. Whether a tree may
remove a protein that the search methods nominated is left to the authors.

Usage: rl_phylo_subtree.py TARGET
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import csv
import os
import sys

ROOT = os.environ.get("PROJECT_ROOT") or (
    (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1"))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "2_candidate_discovery_and_orthology", "search_settings"))
sys.path.insert(0, os.path.join(_REPO_ROOT, "lib"))
# PHYLO_SET names the tip set (for example _boundary or _genecollapsed) and is
# mandatory. An empty value would use the bare 11_phylo_vandam folder, which holds the
# earlier design built from search hits only, whose tip sets are not comparable with
# the boundary or gene-collapsed ones; a bare call would also silently recreate that
# folder if it were deleted. The script therefore refuses to run without it.
_ps = os.environ.get("PHYLO_SET")
if not _ps:
    sys.exit("REFUSED: PHYLO_SET must be set, for example _boundary or _genecollapsed. "
             "Defaulting to the bare 11_phylo_vandam prefix writes into the SUPERSEDED "
             "hits-only track, which is a different DESIGN and not a filter variant.")
OUT = ROOT + "/11_phylo_vandam" + _ps


def rows(path):
    return csv.DictReader((l for l in open(path, errors="replace")
                           if not l.startswith("##")), delimiter="\t")


def support_of(nd):
    """RAxML writes one support value per node. IQ-TREE run with -B and --alrt
    writes "alrt/ufboot", so the last field, the ultrafast bootstrap, is taken.

    The two values are not on the same scale and must never be compared as though
    they were: RAxML's is a standard (Felsenstein) bootstrap proportion, and
    IQ-TREE's ultrafast bootstrap is a different estimator that runs
    systematically higher. The cross-check between the two programs is about
    clade membership, not support, and the report says so."""
    lab = (nd.label or "").strip() if nd.label else ""
    if "/" in lab:
        lab = lab.rsplit("/", 1)[-1].strip()
    try:
        return float(lab)
    except (TypeError, ValueError):
        return None


def fmt(x):
    return ("%.0f" % x) if x is not None else ""


def main():
    import dendropy
    import select_arm_settings as A

    target = sys.argv[1]
    d = os.path.join(OUT, target)

    # The control set must be named (R59_TAG); there is no untagged default
    # anywhere in this pipeline. select_arm_settings reads R59_TAG itself, and a
    # wrong tag here would silently anchor the clade on the wrong controls.
    tag = (os.environ.get("R59_TAG") or "").strip()
    if not tag:
        raise SystemExit("REFUSED: R59_TAG is unset. The clade is anchored on a "
                         "NAMED control set and there is no default.")

    # PHYLO_TREE names an alternative tree (the IQ-TREE cross-check), and
    # PHYLO_SUFFIX keeps its outputs beside the primary rather than overwriting
    # them. Both trees are then judged by the same code, which is the only way the
    # comparison means anything.
    sfx = (os.environ.get("PHYLO_SUFFIX") or "").strip()
    forced = (os.environ.get("PHYLO_TREE") or "").strip()
    if forced:
        tf = forced if os.path.isabs(forced) else os.path.join(d, forced)
        if not os.path.exists(tf):
            raise SystemExit("REFUSED: PHYLO_TREE=%s does not exist" % tf)
    else:
        tf = None
        for cand in ("RAxML_bipartitions.%s" % target, "%s.treefile" % target):
            q = os.path.join(d, cand)
            if os.path.exists(q):
                tf = q
                break
        if tf is None:
            raise SystemExit("no tree found for %s in %s" % (target, d))

    role, spof, screen = {}, {}, {}
    for r in rows(os.path.join(d, "%s_labels.tsv" % target)):
        a = (r.get("PROTEIN_ID") or "").strip()
        role[a] = (r.get("ROLE") or "").strip()
        spof[a] = (r.get("SPECIES") or "").strip()
        screen[a] = (r.get("SCREEN") or "").strip()

    # Tips that the crop dropped because they carry too little of the cropped
    # region. They are not a call and not judged absent; they simply cannot be
    # placed. They are reported so that every candidate in the master table is
    # accounted for.
    dropped = {}
    cop = os.path.join(d, "%s_cropped_out.tsv" % target)
    if os.path.exists(cop):
        for r in rows(cop):
            a = (r.get("PROTEIN_ID") or "").strip()
            if a:
                dropped[a] = (r.get("RESIDUES_IN_CROP") or "").strip()

    t = dendropy.Tree.get(path=tf, schema="newick", preserve_underscores=True)
    tips = {}
    for tip in t.leaf_node_iter():
        lab = (tip.taxon.label if tip.taxon else "") or ""
        tips[tip] = lab.split("|")[-1].strip()
    accs = set(tips.values())

    unlabelled = accs - set(role)
    if unlabelled:
        raise SystemExit("REFUSED: %d tip(s) in the tree carry no label row, e.g. "
                         "%s. The tree and the label table disagree about the tip "
                         "set." % (len(unlabelled), sorted(unlabelled)[:3]))

    ctrl = {a for a in accs if role.get(a) == "CONTROL"}
    held = {a for a in accs if role.get(a) == "CONTROL_HELDOUT"}
    cand = {a for a in accs if role.get(a) == "CANDIDATE"}

    # The purity constraint: declared negatives and decoys only. Withheld controls
    # are explicitly subtracted: they are curated positives and are the test, not
    # the boundary.
    neg = A.negative_controls().get(target, set()) & accs
    dec = A.declared_decoys().get(target, set()) & accs
    bnd = (neg | dec) - ctrl - held
    overlap = held & (neg | dec)

    print("%s [%s]: %d tips = %d CONTROL, %d CONTROL_HELDOUT, %d CANDIDATE"
          % (target, tag, len(tips), len(ctrl), len(held), len(cand)))
    print("  purity constraint: %d boundary markers (%d declared negatives, "
          "%d declared decoys), held-out controls EXCLUDED from it"
          % (len(bnd), len(neg - held), len(dec - held)))
    if overlap:
        print("  WARNING: %d held-out control(s) are also declared negative or "
              "decoy -- a curation conflict, reported not resolved: %s"
              % (len(overlap), ", ".join(sorted(overlap))))
    if dropped:
        print("  %d candidate(s) absent from the tree because the crop left them "
              "under the residue floor; absence says nothing about them"
              % len(dropped))
    if not ctrl:
        raise SystemExit("REFUSED: no declared control in the tree; the criterion "
                         "is anchored on the controls and cannot be applied.")

    below = {}
    for nd in t.postorder_node_iter():
        if nd.is_leaf():
            below[nd] = {tips[nd]}
        else:
            s = set()
            for ch in nd.child_nodes():
                s |= below[ch]
            below[nd] = s

    def splits():
        """Every edge of the unrooted tree, both sides, with that edge's support."""
        for nd in t.postorder_node_iter():
            if nd.parent_node is None:
                continue
            sup = support_of(nd)
            ins = below[nd]
            yield ins, sup
            yield accs - ins, sup

    # ---- CONTROL_CLADE: smallest split side containing every binding control --
    best_cc = None
    for side, sup in splits():
        if ctrl <= side and (best_cc is None or len(side) < len(best_cc[0])):
            best_cc = (side, sup)
    cc_side, cc_sup = best_cc if best_cc else (accs, None)
    cc_whole = best_cc is None or cc_side == accs
    cc_sp = {spof.get(a) for a in cc_side} - {None, ""}
    print("  CONTROL_CLADE: smallest split holding all %d controls has %d tips "
          "(%d candidates, %d species), support %s"
          % (len(ctrl), len(cc_side), len(cc_side & cand), len(cc_sp),
             fmt(cc_sup) if not cc_whole else "n/a (whole tree)"))

    # ---- held-out recovery: the test, reported, never a constraint ------------
    if held:
        rec = held & cc_side
        print("  HELD-OUT RECOVERY: %d/%d held-out controls fall INSIDE the "
              "control clade%s" % (len(rec), len(held),
                                   "" if len(rec) == len(held)
                                   else " -- missed: " + ", ".join(sorted(held - rec))))
    else:
        print("  HELD-OUT RECOVERY: no held-out control for this target")

    # ---- maximum inclusion under the purity constraint ------------------------
    if not bnd:
        mi_status = "UNCONSTRAINED"
        mi_side, mi_sup = accs, None
        print("  MAXIMUM_INCLUSION: UNCONSTRAINED -- no declared negative and no "
              "declared decoy is present in this tree, so the clause 'without "
              "including other known related genes' excludes nothing and the "
              "sub-tree satisfying it is the whole tree. Reported, not worked "
              "around.")
    else:
        best_mi = None
        for side, sup in splits():
            if not (ctrl <= side) or (side & bnd):
                continue
            nsp = len({spof.get(a) for a in side} - {None, ""})
            key = (nsp, len(side))
            if best_mi is None or key > best_mi[0]:
                best_mi = (key, side, sup)
        if best_mi is None:
            # Every split holding all controls also holds a boundary marker.
            # Report how far the anchor can be kept, which localises the problem.
            partial = None
            for side, sup in splits():
                if side & bnd:
                    continue
                k = (len(side & ctrl), len({spof.get(a) for a in side} - {None, ""}))
                if partial is None or k > partial[0]:
                    partial = (k, side, sup)
            mi_status = "NO_SOLUTION"
            mi_side, mi_sup = set(), None
            print("  MAXIMUM_INCLUSION: NO_SOLUTION -- every split containing all "
                  "%d controls also contains a boundary marker. The family and "
                  "the known related genes are NOT separated by this tree, which "
                  "is the result." % len(ctrl))
            if partial:
                print("    largest boundary-free split holds %d/%d controls, "
                      "%d species" % (partial[0][0], len(ctrl), partial[0][1]))
        else:
            mi_status = "SELECTED"
            _k, mi_side, mi_sup = best_mi
            print("  MAXIMUM_INCLUSION: %d tips, %d/%d controls, %d candidates, "
                  "%d species, 0 boundary markers, support %s"
                  % (len(mi_side), len(mi_side & ctrl), len(ctrl),
                     len(mi_side & cand), len({spof.get(a) for a in mi_side} - {None, ""}),
                     fmt(mi_sup)))
            # Did the boundary do any work?
            #
            # A boundary made of remote, low-scoring proteins sits far from the
            # family, so maximum inclusion excludes it and nothing else: the
            # selected clade is exactly "every tip except the boundary", the
            # constraint is satisfied trivially, and MI_STATUS=SELECTED looks like
            # a result while measuring nothing. That is the expected outcome for
            # RICTOR (boundary 5.7 to 26.8 bits) and SIN1 (3.3 to 11.1 bits), whose
            # decoys the crop did not drop, so this check catches it instead.
            #
            # A boundary that is a confusable family interleaves with the target:
            # excluding it also costs some candidates, and that difference is what
            # the tree shows.
            excluded_beyond = (accs - mi_side) - bnd
            mi_trivial = not excluded_beyond
            if mi_trivial:
                print("    *** THE CONSTRAINT DID NO WORK. The selected clade is "
                      "exactly every tip except the %d boundary markers, so "
                      "excluding them cost nothing and bounded nothing. The "
                      "boundary is remote from the family rather than confusable "
                      "with it; treat MI_STATUS=SELECTED here as UNCONSTRAINED in "
                      "substance." % len(bnd))
            else:
                print("    the constraint BINDS: excluding the boundary also "
                      "excludes %d non-boundary tip(s), %d of them candidates -- "
                      "those are what the boundary is telling us about"
                      % (len(excluded_beyond), len(excluded_beyond & cand)))

    # ======================================================================
    # The contaminant call: per tip, by nearest labelled neighbour.
    # ======================================================================
    #
    # Why not clade selection. van Dam selects a sub-tree because in his tree the
    # family is a clade: his tip set held several small-GTPase families and Rheb
    # was one block inside it. Here the tip set is a single family, already
    # separated by the search methods, so the tree is essentially that family's
    # species tree and the controls do not form a clade. Measured: SIN1's controls
    # span the whole topology, and the smallest split holding all 14 is 112 of 137
    # tips.
    #
    # Clade selection on such a tree fails in a way that looks like an answer.
    # Maximising species subject to "free of boundary markers and contains a
    # control" picked, for AMPKg, a 70-tip pocket holding 1 of 33 controls and
    # called 547 of 616 hits contaminants. That is not a contaminant call but an
    # artefact of demanding a clade from a topology that has none.
    #
    # What van Dam's Figure S3 actually encodes is a per-sequence assignment: each
    # sequence belongs to whichever labelled group it sits with. That works whether
    # or not the groups are clades, so it is what is computed here.
    #
    # Rule: from each tip that is not itself a binding control or boundary marker
    # (candidates and withheld controls), walk up until the enclosing clade first
    # contains a labelled tip, a declared CONTROL or a declared BOUNDARY marker. If
    # it holds only controls, the tip sits with the family. If only boundary
    # markers, it sits with the known related genes and is called a contaminant. If
    # both, the call is AMBIGUOUS, and the ambiguity is reported rather than broken
    # by a tie-break.
    call = {}
    if bnd:
        kids = {}
        for nd in t.postorder_node_iter():
            kids[nd] = nd.child_nodes()
        for tip, acc in tips.items():
            if acc in ctrl or acc in bnd:
                continue
            nd = tip
            verdict = "UNRESOLVED"
            while nd.parent_node is not None:
                nd = nd.parent_node
                lab_c = below[nd] & ctrl
                lab_b = below[nd] & bnd
                if lab_c or lab_b:
                    verdict = ("AMBIGUOUS" if (lab_c and lab_b)
                               else "WITH_CONTROLS" if lab_c else "CONTAMINANT")
                    break
            call[acc] = verdict
        # ------------------------------------------------------------------
        # A second, independent signal: patristic distance.
        #
        # The topological call above asks which group a candidate sits with. It
        # cannot rank, and it cannot resolve its own AMBIGUOUS cases. The standard
        # quantitative diagnostic is the patristic distance (branch length summed
        # along the path) to the nearest CONTROL against the nearest BOUNDARY
        # marker; whichever is closer is the group it belongs to.
        #
        # No threshold is introduced: the comparison is d_boundary < d_control.
        # The MARGIN is reported so that a near-tie is visible as a near-tie
        # rather than hidden inside a call.
        #
        # The known failure mode: a fast-evolving sequence sits far from
        # everything, so distance alone can call a genuine but divergent ortholog
        # a contaminant. That is why the terminal branch length is reported too,
        # and why the distance call is a second measure beside the topological
        # call rather than a replacement for it. Where the two disagree, the
        # disagreement is the result.
        pdm = t.phylogenetic_distance_matrix()
        tip_of = {}
        for tp in t.leaf_node_iter():
            lab = (tp.taxon.label if tp.taxon else "") or ""
            tip_of[lab.split("|")[-1].strip()] = tp
        # DendroPy's patristic_distance() takes taxon objects, not nodes.
        ctrl_tax = [tip_of[a].taxon for a in ctrl if a in tip_of]
        bnd_tax = [tip_of[a].taxon for a in bnd if a in tip_of]
        dist = {}
        for a in call:
            tp = tip_of.get(a)
            if tp is None or tp.taxon is None or not ctrl_tax or not bnd_tax:
                continue
            dc = min(pdm.patristic_distance(tp.taxon, x) for x in ctrl_tax)
            db = min(pdm.patristic_distance(tp.taxon, x) for x in bnd_tax)
            bl = tp.edge_length if tp.edge_length is not None else float("nan")
            dist[a] = (dc, db, bl)
        n_dist_ct = sum(1 for v in dist.values() if v[1] < v[0])
        agree_both = sum(1 for a, v in dist.items()
                         if (v[1] < v[0]) == (call.get(a) == "CONTAMINANT"))
        if dist:
            print("  PATRISTIC INSTRUMENT: %d of %d candidates are closer to a boundary "
                  "marker than to any control; the two instruments agree on %d (%.0f%%)"
                  % (n_dist_ct, len(dist), agree_both,
                     100.0 * agree_both / len(dist)))
            amb = [a for a in call if call[a] == "AMBIGUOUS" and a in dist]
            if amb:
                ac = sum(1 for a in amb if dist[a][1] < dist[a][0])
                print("      of the %d AMBIGUOUS by topology, distance puts %d "
                      "nearer the boundary and %d nearer the controls"
                      % (len(amb), ac, len(amb) - ac))

        n_ct = sum(1 for v in call.values() if v == "CONTAMINANT")
        n_cl = sum(1 for v in call.values() if v == "WITH_CONTROLS")
        n_am = sum(1 for v in call.values() if v == "AMBIGUOUS")
        print("  *** CONTAMINANT CALL (nearest labelled neighbour): "
              "%d contaminant, %d with the controls, %d ambiguous, of %d candidates"
              % (n_ct, n_cl, n_am, len(call)))
        agree = sum(1 for a, v in call.items()
                    if v == "CONTAMINANT" and screen.get(a) == "REMOVED")
        print("      of the %d contaminants, the screen already removed %d; "
              "%d are screen-RETAINED and are the disagreement"
              % (n_ct, agree, n_ct - agree))
        with open(os.path.join(d, "%s_contaminants%s.tsv" % (target, sfx)), "w") as fh:
            fh.write("## CONTAMINANT CALL, per tip, by nearest labelled neighbour in "
                     "the tree. A candidate is CONTAMINANT when the smallest clade "
                     "enclosing it that contains any DECLARED tip contains only "
                     "BOUNDARY markers and no CONTROL.\n")
            fh.write("## Clade selection (van Dam's own sub-tree rule) is NOT used: "
                     "our controls do not form a clade, because the tip set is a "
                     "single family already separated upstream, so the tree is that "
                     "family's species tree. Recorded in PHYLO_RUN_LOG.md.\n")
            fh.write("## AMBIGUOUS means the enclosing clade holds BOTH a control and "
                     "a boundary marker; it is reported, never tie-broken.\n")
            fh.write("## Nothing is removed from the master table on this basis.\n")
            fh.write("## TREE_CALL is topological (nearest labelled neighbour). "
                     "DIST_CALL is the independent patristic instrument: whichever of the "
                     "nearest control / nearest boundary marker is closer. "
                     "INSTRUMENTS_AGREE flags where they differ -- that disagreement is "
                     "the result and is not resolved here.\n")
            fh.write("## TERMINAL_BRANCH is reported because a fast-evolving "
                     "sequence sits far from EVERYTHING, which is the known way "
                     "for the distance instrument to mistake a divergent ortholog for a "
                     "contaminant.\n")
            fh.write("TARGET\tPROTEIN_ID\tSPECIES\tROLE\tSCREEN\tTREE_CALL\t"
                     "DIST_TO_CONTROL\tDIST_TO_BOUNDARY\tDIST_CALL\tMARGIN\t"
                     "TERMINAL_BRANCH\tINSTRUMENTS_AGREE\tSCREEN_AGREES\n")
            for a in sorted(call):
                v = call[a]
                sc = screen.get(a, "")
                ag = "yes" if ((v == "CONTAMINANT") == (sc == "REMOVED")) else "no"
                dc, db, bl = dist.get(a, (float("nan"),) * 3)
                dcall = ("" if a not in dist
                         else "CONTAMINANT" if db < dc else "WITH_CONTROLS")
                marg = "" if a not in dist else "%.4f" % abs(dc - db)
                instruments = "" if a not in dist else (
                    "yes" if (dcall == "CONTAMINANT") == (v == "CONTAMINANT")
                    else "no")
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                         % (target, a, spof.get(a, ""), role.get(a, ""), sc, v,
                            "" if a not in dist else "%.4f" % dc,
                            "" if a not in dist else "%.4f" % db,
                            dcall, marg,
                            "" if a not in dist else "%.4f" % bl, instruments, ag))
    else:
        print("  *** CONTAMINANT CALL: not possible -- this tree contains no "
              "declared boundary marker, so nothing can be called a contaminant. "
              "An empty result here means NOT TESTED, not clean.")

    # ======================================================================
    # van Dam's own sub-tree rule (VANDAM_SUBTREE), computed as a diagnostic.
    # ======================================================================
    #
    # van Dam et al. (2011): "From the resulting tree a sub-tree was selected that
    # contained Rheb and as many genes from the searched species as possible
    # without including other known small GTPases." In the paper, sequences
    # outside that sub-tree are the non-orthologs. Here the rule is computed and
    # printed as a diagnostic of whether the family can be separated as a clade
    # at all, and its counts go into the summary table; the contaminant call is
    # the per-tip call above (see the note on the removed clade-based writer
    # below).
    #
    # The anchor is at least one control, not the whole control set. van Dam
    # anchors on Rheb and maximises species. MAXIMUM_INCLUSION above requires
    # every binding control in one clade free of boundary markers, which is
    # strictly stronger than his rule, and on SIN1 it returned NO_SOLUTION for
    # that reason alone, an artefact of the stronger rule rather than a statement
    # about SIN1. Both are therefore reported: the strict one as MI, and this one,
    # faithful to the paper.
    #
    # Ties are broken toward more species first, then more tips, as the sentence
    # reads. No support threshold, as in the paper.
    vd_side, vd_sup, vd_status = set(), None, "UNCONSTRAINED"
    if bnd:
        best_vd = None
        for side, sup in splits():
            if (side & bnd) or not (side & ctrl):
                continue
            nsp = len({spof.get(a) for a in side} - {None, ""})
            key = (nsp, len(side))
            if best_vd is None or key > best_vd[0]:
                best_vd = (key, side, sup)
        if best_vd is None:
            vd_status = "NO_SOLUTION"
            print("  VANDAM_SUBTREE: NO_SOLUTION -- no split containing any control "
                  "is free of boundary markers. The tree does not separate the "
                  "family from the known related genes at all.")
        else:
            _k, vd_side, vd_sup = best_vd
            vd_status = "SELECTED"
    else:
        vd_side = accs

    # Hits (candidates and withheld controls) that fall outside the selected
    # sub-tree; these counts go only into the printed summary and the summary table.
    contam = sorted((cand | held) - vd_side) if vd_status != "NO_SOLUTION" else []
    ctrl_out = sorted(ctrl - vd_side) if vd_status != "NO_SOLUTION" else []

    if vd_status == "SELECTED":
        print("  VANDAM_SUBTREE: %d tips, %d species, %d/%d controls, 0 boundary "
              "markers, support %s"
              % (len(vd_side), len({spof.get(a) for a in vd_side} - {None, ""}),
                 len(vd_side & ctrl), len(ctrl), fmt(vd_sup)))
        print("  *** CONTAMINANTS: %d of %d called hits fall OUTSIDE the sub-tree"
              % (len(contam), len(cand | held)))
        if ctrl_out:
            print("      WARNING: %d declared CONTROL(s) are also outside: %s. A "
                  "curated positive outside the selected sub-tree means the "
                  "criterion and the curation disagree, and that is reported, not "
                  "resolved." % (len(ctrl_out), ", ".join(ctrl_out[:5])))
    elif vd_status == "UNCONSTRAINED":
        print("  VANDAM_SUBTREE: UNCONSTRAINED -- no boundary marker in this tree, "
              "so nothing is excluded and NO contaminant can be called.")

    # The clade-based writer of contaminant calls was removed, not merely
    # superseded. It wrote to the same file name as the per-tip call above and,
    # running after it, silently overwrote that result with its own. It is gone
    # rather than renamed, because a second file of contaminant calls invites the
    # wrong one being quoted. VANDAM_SUBTREE is still computed and printed above
    # as a diagnostic of whether the family can be separated as a clade at all; it
    # is not a call.

    # ---- Why tips fall outside the control clade, measured before the two-by-two
    # table against the screen is read.
    #
    # The smallest split containing every binding control is fixed by which taxa
    # are controls. A lineage that contributes no binding control can be excluded
    # by that split for taxon-sampling reasons alone, which says nothing about the
    # contaminant screen.
    #
    # This is not hypothetical. Measured on SIN1: all 25 tips outside the split
    # are ciliates (Stentor x11, Paramecium x7, Ichthyophthirius, Oxytricha,
    # Pseudocohnilembus, Stylonychia, both Tetrahymena), and none of the 14
    # binding controls is a ciliate. The withheld Tetrahymena control falls
    # outside with its own clade, which was once scored as a recovery miss. On
    # LST8 the outside set is parabasalids, diplomonads, red algae and
    # cryptophytes, including both withheld Giardia controls.
    #
    # An earlier version labelled the outside column "screen may be too
    # permissive" unconditionally. That is a false inference whenever the outside
    # set is simply an unsampled clade, so the interpretation now depends on the
    # measurement below.
    ctrl_sp = {spof.get(a) for a in ctrl} - {None, ""}
    out_tips = accs - cc_side
    out_sp = {spof.get(a) for a in out_tips} - {None, ""}
    shared_sp = out_sp & ctrl_sp

    ret_in = len([a for a in cand if screen.get(a) == "RETAINED" and a in cc_side])
    ret_out = len([a for a in cand if screen.get(a) == "RETAINED" and a not in cc_side])
    rem_in = len([a for a in cand if screen.get(a) == "REMOVED" and a in cc_side])
    rem_out = len([a for a in cand if screen.get(a) == "REMOVED" and a not in cc_side])

    print("  OUTSIDE THE CLADE: %d tips, %d species, of which %d also contribute "
          "a binding control" % (len(out_tips), len(out_sp), len(shared_sp)))
    print("  2x2 vs the screen (candidates only, against CONTROL_CLADE):")
    print("      screen RETAINED  inside %4d | outside %4d" % (ret_in, ret_out))
    print("      screen REMOVED   inside %4d | outside %4d" % (rem_in, rem_out))
    if out_tips and not shared_sp:
        print("      *** THE OUTSIDE COLUMN IS NOT EVIDENCE ABOUT THE SCREEN "
              "HERE. No species outside the clade contributes a binding control, "
              "so the exclusion is fully consistent with control taxon sampling "
              "alone -- an unsampled lineage, not a contaminant call.")
    elif shared_sp:
        print("      %d species fall on BOTH sides, so for those the split is not "
              "explained by taxon sampling and the off-diagonals are "
              "interpretable: outside+RETAINED = the screen kept what the "
              "topology does not place with the controls; inside+REMOVED = the "
              "reverse." % len(shared_sp))

    def cls(a):
        if a in ctrl:
            return "CONTROL"
        if a in held:
            return "CONTROL_HELDOUT"
        if a in bnd:
            return "BOUNDARY"
        return "CANDIDATE"

    with open(os.path.join(d, "%s_subtree_calls%s.tsv" % (target, sfx)), "w") as fh:
        fh.write("## Topology assessment, NOT a filter -- nothing here has been "
                 "subtracted from the master table.\n"
                 "## CLASS: CONTROL anchors the clade; CONTROL_HELDOUT is a TEST "
                 "(inside = recovered) and never constrains the search; BOUNDARY "
                 "is a declared negative or decoy and is the only purity "
                 "constraint.\n"
                 "## IN_CONTROL_CLADE = inside the smallest split containing every "
                 "binding control. IN_MAX_INCLUSION follows van Dam's clause via "
                 "Yang & Smith 2014 maximum inclusion, and is blank when that "
                 "clause is UNCONSTRAINED for this target.\n"
                 "## No bootstrap threshold is applied anywhere. Support is "
                 "reported FOR the clade, never used TO DEFINE it.\n")
        fh.write("TARGET\tPROTEIN_ID\tSPECIES\tROLE\tCLASS\tSCREEN\t"
                 "IN_CONTROL_CLADE\tIN_MAX_INCLUSION\tMI_STATUS\t"
                 "CONTROL_CLADE_SUPPORT\n")
        for a in sorted(accs):
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                     % (target, a, spof.get(a, ""), role.get(a, "?"), cls(a),
                        screen.get(a, ""),
                        "yes" if a in cc_side else "no",
                        ("yes" if a in mi_side else "no") if mi_status == "SELECTED" else "",
                        mi_status, fmt(cc_sup) if not cc_whole else ""))
        for a in sorted(dropped):
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\tNOT_IN_TREE\t\t%s\t\n"
                     % (target, a, spof.get(a, ""), role.get(a, "?"), cls(a),
                        screen.get(a, ""), mi_status))

    with open(os.path.join(d, "%s_subtree_summary%s.tsv" % (target, sfx)), "w") as fh:
        fh.write("TARGET\tCONTROL_SET\tTREE\tTIPS\tCONTROLS\tHELDOUT\tCANDIDATES\t"
                 "BOUNDARY_MARKERS\tNEGATIVES\tDECOYS\tNOT_IN_TREE\t"
                 "CC_TIPS\tCC_CANDIDATES\tCC_SPECIES\tCC_SUPPORT\t"
                 "HELDOUT_RECOVERED\tCAND_OUTSIDE_CC\t"
                 "RETAINED_IN\tRETAINED_OUT\tREMOVED_IN\tREMOVED_OUT\t"
                 "OUT_TIPS\tOUT_SPECIES\tOUT_SPECIES_ALSO_CONTROL\t"
                 "MI_STATUS\tMI_TIPS\tMI_SPECIES\tMI_SUPPORT\t"
                 "MI_EXCLUDED_BEYOND_BOUNDARY\tMI_CONSTRAINT_BINDS\t"
                 "VANDAM_STATUS\tVANDAM_TIPS\tVANDAM_SPECIES\t"
                 "CONTAMINANTS\tCONTAM_SCREEN_AGREES\tCONTROLS_OUTSIDE\n")
        fh.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"
                 "%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                 % (target, tag, os.path.basename(tf), len(tips), len(ctrl),
                    len(held), len(cand), len(bnd), len(neg - held), len(dec - held),
                    len(dropped),
                    len(cc_side), len(cc_side & cand), len(cc_sp),
                    fmt(cc_sup) if not cc_whole else "",
                    ("%d/%d" % (len(held & cc_side), len(held))) if held else "",
                    len(cand) - len(cc_side & cand),
                    ret_in, ret_out, rem_in, rem_out,
                    len(out_tips), len(out_sp), len(shared_sp),
                    mi_status,
                    len(mi_side) if mi_status == "SELECTED" else "",
                    len({spof.get(a) for a in mi_side} - {None, ""}) if mi_status == "SELECTED" else "",
                    fmt(mi_sup) if mi_status == "SELECTED" else "",
                    (len((accs - mi_side) - bnd)) if mi_status == "SELECTED" else "",
                    ("no" if not ((accs - mi_side) - bnd) else "yes")
                    if mi_status == "SELECTED" else "",
                    vd_status,
                    len(vd_side) if vd_status == "SELECTED" else "",
                    len({spof.get(a) for a in vd_side} - {None, ""}) if vd_status == "SELECTED" else "",
                    len(contam) if vd_status == "SELECTED" else "",
                    sum(1 for a in contam if screen.get(a) == "REMOVED") if vd_status == "SELECTED" else "",
                    len(ctrl_out) if vd_status == "SELECTED" else ""))

    print("  wrote %s_subtree_calls%s.tsv and %s_subtree_summary%s.tsv"
          % (target, sfx, target, sfx))


if __name__ == "__main__":
    main()
