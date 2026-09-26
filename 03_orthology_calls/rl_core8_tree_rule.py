#!/usr/bin/env python3
# >>> deposit header
# Purpose: Resolves the candidates of the eight core targets that the final ortholog table leaves unresolved,
#          by the gene-tree call where a tree placed them and by two-arm support where it did not.
# Methods: Orthology calls
# Inputs:  04_final_tables/FINAL_ORTHOLOG_TABLE__active10.tsv; <target>_contaminants_iqtree.tsv and
#          <target>_cropped_out.tsv of the gene-collapsed and boundary gene-tree tracks
# Outputs: 04_final_tables/CORE8_TREE_RULE_RESOLUTION.tsv (written only with --apply; never overwritten)
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
"""Tree rule for the eight core targets: resolve the candidates the final table leaves unresolved.

The final ortholog table of the ten active targets (04_final_tables/FINAL_ORTHOLOG_TABLE__active10.tsv,
whose rows for the eight core targets are identical to the earlier eight-target table) marks two kinds
of candidate IS_ORTHOLOG=unresolved:
  PRESENCE_ONLY   the proteome holds several sequences per gene, so the gene-collapsed tree (one tip per
                  gene) excluded the tip (<T>_excluded_multirecord.tsv); the boundary tree
                  (11_phylo_vandam_boundary) did place most of them;
  NOT_PLACEABLE   fewer than 20 residues fall inside the tree's crop, so no tree places them.
Rows marked USED_AS_BOUNDARY_MARKER are the reference set of the tree, not candidates, and are left as
they are.

One rule applies to every such candidate, the same standard the final table applies to every other
species:
  1. Where a tree placed the candidate, the tree call decides (the gene-collapsed track first,
     otherwise the boundary track), whatever the support from the search arms:
     WITH_CONTROLS -> ORTHOLOG (IS_ORTHOLOG=yes; the call was made from the gene tree);
     CONTAMINANT   -> NON_ORTHOLOG (IS_ORTHOLOG=no; the call was made from the gene tree);
  2. Where the tree abstained (AMBIGUOUS) or no tree could place the candidate, two-arm support decides:
     nominated by both arms -> ORTHOLOG (yes); otherwise NON_ORTHOLOG_ARM_RULE (no).
The tree calls follow the clade criterion of van Dam et al. (2011). The rule was chosen on 16 September
2026, including the choice that the tree decides for candidates that only one search arm nominated.

Writes 04_final_tables/CORE8_TREE_RULE_RESOLUTION.tsv and never overwrites it. The final table itself
is not modified; the builder of the publication tables applies this file to the ortholog table (T1)
and derives the table of retained orthologs (T1b) and the species-state table (T2) from the result.
Dry run by default; --apply writes the file.
"""
import collections, hashlib, os, sys


def _project_root_from_file(start):
    d = os.path.dirname(os.path.abspath(start))
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "05_method_and_rules")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("REFUSED: PROJECT_ROOT unset and no project root found above %s" % start)


ROOT = os.environ.get("PROJECT_ROOT") or _project_root_from_file(__file__)
MASTER = os.path.join(ROOT, "04_final_tables", "FINAL_ORTHOLOG_TABLE__active10.tsv")
OUT = os.path.join(ROOT, "04_final_tables", "CORE8_TREE_RULE_RESOLUTION.tsv")
CORE8 = ["TOR", "RAPTOR", "RICTOR", "LST8", "AMPKa", "AMPKb", "AMPKg", "SIN1"]
TRACKS = ["11_phylo_vandam_genecollapsed", "11_phylo_vandam_boundary"]
APPLY = "--apply" in sys.argv
WORDS = ('"My rule: are they in the trees? If yes, consider, if no, why?"; "Apply the tree rule (Recommended)"; '
         '"Yes, tree decides (Recommended)"')


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def rd(p):
    h, out = None, []
    with open(p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if h is None:
                h = f
                continue
            out.append(dict(zip(h, f)))
    return h, out


if os.path.lexists(OUT):
    sys.exit("REFUSED: %s exists; it is never overwritten" % OUT)
_, master = rd(MASTER)
calls, crop, inputs = {}, {}, [("MASTER", os.path.relpath(MASTER, ROOT), sha(MASTER))]
for t in CORE8:
    for tr in TRACKS:
        f = os.path.join(ROOT, tr, t, "%s_contaminants_iqtree.tsv" % t)
        if not os.path.isfile(f):
            sys.exit("REFUSED: missing tree calls %s" % f)
        inputs.append(("TREE_CALLS", os.path.relpath(f, ROOT), sha(f)))
        for r in rd(f)[1]:
            calls[(tr, t, r["PROTEIN_ID"])] = r["TREE_CALL"]
        c = os.path.join(ROOT, tr, t, "%s_cropped_out.tsv" % t)
        if os.path.isfile(c):
            for r in rd(c)[1]:
                crop[(tr, t, r["PROTEIN_ID"])] = r["RESIDUES_IN_CROP"]

rows, tally = [], collections.Counter()
for r in master:
    if r["TARGET"] not in CORE8 or r["IS_ORTHOLOG"] != "unresolved" or r["DISPOSITION"] == "USED_AS_BOUNDARY_MARKER":
        continue
    if r["DISPOSITION"] not in ("PRESENCE_ONLY", "NOT_PLACEABLE"):
        sys.exit("REFUSED: unexpected unresolved disposition %s for %s" % (r["DISPOSITION"], r["PROTEIN_ID"]))
    t, p = r["TARGET"], r["PROTEIN_ID"]
    track, call = "-", None
    for tr in TRACKS:
        if calls.get((tr, t, p)):
            track, call = tr, calls[(tr, t, p)]
            break
    why_no_tree = ""
    if call is None:
        res = [crop[(tr, t, p)] for tr in TRACKS if (tr, t, p) in crop]
        why_no_tree = ("fewer than 20 residues in the crop (%s)" % "/".join(res)) if res else "not a tip in either tree"
    arms = r["ARM_SUPPORT"]
    if call == "WITH_CONTROLS":
        new = ("ORTHOLOG", "yes", "TREE", "placed with the controls in %s" % track)
    elif call == "CONTAMINANT":
        new = ("NON_ORTHOLOG", "no", "TREE", "placed outside the controls (contaminant) in %s" % track)
    elif call in ("AMBIGUOUS", None):
        lead = ("the tree abstained (AMBIGUOUS) in %s" % track) if call else ("no tree placed it: %s" % why_no_tree)
        new = (("ORTHOLOG", "yes", "R66", "%s; both arms support it" % lead) if arms == "BOTH" else
               ("NON_ORTHOLOG_ARM_RULE", "no", "R66", "%s; two-arm support absent (%s)" % (lead, arms)))
    else:
        sys.exit("REFUSED: unrecognised tree call %r for %s" % (call, p))
    tally[(r["DISPOSITION"], call or "NO_TREE", arms, new[0])] += 1
    rows.append([t, r["SPECIES"], p, r["DISPOSITION"], r["IS_ORTHOLOG"], r["DISPOSITION_SOURCE"], track,
                 call or "-", arms, new[0], new[1], new[2], new[3]])

for k, v in sorted(tally.items()):
    print("%4d  %s" % (v, k))
print("rows: %d (yes %d, no %d)" % (len(rows), sum(1 for x in rows if x[10] == "yes"), sum(1 for x in rows if x[10] == "no")))
if not APPLY:
    print("dry run: nothing written")
    sys.exit(0)
with open(OUT, "x", encoding="utf-8") as fh:
    fh.write("## CORE-8 TREE RULE RESOLUTION. Bruno 2026-09-16, LKB1 (RBX2) session: %s.\n" % WORDS)
    fh.write("## One row per core-8 candidate the master leaves unresolved (PRESENCE_ONLY: multi-record proteome excluded from\n"
             "## the gene-collapsed tree; NOT_PLACEABLE: under 20 residues in the crop). The tree call decides where a tree\n"
             "## placed it (gene-collapsed first, else boundary); where the tree abstained or could not place it, two-arm\n"
             "## support decides (R66). USED_AS_BOUNDARY_MARKER rows are not candidates and are not listed. The master is\n"
             "## not modified; the publication builder applies these rows to T1.\n")
    for kind, rel, h in inputs:
        fh.write("## INPUT %s %s sha256=%s\n" % (kind, rel, h))
    fh.write("\t".join(["TARGET", "SPECIES", "PROTEIN_ID", "OLD_DISPOSITION", "OLD_IS_ORTHOLOG", "OLD_DISPOSITION_SOURCE",
                        "TREE_TRACK_USED", "TREE_CALL_USED", "ARM_SUPPORT", "DISPOSITION", "IS_ORTHOLOG",
                        "DISPOSITION_SOURCE", "BASIS"]) + "\n")
    for x in rows:
        fh.write("\t".join(x) + "\n")
print("WROTE %s sha256 %s" % (OUT, sha(OUT)))
