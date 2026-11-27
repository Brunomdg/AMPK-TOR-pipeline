#!/usr/bin/env python3
# >>> deposit header
# Purpose: Tabulates the number of independent origins of each trait's absence in every reconstruction block, and
#          its range per taxon subset, from the per-block reconstruction files.
# Methods: Ancestral-state reconstructions
# Inputs:  the per-block reconstruction file of each trait under the standard coding, resolved by `rl_rm_q7_lib.py`;
#          `STANDARD_per_species.tsv` and `STANDARD_vs_R37.tsv`
# Outputs: `Q7_ORIGINS.tsv` (one row per trait and block) and `Q7_ORIGINS_summary.tsv` (ranges per trait and subset)
#          in `14_remodeling/Q7_alltargets`, or in the folder named by Q7_OUT
# Tools:   Python 3 standard library; imports `rl_rm_q7_lib.py`; versions: see envs/
# <<< deposit header
"""Origins: the independent origins of every `<T>_absent` trait under the standard coding, read from the
ancestral-state blocks, never recomputed.

What it does. For each trait, one row per block (`TREE`, `ROOT_MODE`, `SUBSET`, `MODEL`) copied from the
blocks file that rl_rm_q7_lib.read_edges resolves (STANDARD_vs_R37.tsv names which reconstruction that is):
ORIGINS is the count of parent-0 to child-1 edges on the hardened reconstruction (what the loss counts report
as ORIGINS_MODEL for ER, SYM and ARD, and as ORIGINS_PARSIMONY for MPR), with N_STATE1, N_STATE0, N_FLAT_PRIOR
and AICC beside it. A second file summarises each trait, coding and subset: the minimum and maximum of ORIGINS
across the blocks of that subset (and separately for parsimony), the number of blocks, the ranges of N_STATE1
and of the tip count, and the number of blocks missing. Nothing is summed across blocks and no block is
preferred; the reading rule of the loss counts applies: compare models within a block by AICC, and read
recurrence off the minimum.

Sign. State 1 of every `<T>_absent` trait is the absence of a called ortholog, so its origins are origins of
absence; between a core component and a replacement, the shared direction is the one compatible with decay,
never reassignment; whether an absence is a loss is decided by the reconstruction, not by the trait name.

Inputs: for each trait, the edges and blocks files under the standard coding. read_edges refuses unless the
edges file cites the per-species table by its md5 checksum, and this script refuses unless the blocks file
lists 112 blocks.
Outputs: 14_remodeling/Q7_alltargets/Q7_ORIGINS.tsv and Q7_ORIGINS_summary.tsv.
Env: PROJECT_ROOT; Q7_OUT overrides the output directory (tests only).
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import collections
import hashlib
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rl_rm_q7_lib as L                                          # noqa: E402

ROOT = os.environ.get("PROJECT_ROOT") or (
    (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1"))
ASRD = os.path.join(ROOT, "14_remodeling", "Q0_replacement", "asr")
OUTD = os.environ.get("Q7_OUT") or os.path.join(ROOT, "14_remodeling", "Q7_alltargets")
PER_SPECIES = os.path.join(ROOT, "14_remodeling", "Q0_replacement", "Q0_REPLACEMENT_per_species.tsv")
CODINGS = (L.CODING,)
SIGN_LINE = ("## SIGN (DECISIONS I7, restated under L43): state 1 of every `<T>_absent` trait is the absence of a called "
             "ortholog; between a core and a replacement the shared direction is the decay-compatible one, never "
             "reassignment; whether an absence is a loss is what the reconstruction says, never the name.")


def refuse(msg):
    sys.exit("REFUSED: %s" % msg)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rows(p):
    import csv
    if not os.path.exists(p):
        refuse("%s is missing" % p)
    with open(p, errors="replace") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")),
                                   delimiter="\t", quoting=csv.QUOTE_NONE))


def traits_from_header():
    return L.traits()


def main():
    traits = traits_from_header()
    src = L.sources()
    out, summ, inputs = [], [], []
    for t in traits:
        keys = {}
        for cod in CODINGS:
            _e, _b, _fe, fb = L.read_edges(t, src)
            data = rows(fb)
            inputs.append(("%s (source %s for %s)" % (os.path.relpath(fb, ROOT), src[t][1], t), md5(fb)))
            keys[cod] = set((r["TREE"], r["ROOT_MODE"], r["SUBSET"], r["MODEL"]) for r in data)
            if len(data) != 112:
                refuse("%s has %d blocks, expected 112" % (fb, len(data)))
            per_subset = collections.defaultdict(list)
            for r in data:
                out.append([t, cod, r["TREE"], r["ROOT_MODE"], r["SUBSET"], r["MODEL"], r["STATUS"],
                            r["N_TIP"], r["N_STATE1"], r["N_STATE0"], r["N_FLAT_PRIOR"], r["AICC"],
                            r["ORIGINS"], r["PRUNED_TREE_MD5"]])
                if r["STATUS"] == "OK" and r["ORIGINS"] not in ("", "NA"):
                    per_subset[r["SUBSET"]].append((int(r["ORIGINS"]), r["MODEL"], int(r["N_STATE1"]), int(r["N_TIP"])))
            for sub, v in sorted(per_subset.items()):
                o = [x[0] for x in v]
                mpr = [x[0] for x in v if x[1] == "MPR"]
                summ.append([t, cod, sub, len(v), min(o), max(o),
                             min(mpr) if mpr else "NA", max(mpr) if mpr else "NA",
                             min(x[2] for x in v), max(x[2] for x in v), min(x[3] for x in v), max(x[3] for x in v),
                             112 // 4 - len(v)])
        if len(keys[L.CODING]) != 112:
            refuse("%s: %d blocks enumerated, expected 112" % (t, len(keys[L.CODING])))
    os.makedirs(OUTD, exist_ok=True)
    prov = ["## PRODUCED BY %s md5=%s" % (os.path.basename(__file__), md5(os.path.abspath(__file__))),
            "## INPUT per_species (STANDARD, L43) = %s md5=%s" % (os.path.relpath(L.STANDARD, ROOT), md5(L.STANDARD)),
            "## INPUT source map = %s md5=%s" % (os.path.relpath(L.VS, ROOT), md5(L.VS))]
    prov += ["## INPUT %s md5=%s" % (p, m) for p, m in inputs]
    p1 = os.path.join(OUTD, "Q7_ORIGINS.tsv")
    with open(p1 + ".tmp", "w") as fh:
        fh.write("## Q7 ORIGINS per block: one row per (TRAIT, CODING=STANDARD, TREE, ROOT_MODE, SUBSET, MODEL), copied from the\n"
                 "## ancestral-state blocks the source map resolves (L43). ORIGINS counts parent-0 to child-1 edges on the hardened states (corHMM marginal\n"
                 "## for ER/SYM/ARD, phangorn MPR for MPR). N_FLAT_PRIOR tips are coded 1&2 (R5), kept in the likelihood and\n"
                 "## excluded from every count. Compare models WITHIN a block by AICC; across blocks the data differ.\n"
                 "## NOTHING IS SUMMED ACROSS BLOCKS AND NO BLOCK IS A VOTE.\n" + SIGN_LINE + "\n" + "\n".join(prov) + "\n")
        fh.write("TRAIT\tCODING\tTREE\tROOT_MODE\tSUBSET\tMODEL\tSTATUS\tN_TIP\tN_STATE1\tN_STATE0\tN_FLAT_PRIOR\tAICC\tORIGINS\tPRUNED_TREE_MD5\n")
        for r in out:
            fh.write("\t".join(str(x) for x in r) + "\n")
    os.replace(p1 + ".tmp", p1)
    p2 = os.path.join(OUTD, "Q7_ORIGINS_summary.tsv")
    with open(p2 + ".tmp", "w") as fh:
        fh.write("## Q7 ORIGINS summary: per (TRAIT, CODING, SUBSET) the minimum and maximum of ORIGINS across the OK blocks of\n"
                 "## that subset (7 TREE x ROOT_MODE combinations x 4 models = 28 blocks per subset), the MPR (parsimony) floor\n"
                 "## separately, and the range of N_STATE1 and N_TIP. Read recurrence off the MINIMUM (Q2's rule); the maximum\n"
                 "## says how far a model or topology can move it. N_BLOCKS_MISSING counts blocks that failed or are absent.\n"
                 "## NEVER a sum, NEVER a mean: a range is reported because the blocks are perturbations, not replicates.\n" + SIGN_LINE + "\n" + "\n".join(prov) + "\n")
        fh.write("TRAIT\tCODING\tSUBSET\tN_BLOCKS\tMIN_ORIGINS\tMAX_ORIGINS\tMIN_ORIGINS_MPR\tMAX_ORIGINS_MPR\tMIN_N_STATE1\tMAX_N_STATE1\tMIN_N_TIP\tMAX_N_TIP\tN_BLOCKS_MISSING\n")
        for r in summ:
            fh.write("\t".join(str(x) for x in r) + "\n")
    os.replace(p2 + ".tmp", p2)
    print("wrote %s (%d rows)\nwrote %s (%d rows)" % (p1, len(out), p2, len(summ)))
    for r in summ:
        if r[2] == "all":
            print("  %-18s all blocks: origins %s to %s (MPR %s to %s), state1 %s to %s" % (r[0], r[4], r[5], r[6], r[7], r[8], r[9]))


if __name__ == "__main__":
    main()
