#!/usr/bin/env python3
# >>> deposit header
# Purpose: Repaired path resolver for the origin and co-transition analyses: finds, for every trait under the
#          standard coding, its per-edge and per-block reconstruction files and checks that they were built
#          from the per-species table they must cite.
# Methods: Ancestral-state reconstructions
# Inputs:  `STANDARD_per_species.tsv`, `STANDARD_vs_R37.tsv` (which reconstruction each trait uses) and
#          `ACTIVE_PUBLICATION_TARGETS.tsv`; the reconstruction files it resolves
# Outputs: none of its own (imported by `rl_rm_q7_origins.py` and by the co-transition count)
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
# Ran as: rl_rm_q7_lib.py
"""Shared path resolver for the origin and co-transition analyses under the standard coding.

The standard coding (STANDARD_per_species.tsv) uses one rule and one `<T>_absent` name for every
trait: 0 where an ortholog was called, 0.5 where the proteome cannot support an absence call, and
1 where no ortholog was called. traits() requires one `<T>_absent_p1` column per active publication
target. For each trait, the SOURCE_ASR column of STANDARD_vs_R37.tsv names the reconstruction to read:
`R37` reuses the reconstruction made for the replacement analysis under its second coding (files
`Q0_replacement/asr/<Q0 name>__R37.*`); `STANDARD` reads the reconstruction made afresh under the
standard coding (files `<T>_absent__STANDARD.*` in Q7_alltargets/asr). In the ten-target build every
trait is reconstructed afresh. In this repaired copy, the folder for the standard coding can be
redirected with the environment variable Q7_STANDARD_ASR. read_edges() checks, for either source,
that the edge file was built from the per-species table it must cite (by md5) and that every block
in it is listed in the blocks file.
"""
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import collections, csv, hashlib, os, sys

ROOT = os.environ.get("PROJECT_ROOT") or (_TA210_SCRATCH + "/ms210_phase2/TOR-AMPK_orthologs_10target_R1")
Q7D = os.path.join(ROOT, "14_remodeling", "Q7_alltargets")
Q0D = os.path.join(ROOT, "14_remodeling", "Q0_replacement")
STANDARD = os.path.join(Q7D, "STANDARD_per_species.tsv")
ACTIVE = os.path.join(ROOT, "05_method_and_rules", "ACTIVE_PUBLICATION_TARGETS.tsv")
VS = os.path.join(Q7D, "STANDARD_vs_R37.tsv")
R37 = os.path.join(Q0D, "Q0_REPLACEMENT_per_species__R37coding.tsv")
CODING = "STANDARD"


def refuse(m):
    sys.exit("REFUSED: %s" % m)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rows(p):
    if not os.path.exists(p):
        refuse("%s is missing" % p)
    with open(p, errors="replace") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t", quoting=csv.QUOTE_NONE))


def banner(p):
    with open(p, errors="replace") as fh:
        return [l.rstrip("\n") for l in fh if l.startswith("##")]


def traits():
    with open(STANDARD, errors="replace") as fh:
        for ln in fh:
            if not ln.startswith("##"):
                t = [c[:-3] for c in ln.rstrip("\n").split("\t") if c.endswith("_p1")]
                break
    n = len([r for r in rows(ACTIVE)]) if os.path.exists(ACTIVE) else 0
    if n == 0 or len(t) != n or any(not x.endswith("_absent") for x in t):
        refuse("STANDARD header must carry one `<T>_absent_p1` column per active publication target (%d), found %s"
               % (n, t))
    return t


def sources():
    """Map each trait to (name in the replacement analysis, source coding, edges path, blocks path,
    per-species table that the reconstruction must cite)."""
    out = {}
    for r in rows(VS):
        t, q0n, src = r["TRAIT"], r["Q0_NAME"], r["SOURCE_ASR"]
        if src == "R37":
            d, cod, ps = os.path.join(Q0D, "asr"), "R37", R37
            name = q0n
        elif src == "STANDARD":
            ## Changed on 18 September 2026, for the standard coding only: Q7_STANDARD_ASR lets this
            ## read the repaired reconstructions in 19_asr_repair/asr; when it is unset, the path is
            ## as before. The branch above, for the second coding, is unchanged: those traits read
            ## Q0_replacement/asr, which the repair did not rebuild.
            d = os.environ.get("Q7_STANDARD_ASR") or os.path.join(Q7D, "asr")
            cod, ps = "STANDARD", STANDARD
            name = t
        else:
            refuse("SOURCE_ASR %r for %s" % (src, t))
        out[t] = (q0n, src, os.path.join(d, "%s__%s.edges.tsv" % (name, cod)), os.path.join(d, "%s__%s.blocks.tsv" % (name, cod)), ps)
    if set(out) != set(traits()):
        refuse("STANDARD_vs_R37.tsv does not list exactly the active standard traits")
    return out


def read_edges(trait, src=None):
    src = src or sources()
    q0n, sc, fe, fb, ps = src[trait]
    for f in (fe, fb):
        if not os.path.exists(f):
            refuse("%s is missing (source %s for %s)" % (f, sc, trait))
    rec = [ln.split("md5=")[-1].strip() for ln in banner(fe) if ln.startswith("## INPUT per_species = ")]
    if not rec or rec[0] != md5(ps):
        refuse("%s was reconstructed from per_species md5 %s; the %s file it must cite has md5 %s" % (fe, rec[0] if rec else "none", sc, md5(ps)))
    blocks = {(r["TREE"], r["ROOT_MODE"], r["SUBSET"], r["MODEL"]): r for r in rows(fb)}
    edges = collections.defaultdict(list)
    for r in rows(fe):
        edges[(r["TREE"], r["ROOT_MODE"], r["SUBSET"], r["MODEL"])].append(r)
    for k in edges:
        if k not in blocks:
            refuse("%s has edges for block %s that %s does not list" % (fe, k, fb))
    return edges, blocks, fe, fb
