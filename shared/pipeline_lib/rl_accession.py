#!/usr/bin/env python3
# >>> deposit header
# Purpose: Provide the one normalisation of protein accessions used to join results across tools, and a
#          join key for FastTree output; run directly, check that normalisation merges no two proteins of a genome.
# Methods: several steps (shared helper)
# Inputs:  when run directly: a folder of proteome FASTA files (default `proteomes_norm/`)
# Outputs: when run directly: counts of rewritten accessions and of genomes with merged accessions, printed
# Tools:   Python 3 standard library; versions: see envs/
# <<< deposit header
"""One normalisation of protein accessions, shared by every script that joins evidence across tools.

Why this exists. Each tool rewrites characters it cannot carry, and each rewrite was once found
separately instead of being handled in one place; the same kind of error cost work three times in
one day:

  FastTree      truncates at ':' (the Newick branch-length separator), so 17 Trypanosoma vivax
                tips, one of them a decoy, disappeared from the FastTree gene tree, which shrinks
                the denominator of a false-positive rate.
  IQ-TREE       replaces ':' and '@' with '_', so tips survived under names that no longer matched
                the sequence table, and rows came back with ROLE='?'.
  OrthoFinder   replaces ':' with '_' in Orthogroups.txt. The proteome set holds
                'TvY486_0200210:mRNA-p1' and the clustering 'TvY486_0200210_mRNA-p1'. All 5,829
                T. vivax proteins are present and none carries the colon, so 7 AMPKα candidates were
                labelled NO_COVERAGE although they had been clustered correctly.

Monocystis agilis accessions also contain '@' and a literal '__'; the latter collides with the
separator of the 'ROLE__SPECIES__PROTEIN' tip-name convention.

The rule. Join accessions through normalise() on both sides. Never compare raw strings across
tools, and never read a failure to match as missing data: an accession absent from a tool's output
means the tool did not report that protein, not that the protein was never searched.

The normalisation is deliberately narrow, touching only the characters these tools are known to
rewrite. Anything broader risks merging two real accessions into one, which is worse than a missed
join; collide() checks that no such merge happens.
"""

# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
_REWRITTEN = ":@"


def normalise(acc):
    """The accession as any of these tools would write it.

    This is the single normalisation. Accessions must be joined through one normalisation on
    both sides; until 4 September 2026 the project had three, so the requirement looked met while
    it was not:

        rl_accession.normalise  ':' and '@' -> '_'
        rl_join.normalise       the same, plus removal of an InterProScan
                                'Species||ACCESSION' prefix
        .split(":")[0]          truncation at ':', written by hand in 16 scripts

    On 'Trypanosoma_vivax||TvY486_0303830:mRNA-p1' these give three different keys, so two tables
    that each followed a normalisation still failed to join. That is how two real KinOrtho MEMBER
    calls were lost from the AMPKα table on 4 September 2026.

    The '||' removal is included here because it is the superset and changes nothing on this
    proteome set: none of the 6,288 accessions in the master call table contains '||' (measured,
    not assumed). rl_join now imports this function instead of defining a second one.

    Truncation is not used: it discards the isoform suffix, so it cannot be mapped back to the
    real accession and it can merge two isoforms of one gene. collide() checks that no merge has
    happened.
    """
    if not acc:
        return acc
    out = acc.split("||")[-1]        # InterProScan  Species||ACCESSION
    for ch in _REWRITTEN:
        out = out.replace(ch, "_")
    return out


def tree_key(acc):
    """Join key for anything that has passed through FastTree. It loses information on purpose.

    Joining through one normalisation on both sides cannot be done with a single function
    (measured on 4 September 2026), because two tools do incompatible things to the same accession:

        FastTree   truncates at ':'   -> AMPKa|Trypanosoma_vivax|TvY486_0303830
        IQ-TREE    replaces ':'       -> TvY486_0303830_mRNA-p1

    Only truncation maps a truncated tip and a substituted tip onto the same key, and truncation
    loses information: it cannot be mapped back to the real accession and it can merge two
    isoforms of one gene.

    Measured on the AMPKα FastTree tree (6,475 tips) against the 786 members:
        truncate    786 tips match      (correct)
        normalise   784 tips match      (silently loses the two T. vivax proteins)

    The rule is therefore one normalisation per join, chosen to match the side that loses more:
    use tree_key() when a FastTree product is on either side and normalise() everywhere else.
    Never write a tree_key out as an accession; it is a key only. collide() checks normalise()
    for merged proteins; it does not test tree_key().
    """
    if not acc:
        return acc
    return acc.split("||")[-1].split(":")[0].replace("@", "_")


def index(accessions):
    """normalised -> original. Raises if two originals collapse to one key."""
    out = {}
    for a in accessions:
        k = normalise(a)
        if k in out and out[k] != a:
            raise ValueError(
                "accession collision: %r and %r both normalise to %r -- "
                "normalisation must never merge two real proteins"
                % (out[k], a, k))
        out[k] = a
    return out


def collide(accessions):
    """Groups of accessions that would merge. Empty is the required result."""
    seen = {}
    bad = {}
    for a in accessions:
        k = normalise(a)
        if k in seen and seen[k] != a:
            bad.setdefault(k, {seen[k]}).add(a)
        seen[k] = a
    return {k: sorted(v) for k, v in bad.items()}


if __name__ == "__main__":
    import os
    import sys
    corpus = (sys.argv[1] if len(sys.argv) > 1
              else (_TA210_SCRATCH + "/proteomes_norm"))
    n_tot = n_rw = 0
    worst = {}
    for fn in sorted(os.listdir(corpus)):
        if not fn.endswith(".faa"):
            continue
        accs = []
        for line in open(os.path.join(corpus, fn), errors="replace"):
            if line.startswith(">"):
                a = line[1:].split()[0]
                accs.append(a)
                n_tot += 1
                if normalise(a) != a:
                    n_rw += 1
        c = collide(accs)
        if c:
            worst[fn[:-4]] = c
    print("corpus accessions          : %d" % n_tot)
    print("rewritten by normalisation : %d (%.3f%%)" % (n_rw, 100.0 * n_rw / n_tot))
    print("within-genome collisions   : %d genomes" % len(worst))
    for sp, c in list(worst.items())[:5]:
        print("   %s: %s" % (sp, list(c.items())[:2]))
    if worst:
        raise SystemExit("REFUSING: normalisation would merge distinct proteins")
    print("safe: no accession pair merges under normalisation")
