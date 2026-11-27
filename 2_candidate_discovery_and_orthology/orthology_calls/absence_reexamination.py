#!/usr/bin/env python3
# >>> deposit header
# Purpose: For one species, rescore every protein against the eight core components from the genome-wide
#          HMMER search without score thresholds, to re-examine absence calls: rank the clean proteins,
#          report both bit-score statistics, and measure how much of the winning profile the top proteins cover.
# Methods: Orthology calls
# Inputs:  `<species>.domtbl.gz` (HMMER domain table of the genome-wide search), `PANEL_DECLARATION_v2.tsv`
#          (target-pole and contaminant-pole profiles), `FINAL_ORTHOLOG_TABLE.tsv`, `07_contaminant_screen.tsv`;
#          the species is named in the environment variable ARS_SPECIES
# Outputs: `_absence_research_parts/<species>.tsv` (the top five clean proteins per target) and
#          `<species>.tsv.gate` (counts of agreement with the contaminant screen)
# Tools:   Python 3 standard library and rl_accession.py; reads HMMER hmmsearch output; versions: see envs/
# <<< deposit header
"""Re-examine every absence against the genome-wide profile search run without score thresholds, and
report both bit-score statistics.

Design principles, from the project's review of the literature on the sensitivity of absence calls:
the detection procedure is the same for all 210 proteomes, and any place where it was not is reported;
a control may not recover itself; accessions are normalised the same way on both sides of every join;
absence of a signal is evidence only where the sensitivity of that signal is known; and a method that
gives no opinion is never a vote against.

This is the second version of the script. Three problems of the first version, found on review on
10 September 2026, are corrected here:

  1. Scores were floored at zero. Version 1 started each side (the target-pole and the contaminant-pole
     profiles) at 0.0 and only raised it, so a genuine negative bit score that arrived after the
     opposite side had created the entry was discarded, and the call depended on the order of lines in
     the domain table (domtbl). The contaminant screen starts from minus infinity and has a NO_SCORES
     state; both disagreements with the screen in version 1 came from this. Each side now starts as
     None, and NO_SCORES exists.
  2. The accession normalisation was wrong. Version 1 used a bare `.split("||")[-1]`. The project uses
     one normalisation, rl_accession.normalise, because on 'Trypanosoma_vivax||TvY486_0303830:mRNA-p1'
     three different hand-written rules give three different keys. Version 1's rule dropped
     Trypanosoma vivax, Monocystis agilis and Heliospora caprellae from the agreement check with the
     screen, and with them 26 of the screen's 2,726 rows; the check then reported 2,698 of 2,700
     against a denominator it had itself narrowed. This version imports rl_accession, and the
     denominator of the check is every row of the screen.
  3. One statistic was wrong for at least one target. Version 1 read only the full-sequence score
     (domtbl field 8). That score sums every repeat hit, and the project had already measured it to be
     useless for the LST8 profile KOG0315: the per-target threshold pilot found 0% recall at a 1%
     false-positive rate from the sequence score alone and no sequence-score cut-off at any
     false-positive rate, and only the best-domain score separates. Both statistics are now written
     for each side, and neither is called the answer.

It also adds what the winner-take-all screen could not see. A rule that reads one protein per species
cannot see a gene split across two accessions: each half scores low, and the absence is confirmed by
the very defect that caused it. Measured case: Toxoplasma XP_018636102.1 (2,320 aa) covers profile
positions 650 to 1131 of PTHR13298 and XP_018636117.1 (2,072 aa) covers positions 74 to 421; together
they are 4,392 aa and 102.4 bits, against a single Neospora protein of 4,471 aa and 105.1 bits. This
blindness works in the direction that confirms absence, which is the direction that favours a finding
that absences hold.

The profile (hmm) coordinates were in every domtbl line but had not been read. This version reads them
and reports how much of the winning profile the best clean protein covers, and how much is covered once
the second, and then the third to fifth, best clean proteins are added. A clean protein is one with a
target-pole hit whose best target-pole full-sequence score is not below its best contaminant-pole
full-sequence score (or that has no contaminant-pole hit).

No threshold decides anything, and the order is the screen's own. A first attempt added proteins
greedily by coverage over every protein that hit the profile, which under `-E 1e5 --max` means
thousands of spurious alignments: for Toxoplasma RICTOR it produced 99.1% coverage of the profile with
a summed score of minus 118 bits. Proteins are therefore taken in the order the screen already ranks
them, by full-sequence score among the clean ones, and the cumulative coverage is reported rank by
rank. No cut-off is applied and no protein is accepted or rejected; the reader sees what the second
protein adds.

What this is not: it is not a more sensitive search, and no search is run here. The script reads the
domain table of the same single `hmmsearch -E 1e5 --domE 1e5 --max` run whose per-sequence table
(tblout) the screen reads: same profiles, same options, same numbers. Sensitivity is set by the
profiles and the options, and neither changed. What changes is the set of proteins and the threshold:
every protein of the proteome instead of the 2,726 nominated ones, without the published gathering
(GA) cut-off that candidate discovery applied at analysis time. That widening is real and useful, and
it produced every recovered candidate here, but it is a change of scope and threshold and is described
as one.
"""
import collections
import csv
import gzip
import hashlib
import os
import sys

def _project_root_from_file(start):
    """Return the project root: the nearest folder above this file that contains 05_method_and_rules.

    A copy of the project (such as the ten-target copy) must never fall back to, or hard-code, the root
    of the original project tree (decided on 16 September 2026). Refusing to guess is the point: a
    hard-coded root made a run against the wrong tree indistinguishable from a correct one, and in its
    unconditional form even setting PROJECT_ROOT did not help.

    The path is made absolute but symbolic links are deliberately not resolved. The ten-target copy
    holds 2,014 symbolic links into the original tree. Walking up from the link finds the tree the
    script is being run as part of; walking up from the link's target would find the tree the file was
    copied from, which for a link into the original tree is exactly the wrong-root error this removes.
    """
    d = os.path.dirname(os.path.abspath(start))
    while d != os.path.dirname(d):
        if os.path.isdir(os.path.join(d, "05_method_and_rules")):
            return d
        d = os.path.dirname(d)
    raise SystemExit("REFUSED: PROJECT_ROOT unset and no project root found above %s" % start)


ROOT = os.environ.get("PROJECT_ROOT") or _project_root_from_file(__file__)
sys.path.insert(0, os.path.join(ROOT, "01_scripts", "lib"))
import rl_accession                                            # noqa: E402

DOMT = os.path.join(ROOT, "03_evidence_by_method/hmmer_genomewide/domtbl")
PANEL = os.path.join(ROOT, "02_reference_data/hmmer_panel/PANEL_DECLARATION_v2.tsv")
FINAL = os.path.join(ROOT, "04_final_tables", "FINAL_ORTHOLOG_TABLE.tsv")
SCREEN = os.path.join(ROOT, "09_per_method_report__setA",
                      "07_contaminant_screen.tsv")
PART = os.path.join(ROOT, "07_findings", "_absence_research_parts")
SELF = os.path.abspath(__file__)

TARGETS = ["TOR", "RAPTOR", "RICTOR", "LST8", "AMPKa", "AMPKb", "AMPKg", "SIN1"]
SELF_CONTAM = {"RAPTOR__RAPTOR_s_own_lineage-specific_in-paralogs": "RAPTOR",
               "KOG0891.faa.final_tree.fa": "TOR",
               "TOR__ATR": "TOR"}
TOPN = 5

# Field positions in the HMMER domain table (domtblout), 0-based, checked against a real line
F_TARGET, F_TLEN, F_MODEL, F_QLEN = 0, 2, 3, 5
F_FULLBITS, F_DOMBITS, F_HMMFROM, F_HMMTO = 7, 13, 15, 16


def rows(p):
    return list(csv.DictReader((l for l in open(p, errors="replace")
                                if not l.startswith("##")), delimiter="\t"))


def poles():
    tgt = collections.defaultdict(set)
    con = collections.defaultdict(set)
    for r in rows(PANEL):
        t = (r.get("FOR_TARGET") or "").strip()
        m = (r.get("MODEL") or "").strip()
        pole = (r.get("POLE") or "").strip()
        if t not in TARGETS or not m:
            continue
        if pole == "target":
            tgt[t].add(m)
        elif pole == "contaminant":
            if SELF_CONTAM.get(m) == t:
                continue
            con[t].add(m)
    return tgt, con


def merge(iv):
    """Union of closed integer intervals, as a sorted disjoint list."""
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def covlen(iv):
    return sum(b - a + 1 for a, b in iv)


def pass1(f, model_t, model_c):
    """Return {protein: {target: [tfull, tdom, tmodel, cfull, cdom, cmodel]}}: for the target-pole and
    the contaminant-pole profiles, the best full-sequence score with the profile that gave it, and the
    best single-domain score. None means that side was never hit, which is not a score of zero."""
    best = collections.defaultdict(dict)
    with gzip.open(f, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) <= F_HMMTO:
                continue
            m = p[F_MODEL]
            ts, cs = model_t.get(m), model_c.get(m)
            if ts is None and cs is None:
                continue
            try:
                fb = float(p[F_FULLBITS])
                db = float(p[F_DOMBITS])
            except ValueError:
                continue
            pid = rl_accession.normalise(p[F_TARGET])
            d = best[pid]
            for pole, names in ((0, ts), (1, cs)):
                if not names:
                    continue
                for t in names:
                    e = d.get(t)
                    if e is None:
                        e = d[t] = [None, None, "", None, None, ""]
                    o = 0 if pole == 0 else 3
                    if e[o] is None or fb > e[o]:
                        e[o], e[o + 2] = fb, m
                    if e[o + 1] is None or db > e[o + 1]:
                        e[o + 1] = db
    return best


def pass2(f, wanted):
    """Second, targeted pass: for each profile in `wanted` ({profile: placeholder}), collect every
    protein's profile intervals (hmm from and to; merged later), its full-sequence score, and the
    profile length."""
    iv = collections.defaultdict(lambda: collections.defaultdict(list))
    sc = collections.defaultdict(dict)
    qlen = {}
    with gzip.open(f, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) <= F_HMMTO:
                continue
            m = p[F_MODEL]
            if m not in wanted:
                continue
            try:
                a, b = int(p[F_HMMFROM]), int(p[F_HMMTO])
                fb = float(p[F_FULLBITS])
            except ValueError:
                continue
            pid = rl_accession.normalise(p[F_TARGET])
            iv[m][pid].append((a, b))
            sc[m][pid] = fb
            qlen[m] = int(p[F_QLEN])
    return iv, sc, qlen


def coverage_series(per, order, qlen_m):
    """Cumulative profile coverage as proteins are added in score order.

    Score order is used rather than a greedy choice by coverage. Greedy selection over every protein
    that hits the profile is meaningless under `-E 1e5 --max`: it picks whichever spurious alignment
    spans the most positions, and for Toxoplasma RICTOR it chose a set summing to -118.4 bits. The
    screen's own ranking is used instead, so no new criterion decides anything.
    """
    cov, out = [], []
    for pid in order:
        v = per.get(pid)
        if v:
            cov = merge(cov + v)
        out.append(covlen(cov) / qlen_m)
    return out


def span(per, pid):
    v = per.get(pid)
    return "" if not v else "%d-%d" % (v[0][0], v[-1][1])


def overlap_frac(per, a, b):
    """Fraction of the smaller of the two proteins' profile coverages that the other protein also
    covers. Zero means the two occupy disjoint regions of the profile, which is what a gene split
    across two accessions looks like."""
    va, vb = per.get(a), per.get(b)
    if not va or not vb:
        return ""
    la, lb = covlen(va), covlen(vb)
    inter = covlen(va) + covlen(vb) - covlen(merge(va + vb))
    return "%.3f" % (inter / min(la, lb)) if min(la, lb) else ""


def main():
    sp = os.environ.get("ARS_SPECIES")
    if not sp:
        sys.exit("REFUSED: ARS_SPECIES must be set to one species name")
    tgt, con = poles()
    for t in TARGETS:
        if not tgt[t] or not con[t]:
            sys.exit("REFUSED: %s has %d target-pole and %d contaminant-pole models"
                     % (t, len(tgt[t]), len(con[t])))
    model_t = collections.defaultdict(list)
    model_c = collections.defaultdict(list)
    for t in TARGETS:
        for m in tgt[t]:
            model_t[m].append(t)
        for m in con[t]:
            model_c[m].append(t)

    f = os.path.join(DOMT, "%s.domtbl.gz" % sp)
    if not os.path.exists(f):
        sys.exit("REFUSED: no domtbl for %s at %s" % (sp, DOMT))
    best = pass1(f, model_t, model_c)

    known = collections.defaultdict(set)
    for r in rows(FINAL):
        if r["SPECIES"] == sp:
            known[r["TARGET"]].add(rl_accession.normalise(r["PROTEIN_ID"]))
    screen = {}
    for r in rows(SCREEN):
        if r["SPECIES"] == sp:
            screen[(r["TARGET"], rl_accession.normalise(r["PROTEIN_ID"]))] = r["CALL"]

    # ---- The call, made as the contaminant screen makes it: a side never hit counts as minus
    # infinity, and a protein with no score on either side is NO_SCORES.
    clean = collections.defaultdict(list)
    gate_n = gate_ok = gate_unjoined = 0
    seen_pairs = set()
    for pid, d in best.items():
        for t, e in d.items():
            tf, td, tm, cf, cd, cm = e
            if tf is None and cf is None:
                call = "NO_SCORES"
            elif cf is not None and (tf is None or cf > tf):
                call = "CONTAMINANT"
            else:
                call = "CLEAN"
            if (t, pid) in screen:
                seen_pairs.add((t, pid))
                gate_n += 1
                gate_ok += (screen[(t, pid)] == call
                            or screen[(t, pid)] == "CONTROL_EXEMPT")
            if call == "CLEAN":
                clean[t].append((tf, td, pid, tm, cf, cd, cm))
    gate_unjoined = len(set(screen) - seen_pairs)

    # ---- Second pass, only for the winning profile of each target (that of its top clean protein)
    wanted = {}
    for t in TARGETS:
        if clean[t]:
            wanted[max(clean[t])[3]] = 1
    iv, sc, qlen = pass2(f, wanted) if wanted else ({}, {}, {})

    out = []
    for t in TARGETS:
        cl = sorted(clean[t], reverse=True)
        novel = [c for c in cl if c[2] not in known[t]]
        ml, cser, sp1, sp2, ov = "", [], "", "", ""
        if cl:
            wm = cl[0][3]
            if wm in iv and wm in qlen and qlen[wm]:
                per = {pid: merge(v) for pid, v in iv[wm].items()}
                order = [c[2] for c in cl[:TOPN]]
                ml = qlen[wm]
                cser = coverage_series(per, order, qlen[wm])
                sp1 = span(per, order[0])
                if len(order) > 1:
                    sp2 = span(per, order[1])
                    ov = overlap_frac(per, order[0], order[1])
        c1 = "%.3f" % cser[0] if cser else ""
        c2 = "%.3f" % cser[1] if len(cser) > 1 else ""
        c5 = "%.3f" % cser[-1] if cser else ""
        for rank, (tf, td, pid, tm, cf, cd, cm) in enumerate(cl[:TOPN], 1):
            out.append([sp, t, rank, pid,
                        "known" if pid in known[t] else "not_in_candidate_pool",
                        "%.1f" % tf, "" if td is None else "%.1f" % td, tm,
                        "" if cf is None else "%.1f" % cf,
                        "" if cd is None else "%.1f" % cd, cm,
                        "" if cf is None else "%.1f" % (tf - cf),
                        len(cl), len(novel),
                        ml, c1, c2, c5, sp1, sp2, ov])
        if not cl:
            out.append([sp, t, "", "", "no_clean_protein"] + [""] * 6
                       + ["", 0, 0] + [""] * 7)

    os.makedirs(PART, exist_ok=True)
    pf = os.path.join(PART, "%s.tsv" % sp)
    with open(pf, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["SPECIES", "TARGET", "RANK", "PROTEIN_ID", "PROVENANCE",
                    "TARGET_BITS_FULLSEQ", "TARGET_BITS_BESTDOM", "TARGET_MODEL",
                    "CONTAM_BITS_FULLSEQ", "CONTAM_BITS_BESTDOM", "CONTAM_MODEL",
                    "MARGIN_FULLSEQ", "N_CLEAN_IN_PROTEOME", "N_CLEAN_NOT_IN_POOL",
                    "MODEL_LEN", "MODEL_COV_RANK1", "MODEL_COV_CUM_RANK2",
                    "MODEL_COV_CUM_TOP5", "RANK1_HMM_SPAN", "RANK2_HMM_SPAN",
                    "RANK1_RANK2_OVERLAP"])
        for r in out:
            w.writerow(r)
    with open(pf + ".gate", "w") as fh:
        fh.write("%s\t%d\t%d\t%d\t%d\n"
                 % (sp, gate_n, gate_ok, gate_unjoined, len(screen)))
    print("%s: %d proteins scored, %d rows, gate %d/%d matched, %d screen rows unjoined"
          % (sp, len(best), len(out), gate_ok, gate_n, gate_unjoined))


if __name__ == "__main__":
    main()
