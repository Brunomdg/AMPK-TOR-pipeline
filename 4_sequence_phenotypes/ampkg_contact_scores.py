#!/usr/bin/env python3
# >>> deposit header
# Purpose: Scores the nucleotide-contact positions of every retained AMPKγ copy: for three sets of ten positions,
#          how many can be scored and how many are identical to the human anchor, under a pairwise (primary) and a
#          joint alignment, with the identity of the non-contact CBS positions as a background.
# Methods: Sequence phenotypes of the present subunits
# Inputs:  `03_ampkg_contacts/AMPKg_orthologs_fulllength.faa` and `RBX1_sequence_manifest.AMPKg.tsv` (from job x1);
#          `AMPKg_contact_reference.tsv` (the three site sets in canonical and anchor numbering); `Homo_sapiens.faa`
#          (anchor NP_001193638.1) and `p54619_canonical.faa`; `gamma_pocket.tsv` of the original pocket analysis
# Outputs: `03_ampkg_contacts/AMPKg_contact_copy.tsv` (one row per copy) and `AMPKg_contact_long.tsv` (one row per
#          copy, alignment, site set and position); the pairwise and joint alignments; in `10_validation`,
#          `AMPKg_numbering_gates.tsv`, `AMPKg_legacy_reproduction.tsv`, `AMPKg_legacy_pairwise_reproduction_delta.tsv`
# Tools:   MAFFT 7.526 [env record] from the orthofinder3 environment (--globalpair --maxiterate 1000); Python 3;
#          versions: see envs/
# <<< deposit header
# Ran as: rbx1_gamma.py
## Job x2 of this chain, AMPKγ branch. The alignment rule for the AMP-contact score was chosen on 14 September
## 2026: the primary method (PAIRWISE_LEGACY_P02) aligns each copy to the anchor pairwise, exactly as the original
## AMPKγ pocket analysis (folder `P02_gamma_pocket`) did, because that keeps this phenotype continuous with the
## earlier pocket-erosion trait (`gamma_eroded`). One joint alignment of the anchor with every copy is an
## alignment-procedure sensitivity declared in advance. The primary site set is site1; site3 and site4 are
## site-definition sensitivities declared in advance. Scores are never pooled, averaged or chosen per protein,
## and a position that cannot be scored is never counted as a lost residue.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "4_sequence_phenotypes", "lib"))
import phenotype_tables_chain1 as L  # noqa: E402
import sequence_scoring as S  # noqa: E402

SELF = os.path.abspath(__file__)
ANCHOR = "NP_001193638.1"          # anchor of the original analysis; its numbering is canonical PRKAG1 plus 9
CANON = "P54619"                   # UniProt canonical PRKAG1
CBS_CANON_FROM, CBS_CANON_TO = 70, 317   # CBS span of the original coverage measure (cbs_coverage), canonical numbering
LOW_COVERAGE = 0.5                 # reliability rule of the original analysis: coverage below 0.5 is LOW
PRIMARY_METHOD = "PAIRWISE_LEGACY_P02"


def load_reference():
    ref = L.rows(os.path.join(L.X, "03_ampkg_contacts/AMPKg_contact_reference.tsv"))
    sites = {}
    for r in ref:
        sites.setdefault(r["SITESET"], []).append((int(r["CANONICAL_POS_PRKAG1"]), r["EXPECTED_AA"],
                                                   int(r["ANCHOR_INDEX"]), r["ANCHOR_AA"]))
    for k in sites:
        sites[k].sort()
        if len(sites[k]) != 10:
            sys.exit("REFUSED: site set %s does not hold ten positions" % k)
    return sites


def score_set(anchor_aln, qry_aln, colmap, positions):
    per, callable_n, exact_n = [], 0, 0
    for canon, exp_aa, anchor_idx, anchor_aa in positions:
        state, exact, aa = S.position_call(anchor_aln, qry_aln, colmap[anchor_idx])
        if state in S.CALLABLE_STATES:
            callable_n += 1
            exact_n += exact
        per.append((canon, exp_aa, anchor_aa, aa, state, exact))
    return per, callable_n, exact_n


def background(anchor_aln, qry_aln, colmap, allpos):
    bg_call = bg_id = cov = 0
    for canon_pos in range(CBS_CANON_FROM, CBS_CANON_TO + 1):
        col = colmap[canon_pos + 9]
        state, exact, _aa = S.position_call(anchor_aln, qry_aln, col)
        if qry_aln[col] != S.GAP:
            cov += 1
        if canon_pos in allpos:
            continue
        if state in S.CALLABLE_STATES:
            bg_call += 1
            bg_id += exact
    return bg_call, bg_id, cov / float(CBS_CANON_TO - CBS_CANON_FROM + 1)


def main():
    L.verify_freezes()
    sites = load_reference()
    allpos = {p for v in sites.values() for p, _e, _a, _aa in v}
    man = {r["ALIGNMENT_KEY"]: r for r in L.rows(os.path.join(L.X, "02_inputs_manifest/RBX1_sequence_manifest.AMPKg.tsv"))
           if r["ALIGNMENT_KEY"] != "-"}
    copies = S.read_fasta(os.path.join(L.X, "03_ampkg_contacts/AMPKg_orthologs_fulllength.faa"))
    human = S.read_fasta(os.path.join(L.PROT, "Homo_sapiens.faa"))
    canon = S.read_fasta(os.path.join(L.CONCORD, "ms210/p54619_canonical.faa"))
    if ANCHOR not in human:
        sys.exit("REFUSED: anchor %s is not in the frozen Homo_sapiens proteome" % ANCHOR)
    anchor_seq, canon_seq = human[ANCHOR], canon[list(canon)[0]]

    gates = []
    badn = [(s, p) for s, v in sites.items() for p, e, ai, aa in v if anchor_seq[ai - 1] != aa]
    gates.append(["anchor residues match the frozen contact record at every anchor_index", "PASS" if not badn else "FAIL", str(badn[:5])])
    badc = [(s, p, e, canon_seq[p - 1]) for s, v in sites.items() for p, e, ai, aa in v if canon_seq[p - 1] != e]
    gates.append(["canonical PRKAG1 (%s) carries the expected residue at every canonical position" % CANON,
                  "PASS" if not badc else "REPORTED_DIFFERENCE", str(badc[:5]) or "30/30"])
    badoff = [(s, p) for s, v in sites.items() for p, e, ai, aa in v if ai - p != 9]
    gates.append(["the frozen +9 anchor offset holds for every position", "PASS" if not badoff else "FAIL", str(badoff[:5])])
    if badn or badoff:
        L.write("10_validation/AMPKg_numbering_gates.tsv", ["GATE", "RESULT", "DETAIL"], gates, script=SELF)
        sys.exit("REFUSED: the frozen contact numbering does not reproduce against the anchor")

    # sensitivity procedure: one joint alignment of the anchor, the canonical reference and every copy
    aln_path = os.path.join(L.X, "03_ampkg_contacts/AMPKg_contact_alignment_joint.faa")
    jcmd, jaln = S.mafft(dict([(ANCHOR, anchor_seq), (CANON, canon_seq)] + sorted(copies.items())), aln_path)
    L.WRITTEN.append(aln_path)
    jcolmap = S.ref_column_map(jaln[ANCHOR])

    pw_dir = os.path.join(L.X, "03_ampkg_contacts/pairwise")
    os.makedirs(pw_dir, exist_ok=True)
    long_rows, copy_rows, per_primary = [], [], {}
    pcmd = ""
    for key in sorted(copies):
        m = man[key]
        sp, acc = key.split("|", 1)
        # primary: the pairwise procedure of the original analysis, each copy against the anchor
        pcmd, pa = S.mafft({ANCHOR: anchor_seq, key: copies[key]},
                           os.path.join(pw_dir, key.replace("|", "__") + ".faa"))
        pcolmap = S.ref_column_map(pa[ANCHOR])
        row = [m["TARGET"], sp, acc, m["SEQUENCE_SHA256"], m["DECISION_SOURCE"], m["MODEL_SENSITIVE"],
               m["SOURCE_TYPE"], m["N_ORTHOLOG_COPIES_IN_SPECIES"], PRIMARY_METHOD, L.PRIMARY_SITESET]
        prim, joint = {}, {}
        for st in L.SITESETS:
            per, cn, en = score_set(pa[ANCHOR], pa[key], pcolmap, sites[st])
            prim[st] = (cn, en, per)
            jper, jcn, jen = score_set(jaln[ANCHOR], jaln[key], jcolmap, sites[st])
            joint[st] = (jcn, jen)
            for tag, rows_ in (("PAIRWISE", per), ("JOINT", jper)):
                for canon_pos, exp_aa, anchor_aa, aa, state, exact in rows_:
                    long_rows.append([m["TARGET"], sp, acc, tag, st, canon_pos, exp_aa, anchor_aa, aa, state,
                                      "yes" if state in S.CALLABLE_STATES else "no", exact])
        per_primary[acc] = prim
        for st in L.SITESETS:
            cn, en, _p = prim[st]
            row += [cn, en, ("%.3f" % (en / float(cn))) if cn else "NA", (en if cn == 10 else "NA")]
        for st in L.SITESETS:
            jcn, jen = joint[st]
            row += [jcn, (jen if jcn == 10 else "NA")]
        for st in L.SITESETS:
            cn, en, _p = prim[st]
            jcn, jen = joint[st]
            ps, js = (en if cn == 10 else None), (jen if jcn == 10 else None)
            row += ["NA" if None in (ps, js) else ps - js,
                    "yes" if (ps is not None and ps == js) else ("no" if None not in (ps, js) else "NOT_COMPARABLE")]
        bg_call, bg_id, coverage = background(pa[ANCHOR], pa[key], pcolmap, allpos)
        s1c, s1e, _ = prim[L.PRIMARY_SITESET]
        bgi = (bg_id / float(bg_call)) if bg_call else None
        cmb = ((s1e / float(s1c)) - bgi) if (s1c and bgi is not None) else None
        row += [bg_call, bg_id, ("%.4f" % bgi) if bgi is not None else "NA",
                ("%.4f" % cmb) if cmb is not None else "NA", "%.3f" % coverage,
                "LOW" if coverage < LOW_COVERAGE else "OK"]
        copy_rows.append(row)

    hdr = ["TARGET", "SPECIES", "PROTEIN_ID", "SEQUENCE_SHA256", "DECISION_SOURCE", "MODEL_SENSITIVE", "SOURCE_TYPE",
           "N_ORTHOLOG_COPIES_IN_SPECIES", "PRIMARY_ALIGNMENT_METHOD", "AMP_CONTACT_SITESET_PRIMARY"]
    for st in L.SITESETS:
        u = st.upper()
        hdr += ["N_CALLABLE_%s_PAIRWISE" % u, "N_EXACT_%s_PAIRWISE" % u, "EXACT_FRACTION_CALLABLE_%s_PAIRWISE" % u,
                "AMP_CONTACT_SCORE_%s_PAIRWISE" % u]
    for st in L.SITESETS:
        u = st.upper()
        hdr += ["N_CALLABLE_%s_JOINT" % u, "AMP_CONTACT_SCORE_%s_JOINT" % u]
    for st in L.SITESETS:
        u = st.upper()
        hdr += ["PAIRWISE_JOINT_SCORE_DELTA_%s" % u, "PAIRWISE_JOINT_IDENTICAL_%s" % u]
    hdr += ["CBS_NONCONTACT_N_CALLABLE", "CBS_NONCONTACT_N_IDENTICAL", "CBS_NONCONTACT_IDENTITY",
            "CONTACT_MINUS_BACKGROUND", "CBS_COVERAGE", "RELIABILITY"]
    L.write("03_ampkg_contacts/AMPKg_contact_copy.tsv", hdr, copy_rows, script=SELF,
            notes=["One row per frozen R1 AMPKg ORTHOLOG copy; no copy is selected by its score.",
                   "PRIMARY is the frozen legacy P02 pairwise procedure; JOINT is the alignment-procedure sensitivity.",
                   "The two are never averaged and never chosen per protein.",
                   "A score is emitted only at ten callable positions; partial information stays in N_CALLABLE.",
                   "Background identity uses every callable non-contact position of canonical CBS span %d-%d," % (CBS_CANON_FROM, CBS_CANON_TO),
                   "excluding the union of the three frozen site sets, on the primary alignment.",
                   "MAFFT primary: " + pcmd, "MAFFT joint: " + jcmd])
    L.write("03_ampkg_contacts/AMPKg_contact_long.tsv",
            ["TARGET", "SPECIES", "PROTEIN_ID", "ALIGNMENT_METHOD", "SITESET", "CANONICAL_POS_PRKAG1", "EXPECTED_AA",
             "ANCHOR_AA", "OBSERVED_AA", "POSITION_STATE", "CALLABLE", "EXACT_MATCH"], long_rows, script=SELF,
            notes=["ALIGNMENT_METHOD PAIRWISE is the primary frozen legacy procedure; JOINT is the sensitivity.",
                   "INTERNAL_DELETION is callable and scores 0; a position beyond the observed sequence span is",
                   "not callable and is never scored as a lost residue."])

    # ---- reproduction of the original site1 scores (gamma_pocket.tsv) under the primary procedure ----
    leg = {r["protein"]: r for r in L.rows(os.path.join(L.P02, "gamma_pocket.tsv"))}
    ip, ij = hdr.index("AMP_CONTACT_SCORE_SITE1_PAIRWISE"), hdr.index("AMP_CONTACT_SCORE_SITE1_JOINT")
    inc = hdr.index("N_CALLABLE_SITE1_PAIRWISE")
    rep, delta, pm, jm, both = [], [], 0, 0, 0
    for r in copy_rows:
        lr = leg.get(r[2])
        if not lr:
            continue
        try:
            ls = int(lr["site1_ident"])
        except ValueError:
            continue
        both += 1
        pm += 1 if r[ip] == ls else 0
        jm += 1 if r[ij] == ls else 0
        rep.append([r[1], r[2], ls, r[ip], r[ij], "yes" if r[ip] == ls else "no", "yes" if r[ij] == ls else "no",
                    lr["cbs_coverage"], lr["reliability"]])
        if r[ip] != ls:
            per = per_primary[r[2]][L.PRIMARY_SITESET][2]
            ncall = r[inc]
            cause = ("CALLABILITY_RULE_DIFFERENCE: the frozen legacy counted every non-matching position, including "
                     "positions beyond the observed sequence span; RBX1 does not score those") if ncall < 10 else \
                    "IMPLEMENTATION_OR_INPUT_DIFFERENCE_UNRESOLVED"
            delta.append([r[1], r[2], ls, r[ip], ncall, r[ij], r[3],
                          ";".join("%d:%s>%s:%s" % (c, e, a, st) for c, e, _aa, a, st, _x in per),
                          "NOT_STORED_IN_FROZEN_RECORD", L.sha(os.path.join(L.PROT, "Homo_sapiens.faa")),
                          pcmd, lr["cbs_coverage"], lr["reliability"], cause,
                          "PASS_WITH_REPRODUCTION_DELTA" if ncall < 10 else "UNRESOLVED_REQUIRES_REVIEW"])
    L.write("10_validation/AMPKg_legacy_reproduction.tsv",
            ["SPECIES", "PROTEIN_ID", "LEGACY_P02_SITE1_IDENT", "RBX1_SITE1_PAIRWISE_PRIMARY", "RBX1_SITE1_JOINT",
             "PAIRWISE_MATCHES_LEGACY", "JOINT_MATCHES_LEGACY", "LEGACY_CBS_COVERAGE", "LEGACY_RELIABILITY"], rep,
            script=SELF, notes=["Reproduction of the frozen legacy phenotype, protein by protein, for every protein",
                                "the legacy also scored. The primary procedure is the one this table shows reproduces it."])
    L.write("10_validation/AMPKg_legacy_pairwise_reproduction_delta.tsv",
            ["SPECIES", "PROTEIN_ID", "LEGACY_SITE1_IDENT", "RBX1_SITE1_PAIRWISE", "RBX1_N_CALLABLE_SITE1",
             "RBX1_SITE1_JOINT", "RBX1_SEQUENCE_SHA256", "RBX1_PER_SITE_CALLS", "LEGACY_PER_SITE_RESIDUES",
             "CURRENT_ANCHOR_PROTEOME_SHA256", "ALIGNMENT_COMMAND", "LEGACY_CBS_COVERAGE", "LEGACY_RELIABILITY",
             "CAUSE", "OUTCOME"], delta, script=SELF,
            notes=["Every protein where the primary pairwise procedure does not reproduce the frozen legacy site1",
                   "value. The frozen record stores counts only, so the legacy per-site residues are not recoverable.",
                   "Nothing is tuned to force agreement; an unresolved row stops the gamma branch at the x7 gate."])
    gates.append(["Homo_sapiens site1 self-score is 10/10 under the primary procedure", "PASS"
                  if any(r[ip] == 10 for r in copy_rows if r[1] == "Homo_sapiens") else "FAIL",
                  str(sorted((r[ip], r[2]) for r in copy_rows if r[1] == "Homo_sapiens"))])
    gates.append(["RBX1 reproduces the frozen legacy P02 site1 scores", "REPORTED",
                  "primary pairwise %d/%d proteins, joint alignment %d/%d" % (pm, both, jm, both)])
    nid = sum(1 for r in copy_rows if r[hdr.index("PAIRWISE_JOINT_IDENTICAL_SITE1")] == "yes")
    gates.append(["the primary and the alignment-sensitivity procedures agree per copy", "REPORTED",
                  "%d of %d copies identical on site1" % (nid, len(copy_rows))])
    gates.append(["every legacy reproduction delta is explained", "PASS" if all(d[-1] == "PASS_WITH_REPRODUCTION_DELTA" for d in delta)
                  else "REPORTED", "%d delta rows: %s" % (len(delta), sorted({d[-2] for d in delta}))])
    L.write("10_validation/AMPKg_numbering_gates.tsv", ["GATE", "RESULT", "DETAIL"], gates, script=SELF)
    L.record_outputs("x2_gamma", SELF)
    n10 = sum(1 for r in copy_rows if r[ip] != "NA")
    print("AMPKg: %d copies, %d with ten callable site1 positions (primary pairwise), legacy reproduction %d/%d"
          % (len(copy_rows), n10, pm, both))
    hard = [g for g in gates if g[1] not in ("PASS", "REPORTED", "REPORTED_DIFFERENCE")]
    if hard:
        sys.exit("REFUSED: a numbering gate failed: %s" % hard)


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
