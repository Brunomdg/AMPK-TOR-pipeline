#!/usr/bin/env python3
# >>> deposit header
# Purpose: Makes the per-copy and per-species CAMKK subfamily calls from the reference-masked tree, compares them with
#          the unmasked tree and the earlier membership-screen tree, and tests reference recovery by leave-one-out.
# Methods: Orthology calls
# Inputs:  `04_camkk_subfamilies/CAMKK_kinase_domain_manifest.tsv`, `CAMKK_tree_calls_masked.tsv`,
#          `CAMKK_tree_calls_unmasked.tsv`, `02_inputs/tree_qc.tsv`, `02_inputs/alignment_qc.tsv`,
#          `CAMKK_alignment_mask_accounting.tsv`, `CAMKK/CAMKK_contaminants_iqtree.tsv`, `T1_ortholog_table.tsv`,
#          the per-species domain scans `<species>.csv.gz`
# Outputs: `04_camkk_subfamilies/CAMKK_constituent_calls.tsv`, `CAMKK_constituent_species.tsv`,
#          `CAMKK_constituent_class_validation.tsv`, `CAMKK_reference_loo_detail.tsv`, `CAMKK_constituent_counts.tsv`
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
## The primary CAMKK tree is the reference-masked one and the unmasked tree is the sensitivity analysis. This order
## was fixed by the alignment check (rbx3_aln_qc.py), whose occupancy threshold was set before any topology was read
## and which the unmasked CAMKK alignment failed; it was not chosen by which tree gives cleaner subfamilies. A copy
## that the two trees place, with support, in different first-level classes is ALIGNMENT_SENSITIVE and supports no
## subfamily-specific evolutionary statement; its CAMKK presence in the species-state table is unchanged either
## way. Subfamily uncertainty never touches CAMKK presence or absence: subfamily labels say what kind of upstream
## kinase is retained, not whether one is.
import gzip, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "2_candidate_discovery_and_orthology", "upstream_kinases", "lib"))
import upstream_kinase_tables as C  # noqa: E402

SELF = os.path.abspath(__file__)
SUB = "04_camkk_subfamilies"
PROSITE = ["PS00107", "PS00108"]
CDD_L1 = {"cd14200": "CAMKK_METAZOAN_LIKE", "cd14199": "CAMKK_METAZOAN_LIKE",
          "cd14118": "CAMKK_METAZOAN_LIKE", "cd14008": "SHARED_LKB1_CAMKK_NODE_NOT_DIAGNOSTIC"}
NAMED_L2 = ("SAK1_TOS3_ELM1_LIKE", "GRIK_SNAK_LIKE")


def annotations(species, acc):
    p = os.path.join(C.COMBINED, "%s.csv.gz" % species)
    cdd, funfam, ps = [], [], []
    if os.path.exists(p):
        with gzip.open(p, "rt", errors="replace") as fh:
            for ln in fh:
                f = ln.rstrip("\n").split(",")
                if len(f) < 11 or f[1] != acc or f[2] != "fwd":
                    continue
                if f[3] == "CDD" and f[4] in CDD_L1:
                    cdd.append(f[4])
                elif f[3] == "IPS:FunFam" and f[5]:
                    funfam.append(f[5][:60])
                elif f[3] == "IPS:ProSitePatterns" and f[4] in PROSITE:
                    ps.append(f[4])
    ps = sorted(set(ps))
    st = "CANONICAL_BOTH" if len(ps) == 2 else ("PARTIAL_" + ps[0] if len(ps) == 1 else "NONE_DETECTED")
    return ";".join(sorted(set(cdd))) or "-", ";".join(sorted(set(funfam))[:2]) or "-", ";".join(ps) or "-", st


def _qc_desc():
    """Alignment numbers quoted in the table notes, read from this run's alignment-check tables, not typed in."""
    qc = {r["TARGET"]: r for r in C.rows(os.path.join(C.X3, "02_inputs/alignment_qc.tsv"))}["CAMKK"]
    mk = {r["QUANTITY"]: r["VALUE"] for r in C.rows(os.path.join(C.X3, SUB, "CAMKK_alignment_mask_accounting.tsv"))}
    occ = "%.1f" % (100.0 * float(qc["FRACTION_COLUMNS_OCCUPIED_BY_ANY_POSITIVE_REFERENCE"]))
    return ((mk["COLUMNS_AFTER"], mk["MEDIAN_GAP_FRACTION_AFTER"]),
            (qc["ALIGNMENT_COLUMNS"], qc["MEDIAN_GAP_FRACTION_PER_TIP"], occ,
             "failed" if qc["QC_VERDICT"] != "OK" else "passed"))


def main():
    C.verify_freezes()
    man = C.rows(os.path.join(C.X3, SUB, "CAMKK_kinase_domain_manifest.tsv"))
    M = {r["TIP"]: r for r in C.rows(os.path.join(C.X3, SUB, "CAMKK_tree_calls_masked.tsv"))}
    U = {r["TIP"]: r for r in C.rows(os.path.join(C.X3, SUB, "CAMKK_tree_calls_unmasked.tsv"), required=False)}
    if not M:
        sys.exit("REFUSED: the masked tree calls are primary and are missing")
    if annotations("Homo_sapiens", "NP_757343.2")[3] != "CANONICAL_BOTH":
        sys.exit("REFUSED: the CAMKK anchor fails its own signature instrument")
    t1 = {r["PROTEIN_ID"]: r for r in C.rows(os.path.join(C.PUB, "T1_ortholog_table_all20.tsv"))
          if r["TARGET"] == "CAMKK"}
    ## The comparison with the earlier CAMKK membership-screen tree (`CAMKK_contaminants_iqtree.tsv`) is made for each
    ## copy, before the species summary. A copy named here with support of at least 95 on the masked tree that the
    ## screen tree called a contaminant is demoted to unresolved; the species becomes unresolved only when no other
    ## copy provides a resolved, non-conflicting call, which the species rule below then handles.
    frozen_call = {}
    for r in C.rows(os.path.join(C.PHYLO, "CAMKK/CAMKK_contaminants_iqtree.tsv"), required=False):
        if r.get("PROTEIN_ID"):
            frozen_call[r["PROTEIN_ID"].split("||")[-1]] = r.get("TREE_CALL", "-")
    lens = [int(r["DOMAIN_LENGTH"]) for r in man
            if r["DOMAIN_LENGTH"] not in ("NA", "") and "REFERENCE" in r["ROLE"]]
    lo_ok, hi_ok = (min(lens) * 0.6, max(lens) * 1.6) if lens else (0, 1e9)

    out = []
    for r in sorted(man, key=lambda x: (x["SPECIES"], x["PROTEIN_ID"])):
        if r["ROLE"] == "NEGATIVE":
            continue
        tip = "%s__%s" % (r["SPECIES"], r["PROTEIN_ID"].replace(".", "_"))
        m = C.tip_in(M, tip)
        if m is None:
            sys.exit("REFUSED: no masked-tree call for %s; every copy must be accounted for" % tip)
        u = C.tip_in(U, tip) or {}
        mL1, mL2, mC, mS = m["CLASS_L1"], m["CLASS_L2"], m["CALL_CONFIDENCE"], m["TREE_SUPPORT"]
        uL1 = u.get("CLASS_L1", "NOT_AVAILABLE")
        uL2, uC, uS = u.get("CLASS_L2", "-"), u.get("CALL_CONFIDENCE", "NOT_AVAILABLE"), u.get("TREE_SUPPORT", "NA")
        # A conflict counts only when both trees name a first-level class with support of at least 70.
        both_named = mL1 != "SUBFAMILY_UNRESOLVED" and uL1 not in ("SUBFAMILY_UNRESOLVED", "NOT_AVAILABLE")
        supported_both = mC in ("SUPPORTED", "PROVISIONAL_PLACEMENT") and uC in ("SUPPORTED",
                                                                                 "PROVISIONAL_PLACEMENT")
        sens1 = "yes" if (both_named and supported_both and mL1 != uL1) else "no"
        m_named2, u_named2 = mL2 in NAMED_L2, uL2 in NAMED_L2
        sens2 = "yes" if (m_named2 and u_named2 and mL2 != uL2) else "no"
        call_l1 = "ALIGNMENT_SENSITIVE" if sens1 == "yes" else mL1
        call_l2 = "L2_UNRESOLVED" if (sens1 == "yes" or sens2 == "yes") else mL2
        fcall = frozen_call.get(r["PROTEIN_ID"], "NOT_IN_FROZEN_TREE")
        named_now = call_l1 not in ("SUBFAMILY_UNRESOLVED", "ALIGNMENT_SENSITIVE", "-")
        cconf = "yes" if (named_now and mC == "SUPPORTED" and fcall == "CONTAMINANT") else "no"
        if cconf == "yes":
            call_l1, call_l2 = "SUBFAMILY_UNRESOLVED", "L2_UNRESOLVED"
        dl = int(r["DOMAIN_LENGTH"]) if r["DOMAIN_LENGTH"] not in ("NA", "") else None
        cdd, funfam, ps, mstat = annotations(r["SPECIES"], r["PROTEIN_ID"])
        ann = {CDD_L1[x] for x in cdd.split(";") if x in CDD_L1} - {"SHARED_LKB1_CAMKK_NODE_NOT_DIAGNOSTIC"}
        conc = "NOT_DIAGNOSTIC" if not ann else ("AGREES" if call_l1 in ann else "DISAGREES")
        out.append([r["SPECIES"], r["PROTEIN_ID"], t1.get(r["PROTEIN_ID"], {}).get("FINAL_CLASS", "-"),
                    r["DECISION_SOURCE"], m["ROLE"], m["VALIDATION_TYPE"],
                    call_l1, call_l2, mL1, mS, mC, mL2, uL1, uS, uC, uL2, sens1, sens2, fcall, cconf,
                    m["REFERENCE_IDS"], m["NEGATIVE_TIPS"], m["CALL_BASIS"], mstat, ps,
                    dl if dl else "NA", "compatible" if (dl and lo_ok <= dl <= hi_ok) else "atypical",
                    cdd, funfam, conc, r["MODEL_SENSITIVE"]])
    C.write("04_camkk_subfamilies/CAMKK_constituent_calls.tsv",
            ["SPECIES", "PROTEIN_ID", "R1_FINAL_CLASS", "DECISION_SOURCE", "ROLE_IN_TREE", "VALIDATION_TYPE",
             "CAMKK_L1_CLASS", "CAMKK_L2_CLASS", "MASKED_L1_CALL", "MASKED_L1_SUPPORT", "MASKED_L1_CONFIDENCE",
             "MASKED_L2_CALL", "UNMASKED_L1_CALL", "UNMASKED_L1_SUPPORT", "UNMASKED_L1_CONFIDENCE",
             "UNMASKED_L2_CALL", "ALIGNMENT_SENSITIVE_L1", "ALIGNMENT_SENSITIVE_L2", "FROZEN_TREE_CALL",
             "FROZEN_TREE_CONTINUITY_CONFLICT",
             "REFERENCE_IDS_IN_ASSIGNING_SPLIT", "NEGATIVE_IDS_IN_ASSIGNING_SPLIT", "CALL_BASIS",
             "CATALYTIC_MOTIF_STATUS", "PROSITE_SIGNATURES", "DOMAIN_LENGTH", "DOMAIN_ARCHITECTURE",
             "CDD_ANNOTATION", "FUNFAM_ANNOTATION", "ANNOTATION_CONCORDANCE", "MODEL_SENSITIVE"], out,
            script=SELF,
            notes=["PRIMARY_CAMKK_TREE = REFERENCE_OCCUPANCY_MASKED (%s columns, median gap %s);" % _qc_desc()[0],
                   "SENSITIVITY_CAMKK_TREE = UNMASKED (%s columns, median gap %s, %s per cent reference" % _qc_desc()[1][:3],
                   "occupancy, which %s the trigger declared before any topology was read). CAMKK_L1_CLASS" % _qc_desc()[1][3],
                   "is the masked call unless the two trees give supported conflicting placements, in which",
                   "case it reads ALIGNMENT_SENSITIVE and the copy contributes to no constituent-specific",
                   "evolutionary claim. Level-2 disagreement resolves to L2_UNRESOLVED rather than forcing a",
                   "subtype, because the Level-2 anchor depth is shallow (metazoan 2, fungal 3, plant 2).",
                   "None of this changes the frozen R1 CAMKK presence or absence call: constituent labels say",
                   "what kind of upstream kinase is retained, never whether one is retained.",
                   "A frozen-tree continuity conflict (a confident constituent call on a copy the frozen",
                   "CAMKK tree grouped with boundary markers only) demotes that copy to SUBFAMILY_UNRESOLVED.",
                   "It makes the species unresolved only if no other retained copy supplies a nonconflicting",
                   "resolved assignment; a difficult paralog never discards a concordant sibling copy.",
                   "CDD and FunFam are CORROBORATION_ONLY and sit beside the decision, never inside it.",
                   "The evidence stack is phylogenetic placement plus independently computed",
                   "sequence-signature and domain corroboration, not independent biological validation."])

    by_sp, clade = {}, C.clade_of()
    for r in out:
        by_sp.setdefault(r[0], []).append(r)
    sp_rows = []
    for sp, g in sorted(by_sp.items()):
        usable = [r for r in g if r[6] not in ("SUBFAMILY_UNRESOLVED", "ALIGNMENT_SENSITIVE")]
        unres = [r for r in g if r[6] == "SUBFAMILY_UNRESOLVED"]
        alsen = [r for r in g if r[6] == "ALIGNMENT_SENSITIVE"]
        l1s = sorted({r[6] for r in usable})
        if not usable:
            cls, basis = "SUBFAMILY_UNRESOLVED", ("all copies alignment-sensitive" if alsen
                                                  else "no copy resolved")
        elif len(l1s) > 1:
            cls, basis = "MULTIFAMILY", "resolved copies in more than one constituent family"
        else:
            cls, basis = l1s[0], "all resolved copies agree"
        l2s = sorted({r[7] for r in usable if r[7] in NAMED_L2})
        l2 = l2s[0] if len(l2s) == 1 else ("MULTIFAMILY_L2" if len(l2s) > 1 else
             ("SNF1_ACTIVATOR_SUBTYPE_UNRESOLVED" if cls == "SNF1_SNRK1_ACTIVATOR_LIKE" else "-"))
        ccf = [r for r in g if r[19] == "yes"]   # index 19 is `FROZEN_TREE_CONTINUITY_CONFLICT`
        sp_rows.append([sp, clade.get(sp, "UNKNOWN"), len(g), len(usable), len(unres), len(alsen), len(ccf),
                        cls, l2, "yes" if unres and usable else "no", "yes" if alsen else "no",
                        "yes" if ccf else "no", basis])
    C.write("04_camkk_subfamilies/CAMKK_constituent_species.tsv",
            ["SPECIES", "CLADE", "N_COPIES", "N_RESOLVED_COPIES", "N_UNRESOLVED_COPIES",
             "N_ALIGNMENT_SENSITIVE_COPIES", "N_CONTINUITY_CONFLICT_COPIES", "CAMKK_L1_CLASS",
             "CAMKK_L2_CLASS", "UNRESOLVED_COPY_PRESENT", "ALIGNMENT_SENSITIVE_COPY_PRESENT",
             "CONTINUITY_CONFLICT_COPY_PRESENT", "CALL_BASIS"], sp_rows, script=SELF,
            notes=["Frozen duplicate rule. Copies are never counted as independent species. An",
                   "alignment-sensitive copy is excluded from the species constituent call and recorded."])

    cand = [r for r in out if r[4] == "CANDIDATE"]
    loo = [r for r in out if r[4] != "CANDIDATE"]
    n = lambda rs, i, v: sum(1 for r in rs if r[i] == v)
    sp_cls = lambda v: sum(1 for r in sp_rows if r[7] == v)      # index 7 is the species first-level class
    sp_l2 = lambda v: sum(1 for r in sp_rows if r[8] == v)       # index 8 is the species second-level subtype
    ident = sum(1 for r in cand if r[8] == r[12])
    ## A reference counts as recovered only when its leave-one-out call equals its declared class, not merely
    ## when it is resolved.
    expect_l1 = {acc: l1 for _rid, _sp, acc, _g, l1, _l2 in C.CAMKK_REFS}
    expect_l2 = {acc: l2 for _rid, _sp, acc, _g, _l1, l2 in C.CAMKK_REFS}
    rec = [(r[1], expect_l1.get(r[1], "-"), r[6], expect_l2.get(r[1], "-"), r[7], r[9], r[10])
           for r in loo]
    n_rec = sum(1 for a in rec if a[1] == a[2])
    ufb = {r["TREE_ROLE"]: r["UFBOOT_CONVERGENCE"] for r in
           C.rows(os.path.join(C.X3, "02_inputs/tree_qc.tsv"), required=False)}
    uf = ufb.get("CAMKK_MASKED_PRIMARY", "NOT_RECORDED")
    cls_rows = []
    for l1, anchors in (("CAMKK_METAZOAN_LIKE", 2), ("SNF1_SNRK1_ACTIVATOR_LIKE", 5)):
        mine = [a for a in rec if a[1] == l1]
        ok = sum(1 for a in mine if a[2] == l1)
        status = ("REFERENCE_LOO_SUPPORTED" if mine and ok == len(mine)
                  else "ANCHOR_LIMITED_PROVISIONAL")
        cls_rows.append([l1, anchors, len(mine), ok, len(mine) - ok, status, uf,
                         n(out, 6, l1), sp_cls(l1),
                         "every declared reference of this class returned to it under fixed-tree"
                         " reference-label leave-one-out" if status == "REFERENCE_LOO_SUPPORTED" else
                         "at least one declared reference of this class did not return to it; with only %d"
                         " anchors, holding one out leaves %d same-class anchor(s), so constituent calls in"
                         " this class are provisional and must not be presented as equal in weight to a"
                         " reference-LOO-supported class" % (anchors, anchors - 1)])
    for l2, anchors in (("SAK1_TOS3_ELM1_LIKE", 3), ("GRIK_SNAK_LIKE", 2)):
        mine = [a for a in rec if a[3] == l2]
        ok = sum(1 for a in mine if a[4] == l2)
        cls_rows.append([l2, anchors, len(mine), ok, len(mine) - ok,
                         "REFERENCE_LOO_SUPPORTED" if mine and ok == len(mine)
                         else "ANCHOR_LIMITED_PROVISIONAL", uf, n(out, 7, l2), sp_l2(l2),
                         "Level 2; shallow anchor depth means a successful leave-one-out is useful but is"
                         " not deep Level-2 validation"])
    C.write("04_camkk_subfamilies/CAMKK_constituent_class_validation.tsv",
            ["CONSTITUENT_CLASS", "N_DECLARED_ANCHORS", "N_REFERENCES_TESTED", "N_RETURNED_TO_OWN_CLASS",
             "N_FAILED_TO_RETURN", "CLASS_VALIDATION_STATUS", "PRIMARY_TREE_UFBOOT_CONVERGENCE",
             "N_COPIES_CALLED", "N_SPECIES_CALLED", "BASIS"], cls_rows, script=SELF,
            notes=["Class-level validation, machine-readable, so a constituent annotation can never be",
                   "presented as stronger than its anchor set supports. Recovery is an exact match to the",
                   "reference's own declared class, not merely a non-unresolved call. The procedure is",
                   "reference-label leave-one-out on a FIXED inferred tree: the held-out sequence remains in",
                   "the alignment and topology and only its reference role is withheld, so this is",
                   "self-consistency and not full phylogenetic leave-one-out. The primary tree's ultrafast",
                   "bootstrap convergence state is carried here for the same reason."])
    C.write("04_camkk_subfamilies/CAMKK_reference_loo_detail.tsv",
            ["PROTEIN_ID", "EXPECTED_L1", "OBSERVED_L1", "EXPECTED_L2", "OBSERVED_L2", "MASKED_L1_SUPPORT",
             "MASKED_L1_CONFIDENCE"], [list(a) for a in rec], script=SELF,
            notes=["Per-reference fixed-tree leave-one-out outcome against the declared expected class."])
    C.write("04_camkk_subfamilies/CAMKK_constituent_counts.tsv", ["QUANTITY", "COPIES", "SPECIES", "NOTE"],
            [["N_CAMKK_ORTHOLOG_COPIES", len(out), len(sp_rows), "references included"],
             ["N_CANDIDATES", len(cand), "-", "copies that are not declared references"],
             ["N_REFERENCE_LOO", len(loo), "-", "declared references, leave-one-out self-consistency"],
             ["CAMKK_METAZOAN_LIKE", n(out, 6, "CAMKK_METAZOAN_LIKE"), sp_cls("CAMKK_METAZOAN_LIKE"),
              "copies over all 194, candidates plus declared references, so COPIES and SPECIES share a"
              " denominator"],
             ["SNF1_SNRK1_ACTIVATOR_LIKE", n(out, 6, "SNF1_SNRK1_ACTIVATOR_LIKE"),
              sp_cls("SNF1_SNRK1_ACTIVATOR_LIKE"), "includes the 5 declared Snf1-activator references"],
             ["OTHER_CAMKK_CATEGORY", n(out, 6, "OTHER_CAMKK_CATEGORY"), sp_cls("OTHER_CAMKK_CATEGORY"), "-"],
             ["SUBFAMILY_UNRESOLVED", n(out, 6, "SUBFAMILY_UNRESOLVED"), sp_cls("SUBFAMILY_UNRESOLVED"),
              "expected to be substantial and not optimised downward"],
             ["CAMKK_METAZOAN_LIKE_CANDIDATES_ONLY", n(cand, 6, "CAMKK_METAZOAN_LIKE"), "-",
              "excluding the declared references, for comparison"],
             ["SNF1_SNRK1_ACTIVATOR_LIKE_CANDIDATES_ONLY", n(cand, 6, "SNF1_SNRK1_ACTIVATOR_LIKE"), "-",
              "excluding the declared references, for comparison"],
             ["ALIGNMENT_SENSITIVE_L1", n(cand, 16, "yes"), sum(1 for r in sp_rows if r[10] == "yes"),
              "supported conflicting placement between the masked and unmasked trees"],
             ["FROZEN_TREE_CONTINUITY_CONFLICT", n(cand, 19, "yes"),
              sum(1 for r in sp_rows if r[11] == "yes"),
              "confident call on a copy the frozen tree grouped with boundary markers only; copy demoted"],
             ["L1_SUPPORTED", n(cand, 10, "SUPPORTED"), "-", "masked tree UFBoot >= 95"],
             ["L1_PROVISIONAL", n(cand, 10, "PROVISIONAL_PLACEMENT"), "-", "masked tree UFBoot 70 to 94.999"],
             ["L1_UNRESOLVED", n(cand, 10, "UNRESOLVED"), "-", "below 70, mixed, or negative in the split"],
             ["L2_SAK1_TOS3_ELM1_LIKE", n(out, 7, "SAK1_TOS3_ELM1_LIKE"), sp_l2("SAK1_TOS3_ELM1_LIKE"),
              "copies over all 194; includes SAK1, TOS3 and ELM1 themselves"],
             ["L2_GRIK_SNAK_LIKE", n(out, 7, "GRIK_SNAK_LIKE"), sp_l2("GRIK_SNAK_LIKE"),
              "copies over all 194; includes GRIK1 and GRIK2 themselves"],
             ["L2_UNRESOLVED", sum(1 for r in out if r[7] not in NAMED_L2), "-",
              "a high Level-2 unresolved fraction is an expected and valid result"],
             ["MASKED_UNMASKED_L1_IDENTICAL", ident, "-",
              "candidates whose Level-1 call is the same string under both trees"],
             ["MASKED_UNMASKED_L1_DIFFER", len(cand) - ident, "-",
              "any difference, including one tree resolving where the other does not"],
             ["DIFFER_MASKED_NAMED_UNMASKED_UNRESOLVED",
              sum(1 for r in cand if r[8] != "SUBFAMILY_UNRESOLVED" and r[12] == "SUBFAMILY_UNRESOLVED"),
              "-", "primary resolves, sensitivity does not; not a contradiction"],
             ["DIFFER_UNMASKED_NAMED_MASKED_UNRESOLVED",
              sum(1 for r in cand if r[8] == "SUBFAMILY_UNRESOLVED" and r[12] != "SUBFAMILY_UNRESOLVED"),
              "-", "sensitivity resolves, primary does not; the primary call stands and is unresolved"],
             ["DIFFER_BOTH_NAMED_AND_CONTRADICTORY",
              sum(1 for r in cand if r[8] != "SUBFAMILY_UNRESOLVED" and r[12] not in
                  ("SUBFAMILY_UNRESOLVED", "NOT_AVAILABLE") and r[8] != r[12]), "-",
              "the only kind of difference that counts as alignment sensitivity"],
             ["REFERENCE_LOO_RECOVERED_TO_OWN_CLASS", n_rec, "-",
              "declared references returning to their own constituent class under leave-one-out on the "
              "primary tree; the metazoan arm has only two anchors, so holding one out leaves a single "
              "same-class anchor and recovery depends on the support of that one sister split"],
             ["L2_ANCHOR_DEPTH", "metazoan 2, fungal 3, plant 2", "-",
              "holding out one of a two-reference class leaves a single same-class anchor, so a successful "
              "leave-one-out is useful but is not deep Level-2 validation"]], script=SELF)
    C.record_outputs("r3_CAMKK", SELF)
    print("CAMKK: %d copies (%d candidates, %d LOO refs), %d species; alignment-sensitive %d, "
          "continuity-conflict %d" % (len(out), len(cand), len(loo), len(sp_rows), n(cand, 16, "yes"),
                                      n(cand, 19, "yes")))


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
