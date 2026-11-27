#!/usr/bin/env python3
# >>> deposit header
# Purpose: Combines the LKB1 tree placements with ProSite kinase signatures and domain length into per-copy
#          sequence-support statuses, and checks that every declared positive returns to LKB1 under leave-one-out.
# Methods: Orthology calls
# Inputs:  `10_validation/LKB1_kinase_domain_manifest.tsv`, `LKB1_tree_calls_masked.tsv`,
#          `LKB1_kinase_domain.masked.aln.faa`, `02_inputs/tree_qc.tsv`, the ortholog table `T1_ortholog_table.tsv`,
#          the per-species domain scans `<species>.csv.gz`
# Outputs: `10_validation/LKB1_sequence_validation.tsv`, `10_validation/LKB1_validation_counts.tsv` and, when a
#          copy needs review, `10_validation/R2_CANDIDATE_CHANGE_REQUIRED.tsv`
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
# Ran as: rbx3_calls.py
## Turns the tree-based assignments of rbx3_tree_calls.R into per-copy LKB1 sequence-support statuses (strongly
## supported, compatible or unresolved), which combine the tree placement with the ProSite kinase signatures and
## the kinase-domain length. The script refuses CAMKK, which rbx3_camkk_calls.py handles, so its CAMKK branch below
## is never reached. Domain annotations (CDD subfamily models, CATH FunFam names) are reported beside the call for
## comparison only and never decide it: those names were themselves assigned by homology, so letting them decide
## would be circular. The orthology call of the ortholog table (FINAL_CLASS) is carried through unchanged; a copy
## left unresolved in a split that holds negative-family sequences is listed in `R2_CANDIDATE_CHANGE_REQUIRED.tsv`
## for review, and the publication tables are not modified.
import gzip, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "2_candidate_discovery_and_orthology", "upstream_kinases", "lib"))
import upstream_kinase_tables as C  # noqa: E402

SELF = os.path.abspath(__file__)
TARGET = sys.argv[1] if len(sys.argv) > 1 else "LKB1"
if TARGET == "CAMKK":
    sys.exit("REFUSED: CAMKK is written only by rbx3_camkk_calls.py from the masked-primary tree. This "
             "generic writer would emit the same filenames from a single unmasked tree and could silently "
             "replace the authoritative primary products. Use r3b_camkk.sbatch.")
SUB = "04_camkk_subfamilies" if TARGET == "CAMKK" else "10_validation"
TREETAG = sys.argv[2] if len(sys.argv) > 2 else ""
## Integrity of the catalytic motifs is taken from the per-protein ProSite matches of the domain scans, not from
## regular expressions written here. PS00107 is the protein-kinase ATP-binding region signature and PS00108 the
## Ser/Thr kinase active-site signature. Fixed motif patterns would fail here: LKB1 is an atypical kinase whose
## beta-3 lysine sits in AVK rather than V..K, and human STK11 matches none of VAIK, HRD or DFG literally, so a
## pattern panel would call the whole family divergent, including the human reference itself. The check in main()
## stops the run if the human reference does not itself carry both signatures.
PROSITE = ["PS00107", "PS00108"]
## CDD subfamily models exist only for the metazoan CAMKKs; SAK1, TOS3, GRIK1 and GRIK2 fall to the shared
## cd14008 model, and ELM1 has no specific CDD model and no FunFam entry. Annotation can therefore support a
## metazoan placement but almost never a second-level subtype; this asymmetry is reported, not corrected for.
CDD_L1 = {"cd14200": "CAMKK_METAZOAN_LIKE", "cd14199": "CAMKK_METAZOAN_LIKE",
          "cd14118": "CAMKK_METAZOAN_LIKE", "cd14008": "SHARED_LKB1_CAMKK_NODE_NOT_DIAGNOSTIC"}


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
    status = ("CANONICAL_BOTH" if len(ps) == 2 else
              ("PARTIAL_" + ps[0] if len(ps) == 1 else "NONE_DETECTED"))
    return (";".join(sorted(set(cdd))) or "-", ";".join(sorted(set(funfam))[:2]) or "-",
            ";".join(ps) or "-", status)


def ref_columns(aln, ref_tip):
    s, cols, pos = aln[ref_tip], {}, 0
    for i, ch in enumerate(s):
        if ch != "-":
            pos += 1
            cols[pos] = i
    return cols, s.replace("-", "")


def _src_desc(out):
    """Counts the non-reference copies by the kind of evidence behind their orthology call, for the table note."""
    n = {}
    for x in out:
        if x[4] == "CANDIDATE":
            n[x[3]] = n.get(x[3], 0) + 1
    return ", ".join("%s %d" % (k, v) for k, v in sorted(n.items())) or "none"


def main():
    C.verify_freezes()
    man = C.rows(os.path.join(C.X3, SUB, "%s_kinase_domain_manifest.tsv" % TARGET))
    ## The LKB1 calls rest on the reference-masked tree, so the driver passes "masked". The file name is chosen
    ## explicitly rather than by default, so a run cannot quietly read whichever call table happens to be on disk.
    ## Without "masked" the script reads `<target>_tree_calls.tsv`, a name the current drivers do not write.
    cf = "%s_tree_calls_masked.tsv" % TARGET if TREETAG == "masked" else "%s_tree_calls.tsv" % TARGET
    calls = {r["TIP"]: r for r in C.rows(os.path.join(C.X3, SUB, cf))}
    alnf = "%s_kinase_domain%s.aln.faa" % (TARGET, ".masked" if TREETAG == "masked" else "")
    aln = C.read_fasta(os.path.join(C.X3, SUB, alnf))
    by_tip = {"%s__%s" % (r["SPECIES"], r["PROTEIN_ID"].replace(".", "_")): r for r in man}
    anchor_acc = "NP_757343.2" if TARGET == "CAMKK" else "NP_000446.1"
    anchor_status = annotations("Homo_sapiens", anchor_acc)[3]
    if anchor_status != "CANONICAL_BOTH":
        sys.exit("REFUSED: anchor %s does not itself carry both ProSite kinase signatures (%s); a motif "
                 "instrument that fails on its own anchor would mark every copy divergent"
                 % (anchor_acc, anchor_status))
    t1 = {(r["TARGET"], r["PROTEIN_ID"]): r for r in C.rows(os.path.join(C.PUB, "T1_ortholog_table_all20.tsv"))}
    ## The strongly supported versus compatible distinction depends partly on ultrafast bootstrap support, so the
    ## bootstrap convergence state of the tree is written beside every result.
    ufb = {r["TREE_ROLE"]: r["UFBOOT_CONVERGENCE"] for r in
           C.rows(os.path.join(C.X3, "02_inputs/tree_qc.tsv"), required=False)}
    UFC = ufb.get(("LKB1_MASKED_PRIMARY" if TREETAG == "masked" else "LKB1_UNMASKED_SENSITIVITY")
                  if TARGET == "LKB1" else "CAMKK_MASKED_PRIMARY", "NOT_RECORDED")
    lens = [int(r["DOMAIN_LENGTH"]) for r in man if r["DOMAIN_LENGTH"] not in ("NA", "")
            and "REFERENCE" in r["ROLE"]]
    lo_ok, hi_ok = (min(lens) * 0.6, max(lens) * 1.6) if lens else (0, 1e9)

    out = []
    for tip, r in sorted(by_tip.items()):
        if r["ROLE"] == "NEGATIVE":
            continue
        c = C.tip_in(calls, tip)
        if c is None:
            sys.exit("REFUSED: no tree call for %s; every copy must be accounted for" % tip)
        dl = int(r["DOMAIN_LENGTH"]) if r["DOMAIN_LENGTH"] not in ("NA", "") else None
        arch_ok = dl is not None and lo_ok <= dl <= hi_ok
        cdd, funfam, ps, motif_status = annotations(r["SPECIES"], r["PROTEIN_ID"])
        ann_l1 = {CDD_L1[x] for x in cdd.split(";") if x in CDD_L1} - {"SHARED_LKB1_CAMKK_NODE_NOT_DIAGNOSTIC"}
        conc = ("NOT_DIAGNOSTIC" if not ann_l1 else
                "AGREES" if c["CLASS_L1"] in ann_l1 else "DISAGREES")
        out.append([r["SPECIES"], r["PROTEIN_ID"], t1.get((TARGET, r["PROTEIN_ID"]), {}).get("FINAL_CLASS", "-"),
                    r["DECISION_SOURCE"], c["ROLE"], c["VALIDATION_TYPE"], c["CLASS_L1"], c["CLASS_L2"],
                    c["TREE_SUPPORT"], c["CALL_CONFIDENCE"], c["REFERENCE_IDS"], c["NEGATIVE_TIPS"],
                    c["CALL_BASIS"], motif_status, ps,
                    dl if dl else "NA", "compatible" if arch_ok else "atypical", cdd, funfam, conc,
                    r["DOMAIN_SOURCE"], r["MODEL_SENSITIVE"]])

    HDR = ["SPECIES", "PROTEIN_ID", "R1_FINAL_CLASS", "DECISION_SOURCE", "ROLE_IN_TREE", "VALIDATION_TYPE",
           "TREE_CONSTITUENT_CALL", "CLASS_L2", "TREE_SUPPORT", "CALL_CONFIDENCE",
           "REFERENCE_IDS_IN_ASSIGNING_SPLIT", "NEGATIVE_IDS_IN_ASSIGNING_SPLIT", "CALL_BASIS",
           "CATALYTIC_MOTIF_STATUS", "PROSITE_SIGNATURES", "DOMAIN_LENGTH", "DOMAIN_ARCHITECTURE", "CDD_ANNOTATION",
           "FUNFAM_ANNOTATION", "ANNOTATION_CONCORDANCE", "DOMAIN_SOURCE", "MODEL_SENSITIVE"]
    ANN_NOTE = ["CDD_ANNOTATION and FUNFAM_ANNOTATION are CORROBORATION_ONLY and sit beside the decision,",
                "never inside it: names such as TOS3 or GRIK2 are themselves transferred from homology, so",
                "letting them assign a family would manufacture its own confirmation. CDD subfamily models",
                "exist only for the metazoan arm, so ANNOTATION_CONCORDANCE reads NOT_DIAGNOSTIC for most",
                "Snf1-activator copies; that asymmetry is recorded rather than worked around."]

    if TARGET == "CAMKK":
        C.write("04_camkk_subfamilies/CAMKK_constituent_calls.tsv", HDR, out, script=SELF,
                notes=["A constituent label is a post-freeze annotation and never an orthology vote: R1",
                       "FINAL_CLASS is carried through unchanged. Taxonomy never assigns a family."] + ANN_NOTE)
        by_sp = {}
        for r in out:
            by_sp.setdefault(r[0], []).append(r)
        sp_rows = []
        clade = C.clade_of()
        for sp, g in sorted(by_sp.items()):
            res = [r for r in g if r[6] != "SUBFAMILY_UNRESOLVED"]
            unres = [r for r in g if r[6] == "SUBFAMILY_UNRESOLVED"]
            l1s = sorted({r[6] for r in res})
            if not res:
                cls, basis = "SUBFAMILY_UNRESOLVED", "no copy resolved"
            elif len(l1s) > 1:
                cls, basis = "MULTIFAMILY", "resolved copies in more than one constituent family"
            else:
                cls, basis = l1s[0], "all resolved copies agree"
            l2s = sorted({r[7] for r in res if r[7] not in ("-", "", "SNF1_ACTIVATOR_SUBTYPE_UNRESOLVED")})
            l2 = l2s[0] if len(l2s) == 1 else ("MULTIFAMILY_L2" if len(l2s) > 1 else
                 ("SNF1_ACTIVATOR_SUBTYPE_UNRESOLVED" if cls == "SNF1_SNRK1_ACTIVATOR_LIKE" else "-"))
            sp_rows.append([sp, clade.get(sp, "UNKNOWN"), len(g), len(res), len(unres), cls, l2,
                            "yes" if unres and res else "no",
                            ";".join(sorted({r[9] for r in res})) or "-", basis])
        C.write("04_camkk_subfamilies/CAMKK_constituent_species.tsv",
                ["SPECIES", "CLADE", "N_COPIES", "N_RESOLVED_COPIES", "N_UNRESOLVED_COPIES", "CAMKK_L1_CLASS",
                 "CAMKK_L2_CLASS", "UNRESOLVED_COPY_PRESENT", "COPY_CONFIDENCES", "CALL_BASIS"], sp_rows,
                script=SELF,
                notes=["Frozen duplicate rule: all resolved copies agreeing gives that class; resolved copies",
                       "in different classes gives MULTIFAMILY; resolved plus unresolved keeps the resolved",
                       "class with UNRESOLVED_COPY_PRESENT=yes; only unresolved gives SUBFAMILY_UNRESOLVED.",
                       "No copy is ever selected because it gives a cleaner evolutionary pattern."])
        cand = [r for r in out if r[4] == "CANDIDATE"]
        loo = [r for r in out if r[4] != "CANDIDATE"]
        cnt = lambda rs, i, v: sum(1 for r in rs if r[i] == v)
        C.write("04_camkk_subfamilies/CAMKK_constituent_counts.tsv",
                ["QUANTITY", "VALUE", "NOTE"],
                [["N_CAMKK_ORTHOLOG_COPIES", len(out), "all frozen ORTHOLOG copies, references included"],
                 ["N_CANDIDATES", len(cand), "copies that are not declared references"],
                 ["N_REFERENCE_LOO", len(loo), "declared references, classified leave-one-out"],
                 ["N_L1_SUPPORTED", cnt(cand, 9, "SUPPORTED"), "UFBoot >= 95, single family, no negative"],
                 ["N_L1_PROVISIONAL", cnt(cand, 9, "PROVISIONAL_PLACEMENT"), "UFBoot 70 to 94.999"],
                 ["N_L1_UNRESOLVED", cnt(cand, 9, "UNRESOLVED"), "below 70, mixed, or negative in the split"],
                 ["N_L2_SUPPORTED", sum(1 for r in cand if r[7] in ("SAK1_TOS3_ELM1_LIKE", "GRIK_SNAK_LIKE")
                                        and r[9] == "SUPPORTED"), "a named Level-2 subtype at UFBoot >= 95"],
                 ["N_L2_PROVISIONAL", sum(1 for r in cand if r[7] in ("SAK1_TOS3_ELM1_LIKE", "GRIK_SNAK_LIKE")
                                          and r[9] == "PROVISIONAL_PLACEMENT"), "named subtype, 70 to 94.999"],
                 ["N_L2_UNRESOLVED", sum(1 for r in cand if r[7] not in ("SAK1_TOS3_ELM1_LIKE",
                                                                         "GRIK_SNAK_LIKE")), "no named subtype"],
                 ["L2_ANCHOR_DEPTH", "metazoan 2, fungal 3, plant 2",
                  "holding out one of a two-reference class leaves a single same-class anchor, so a successful "
                  "leave-one-out is useful but is not deep Level-2 validation; a high Level-2 unresolved "
                  "fraction is an expected and valid result"]], script=SELF)
    else:
        val, flags = [], []
        for r in out:
            l1, conf, negs = r[6], r[9], r[11]
            placed = l1 == "LKB1_REFERENCE"
            if placed and conf == "SUPPORTED" and negs == "-" and r[13] == "CANONICAL_BOTH" and r[16] == "compatible":
                st = "LKB1_STRONGLY_SUPPORTED"
            elif (placed and negs == "-" and conf != "UNRESOLVED"
                  and r[13] == "CANONICAL_BOTH" and r[16] == "compatible"):
                ## Compatible is a combined status, not only a bootstrap band: it requires placement with the LKB1
                ## references at support of at least 70, no negative-family sequence in the assigning split, both
                ## ProSite kinase signatures and a compatible domain length. Without the signature condition,
                ## copies would enter sensitivity tier A on placement alone.
                st = "LKB1_COMPATIBLE"
            else:
                st = "LKB1_SEQUENCE_UNRESOLVED"
            basis = "%s; support %s; kinase signatures %s; architecture %s" % (r[12], r[8], r[13], r[16])
            val.append([r[0], r[1], r[2], r[3], r[4], r[5], st, r[8], r[9], UFC, r[10], negs, r[13], r[14],
                        r[15], r[16], r[17], r[18], r[19], basis])
            if negs != "-" and st == "LKB1_SEQUENCE_UNRESOLVED":
                flags.append([r[0], r[1], r[3], "placed in a supported split containing negative-family "
                              "references", "R2_CANDIDATE_CHANGE_REQUIRED_CANDIDATE", "frozen R1 unchanged"])
        C.write("10_validation/LKB1_sequence_validation.tsv",
                ["SPECIES", "PROTEIN_ID", "R1_FINAL_CLASS", "DECISION_SOURCE", "ROLE_IN_TREE",
                 "VALIDATION_TYPE", "LKB1_SEQUENCE_SUPPORT", "TREE_SUPPORT", "PLACEMENT_CONFIDENCE",
                 "TREE_UFBOOT_CONVERGENCE", "REFERENCE_IDS_IN_ASSIGNING_SPLIT", "NEGATIVE_IDS_IN_ASSIGNING_SPLIT", "CATALYTIC_MOTIF_STATUS",
                 "PROSITE_SIGNATURES", "DOMAIN_LENGTH", "DOMAIN_ARCHITECTURE", "CDD_ANNOTATION",
                 "FUNFAM_ANNOTATION", "ANNOTATION_CONCORDANCE", "LKB1_SEQUENCE_SUPPORT_BASIS"], val,
                script=SELF,
                notes=["Post-freeze annotation of the frozen LKB1 ORTHOLOG copies. LKB1_COMPATIBLE is an",
                       "integrated status, not merely bootstrap 70 to 94: it requires consistent LKB1",
                       "placement, no conflicting negative-family placement, and compatible motifs and domain",
                       "architecture. No biological recoding happens here and there is no NON_LKB1 status.",
                       "Context (counts read from this run): the retained non-reference LKB1 copies were",
                       "admitted by %s. The %d" % (_src_desc(out), sum(1 for x in out if x[4] != "CANDIDATE")),
                       "declared positives establish the family boundary used to judge the other %d" % sum(1 for x in out if x[4] == "CANDIDATE"),
                       "proteins, so their leave-one-out result is self-consistency, not independent",
                       "validation. It is reference-label leave-one-out on a FIXED inferred tree: the held-out",
                       "sequence stays in the alignment and the topology, and only its reference role is",
                       "withheld during assignment. It is not full phylogenetic leave-one-out.",
                       "validation. The evidence stack is phylogenetic placement plus independently computed",
                       "sequence-signature and domain corroboration; the ProSite signatures are computed from",
                       "the same protein sequence and are therefore not independent biological validation."] + ANN_NOTE)
        if flags:
            C.write("10_validation/R2_CANDIDATE_CHANGE_REQUIRED.tsv",
                    ["SPECIES", "PROTEIN_ID", "DECISION_SOURCE", "OBSERVATION", "DISPOSITION", "R1_EFFECT"],
                    flags, script=SELF)
        loo = [v for v in val if v[4] != "CANDIDATE"]
        r66 = [v for v in val if v[4] == "CANDIDATE"]
        cnt = lambda rs, v: sum(1 for r in rs if r[6] == v)
        ## The check, fixed in advance, is that every declared positive control, held out in turn, is placed back
        ## with the LKB1 references (strongly supported or compatible, with no negative-family sequence in its
        ## split); this shows that the family boundary is self-consistent. It does not require ultrafast bootstrap
        ## of at least 95: several declared positives are the only representative of their lineage, so holding one
        ## out leaves few references of the same class and its support falls into the provisional band. That
        ## effect is reported beside the check.
        ## Index 11 of each row is NEGATIVE_IDS_IN_ASSIGNING_SPLIT (index 9 holds TREE_UFBOOT_CONVERGENCE).
        returned = [v for v in loo if v[6] != "LKB1_SEQUENCE_UNRESOLVED" and v[11] == "-"]
        gate = "PASS" if len(returned) == len(loo) else "FAIL"
        n_strong_loo = cnt(loo, "LKB1_STRONGLY_SUPPORTED")
        C.write("10_validation/LKB1_validation_counts.tsv",
                ["GROUP", "N", "LKB1_STRONGLY_SUPPORTED", "LKB1_COMPATIBLE", "LKB1_SEQUENCE_UNRESOLVED",
                 "TREE_UFBOOT_CONVERGENCE", "GATE"],
                [["DECLARED_POSITIVE_LEAVE_ONE_OUT", len(loo), cnt(loo, "LKB1_STRONGLY_SUPPORTED"),
                  cnt(loo, "LKB1_COMPATIBLE"), cnt(loo, "LKB1_SEQUENCE_UNRESOLVED"), UFC,
                  "ALL_RETURN_TO_LKB1_GROUP=%s (%d of %d also reach UFBoot>=95)"
                  % (gate, n_strong_loo, len(loo))],
                 ["R66_RESOLVED_CANDIDATES", len(r66), cnt(r66, "LKB1_STRONGLY_SUPPORTED"),
                  cnt(r66, "LKB1_COMPATIBLE"), cnt(r66, "LKB1_SEQUENCE_UNRESOLVED"), UFC, "-"]], script=SELF,
                notes=["TREE_UFBOOT_CONVERGENCE qualifies every support value here (this LKB1 tree: %s), so the" % UFC,
                       "STRONG versus COMPATIBLE split, which is",
                       "defined partly by UFBoot support, must never be reported without that limitation.",
                       "The gate: every declared positive must return to the LKB1 reference group under",
                       "leave-one-out, or candidate classification stops. That is self-consistency of the",
                       "family boundary, not independent validation of it: these %d define the boundary being" % len(loo),
                       "used to judge the other %d proteins. Reaching UFBoot >= 95 under leave-one-out is" % len(r66),
                       "reported separately and is not the gate, because a declared positive that is the sole",
                       "representative of its lineage loses its same-class anchors when held out and falls",
                       "into the provisional band for that reason alone."])
        print("LKB1 declared-positive LOO gate (return to LKB1 group): %s; %d of %d also reach UFBoot>=95"
              % (gate, n_strong_loo, len(loo)))
        if gate != "PASS":
            sys.exit("REFUSED: not every declared LKB1 positive returns to the LKB1 group under "
                     "leave-one-out; candidate classification stops")
    C.record_outputs("r3_%s" % TARGET, SELF)
    print("%s: %d copies annotated" % (TARGET, len(out)))


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
