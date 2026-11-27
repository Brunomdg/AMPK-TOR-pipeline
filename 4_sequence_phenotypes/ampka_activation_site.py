#!/usr/bin/env python3
# >>> deposit header
# Purpose: Codes the upstream-kinase repertoire (LKB1 and the CAMKK/SAK1/GRIK category) of every species with
#          AMPKα, and reads the residue at the AMPKα activation-loop site (the Thr172 site) in every AMPKα copy,
#          with the identity of the surrounding window and of the kinase domain, per copy and per species.
# Methods: Phylogenetic regressions
# Inputs:  `AMPKa_activation_reference.tsv` (from job y0); the ten-target three-state table; the AMPKα alignment
#          `AMPKa_CTD_alignment.faa` and `RBX1_sequence_manifest.AMPKa.tsv` of the first chain
# Outputs: `04_ampka_activation/upstream_kinase_states.tsv`, `AMPKa_activation_copy.tsv` and
#          `AMPKa_activation_species.tsv`; `11_validation/AMPKa_toxoplasma_activation_check.tsv` and
#          `AMPKa_activation_gates.tsv`
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
## Job y3 of the second chain. Upstream-kinase repertoire states and the AMPKα activation-loop site, scored through
## the AMPKα alignment of the first chain against the motif-anchored reference coordinates written by job y0.
## "Thr172" is only the conventional name; the sequence index comes from the reference table. NEITHER_CAPTURED
## means that neither LKB1 nor the CAMKK/SAK1/GRIK category (as defined for this project) was detected. It never
## means that the species has no AMPK kinase.
# deposit path defaults: begin
_TA210_SCRATCH = __import__("os").environ.get("TA210_SCRATCH") or "/gpfs2/scratch/bdigenov"
# deposit path defaults: end
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phenotype_tables_chain2 as C  # noqa: E402
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "4_sequence_phenotypes", "lib"))
import sequence_scoring as S  # noqa: E402

SELF = os.path.abspath(__file__)
ALN = os.path.join(C.X1, "04_ampka_ctd/AMPKa_CTD_alignment.faa")
WINDOW = 5
KIN_FROM, KIN_TO = 27, 294        # Pkinase envelope on the anchor, from the genome-wide HMMER domain table
TOXO = "Toxoplasma_gondii"
TOXO_RECORDED_SITE = 221          # activation-loop threonine of Toxoplasma gondii AMPKα in the literature (PMID 36702847)


def ref_row():
    r = C.rows(os.path.join(C.X, "04_ampka_activation/AMPKa_activation_reference.tsv"))
    if len(r) != 1:
        sys.exit("REFUSED: the activation reference must hold exactly one row")
    return r[0]


def upstream_states(t2):
    rows, state = [], {}
    for sp in sorted({s for (_t, s) in t2}):
        a = t2[("AMPKa", sp)]
        if a["STATE"] != "PRESENT":
            continue
        v = {k: C.binary_absence_state(t2[(k, sp)]["STATE"], t2[(k, sp)]["ABSENCE_CALLABLE"]) for k in C.UPSTREAM}
        l, c = v["LKB1"], v["CAMKK"]
        if l == 0 and c == 0:
            s = "BOTH"
        elif l == 0 and c == 1:
            s = "LKB1_ONLY"
        elif l == 1 and c == 0:
            s = "CAMKK_CATEGORY_ONLY"
        elif l == 1 and c == 1:
            s = "NEITHER_CAPTURED"
        else:
            s = "UNKNOWN"
        any_cap = "YES" if 0 in (l, c) else ("NO" if (l == 1 and c == 1) else "UNKNOWN")
        state[sp] = dict(STATE=s, ANY=any_cap)
        rows.append([sp, t2[("TOR", sp)]["SUPERGROUP"], t2[("LKB1", sp)]["STATE"], t2[("LKB1", sp)]["ABSENCE_CALLABLE"],
                     t2[("CAMKK", sp)]["STATE"], t2[("CAMKK", sp)]["ABSENCE_CALLABLE"], s, any_cap,
                     "NONE_CAPTURED" if any_cap == "NO" else ("ANY_CAPTURED" if any_cap == "YES" else "NA")])
    C.write("04_ampka_activation/upstream_kinase_states.tsv",
            ["SPECIES", "SUPERGROUP", "LKB1_STATE", "LKB1_CALLABLE", "CAMKK_STATE", "CAMKK_CALLABLE",
             "UPSTREAM_KINASE_STATE", "ANY_CAPTURED_UPSTREAM", "NONE_CAPTURED_BINARY"], rows, script=SELF,
            notes=["Among species where AMPKa is PRESENT. CAMKK keeps its frozen project definition, the",
                   "CAMKK/SAK1/GRIK functional category, and is never read as strict mammalian CAMKK1/2.",
                   "NEITHER_CAPTURED means neither captured category was detected; other regulatory kinases may",
                   "exist, and this state is never described as the absence of an AMPK kinase.",
                   "A state is UNKNOWN whenever either category is ambiguous or its absence is not callable."])
    return state


def main():
    C.verify_freezes()
    ref = ref_row()
    anchor_id, anchor_pos = ref["ANCHOR_ID"], int(ref["ANCHOR_SITE_POS"])
    motif, mi = ref["MOTIF"], int(ref["MOTIF_SITE_INDEX"])
    t2 = C.read_frozen_t2()
    ustate = upstream_states(t2)

    aln = S.read_fasta(ALN)
    if anchor_id not in aln:
        sys.exit("REFUSED: the frozen RBX1 alignment lacks the anchor %s" % anchor_id)
    anc = aln[anchor_id]
    colmap = S.ref_column_map(anc)
    col = colmap[anchor_pos]
    win = [colmap[p] for p in range(anchor_pos - WINDOW, anchor_pos + WINDOW + 1) if p in colmap]
    bg = [colmap[p] for p in range(KIN_FROM, KIN_TO + 1) if p in colmap and p not in
          range(anchor_pos - WINDOW, anchor_pos + WINDOW + 1)]
    gates = [["the frozen alignment carries the anchor and the recorded site",
              "PASS" if anc[col] != S.GAP and anc.replace(S.GAP, "")[anchor_pos - 1] == "T" else "FAIL",
              "%s %d" % (anchor_id, anchor_pos)]]
    canon_key = [k for k in aln if k.startswith("Q13131")]
    if canon_key:
        q = aln[canon_key[0]]
        st, ex, aa = S.position_call(anc, q, col)
        gates.append(["the frozen alignment maps the anchor site onto the canonical site",
                      "PASS" if aa == "T" else "FAIL", "canonical residue at the anchor column: %s (%s)" % (aa, st)])

    man = {r["ALIGNMENT_KEY"]: r for r in
           C.rows(os.path.join(C.X1, "02_inputs_manifest/RBX1_sequence_manifest.AMPKa.tsv")) if r["ALIGNMENT_KEY"] != "-"}
    copy_rows = []
    for key in sorted(k for k in aln if "|" in k):
        sp, acc = key.split("|", 1)
        q = aln[key]
        st, ex, aa = S.position_call(anc, q, col)
        callable_site = st in S.CALLABLE_STATES
        qpos = S.query_pos_at_column(q, col)
        nw = ne = 0
        for c2 in win:
            s2, e2, _a2 = S.position_call(anc, q, c2)
            if s2 in S.CALLABLE_STATES:
                nw += 1
                ne += e2
        nb = nbi = 0
        for c2 in bg:
            s2, e2, _a2 = S.position_call(anc, q, c2)
            if s2 in S.CALLABLE_STATES:
                nb += 1
                nbi += e2
        wf = (ne / float(nw)) if nw else None
        bi = (nbi / float(nb)) if nb else None
        m = man.get(key, {})
        copy_rows.append([sp, acc, m.get("SEQUENCE_SHA256", "-"), m.get("DECISION_SOURCE", "-"),
                          m.get("MODEL_SENSITIVE", "-"), m.get("SOURCE_TYPE", "-"),
                          m.get("N_ORTHOLOG_COPIES_IN_SPECIES", "-"), col, qpos if qpos else "NA", aa, st,
                          "yes" if callable_site else "no",
                          (1 if aa == "T" else 0) if callable_site else "NA",
                          (1 if aa in ("T", "S") else 0) if callable_site else "NA",
                          nw, ne, "%.4f" % wf if wf is not None else "NA", nb, nbi,
                          "%.4f" % bi if bi is not None else "NA",
                          "%.4f" % (wf - bi) if None not in (wf, bi) else "NA"])
    hdr = ["SPECIES", "PROTEIN_ID", "SEQUENCE_SHA256", "DECISION_SOURCE", "MODEL_SENSITIVE", "SOURCE_TYPE",
           "N_ORTHOLOG_COPIES_IN_SPECIES", "ACTIVATION_SITE_ALIGNMENT_COL", "ACTIVATION_SITE_QUERY_POS",
           "ACTIVATION_SITE_AA", "ACTIVATION_SITE_STATE", "ACTIVATION_SITE_CALLABLE", "ACT_THR_EXACT",
           "ACT_PHOSPHO_ACCEPTOR", "ACT_LOOP_WINDOW_N_CALLABLE", "ACT_LOOP_WINDOW_N_EXACT",
           "ACT_LOOP_WINDOW_EXACT_FRACTION", "KINASE_BACKGROUND_N_CALLABLE", "KINASE_BACKGROUND_N_IDENTICAL",
           "KINASE_BACKGROUND_IDENTITY", "ACT_LOOP_MINUS_BACKGROUND"]
    C.write("04_ampka_activation/AMPKa_activation_copy.tsv", hdr, copy_rows, script=SELF,
            notes=["ACTIVATION_SITE_LABEL = AMPKa activation-loop Thr172-equivalent.",
                   "REFERENCE_COORDINATE_SYSTEM = %s, REFERENCE_SITE_POSITION = %s, REFERENCE_SITE_AA = %s."
                   % (ref["CANONICAL_ID"], ref["CANONICAL_SITE_POS"], ref["CANONICAL_SITE_AA"]),
                   "Mapping runs through the frozen RBX1 alignment and the frozen motif %s at index %d, never" % (motif, mi),
                   "through a universal numeric offset. A position beyond the observed span is not callable and",
                   "is never scored as a lost residue. The background is the frozen anchor Pkinase envelope",
                   "%d-%d excluding the plus or minus %d activation-loop window." % (KIN_FROM, KIN_TO, WINDOW)])

    # Toxoplasma check against the published site (T221); the comparison never redefines the mapping
    tox = [r for r in copy_rows if r[0] == TOXO]
    tg = [[r[1], r[8], r[9], r[10], TOXO_RECORDED_SITE,
           "MATCH" if str(r[8]) == str(TOXO_RECORDED_SITE) and r[9] == "T" else "DISCREPANCY"] for r in tox]
    C.write("11_validation/AMPKa_toxoplasma_activation_check.tsv",
            ["PROTEIN_ID", "MAPPED_QUERY_POSITION", "MAPPED_RESIDUE", "POSITION_STATE", "PROJECT_RECORDED_SITE",
             "RESULT"], tg, script=SELF,
            notes=["The project record (05_method_and_rules/PP2C_LITERATURE_BASIS.md, from PMID 36702847) states",
                   "TgAMPKalpha carries a mapped activation-loop threonine T221. RBX2 maps the frozen site through",
                   "the alignment and compares; it never forces the mapping onto 221."])
    gates.append(["Toxoplasma activation-site mapping reconciles with the project T221 record",
                  "PASS" if any(g[-1] == "MATCH" for g in tg) else "REPORTED_DISCREPANCY", str(tg)])

    by = {}
    for r in copy_rows:
        by.setdefault(r[0], []).append(r)
    sp_rows = []
    ix = {k: i for i, k in enumerate(hdr)}
    for sp in sorted(by):
        cs = by[sp]
        call = [r for r in cs if r[ix["ACTIVATION_SITE_CALLABLE"]] == "yes"]
        thr = {r[ix["ACT_THR_EXACT"]] for r in call}
        w = C.species_copy_summary(cs, lambda r: float(r[ix["ACT_LOOP_WINDOW_EXACT_FRACTION"]])
                                   if r[ix["ACT_LOOP_WINDOW_EXACT_FRACTION"]] != "NA" else None)
        bgm = C.species_copy_summary(cs, lambda r: float(r[ix["KINASE_BACKGROUND_IDENTITY"]])
                                     if r[ix["KINASE_BACKGROUND_IDENTITY"]] != "NA" else None)
        u = ustate.get(sp, {})
        sp_rows.append([sp, u.get("STATE", "AMPKa_NOT_PRESENT"), u.get("ANY", "NA"), len(cs), len(call),
                        "NA" if not thr else ("1" if thr == {1} else ("0" if thr == {0} else "DUPLICATE_DISCORDANT")),
                        "NA" if not thr else ("1" if 1 in thr else "0"),
                        "yes" if len(thr) > 1 else "no",
                        "%.4f" % w["MEDIAN"] if w["MEDIAN"] is not None else "NA",
                        "%.4f" % w["MIN"] if w["MIN"] is not None else "NA",
                        "%.4f" % w["MAX"] if w["MAX"] is not None else "NA",
                        "%.4f" % bgm["MEDIAN"] if bgm["MEDIAN"] is not None else "NA"])
    C.write("04_ampka_activation/AMPKa_activation_species.tsv",
            ["SPECIES", "UPSTREAM_KINASE_STATE", "ANY_CAPTURED_UPSTREAM", "N_AMPKa_COPIES", "N_ACT_SITE_CALLABLE",
             "ACT_THR_ALL_COPIES", "ACT_THR_ANY_COPY", "COPY_DISCORDANT", "ACT_LOOP_WINDOW_MEDIAN",
             "ACT_LOOP_WINDOW_MIN", "ACT_LOOP_WINDOW_MAX", "KINASE_BACKGROUND_IDENTITY_MEDIAN"], sp_rows, script=SELF,
            notes=["Every copy is kept and none is chosen by its activation-loop score. Where callable copies",
                   "disagree the species reads DUPLICATE_DISCORDANT rather than being forced to a state."])
    C.write("11_validation/AMPKa_activation_gates.tsv", ["GATE", "RESULT", "DETAIL"], gates, script=SELF)
    C.record_outputs("y3", SELF)
    nt = sum(1 for r in sp_rows if r[5] == "1")
    n0 = sum(1 for r in sp_rows if r[5] == "0")
    print("y3: %d copies, %d species; activation threonine retained in %d species, absent in %d, discordant %d"
          % (len(copy_rows), len(sp_rows), nt, n0, sum(1 for r in sp_rows if r[5] == "DUPLICATE_DISCORDANT")))
    if any(g[1] == "FAIL" for g in gates):
        sys.exit("REFUSED: an activation-site gate failed: %s" % [g for g in gates if g[1] == "FAIL"])


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
