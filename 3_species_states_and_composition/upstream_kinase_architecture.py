#!/usr/bin/env python3
# >>> deposit header
# Purpose: Gives each of the 210 species an upstream-kinase state (LKB1 and the CAMKK category both found, one of
#          them, neither, or unknown), with LKB1 sequence-support tiers, and writes descriptive breakdowns.
# Methods: Phylogenetic regressions
# Inputs:  the species-state table `T2_presence_long.tsv` and ortholog table `T1_ortholog_table.tsv`;
#          `10_validation/LKB1_sequence_validation.tsv`, `04_camkk_subfamilies/CAMKK_constituent_species.tsv`;
#          the rebuilt `TORC_core_states.tsv`, `AMPKa_activation_species.tsv` and `upstream_kinase_states.tsv`
# Outputs: `05_upstream_architecture/CAMKK_LKB1_architecture.tsv`, `RBX2_ARCHITECTURE_REPRODUCTION.tsv`,
#          `architecture_by_AMPKa_state.tsv`, `TOR_upstream_decomposition.tsv`, `CAMKK_constituent_by_clade.tsv`
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
## Assigns each of the 210 species an upstream-kinase state from its LKB1 and CAMKK-category states and writes
## descriptive breakdowns. The state describes which upstream kinases were found and does not require AMPKα; AMPKα
## presence is a column, so the subset of species that retain AMPKα can be taken later. The TORC1 and TORC2 core
## states are derived from the species-state table and checked against the rebuilt core-state table. Every row
## gives the kind of evidence behind its LKB1 call, so the table shows how many LKB1-positive species rest on
## candidate copies rather than on declared positive controls. An LKB1 presence without sequence support becomes
## UNKNOWN_LKB1_SUPPORT in the sequence-supported state and is never turned into an absence.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "2_candidate_discovery_and_orthology", "upstream_kinases", "lib"))
import upstream_kinase_tables as C  # noqa: E402

SELF = os.path.abspath(__file__)


def arch(lk_state, lk_call, ck_state, ck_call):
    lk_abs = lk_state == "ABSENT" and lk_call == "yes"
    ck_abs = ck_state == "ABSENT" and ck_call == "yes"
    if lk_state == "PRESENT" and ck_state == "PRESENT":
        return "BOTH_CAPTURED"
    if lk_state == "PRESENT" and ck_abs:
        return "LKB1_ONLY"
    if lk_abs and ck_state == "PRESENT":
        return "CAMKK_CATEGORY_ONLY"
    if lk_abs and ck_abs:
        return "NEITHER_CAPTURED"
    return "UNKNOWN"


def main():
    C.verify_freezes()
    t2 = C.read_frozen_t2()
    species = sorted({k[1] for k in t2})
    clade, sgrp = C.clade_of(), C.supergroup_of()
    lkv = {r["PROTEIN_ID"]: r for r in
           C.rows(os.path.join(C.X3, "10_validation/LKB1_sequence_validation.tsv"))}
    lk_by_sp = {}
    for r in lkv.values():
        lk_by_sp.setdefault(r["SPECIES"], []).append(r)
    ck = {r["SPECIES"]: r for r in
          C.rows(os.path.join(C.X3, "04_camkk_subfamilies/CAMKK_constituent_species.tsv"))}
    act = {r["SPECIES"]: r for r in
           C.rows(os.path.join(C.X2, "04_ampka_activation/AMPKa_activation_species.tsv"))}
    frozen_core = {r["SPECIES"]: r for r in
                   C.rows(os.path.join(C.X2, "03_torc_core/TORC_core_states.tsv"))}
    t1 = {}
    for r in C.rows(os.path.join(C.PUB, "T1_ortholog_table_all20.tsv")):
        if r["TARGET"] == "LKB1" and r["FINAL_CLASS"] == "ORTHOLOG":
            t1.setdefault(r["SPECIES"], set()).add(r["DECISION_SOURCE"])

    rows, core_mismatch = [], 0
    for sp in species:
        g = lambda tgt: t2[(tgt, sp)]
        c1 = C.core_state([g(x) for x in C.TORC1_CORE])
        c2 = C.core_state([g(x) for x in C.TORC2_CORE])
        f = frozen_core.get(sp)
        if f and (f.get("TORC1_CORE_STATE") != c1 or f.get("TORC2_CORE_STATE") != c2):
            core_mismatch += 1
        lk, cka = g("LKB1"), g("CAMKK")
        a = arch(lk["STATE"], lk["ABSENCE_CALLABLE"], cka["STATE"], cka["ABSENCE_CALLABLE"])
        ## Decided per species, not per copy, because copies are not independent species. A declared positive
        ## control is LKB1 by external definition and is never demoted by the leave-one-out test, which removes
        ## references on purpose, so species with a declared positive are usable in both sensitivity tiers.
        copies = lk_by_sp.get(sp, [])
        decl = [r for r in copies if r["DECISION_SOURCE"] == "DECLARED_POSITIVE"]
        r66 = [r for r in copies if r["DECISION_SOURCE"] != "DECLARED_POSITIVE"]
        r66sup = ("STRONG" if any(r["LKB1_SEQUENCE_SUPPORT"] == "LKB1_STRONGLY_SUPPORTED" for r in r66)
                  else "COMPATIBLE" if any(r["LKB1_SEQUENCE_SUPPORT"] == "LKB1_COMPATIBLE" for r in r66)
                  else "UNRESOLVED" if r66 else "-")
        if not copies:
            lk_sup, basis = ("-", "no LKB1 ortholog copy") if lk["STATE"] != "PRESENT" else \
                            ("LKB1_SEQUENCE_UNRESOLVED", "present but no validated copy")
        elif decl:
            lk_sup, basis = "LKB1_DECLARED_POSITIVE", "declared positive; externally declared LKB1 truth"
        elif r66sup == "STRONG":
            lk_sup, basis = "LKB1_STRONGLY_SUPPORTED", "at least one strongly supported R66 copy"
        elif r66sup == "COMPATIBLE":
            lk_sup, basis = "LKB1_COMPATIBLE", "at least one compatible R66 copy, none strongly supported"
        else:
            lk_sup, basis = "LKB1_SEQUENCE_UNRESOLVED", "every retained R66 copy unresolved"
        present = lk["STATE"] == "PRESENT"
        useA = "yes" if (not present or decl or r66sup in ("STRONG", "COMPATIBLE")) else "NA"
        useB = "yes" if (not present or decl or r66sup == "STRONG") else "NA"
        ## An LKB1-positive species whose copies do not reach tier A becomes unknown here, never absent.
        if present and useA == "NA":
            a_seq = "UNKNOWN_LKB1_SUPPORT"
        else:
            a_seq = a
        ckr = ck.get(sp, {})
        rows.append([sp, clade.get(sp, "UNKNOWN"), sgrp.get(sp, "UNKNOWN"),
                     g("AMPKa")["STATE"], g("AMPKa")["ABSENCE_CALLABLE"],
                     "yes" if g("AMPKa")["STATE"] == "PRESENT" else "no",
                     lk["STATE"], lk["ABSENCE_CALLABLE"],
                     ",".join(sorted(t1.get(sp, []))) or "-", lk_sup, basis, r66sup, useA, useB,
                     cka["STATE"], cka["ABSENCE_CALLABLE"],
                     ckr.get("CAMKK_L1_CLASS", "-"), ckr.get("CAMKK_L2_CLASS", "-"),
                     ckr.get("N_COPIES", "0"), ckr.get("UNRESOLVED_COPY_PRESENT", "-"),
                     a, a_seq, g("TOR")["STATE"], c1, c2,
                     act.get(sp, {}).get("ACT_THR_ALL_COPIES", "-"), lk["SOURCE_TYPE"]])
    if core_mismatch:
        sys.exit("REFUSED: %d core states derived from T2 disagree with the frozen RBX2 table" % core_mismatch)
    C.write("05_upstream_architecture/CAMKK_LKB1_architecture.tsv",
            ["SPECIES", "CLADE", "SUPERGROUP", "AMPKa_STATE", "AMPKa_ABSENCE_CALLABLE",
             "AMPKa_PRESENT_FOR_MECHANISTIC_SUBSET", "LKB1_STATE_FROZEN", "LKB1_CALLABILITY",
             "LKB1_DECISION_SOURCE", "LKB1_SEQUENCE_SUPPORT", "LKB1_SEQUENCE_SUPPORT_BASIS",
             "LKB1_R66_SEQUENCE_SUPPORT", "LKB1_USABLE_SENSITIVITY_A", "LKB1_USABLE_SENSITIVITY_B",
             "CAMKK_STATE",
             "CAMKK_CALLABILITY", "CAMKK_L1_CLASS", "CAMKK_L2_CLASS", "CAMKK_N_COPIES", "CAMKK_UNRESOLVED_COPY_PRESENT",
             "UPSTREAM_ARCHITECTURE_FROZEN", "UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED", "TOR_STATE",
             "TORC1_CORE_STATE", "TORC2_CORE_STATE", "ACTIVATION_THR_STATE", "SOURCE_TYPE"], rows, script=SELF,
            notes=["All 210 species. The architecture is a state of the captured upstream-kinase repertoire and",
                   "is assigned from callable LKB1 and CAMKK states alone, regardless of whether AMPKa is",
                   "retained; conditioning it on AMPKa would condition on a downstream variable.",
                   "UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED reads UNKNOWN_LKB1_SUPPORT where a frozen LKB1-present",
                   "cell has no sequence-supported copy, and never converts such a cell to LKB1_ABSENT.",
                   "Core states are derived from frozen T2 by core_state() and asserted equal to the frozen",
                   "RBX2 TORC_core_states table.",
                   "LKB1 sensitivity tiers, species level: A admits declared positives plus R66 species with",
                   "a strongly supported or compatible copy; B admits declared positives plus strongly",
                   "supported R66 species only. A frozen LKB1-positive species that does not reach the tier",
                   "becomes NA and is never turned into an absence. Declared positives are externally declared",
                   "truth and are never demoted by the anchor-depleted leave-one-out exercise."])

    # ---- check against the rebuilt upstream-kinase table, on the species that retain AMPKα only
    sub = [r for r in rows if r[5] == "yes"]
    counts = {}
    for r in sub:
        counts[r[20]] = counts.get(r[20], 0) + 1
    ## The expected counts are read from the rebuilt upstream-kinase table, not typed in. That table labels the
    ## state with both kinases found `BOTH` and this script `BOTH_CAPTURED`; only the label is mapped.
    _up = C.rows(os.path.join(C.X2, "04_ampka_activation/upstream_kinase_states.tsv"))
    expect = {}
    for _r in _up:
        _k = "BOTH_CAPTURED" if _r["UPSTREAM_KINASE_STATE"] == "BOTH" else _r["UPSTREAM_KINASE_STATE"]
        expect[_k] = expect.get(_k, 0) + 1
    gate = [["N_AMPKa_PRESENT", str(len(sub)), str(len(_up)), "PASS" if len(sub) == len(_up) else "FAIL"]] + \
           [[k, str(counts.get(k, 0)), str(expect.get(k, 0)), "PASS" if counts.get(k, 0) == expect.get(k, 0) else "FAIL"]
            for k in sorted(set(expect) | set(counts))]
    C.write("05_upstream_architecture/RBX2_ARCHITECTURE_REPRODUCTION.tsv",
            ["QUANTITY", "RBX3_ALL210_RESTRICTED_TO_AMPKa_PRESENT", "REBUILT_RBX2", "RESULT"], gate, script=SELF,
            notes=["The 210-species construction must not have altered the already frozen RBX2 mechanistic",
                   "subset. This gate applies only on the AMPKa-present restriction, which is the scope of the",
                   "frozen upstream_kinase_states table."])
    if any(r[3] != "PASS" for r in gate):
        sys.exit("REFUSED: the AMPKa-present restriction does not reproduce the frozen RBX2 architecture")

    # ---- descriptive decompositions
    d1 = {}
    for r in rows:
        d1.setdefault((r[3], r[20]), 0)
        d1[(r[3], r[20])] += 1
    C.write("05_upstream_architecture/architecture_by_AMPKa_state.tsv",
            ["AMPKa_STATE", "UPSTREAM_ARCHITECTURE", "N"],
            [[k[0], k[1], v] for k, v in sorted(d1.items())], script=SELF,
            notes=["Descriptive: whether CAMKK and LKB1 disappearance accompanies AMPKa loss itself.",
                   "No inferential family is attached to this table."])
    tor = {}
    for r in sub:
        st = r[22]
        e = tor.setdefault(st, dict(n=0, both=0, ck=0, lk=0, neither=0, ck_abs=0, ck_den=0, lk_abs=0,
                                    lk_den=0, nc=0, nc_den=0))
        e["n"] += 1
        e[{"BOTH_CAPTURED": "both", "CAMKK_CATEGORY_ONLY": "ck", "LKB1_ONLY": "lk",
           "NEITHER_CAPTURED": "neither"}.get(r[20], "n")] += 1 if r[20] in (
            "BOTH_CAPTURED", "CAMKK_CATEGORY_ONLY", "LKB1_ONLY", "NEITHER_CAPTURED") else 0
        if C.binary_absence_state(r[14], r[15]) is not None:
            e["ck_den"] += 1
            e["ck_abs"] += C.binary_absence_state(r[14], r[15])
        if C.binary_absence_state(r[6], r[7]) is not None:
            e["lk_den"] += 1
            e["lk_abs"] += C.binary_absence_state(r[6], r[7])
        if r[20] != "UNKNOWN":
            e["nc_den"] += 1
            e["nc"] += 1 if r[20] == "NEITHER_CAPTURED" else 0
    C.write("05_upstream_architecture/TOR_upstream_decomposition.tsv",
            ["TOR_STATE", "N_USABLE", "N_BOTH_CAPTURED", "N_CAMKK_ONLY", "N_LKB1_ONLY", "N_NEITHER_CAPTURED",
             "CAMKK_ABSENT_N", "CAMKK_ABSENT_DENOM", "LKB1_ABSENT_N", "LKB1_ABSENT_DENOM",
             "NEITHER_CAPTURED_N", "NEITHER_CAPTURED_DENOM"],
            [[k, v["n"], v["both"], v["ck"], v["lk"], v["neither"], v["ck_abs"], v["ck_den"], v["lk_abs"],
              v["lk_den"], v["nc"], v["nc_den"]] for k, v in sorted(tor.items())], script=SELF,
            notes=["Callable AMPKa-present species, so the reader can see whether the frozen RBX2 association",
                   "between NEITHER_CAPTURED and TOR loss arises because both upstream systems disappear",
                   "together or because one dominates. Descriptive; no test is attached."])
    cl = {}
    for r in rows:
        if r[14] != "PRESENT":
            continue
        e = cl.setdefault(r[1], dict(n=0, met=0, snf=0, sak=0, grik=0, oth=0, unr=0))
        e["n"] += 1
        e["met"] += r[16] == "CAMKK_METAZOAN_LIKE"
        e["snf"] += r[16] == "SNF1_SNRK1_ACTIVATOR_LIKE"
        e["sak"] += r[17] == "SAK1_TOS3_ELM1_LIKE"
        e["grik"] += r[17] == "GRIK_SNAK_LIKE"
        e["oth"] += r[16] in ("OTHER_CAMKK_CATEGORY", "MULTIFAMILY", "ALIGNMENT_SENSITIVE")
        e["unr"] += r[16] in ("SUBFAMILY_UNRESOLVED", "-")
    C.write("05_upstream_architecture/CAMKK_constituent_by_clade.tsv",
            ["CLADE", "N_CAMKK_PRESENT", "N_CAMKK_METAZOAN_LIKE", "N_SNF1_SNRK1_ACTIVATOR_LIKE",
             "N_SAK1_TOS3_ELM1_LIKE", "N_GRIK_SNAK_LIKE", "N_OTHER", "N_UNRESOLVED"],
            [[k, v["n"], v["met"], v["snf"], v["sak"], v["grik"], v["oth"], v["unr"]]
             for k, v in sorted(cl.items())], script=SELF,
            notes=["Descriptive distribution. A clade is never given a significance test automatically and",
                   "constituent-family distributions are not treated as independent lineage replicates."])
    C.record_outputs("r4", SELF)
    print("r4: %d species, AMPKa-present subset %d, RBX2 architecture reproduction PASS" % (len(rows), len(sub)))


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
