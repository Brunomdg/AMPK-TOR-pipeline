#!/usr/bin/env python3
# >>> deposit header
# Purpose: Derives the canonical TORC1 (TOR, RAPTOR, LST8) and TORC2 (TOR, RICTOR, SIN1, LST8) core states of every
#          species, the four combined RAPTOR and RICTOR states, and species counts of core state against AMPK
#          subunit loss with odds ratios; descriptive only.
# Methods: Phylogenetic regressions
# Inputs:  the ten-target three-state table `T2_presence_long.tsv` and the retained ortholog copies (through
#          `phenotype_tables_chain2.py`)
# Outputs: `03_torc_core/TORC_core_states.tsv`, `RAPTOR_RICTOR_four_state.tsv` and `core_AMPK_descriptive.tsv`
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
# Ran as: rbx2_torc_core.py
## Job y1 of the second chain. Canonical TORC1 and TORC2 core-complement states derived by fixed rules from the
## ten-target three-state table, the descriptive four combined RAPTOR and RICTOR states, and the six two-by-two
## tables of core state against AMPK subunit loss. One absence in a proteome where absence can be assessed disrupts
## a canonical core; an ambiguous or not-assessable component never does on its own, and never cancels an
## assessable absence of another component.
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "4_sequence_phenotypes", "lib"))
import phenotype_tables_chain2 as C  # noqa: E402

SELF = os.path.abspath(__file__)
Z = 1.959964


def odds_ratio(a, b, c, d):
    corr = "none"
    if min(a, b, c, d) == 0:
        a, b, c, d, corr = a + 0.5, b + 0.5, c + 0.5, d + 0.5, "haldane_0.5"
    orv = (a * d) / float(b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return orv, math.exp(math.log(orv) - Z * se), math.exp(math.log(orv) + Z * se), corr


def main():
    C.verify_freezes()
    t2 = C.read_frozen_t2()
    species = sorted({s for (_t, s) in t2})
    C.assert_unique_keys(species, lambda s: s, "species list")
    prov = {}
    for tgt in C.TORC1_CORE + C.TORC2_CORE + C.SUBUNITS:
        for r in C.t1_orthologs(tgt):
            prov.setdefault((tgt, r["SPECIES"]), []).append(r)

    core_rows, four_rows = [], []
    state = {}
    for sp in species:
        row = [sp, t2[("TOR", sp)]["SUPERGROUP"], t2[("TOR", sp)]["CLADE"], t2[("TOR", sp)]["DATA_GRADE"]]
        cells = {}
        for tgt in ("TOR", "RAPTOR", "RICTOR", "SIN1", "LST8"):
            c = t2[(tgt, sp)]
            cells[tgt] = c
            row += [c["STATE"], c["ABSENCE_CALLABLE"]]
        c1 = [cells[t] for t in C.TORC1_CORE]
        c2 = [cells[t] for t in C.TORC2_CORE]
        s1, s2 = C.core_state(c1), C.core_state(c2)
        row += [s1, C.core_basis(c1), s2, C.core_basis(c2)]
        # descriptive RAPTOR and RICTOR state; only an assessable absence gives a minus
        v = {t: C.binary_absence_state(cells[t]["STATE"], cells[t]["ABSENCE_CALLABLE"]) for t in ("RAPTOR", "RICTOR")}
        scaf = "UNKNOWN" if None in v.values() else "%s/%s" % ("RAPTOR-" if v["RAPTOR"] else "RAPTOR+",
                                                              "RICTOR-" if v["RICTOR"] else "RICTOR+")
        row.append(scaf)
        subs = {}
        for sub in C.SUBUNITS:
            cell = t2[(sub, sp)]
            b = C.binary_absence_state(cell["STATE"], cell["ABSENCE_CALLABLE"])
            subs[sub] = b
            cs = prov.get((sub, sp), [])
            row += [cell["STATE"], cell["ABSENCE_CALLABLE"], "NA" if b is None else b,
                    ("yes" if all(x["DECISION_SOURCE"] not in C.HIGH_PROV for x in cs) else "no") if cs else "-",
                    ("yes" if all(x["MODEL_SENSITIVE"] == "yes" for x in cs) else "no") if cs else "-"]
        core_rows.append(row)
        state[sp] = dict(TORC1=s1, TORC2=s2, SCAFFOLD=scaf, **subs)

    hdr = ["SPECIES", "SUPERGROUP", "CLADE", "DATA_GRADE"]
    for tgt in ("TOR", "RAPTOR", "RICTOR", "SIN1", "LST8"):
        hdr += ["%s_STATE" % tgt, "%s_CALLABLE" % tgt]
    hdr += ["TORC1_CORE_STATE", "TORC1_CORE_BASIS", "TORC2_CORE_STATE", "TORC2_CORE_BASIS", "SCAFFOLD_STATE"]
    for sub in C.SUBUNITS:
        hdr += ["%s_STATE" % sub, "%s_CALLABLE" % sub, "%s_absent" % sub, "%s_R66_ONLY" % sub,
                "%s_ALL_MODEL_SENSITIVE" % sub]
    C.write("03_torc_core/TORC_core_states.tsv", hdr, core_rows, script=SELF,
            notes=["TORC1_CORE = %s; TORC2_CORE = %s. These are canonical core-complement states, not a claim"
                   % (" + ".join(C.TORC1_CORE), " + ".join(C.TORC2_CORE)),
                   "that every eukaryote assembles a mammalian-identical complex.",
                   "INTACT requires every component PRESENT. DISRUPTED requires at least one callable ABSENT.",
                   "UNKNOWN covers the rest: an ambiguous or non-callable component never disrupts a core by",
                   "itself, and it never erases a proven callable loss of another required component."])

    for s in ("RAPTOR+/RICTOR+", "RAPTOR-/RICTOR+", "RAPTOR+/RICTOR-", "RAPTOR-/RICTOR-", "UNKNOWN"):
        sps = [sp for sp in species if state[sp]["SCAFFOLD"] == s]
        row = [s, len(sps)]
        for sub in C.SUBUNITS:
            u = [sp for sp in sps if state[sp][sub] is not None]
            k = sum(1 for sp in u if state[sp][sub] == 1)
            row += [len(u), k, "%.4f" % (k / float(len(u))) if u else "NA"]
        row.append(";".join(sorted({t2[("TOR", sp)]["SUPERGROUP"] for sp in sps})) or "-")
        four_rows.append(row)
    fh = ["SCAFFOLD_STATE", "N_SPECIES"]
    for sub in C.SUBUNITS:
        fh += ["N_CALLABLE_%s" % sub, "N_LOST_%s" % sub, "LOSS_FRACTION_%s" % sub]
    fh += ["SUPERGROUPS"]
    C.write("03_torc_core/RAPTOR_RICTOR_four_state.tsv", fh, four_rows, script=SELF,
            notes=["Descriptive. A minus is called only from a callable ABSENT. This descriptor never replaces",
                   "the formal TORC1 and TORC2 core definitions and is not interpreted without the models."])

    desc = []
    for core in ("TORC1", "TORC2"):
        for sub in C.SUBUNITS:
            use = [sp for sp in species if state[sp][core] in ("INTACT", "DISRUPTED") and state[sp][sub] is not None]
            a = sum(1 for sp in use if state[sp][core] == "DISRUPTED" and state[sp][sub] == 1)
            b = sum(1 for sp in use if state[sp][core] == "DISRUPTED" and state[sp][sub] == 0)
            c = sum(1 for sp in use if state[sp][core] == "INTACT" and state[sp][sub] == 1)
            d = sum(1 for sp in use if state[sp][core] == "INTACT" and state[sp][sub] == 0)
            fi = c / float(c + d) if c + d else None
            fd = a / float(a + b) if a + b else None
            orv, lo, hi, corr = odds_ratio(a, b, c, d) if (a + b) and (c + d) else (None, None, None, "not_estimable")
            desc.append([core, sub, len(species), len(use),
                         sum(1 for sp in species if state[sp][core] == "INTACT"),
                         sum(1 for sp in species if state[sp][core] == "DISRUPTED"),
                         sum(1 for sp in species if state[sp][core] == "UNKNOWN"),
                         c, c + d, a, a + b,
                         "%.4f" % fi if fi is not None else "NA", "%.4f" % fd if fd is not None else "NA",
                         "%.4f" % (fd - fi) if None not in (fi, fd) else "NA",
                         "%.3f" % orv if orv else "NA", "%.3f" % lo if lo else "NA", "%.3f" % hi if hi else "NA", corr])
    C.write("03_torc_core/core_AMPK_descriptive.tsv",
            ["CORE", "SUBUNIT", "N_TOTAL", "N_USABLE", "N_CORE_INTACT", "N_CORE_DISRUPTED", "N_CORE_UNKNOWN",
             "AMPK_LOSS_INTACT_N", "AMPK_LOSS_INTACT_DENOM", "AMPK_LOSS_DISRUPTED_N", "AMPK_LOSS_DISRUPTED_DENOM",
             "LOSS_FREQ_INTACT", "LOSS_FREQ_DISRUPTED", "ABSOLUTE_DIFFERENCE", "ODDS_RATIO", "CI_LOW", "CI_HIGH",
             "ZERO_CELL_CORRECTION"], desc, script=SELF,
            notes=["Counts of species, not of independent evolutionary events. Species with an UNKNOWN core or a",
                   "non-callable subunit are excluded from the usable denominator and counted separately."])
    C.record_outputs("y1", SELF)
    n1 = sum(1 for sp in species if state[sp]["TORC1"] == "DISRUPTED")
    n2 = sum(1 for sp in species if state[sp]["TORC2"] == "DISRUPTED")
    print("y1: %d species, TORC1 disrupted %d intact %d, TORC2 disrupted %d intact %d"
          % (len(species), n1, sum(1 for sp in species if state[sp]["TORC1"] == "INTACT"),
             n2, sum(1 for sp in species if state[sp]["TORC2"] == "INTACT")))


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
