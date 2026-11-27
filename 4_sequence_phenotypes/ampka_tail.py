#!/usr/bin/env python3
# >>> deposit header
# Purpose: Measures the C-terminal region of every retained AMPKα copy: its length after the residue homologous to
#          human PRKAA1 residue 290 and its predicted disorder, with the length from the kinase-domain end as a
#          sensitivity.
# Methods: Sequence phenotypes of the retained subunits
# Inputs:  `04_ampka_ctd/AMPKa_orthologs_fulllength.faa` and `RBX1_sequence_manifest.AMPKa.tsv` (from job x1); the
#          human proteome `Homo_sapiens.faa` (anchor NP_996790.3) and `q13131_canonical.faa`; the genome-wide
#          HMMER domain tables `<species>.domtbl.gz`
# Outputs: `04_ampka_ctd/AMPKa_CTD_alignment.faa`, `AMPKa_CTD_copy.tsv` (one row per copy) and
#          `AMPKa_CTD_disorder_long.tsv` (one row per scored residue); `10_validation/AMPKa_boundary_gates.tsv`
#          (checks of the boundary residue and of the human C-terminal length)
# Tools:   metapredict 3.0.2 (version checked at start); MAFFT 7.526 [env record] from the orthofinder3
#          environment (--globalpair --maxiterate 1000); Python 3; versions: see envs/
# <<< deposit header
## Job x2 of this chain, AMPKα branch. Repeats the original AMPKα C-terminal measurement (folder `P03_actd`) on
## every accepted AMPKα ortholog copy: the length C-terminal of the residue homologous to canonical PRKAA1 residue
## 290 (residue 305 of the anchor NP_996790.3), and the metapredict disorder over that region. Runs in the
## environment that holds metapredict 3.0.2. The length measured from the end of the kinase domain is a separate
## sensitivity column and never replaces the primary length.
import gzip, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "4_sequence_phenotypes", "lib"))
import phenotype_tables_chain1 as L  # noqa: E402
import sequence_scoring as S  # noqa: E402

SELF = os.path.abspath(__file__)
ANCHOR = "NP_996790.3"
CANON = "Q13131"
ANCHOR_BOUNDARY = 305            # canonical 290 (E) = anchor 305 (E), fixed in advance in `p03_boundary_freeze.tsv`
COV_FROM, COV_TO = 285, 325      # coverage window of the original analysis, in anchor numbering
COV_MIN = 10                     # original rule: fewer than ten non-gap query residues in the window leave the anchor unplaced
HUMAN_CTD_EXPECTED = 269         # C-terminal length of both human references (checks F1 and F2 of the original analysis)
LEGACY_EXPANSION_LEN = 300       # expansion length of the original analysis, used here only to label an observed extension
SANITISE = {"U": "C", "B": "N", "Z": "Q", "J": "L", "O": "K"}
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def sanitise(seq):
    """Substitution of the original analysis: the disorder predictor never sees a non-standard residue (U to C,
    B to N, Z to Q, J to L, O to K, anything else to G); lengths use the unmodified sequence."""
    out, n = [], 0
    for c in seq.upper():
        if c in STANDARD:
            out.append(c)
        else:
            out.append(SANITISE.get(c, "G"))
            n += 1
    return "".join(out), n


def kinase_domain_ends(species_to_accs):
    """End of the best-scoring Pkinase domain envelope (lowest independent E-value) per ortholog copy, from the
    genome-wide HMMER domain tables (`hmmer_genomewide/domtbl`)."""
    d = os.path.join(L.R1, "03_evidence_by_method/hmmer_genomewide/domtbl")
    out = {}
    for sp, accs in species_to_accs.items():
        p = os.path.join(d, sp + ".domtbl.gz")
        if not os.path.exists(p):
            continue
        want = {"%s||%s" % (sp, a): a for a in accs}
        with gzip.open(p, "rt", errors="replace") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                f = ln.split()
                if len(f) < 21 or f[0] not in want or not f[3].startswith("Pkinase"):
                    continue
                iev, env_to = float(f[12]), int(f[20])
                k = (sp, want[f[0]])
                if k not in out or iev < out[k][0]:
                    out[k] = (iev, env_to, f[3])
    return out


def main():
    L.verify_freezes()
    import metapredict as meta
    if meta.__version__ != L.METAPREDICT_VERSION:
        sys.exit("REFUSED: metapredict %s is not the pinned %s" % (meta.__version__, L.METAPREDICT_VERSION))
    man = {r["ALIGNMENT_KEY"]: r for r in L.rows(os.path.join(L.X, "02_inputs_manifest/RBX1_sequence_manifest.AMPKa.tsv"))
           if r["ALIGNMENT_KEY"] != "-"}
    copies = S.read_fasta(os.path.join(L.X, "04_ampka_ctd/AMPKa_orthologs_fulllength.faa"))
    human = S.read_fasta(os.path.join(L.PROT, "Homo_sapiens.faa"))
    canon = S.read_fasta(os.path.join(L.CONCORD, "ms210/q13131_canonical.faa"))
    ckey = list(canon)[0]
    if ANCHOR not in human:
        sys.exit("REFUSED: anchor %s is not in the frozen Homo_sapiens proteome" % ANCHOR)
    anchor_seq, canon_seq = human[ANCHOR], canon[ckey]
    gates = [["anchor %s residue %d is E (frozen boundary)" % (ANCHOR, ANCHOR_BOUNDARY),
              "PASS" if anchor_seq[ANCHOR_BOUNDARY - 1] == "E" else "FAIL", anchor_seq[ANCHOR_BOUNDARY - 1]],
             ["canonical %s residue %d is E" % (CANON, L.CANON_CTD_BOUNDARY),
              "PASS" if canon_seq[L.CANON_CTD_BOUNDARY - 1] == "E" else "FAIL", canon_seq[L.CANON_CTD_BOUNDARY - 1]]]
    if gates[0][1] == "FAIL" or gates[1][1] == "FAIL":
        sys.exit("REFUSED: the frozen CTD boundary does not reproduce")

    seqs = {ANCHOR: anchor_seq, CANON: canon_seq}
    seqs.update(copies)
    aln_path = os.path.join(L.X, "04_ampka_ctd/AMPKa_CTD_alignment.faa")
    cmd, aln = S.mafft(seqs, aln_path)
    L.WRITTEN.append(aln_path)
    anchor_aln = aln[ANCHOR]
    colmap = S.ref_column_map(anchor_aln)
    bcol = colmap[ANCHOR_BOUNDARY]
    cov_cols = [colmap[i] for i in range(COV_FROM, COV_TO + 1)]

    sp_accs = {}
    for k in copies:
        sp, acc = k.split("|", 1)
        sp_accs.setdefault(sp, set()).add(acc)
    dom = kinase_domain_ends(sp_accs)

    rows, long_rows = [], []
    for key in sorted(list(copies) + [ANCHOR, CANON]):
        raw = seqs[key]
        q = aln[key]
        cov = sum(1 for c in cov_cols if q[c] != S.GAP)
        at_col = S.query_pos_at_column(q, bcol)
        anchor_state = "MATCHED_COLUMN" if at_col else "GAP_AT_BOUNDARY_COLUMN"
        pos = at_col or S.last_query_pos_before(q, bcol)
        callable_anchor = "yes" if cov >= COV_MIN and pos > 0 else "no"
        ctd = raw[pos:] if callable_anchor == "yes" else ""
        clean, nns = sanitise(raw)
        d = meta.predict_disorder(clean)
        d = list(getattr(d, "disorder", d))
        ds = d[pos:] if callable_anchor == "yes" else []
        run = best = 0
        for v in ds:
            run = run + 1 if v > 0.5 else 0
            best = max(best, run)
        m = man.get(key, {})
        src = m.get("SOURCE_TYPE", "genome" if key in (ANCHOR, CANON) else "-")
        grade = m.get("DATA_GRADE", "G" if key in (ANCHOR, CANON) else "-")
        if grade == "G":
            cstat = "CALLABLE" if callable_anchor == "yes" else "ANCHOR_NOT_CALLABLE"
        elif callable_anchor == "yes" and len(ctd) >= LEGACY_EXPANSION_LEN:
            cstat = "POSITIVE_EXTENSION_ONLY"
        else:
            cstat = "C_TERMINUS_UNRESOLVED"
        dk = dom.get(tuple(key.split("|", 1))) if "|" in key else None
        legacy = "NA"
        if callable_anchor == "yes" and ds:
            legacy = "EXPANDED" if (len(ctd) >= LEGACY_EXPANSION_LEN and statistics.mean(ds) >= 0.50) else "NOT_EXPANDED"
        rows.append([m.get("TARGET", "REFERENCE"), key.split("|")[0], key.split("|")[-1], m.get("SEQUENCE_SHA256", "-"),
                     m.get("DECISION_SOURCE", "-"), m.get("MODEL_SENSITIVE", "-"), src, grade,
                     m.get("N_ORTHOLOG_COPIES_IN_SPECIES", "-"), len(raw), pos, anchor_state, cov, callable_anchor,
                     len(ctd) if callable_anchor == "yes" else "NA",
                     ("%.4f" % statistics.mean(ds)) if ds else "NA",
                     ("%.4f" % statistics.median(ds)) if ds else "NA",
                     ("%.4f" % (sum(1 for v in ds if v > 0.5) / float(len(ds)))) if ds else "NA",
                     best if ds else "NA", len(ds), nns, cstat, legacy,
                     (len(raw) - dk[1]) if dk else "NA", ("%.2g" % dk[0]) if dk else "NA", dk[2] if dk else "-"])
        if callable_anchor == "yes" and key in copies:
            for i, v in enumerate(ds):
                long_rows.append([m.get("TARGET", "-"), key.split("|")[0], key.split("|")[1], pos + i + 1,
                                  ctd[i] if i < len(ctd) else "-", "%.4f" % v])

    hdr = ["TARGET", "SPECIES", "PROTEIN_ID", "SEQUENCE_SHA256", "DECISION_SOURCE", "MODEL_SENSITIVE", "SOURCE_TYPE",
           "DATA_GRADE", "N_ORTHOLOG_COPIES_IN_SPECIES", "PROTEIN_LENGTH", "ANCHOR_QUERY_POSITION", "ANCHOR_STATE",
           "ANCHOR_WINDOW_COVERAGE", "ANCHOR_CALLABLE", "CTD_LENGTH", "CTD_MEAN_DISORDER", "CTD_MEDIAN_DISORDER",
           "CTD_FRACTION_DISORDER_GT_0_5", "CTD_LONGEST_DISORDER_RUN_GT_0_5", "N_CTD_RESIDUES_SCORED",
           "N_NONSTANDARD_SUBSTITUTED", "CTD_SEQUENCE_STATUS", "LEGACY_EXPANSION_CALL", "CTD_LENGTH_DOMAIN_END",
           "PKINASE_I_EVALUE", "PKINASE_MODEL"]
    L.write("04_ampka_ctd/AMPKa_CTD_copy.tsv", hdr, rows, script=SELF,
            notes=["Primary CTD is residues C-terminal of the homolog of canonical PRKAA1-290 (anchor %s residue %d)." % (ANCHOR, ANCHOR_BOUNDARY),
                   "CTD_LENGTH_DOMAIN_END is the sensitivity definition from the frozen Pkinase envelope end and never replaces it.",
                   "CTD_SEQUENCE_STATUS: genome-derived records read CALLABLE; transcriptome-derived records cannot be",
                   "FULL_LENGTH_CONFIRMED because the frozen source metadata carry no per-protein completeness, so they read",
                   "POSITIVE_EXTENSION_ONLY (an extension of at least %d residues is observed) or C_TERMINUS_UNRESOLVED," % LEGACY_EXPANSION_LEN,
                   "and only genome-derived records enter the primary CTD length comparisons.",
                   "LEGACY_EXPANSION_CALL replicates the frozen rule and is a secondary replication variable only.",
                   "The last two rows are the human anchor and canonical references, kept as the numbering gate.",
                   "MAFFT: " + cmd + " | metapredict " + meta.__version__])
    L.write("04_ampka_ctd/AMPKa_CTD_disorder_long.tsv",
            ["TARGET", "SPECIES", "PROTEIN_ID", "RESIDUE_POSITION", "RESIDUE", "DISORDER"], long_rows, script=SELF)
    ref = {r[2]: r for r in rows if r[0] == "REFERENCE"}
    for k in (ANCHOR, CANON):
        got = ref.get(k, [None] * 15)[14]
        gates.append(["%s CTD length reproduces the frozen P03 value %d" % (k, HUMAN_CTD_EXPECTED),
                      "PASS" if got == HUMAN_CTD_EXPECTED else "FAIL", str(got)])
    L.write("10_validation/AMPKa_boundary_gates.tsv", ["GATE", "RESULT", "DETAIL"], gates, script=SELF)
    L.record_outputs("x2_alpha", SELF)
    nc = sum(1 for r in rows if r[13] == "yes" and r[0] != "REFERENCE")
    print("AMPKa: %d copies, %d with a callable anchor, %d CTD residues scored"
          % (len(rows) - 2, nc, len(long_rows)))
    if any(g[1] == "FAIL" for g in gates):
        sys.exit("REFUSED: a boundary gate failed: %s" % [g for g in gates if g[1] == "FAIL"])


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
