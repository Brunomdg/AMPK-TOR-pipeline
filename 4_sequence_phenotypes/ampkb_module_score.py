#!/usr/bin/env python3
# >>> deposit header
# Purpose: Scores the AMPKβ carbohydrate-binding module of every retained AMPKβ copy at five binding positions and
#          seven structural positions (identity to human PRKAB1), under a pairwise (primary) and a joint alignment,
#          with the identity of the other module positions as background, per copy and per species.
# Methods: Sequence phenotypes of the retained subunits
# Inputs:  `AMPKb_CBM_reference.tsv` (from job y0); `CLEAN_literature_controls.tsv` (human PRKAB1 NP_006244.2); the
#          retained AMPKβ copies and species states; the proteome FASTA files
# Outputs: in `05_ampkb_cbm`: `AMPKb_sequence_manifest.tsv` (one row per copy), the pairwise and joint alignments,
#          `AMPKb_CBM_copy.tsv`, `AMPKb_contact_long.tsv` and `AMPKb_CBM_species.tsv`
# Tools:   MAFFT 7.526 [env record] from the orthofinder3 environment (--globalpair --maxiterate 1000); Python 3;
#          versions: see envs/
# <<< deposit header
## Job y4 of the second chain. The AMPKβ carbohydrate-binding module, scored on every accepted AMPKβ ortholog copy
## against the human PRKAB1 reference (NP_006244.2), which carries the numbering of the rat structure 1Z0N at
## offset zero. The primary alignment is each copy against the reference, pairwise; one joint alignment of all
## copies is a declared sensitivity. The score measures conservation of the module's binding positions, not a
## general glycogen-sensing phenotype.
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
REF_ID = "NP_006244.2"
CBM_FROM, CBM_TO = 78, 161        # AMPK1_CBM domain envelope on the human reference, from the genome-wide HMMER domain table
HIGH_PROV = C.HIGH_PROV


def norm(acc):
    return acc.replace(":", "_").replace("@", "_")


def reference():
    p = os.path.join(C.R1, "02_reference_data/CLEAN_literature_controls.tsv")
    for r in C.rows(p):
        if r.get("PROTEIN_ID") == REF_ID and r.get("TARGET") == "AMPKb":
            return (r.get("SEQUENCE") or "").strip().rstrip("*")
    sys.exit("REFUSED: the frozen literature controls carry no %s sequence" % REF_ID)


def sites():
    out = {}
    for r in C.rows(os.path.join(C.X, "05_ampkb_cbm/AMPKb_CBM_reference.tsv")):
        out.setdefault(r["SITESET"], []).append((r["SITE_ID"], int(r["REF_POSITION"]), r["REF_AA"]))
    for k in out:
        out[k].sort(key=lambda x: x[1])
    if len(out.get("CBM_CORE5", [])) != 5 or len(out.get("CBM_STRUCTURAL7", [])) != 7:
        sys.exit("REFUSED: the frozen CBM panel is not 5 and 7 positions")
    return out


def extract():
    t2 = C.read_frozen_t2()
    copies, man = {}, []
    by = {}
    for r in C.t1_orthologs("AMPKb"):
        by.setdefault(r["SPECIES"], []).append(r)
    cache = {}
    for sp in sorted(by):
        path = os.path.join(C.PROT, sp + ".faa")
        if path not in cache:
            cache.clear()
            pool = S.read_fasta(path) if os.path.exists(path) else {}
            cache[path] = (pool, {norm(a): a for a in pool})
        pool, nidx = cache[path]
        for r in sorted(by[sp], key=lambda x: x["PROTEIN_ID"]):
            acc = r["PROTEIN_ID"]
            seq, how = pool.get(acc), "EXACT"
            if seq is None and norm(acc) in nidx:
                seq, how = pool[nidx[norm(acc)]], "NORMALISED_LABEL_PUNCTUATION"
            cell = t2.get(("AMPKb", sp), {})
            if seq is None:
                man.append([sp, acc, "SEQUENCE_UNAVAILABLE", 0, cell.get("SOURCE_TYPE", "-"), r["DECISION_SOURCE"],
                            r["TREE_CALL"], r["MODEL_SENSITIVE"], len(by[sp]), "NONE"])
                continue
            key = "%s|%s" % (sp, acc)
            copies[key] = seq
            man.append([sp, acc, S.seq_sha(seq), len(seq), cell.get("SOURCE_TYPE", "-"), r["DECISION_SOURCE"],
                        r["TREE_CALL"], r["MODEL_SENSITIVE"], len(by[sp]), how])
    C.write("05_ampkb_cbm/AMPKb_sequence_manifest.tsv",
            ["SPECIES", "PROTEIN_ID", "SEQUENCE_SHA256", "SEQUENCE_LENGTH", "SOURCE_TYPE", "DECISION_SOURCE",
             "TREE_CALL", "MODEL_SENSITIVE", "N_AMPKb_COPIES_SPECIES", "ID_MATCH"], man, script=SELF,
            notes=["Every frozen R1 AMPKb ORTHOLOG copy. No raw hit enters. A copy that resolves to no sequence",
                   "reads SEQUENCE_UNAVAILABLE and fails the freeze gate rather than disappearing."])
    return copies, {r[0] + "|" + r[1]: r for r in man}


def score_set(ref_aln, qry_aln, colmap, panel):
    per, n_call, n_exact = [], 0, 0
    for sid, pos, aa in panel:
        st, ex, obs = S.position_call(ref_aln, qry_aln, colmap[pos])
        if st in S.CALLABLE_STATES:
            n_call += 1
            n_exact += ex
        per.append((sid, pos, aa, obs, st, ex))
    return per, n_call, n_exact


def main():
    C.verify_freezes()
    ref = reference()
    panel = sites()
    bad = [(sid, pos, aa, ref[pos - 1]) for sid, pos, aa in panel["CBM_STRUCTURAL7"] if ref[pos - 1] != aa]
    if bad:
        sys.exit("REFUSED: the CBM numbering gate failed on the frozen reference: %s" % bad)
    copies, man = extract()
    allpos = {p for _s, p, _a in panel["CBM_STRUCTURAL7"]}
    dom = C.rows(os.path.join(C.R1, "03_evidence_by_method/hmmer_genomewide/HMMER_EVIDENCE_ROOTS.tsv"), required=False)

    jpath = os.path.join(C.X, "05_ampkb_cbm/AMPKb_CBM_alignment_joint.faa")
    jcmd, jaln = S.mafft(dict([(REF_ID, ref)] + sorted(copies.items())), jpath)
    C.WRITTEN.append(jpath)
    jcol = S.ref_column_map(jaln[REF_ID])
    pw_dir = os.path.join(C.X, "05_ampkb_cbm/pairwise")
    os.makedirs(pw_dir, exist_ok=True)

    long_rows, copy_rows, pcmd = [], [], ""
    for key in sorted(copies):
        sp, acc = key.split("|", 1)
        m = man[key]
        pcmd, pa = S.mafft({REF_ID: ref, key: copies[key]}, os.path.join(pw_dir, key.replace("|", "__") + ".faa"))
        pcol = S.ref_column_map(pa[REF_ID])
        row = [sp, acc, m[2], m[5], m[7], m[4], m[8], "PAIRWISE_HUMAN_PRKAB1"]
        scores = {}
        for tag, rref, rq, cmapx in (("PAIRWISE", pa[REF_ID], pa[key], pcol), ("JOINT", jaln[REF_ID], jaln[key], jcol)):
            for setname in ("CBM_CORE5", "CBM_STRUCTURAL7"):
                per, nc, ne = score_set(rref, rq, cmapx, panel[setname])
                scores[(tag, setname)] = (nc, ne)
                if setname == "CBM_STRUCTURAL7":
                    for sid, pos, aa, obs, st, ex in per:
                        long_rows.append([sp, acc, tag, sid, pos, aa, obs, st,
                                          "yes" if st in S.CALLABLE_STATES else "no", ex])
        for tag in ("PAIRWISE", "JOINT"):
            nc, ne = scores[(tag, "CBM_CORE5")]
            row += [nc, ne, (ne if nc == 5 else "NA")]
            nc7, ne7 = scores[(tag, "CBM_STRUCTURAL7")]
            row += [nc7, ne7, (ne7 if nc7 == 7 else "NA")]
        p5, j5 = row[8 + 2], row[8 + 8]     # CORE5 scores under each procedure
        row += ["NA" if "NA" in (str(p5), str(j5)) else p5 - j5,
                "yes" if (p5 != "NA" and p5 == j5) else ("no" if "NA" not in (str(p5), str(j5)) else "NOT_COMPARABLE")]
        nb = nbi = ncov = 0
        for pos in range(CBM_FROM, CBM_TO + 1):
            st, ex, _a = S.position_call(pa[REF_ID], pa[key], pcol[pos])
            if pa[key][pcol[pos]] != S.GAP:
                ncov += 1
            if pos in allpos:
                continue
            if st in S.CALLABLE_STATES:
                nb += 1
                nbi += ex
        bgi = (nbi / float(nb)) if nb else None
        c5 = scores[("PAIRWISE", "CBM_CORE5")]
        frac = (c5[1] / float(c5[0])) if c5[0] else None
        row += [nb, nbi, "%.4f" % bgi if bgi is not None else "NA",
                "%.4f" % frac if frac is not None else "NA",
                "%.4f" % (frac - bgi) if None not in (frac, bgi) else "NA",
                "%.3f" % (ncov / float(CBM_TO - CBM_FROM + 1)),
                "DOMAIN_ANNOTATED" if ncov >= 10 else "DOMAIN_REGION_UNRESOLVED"]
        copy_rows.append(row)

    hdr = ["SPECIES", "PROTEIN_ID", "SEQUENCE_SHA256", "DECISION_SOURCE", "MODEL_SENSITIVE", "SOURCE_TYPE",
           "N_AMPKb_COPIES_SPECIES", "PRIMARY_ALIGNMENT_METHOD"]
    for tag in ("PAIRWISE", "JOINT"):
        hdr += ["CBM_CORE5_N_CALLABLE_%s" % tag, "CBM_CORE5_N_EXACT_%s" % tag, "CBM_CORE5_SCORE_0_5_%s" % tag,
                "CBM_STRUCTURAL7_N_CALLABLE_%s" % tag, "CBM_STRUCTURAL7_N_EXACT_%s" % tag,
                "CBM_STRUCTURAL7_SCORE_0_7_%s" % tag]
    hdr += ["CBM_PAIRWISE_JOINT_DELTA", "ALIGNMENT_SENSITIVE", "CBM_NONCONTACT_N_CALLABLE",
            "CBM_NONCONTACT_N_IDENTICAL", "CBM_NONCONTACT_IDENTITY", "CBM_CORE5_FRACTION",
            "CBM_CONTACT_MINUS_BACKGROUND", "CBM_DOMAIN_SPAN_COVERAGE", "CBM_DOMAIN_ANNOTATION"]
    C.write("05_ampkb_cbm/AMPKb_CBM_copy.tsv", hdr, copy_rows, script=SELF,
            notes=["One row per frozen AMPKb ORTHOLOG copy; no copy is selected by its score. The primary score is",
                   "the pairwise mapping to the frozen human PRKAB1 reference and the joint alignment is the",
                   "declared alignment-procedure sensitivity; the two are never averaged.",
                   "A score is emitted only when every position of the panel is callable; a position beyond the",
                   "observed span is not callable and is never scored as functional erosion.",
                   "Background identity uses every callable non-contact position of the frozen CBM envelope %d-%d." % (CBM_FROM, CBM_TO),
                   "Domain annotation is annotation only and is never an orthology instrument.",
                   "MAFFT primary: " + pcmd, "MAFFT joint: " + jcmd])
    C.write("05_ampkb_cbm/AMPKb_contact_long.tsv",
            ["SPECIES", "PROTEIN_ID", "ALIGNMENT_METHOD", "SITE_ID", "REF_POSITION", "REF_AA", "OBSERVED_AA",
             "POSITION_STATE", "CALLABLE", "EXACT_MATCH"], long_rows, script=SELF)

    ix = {k: i for i, k in enumerate(hdr)}
    by = {}
    for r in copy_rows:
        by.setdefault(r[0], []).append(r)
    sp_rows = []
    for sp in sorted(by):
        cs = by[sp]
        num = lambda r, k: (float(r[ix[k]]) if r[ix[k]] not in ("NA", "") else None)
        a = C.species_copy_summary(cs, lambda r: num(r, "CBM_CORE5_SCORE_0_5_PAIRWISE"))
        j = C.species_copy_summary(cs, lambda r: num(r, "CBM_CORE5_SCORE_0_5_JOINT"))
        s7 = C.species_copy_summary(cs, lambda r: num(r, "CBM_STRUCTURAL7_SCORE_0_7_PAIRWISE"))
        hp = C.species_copy_summary(cs, lambda r: num(r, "CBM_CORE5_SCORE_0_5_PAIRWISE"),
                                    keep=lambda r: r[ix["DECISION_SOURCE"]] in HIGH_PROV)
        nm = C.species_copy_summary(cs, lambda r: num(r, "CBM_CORE5_SCORE_0_5_PAIRWISE"),
                                    keep=lambda r: r[ix["MODEL_SENSITIVE"]] != "yes")
        bg = C.species_copy_summary(cs, lambda r: num(r, "CBM_NONCONTACT_IDENTITY"))
        cmb = C.species_copy_summary(cs, lambda r: num(r, "CBM_CONTACT_MINUS_BACKGROUND"))
        f = lambda v, nd=4: "NA" if v is None else "%.*f" % (nd, v)
        sp_rows.append([sp, len(cs), a["N_CALLABLE"], f(a["MEDIAN"], 1), f(a["MIN"], 1), f(a["MAX"], 1),
                        a["N_UNIQUE"], a["DISCORDANT"], f(j["MEDIAN"], 1), f(s7["MEDIAN"], 1), f(hp["MEDIAN"], 1),
                        f(nm["MEDIAN"], 1), f(bg["MEDIAN"]), f(cmb["MEDIAN"]),
                        "yes" if all(r[ix["MODEL_SENSITIVE"]] == "yes" for r in cs) else "no",
                        "yes" if all(r[ix["DECISION_SOURCE"]] not in HIGH_PROV for r in cs) else "no"])
    C.write("05_ampkb_cbm/AMPKb_CBM_species.tsv",
            ["SPECIES", "N_AMPKb_COPIES", "N_CALLABLE_COPIES", "CBM_CORE5_MEDIAN", "CBM_CORE5_MIN", "CBM_CORE5_MAX",
             "N_UNIQUE", "COPY_DISCORDANT", "CBM_CORE5_MEDIAN_JOINT", "CBM_STRUCTURAL7_MEDIAN",
             "CBM_CORE5_MEDIAN_HIGHPROV", "CBM_CORE5_MEDIAN_NO_MODEL_SENSITIVE", "CBM_NONCONTACT_IDENTITY_MEDIAN",
             "CBM_CONTACT_MINUS_BACKGROUND_MEDIAN", "ALL_COPIES_MODEL_SENSITIVE", "R66_ONLY_SPECIES"], sp_rows,
            script=SELF, notes=["Median over callable copies is the primary species value; min and max define the",
                                "declared duplicate sensitivities and no copy is chosen to strengthen a relationship."])
    C.record_outputs("y4", SELF)
    nal = sum(1 for r in copy_rows if r[ix["ALIGNMENT_SENSITIVE"]] == "no")
    print("y4: %d AMPKb copies in %d species, %d with a full CORE5 score, %d copies alignment-concordant"
          % (len(copy_rows), len(sp_rows), sum(1 for r in copy_rows if r[ix["CBM_CORE5_SCORE_0_5_PAIRWISE"]] != "NA"), nal))


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
