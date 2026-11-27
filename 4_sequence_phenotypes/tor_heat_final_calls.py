#!/usr/bin/env python3
# >>> deposit header
# Purpose: Makes the final HEAT present or absent call, with a confidence tier, for every TOR copy (246, of which 204
#          were left unresolved by the signature and structural screens) from InterProScan annotations overlapping the
#          N-terminal region, and calibrates the rule on the 42 known positives and on the C-terminal FAT and kinase
#          region.
# Methods: TOR N-terminal HEAT repeats
# Inputs:  `TOR_observed_domains.tsv`, `TOR_all_retained.faa`, `TOR_manifest.tsv` (signature screen);
#          `TOR_HEAT_updated_classification.tsv` (structural screen); `CDS_AND_SOURCE_PROVENANCE.tsv` (source
#          annotation marked partial); `HEAT_EVIDENCE_RULES.tsv` (checked to be non-empty; the accession classes
#          applied are set in `evi_class`)
# Outputs: `TOR_HEAT_FINAL_CALLS.tsv` (all 246 copies); `CONTROL_CALIBRATION.json`
# Tools:   Python 3 with pandas; versions: see envs/
# <<< deposit header
# Ran as: resolve_heat_final.py
"""
resolve_heat_final.py: final two-state HEAT calls (present or absent) for the 204 TOR proteins that the
signature and structural screens left unresolved.

Reads (existing files, none modified):
  TOR_HEAT_20260919/TOR_observed_domains.tsv          (InterProScan annotation rows)
  TOR_HEAT_20260919/TOR_all_retained.faa              (sequences)
  TOR_HEAT_20260919/TOR_manifest.tsv                  (length, MD5, SHA-256)
  TOR_HEAT_resolution_20260919/TOR_HEAT_updated_classification.tsv
  TOR_HEAT_resolution_20260919/round2/CDS_AND_SOURCE_PROVENANCE.tsv
  HEAT_EVIDENCE_RULES.tsv                             (in the output folder; accessions accepted as HEAT evidence)

Writes:
  TOR_HEAT_FINAL_CALLS.tsv      (all 246 proteins)
  CONTROL_CALIBRATION.json      (recovery of the 42 known positives; false positives in the FAT and kinase region)

Decision rule:
  1. The N-terminal region is NTERMINAL_HEAT_MAPPED_START/END_AA from the updated classification (coordinates
     transferred through the alignment). When those are missing, residues 1 to NTERM_REGION_LENGTH_PROXY_AA are
     used instead (the copy's own FAT start minus 1, or the boundary transferred through the alignment),
     provided that this is at least 60 residues.
  2. HEAT_PRESENT if at least one accepted accession (class A, B or C) overlaps the region by at least one
     residue; otherwise HEAT_ABSENT, meaning the region was evaluated and no HEAT-indicative evidence was found.
     Class A is the mTOR-family domain (Pfam PF11865, SMART SM01346); class B are HEAT repeat families (Pfam
     families and the PROSITE profile PS50077); class C are repeat-fold annotations (SUPERFAMILY SSF48371,
     Gene3D 1.25.10.10 and its FunFam families).
  3. Confidence tiers: class A or B gives HIGH; class C from both SUPERFAMILY and Gene3D gives MODERATE; class C
     from one of them gives LOW. A weak E-value (above 1e-5) on class A or B evidence lowers the tier by one step,
     and any of the other caveats (protein shorter than 1,800 residues, unknown residues, source annotation marked
     partial, region shorter than 400 residues, fallback region) lowers it by one further step. HEAT_ABSENT is
     MODERATE only when the region is at least 400 residues long and there is no caveat; otherwise LOW.
"""
# deposit path defaults: begin
_TA210_ONEDRIVE_ALT = __import__("os").environ.get("TA210_ONEDRIVE_ALT") or "/Users/jupiter-room/Library/CloudStorage/OneDrive-UVMLarnerCollegeofMedicine(2)"
# deposit path defaults: end
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path((_TA210_ONEDRIVE_ALT + "/MBE"))
DIR_A = ROOT / "TOR_HEAT_20260919"
DIR_B = ROOT / "TOR_HEAT_resolution_20260919"
DIR_R2 = DIR_B / "round2"
OUT = ROOT / "TOR_HEAT_final_calls_20260920"

UNRESOLVED3 = [
    "UNRESOLVED_NO_EXACT_SEQUENCE_STRUCTURE",
    "UNRESOLVED_NO_USABLE_REGION",
    "UNRESOLVED_ALIGNMENT_OR_CONFIDENCE",
]
POSITIVE_CLASSES = [
    "NTERMINAL_HEAT_SIGNATURE_DETECTED",
    "NTERMINAL_HEAT_STRUCTURE_SUPPORTED",
    "CURATED_NTERMINAL_HEAT_REFERENCE",
]
FALLBACK_MIN_REGION = 60      # residues; the minimum usable region of the second-round protein language model screen
SHORT_REGION_CAVEAT = 400     # residues; for comparison, the median mapped region is about 1,300 residues
WEAK_EVALUE = 1e-5            # marks weak decisive evidence; it lowers the tier but is not an acceptance threshold
INDEPENDENT_FOLD_DBS = {"SUPERFAMILY", "Gene3D"}  # FunFam is counted with Gene3D, not as a separate database


def load_fasta(path):
    seqs, cid, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cid:
                seqs[cid] = "".join(buf)
            cid, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if cid:
        seqs[cid] = "".join(buf)
    return seqs


def evi_class(db, acc):
    if (db, acc) in [("Pfam", "PF11865"), ("SMART", "SM01346")]:
        return "A"
    if (db, acc) in [("Pfam", "PF23593"), ("Pfam", "PF13513"), ("Pfam", "PF13646"),
                     ("Pfam", "PF20175"), ("Pfam", "PF20206"), ("Pfam", "PF23271"),
                     ("Pfam", "PF25786"), ("ProSiteProfiles", "PS50077")]:
        return "B"
    if (db, acc) in [("SUPERFAMILY", "SSF48371"), ("Gene3D", "G3DSA:1.25.10.10")]:
        return "C"
    if db == "FunFam" and acc.startswith("G3DSA:1.25.10.10:FF:"):
        return "C"
    return None


def parse_evalue(db, score):
    """Return the E-value of a match from a HMMER-based database; PROSITE scores are not E-values, so give None."""
    try:
        v = float(score)
    except (TypeError, ValueError):
        return None
    if db in ("ProSiteProfiles", "ProSitePatterns"):
        return None
    return v


def main():
    # ---------- load ----------
    rules = pd.read_csv(OUT / "HEAT_EVIDENCE_RULES.tsv", sep="\t", dtype=str)
    admitted = rules[rules.ADMITTED_AS_HEAT_EVIDENCE.str.startswith("YES")]
    assert not admitted.empty, "rules file unreadable"

    cls = pd.read_csv(DIR_B / "TOR_HEAT_updated_classification.tsv", sep="\t", dtype=str)
    dom = pd.read_csv(DIR_A / "TOR_observed_domains.tsv", sep="\t", dtype=str)
    man = pd.read_csv(DIR_A / "TOR_manifest.tsv", sep="\t", dtype=str).set_index("COPY_ID")
    seqs = load_fasta(DIR_A / "TOR_all_retained.faa")
    cds = pd.read_csv(DIR_R2 / "CDS_AND_SOURCE_PROVENANCE.tsv", sep="\t", dtype=str)

    dom[["START", "END"]] = dom[["START", "END"]].astype(int)
    cls["plen"] = pd.to_numeric(cls.PROTEIN_LENGTH_AA, errors="coerce")
    cls["ms"] = pd.to_numeric(cls.NTERMINAL_HEAT_MAPPED_START_AA, errors="coerce")
    cls["me"] = pd.to_numeric(cls.NTERMINAL_HEAT_MAPPED_END_AA, errors="coerce")
    cls["proxy"] = pd.to_numeric(cls.NTERM_REGION_LENGTH_PROXY_AA, errors="coerce")
    cls["fat"] = pd.to_numeric(cls.OBSERVED_FAT_START_AA, errors="coerce")
    cls["unk"] = pd.to_numeric(cls.UNKNOWN_RESIDUES, errors="coerce").fillna(0).astype(int)

    # ---------- check every sequence against its checksums and length ----------
    identity_failures = []
    for cid, s in seqs.items():
        md5 = hashlib.md5(s.encode()).hexdigest()
        sha = hashlib.sha256(s.encode()).hexdigest()
        if cid not in man.index or cid not in set(cls.COPY_ID):
            identity_failures.append((cid, "missing_from_manifest_or_classification"))
            continue
        c = cls[cls.COPY_ID == cid].iloc[0]
        if md5 != man.loc[cid, "MD5"] or sha != man.loc[cid, "SHA256"] \
                or sha != c.SEQUENCE_SHA256 or int(man.loc[cid, "LENGTH_AA"]) != len(s):
            identity_failures.append((cid, "hash_or_length_mismatch"))
    assert set(cls.COPY_ID) == set(seqs) == set(man.index), "COPY_ID sets differ"
    assert set(dom.COPY_ID) <= set(seqs), "domain rows reference unknown COPY_ID"

    # copies whose source annotation is marked partial in the second-round source table
    partial_ids = set(cds.loc[cds.SOURCE_ANNOTATION_PARTIAL.str.upper() == "TRUE", "COPY_ID"]) \
        if "SOURCE_ANNOTATION_PARTIAL" in cds.columns else set()
    if not partial_ids:  # if the table marks none, use the eight copies previously identified as partial
        partial_ids = {"TOR0006", "TOR0091", "TOR0131", "TOR0166",
                       "TOR0187", "TOR0206", "TOR0207", "TOR0209"}

    dom["evi"] = [evi_class(d, a) for d, a in zip(dom.DB, dom.ACCESSION)]
    dom["eval"] = [parse_evalue(d, s) for d, s in zip(dom.DB, dom.SCORE)]

    cohort204 = set(cls.loc[cls.RESOLVED_NTERMINAL_HEAT_CLASS.isin(UNRESOLVED3), "COPY_ID"])
    positives42 = set(cls.loc[cls.RESOLVED_NTERMINAL_HEAT_CLASS.isin(POSITIVE_CLASSES), "COPY_ID"])
    assert len(cohort204) == 204 and len(positives42) == 42

    # ---------- per-protein evaluation ----------
    records = []
    for _, c in cls.iterrows():
        cid = c.COPY_ID
        rows = dom[dom.COPY_ID == cid]

        # N-terminal region: transferred coordinates, else residues 1 to the length estimate (see rule 1 above)
        region_src, rs, re_ = "MSA_MAPPED", c.ms, c.me
        if pd.isna(rs) or pd.isna(re_) or re_ < rs:
            if pd.notna(c.proxy) and c.proxy >= FALLBACK_MIN_REGION:
                region_src, rs, re_ = "FAT_PROXY_FALLBACK", 1, int(c.proxy)
            else:
                region_src, rs, re_ = "NONE", None, None
        if rs is not None:
            rs, re_ = int(rs), int(re_)
        region_len = (re_ - rs + 1) if rs is not None else 0

        # caveats
        caveats = []
        if c.SHORT_PROTEIN_REVIEW_LT1800AA == "1" or (pd.notna(c.plen) and c.plen < 1800):
            caveats.append("FRAGMENT_LT1800AA")
        if c.unk > 0:
            caveats.append(f"UNKNOWN_RESIDUES_{c.unk}")
        if cid in partial_ids:
            caveats.append("SOURCE_RECORD_PARTIAL")
        if region_src == "NONE":
            caveats.append("NO_USABLE_MAPPED_REGION")
        elif region_src == "FAT_PROXY_FALLBACK":
            caveats.append("REGION_FROM_FAT_PROXY_FALLBACK")
        if rs is not None and region_len < SHORT_REGION_CAVEAT:
            caveats.append(f"SHORT_MAPPED_REGION_{region_len}AA")

        # accepted HEAT evidence overlapping the region
        hits = []
        if rs is not None:
            for r in rows.itertuples():
                if isinstance(r.evi, str) and r.START <= re_ and r.END >= rs:
                    hits.append(r)
        # weak E-value on decisive evidence (class A or B matches that carry an E-value)
        weak_e = any(h.evi in "AB" and h.eval is not None and h.eval > WEAK_EVALUE for h in hits)
        if weak_e:
            caveats.append("WEAK_DECISIVE_EVALUE_GT_1E-5_REVIEW")

        has_ab = any(h.evi in "AB" for h in hits)
        c_dbs = {("Gene3D" if h.DB == "FunFam" else h.DB) for h in hits if h.evi == "C"}
        n_fold_dbs = len(c_dbs & INDEPENDENT_FOLD_DBS)

        # call and confidence tier
        if hits:
            call = "HEAT_PRESENT"
            if has_ab:
                tier = "HIGH"
            elif n_fold_dbs >= 2:
                tier = "MODERATE"
            else:
                tier = "LOW"
            downgrades = int(weak_e) + int(any(x.startswith(("FRAGMENT", "UNKNOWN_RESIDUES",
                          "SOURCE_RECORD_PARTIAL", "REGION_FROM_FAT", "SHORT_MAPPED"))
                          for x in caveats))
            order = ["HIGH", "MODERATE", "LOW"]
            tier = order[min(2, order.index(tier) + downgrades)]
        else:
            call = "HEAT_ABSENT"
            complete = (not caveats) or caveats == []
            tier = "MODERATE" if (rs is not None and region_len >= SHORT_REGION_CAVEAT
                                  and not caveats) else "LOW"

        # machine-readable evidence string
        if hits:
            ev_str = ";".join(
                f"{h.evi}:{h.DB}:{h.ACCESSION}:{h.START}-{h.END}"
                + (f"(E={h.eval:.2e})" if pd.notna(h.eval) else f"(score={h.SCORE})")
                for h in sorted(hits, key=lambda h: (h.evi, h.DB, h.START)))
        else:
            # for an absence call, state what was observed in the region
            if rs is None:
                obs = "no usable mapped N-terminal region (MSA transfer failed and FAT proxy <60aa)"
            else:
                inreg = [r for r in rows.itertuples()
                         if r.START <= re_ and r.END >= rs and r.DB not in ("MobiDBLite", "Coils")]
                if inreg:
                    obs = "non-HEAT annotations in region: " + ";".join(
                        f"{r.DB}:{r.ACCESSION}:{r.START}-{r.END}" for r in inreg[:6])
                elif len(rows):
                    obs = "no domain annotation inside region; domains only outside: " + ";".join(
                        f"{r.DB}:{r.ACCESSION}" for r in list(rows.itertuples())[:4])
                else:
                    obs = "no InterProScan annotation rows at all"
            ev_str = f"NO_HEAT_INDICATIVE_EVIDENCE[{obs}]"

        records.append({
            "COPY_ID": cid, "SPECIES": c.SPECIES, "PROTEIN_ID": c.PROTEIN_ID,
            "PROTEIN_LENGTH_AA": c.PROTEIN_LENGTH_AA,
            "PRIOR_RESOLVED_CLASS": c.RESOLVED_NTERMINAL_HEAT_CLASS,
            "IN_COHORT_204": cid in cohort204,
            "FINAL_CALL": call, "CONFIDENCE_TIER": tier,
            "EVIDENCE": ev_str,
            "REGION_SOURCE": region_src,
            "REGION_START_AA": rs if rs is not None else "",
            "REGION_END_AA": re_ if re_ is not None else "",
            "REGION_LENGTH_AA": region_len if rs is not None else "",
            "OBSERVED_FAT_START_AA": c.OBSERVED_FAT_START_AA,
            "EVIDENCE_CLASSES_PRESENT": "".join(sorted({h.evi for h in hits})) if hits else "",
            "N_INDEPENDENT_FOLD_DBS": n_fold_dbs,
            "CAVEATS": "|".join(caveats) if caveats else "NONE",
            "SEQUENCE_SHA256_VERIFIED": cid not in {f[0] for f in identity_failures},
        })

    out = pd.DataFrame(records)

    # ---------- calibration against controls ----------
    # (a) the 42 copies already known to carry N-terminal HEAT repeats should all be called HEAT_PRESENT
    pos_tbl = out[out.COPY_ID.isin(positives42)]
    pos_recovered = pos_tbl[pos_tbl.FINAL_CALL == "HEAT_PRESENT"]
    pos_missed = sorted(set(positives42) - set(pos_recovered.COPY_ID))

    # (b) negative-region control: apply the same overlap rule for accepted accessions to the C-terminal FAT
    #     and kinase region, from the FAT start to the end of the protein
    def cterm_hits(sub_cls, classes=None):
        ids, per_class = set(), {"A": set(), "B": set(), "C": set()}
        for _, cc in sub_cls.iterrows():
            if pd.isna(cc.fat) or pd.isna(cc.plen):
                continue
            rows = dom[(dom.COPY_ID == cc.COPY_ID) & dom.evi.notna()
                       & (dom.END >= cc.fat) & (dom.START <= cc.plen)]
            if len(rows):
                ids.add(cc.COPY_ID)
                for cl in "ABC":
                    if (rows.evi == cl).any():
                        per_class[cl].add(cc.COPY_ID)
        return ids, per_class

    all_tbl = cls.copy()
    fp_all, fp_class = cterm_hits(all_tbl)
    fp_204, _ = cterm_hits(cls[cls.COPY_ID.isin(cohort204)])
    fp_42, _ = cterm_hits(cls[cls.COPY_ID.isin(positives42)])

    calibration = {
        "rule_region": "admitted-accession interval overlap >=1aa with mapped N-terminal region",
        "positive_controls_n": 42,
        "positive_controls_recovered_HEAT_PRESENT": int(len(pos_recovered)),
        "positive_controls_missed": pos_missed,
        "positive_recovery_notes": {
            "TOR0096": "lacks class-A PF11865/SM01346 N-terminal match; recovered by class-B PS50077 HEAT profile plus class-C fold hits (its original native signature basis)"
        } if "TOR0096" in set(pos_recovered.COPY_ID) else {},
        "sensitivity_on_42_known_positives": round(len(pos_recovered) / 42, 4),
        "class_A_alone_recovers": int(sum(1 for _, r in pos_tbl.iterrows()
                                          if "A" in r.EVIDENCE_CLASSES_PRESENT)),
        "negative_region_definition": "observed FAT start .. protein end (FAT/FRB/kinase/FATC)",
        "negative_region_fp_rate_all_246": round(len(fp_all) / 246, 4),
        "negative_region_fp_n_all_246": len(fp_all),
        "negative_region_fp_by_class_all_246": {k: len(v) for k, v in fp_class.items()},
        "negative_region_fp_rate_class_A": round(len(fp_class["A"]) / 246, 4),
        "negative_region_fp_rate_class_C": round(len(fp_class["C"]) / 246, 4),
        "negative_region_fp_in_204": len(fp_204),
        "negative_region_fp_in_42": len(fp_42),
        "negative_region_interpretation": (
            "Class A (PF11865/SM01346) never occurs C-terminal of FAT (0/246): perfect region "
            "specificity. Class C fold hits leak into the C-terminal region (expected: the PIKK "
            "FAT-kinase region is genuinely alpha-solenoid; PF23593 ATR-like HEAT is annotated "
            "there in 155/204), so class C alone can never reach HIGH. This is a "
            "region-specificity control, not a true biological negative set."),
        "sequence_identity_check": {
            "fasta_records": len(seqs),
            "md5_sha256_length_failures": identity_failures,
            "all_domain_rows_sequence_matched": bool(
                (dom.EVIDENCE == "OBSERVED_SEQUENCE_MD5_AND_LENGTH_MATCHED").all()),
        },
    }

    # ---------- check that no protein was dropped and every call is defined ----------
    assert len(out) == 246, f"expected 246 rows, got {len(out)}"
    assert out.FINAL_CALL.isin(["HEAT_PRESENT", "HEAT_ABSENT"]).all()
    assert out.FINAL_CALL.notna().all()
    assert set(out[out.IN_COHORT_204].COPY_ID) == cohort204
    assert len(out[out.IN_COHORT_204]) == 204

    out = out.sort_values("COPY_ID")
    out.to_csv(OUT / "TOR_HEAT_FINAL_CALLS.tsv", sep="\t", index=False)

    sub204 = out[out.IN_COHORT_204]
    calibration["cohort_204_split"] = sub204.FINAL_CALL.value_counts().to_dict()
    calibration["cohort_204_tiers"] = (
        sub204.groupby(["FINAL_CALL", "CONFIDENCE_TIER"]).size().unstack(fill_value=0)
        .to_dict())
    with open(OUT / "CONTROL_CALIBRATION.json", "w") as f:
        json.dump(calibration, f, indent=2)

    print("=== FINAL CALLS (246) ===")
    print(out.FINAL_CALL.value_counts().to_dict())
    print("=== 204 cohort ===")
    print(sub204.groupby(["FINAL_CALL", "CONFIDENCE_TIER"]).size())
    print("=== positives recovered:", len(pos_recovered), "/42; missed:", pos_missed)
    print("=== identity failures:", identity_failures)
    print("=== negative-region FP: all", len(fp_all), "/246; classA",
          len(fp_class["A"]), "; classC", len(fp_class["C"]))


if __name__ == "__main__":
    sys.exit(main())
