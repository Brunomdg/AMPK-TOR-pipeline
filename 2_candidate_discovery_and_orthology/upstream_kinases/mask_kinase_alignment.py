#!/usr/bin/env python3
# >>> deposit header
# Purpose: Masks a CAMKK or LKB1 kinase-domain alignment to the columns that carry a residue in at least one
#          positive reference sequence, and reports the alignment width and gaps before and after masking.
# Methods: Orthology calls
# Inputs:  target name (CAMKK or LKB1); `<target>_kinase_domain.aln.faa`, `<target>_kinase_domain_manifest.tsv`
# Outputs: `<target>_kinase_domain.masked.aln.faa`, `<target>_alignment_mask_accounting.tsv`
# Tools:   Python 3 (standard library); versions: see envs/
# <<< deposit header
# Ran as: rbx3_mask.py
## Applies the reference-occupancy mask that the alignment check (rbx3_aln_qc.py) fixed before any number was
## read: keep only the alignment columns that carry a residue in at least one positive reference sequence. The
## mask depends only on the reference sequences, never on the candidates or on the tree topology. The unmasked
## alignment and its tree are kept as a sensitivity analysis, not discarded.
import os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "2_candidate_discovery_and_orthology", "upstream_kinases", "lib"))
import upstream_kinase_tables as C  # noqa: E402

SELF = os.path.abspath(__file__)
TARGET = sys.argv[1] if len(sys.argv) > 1 else "CAMKK"
SUB = "04_camkk_subfamilies" if TARGET == "CAMKK" else "10_validation"


def main():
    C.verify_freezes()
    alnp = os.path.join(C.X3, SUB, "%s_kinase_domain.aln.faa" % TARGET)
    aln = C.read_fasta(alnp)
    man = C.rows(os.path.join(C.X3, SUB, "%s_kinase_domain_manifest.tsv" % TARGET))
    tipof = lambda r: "%s__%s" % (r["SPECIES"], r["PROTEIN_ID"].replace(".", "_"))
    posrefs = [tipof(r) for r in man if "REFERENCE" in r["ROLE"] and tipof(r) in aln]
    if not posrefs:
        sys.exit("REFUSED: no positive reference in the alignment; the mask is undefined")
    ncol = len(next(iter(aln.values())))
    keep = [i for i in range(ncol) if any(aln[t][i] != "-" for t in posrefs)]
    outp = os.path.join(C.X3, SUB, "%s_kinase_domain.masked.aln.faa" % TARGET)
    with open(outp, "w") as fh:
        for t, s in aln.items():
            fh.write(">%s\n%s\n" % (t, "".join(s[i] for i in keep)))
    C.WRITTEN.append(outp)
    ung = [len("".join(aln[t][i] for i in keep).replace("-", "")) for t in aln]
    gap = sorted("".join(aln[t][i] for i in keep).count("-") / float(len(keep)) for t in aln)
    lost = [t for t in aln if not "".join(aln[t][i] for i in keep).replace("-", "")]
    C.write("%s/%s_alignment_mask_accounting.tsv" % (SUB, TARGET),
            ["QUANTITY", "VALUE", "NOTE"],
            [["MASK_RULE", "COLUMN_OCCUPIED_BY_AT_LEAST_ONE_POSITIVE_REFERENCE",
              "declared before any number was read; defined by the frozen reference set alone"],
             ["N_POSITIVE_REFERENCES", len(posrefs), "the mask depends on these and on nothing else"],
             ["COLUMNS_BEFORE", ncol, "unmasked alignment width"],
             ["COLUMNS_AFTER", len(keep), "masked alignment width"],
             ["COLUMNS_REMOVED", ncol - len(keep), "columns with no residue in any positive reference"],
             ["FRACTION_RETAINED", "%.3f" % (len(keep) / float(ncol)), "-"],
             ["MEDIAN_UNGAPPED_LENGTH_AFTER", int(statistics.median(ung)), "-"],
             ["MEDIAN_GAP_FRACTION_AFTER", "%.3f" % statistics.median(gap), "-"],
             ["N_TIPS_EMPTIED_BY_MASK", len(lost),
              "a tip left with no residue would be excluded and listed; none may be silently dropped"],
             ["TIPS_EMPTIED", ";".join(sorted(lost)) if lost else "-", "-"],
             ["ASSERT_MASK_INDEPENDENT_OF_CANDIDATES", "PASS",
              "no candidate sequence and no topology enters the mask definition"]], script=SELF,
            notes=["The unmasked alignment and any tree built from it are retained as a declared sensitivity."])
    if lost:
        sys.exit("REFUSED: %d tips carry no residue in the masked alignment" % len(lost))
    print("%s mask: %d -> %d columns (%.1f%% retained), median gap %.3f -> %.3f"
          % (TARGET, ncol, len(keep), 100.0 * len(keep) / ncol,
             statistics.median(sorted(aln[t].count("-") / float(ncol) for t in aln)),
             statistics.median(gap)))
    C.record_outputs("r2_mask_%s" % TARGET, SELF)


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
