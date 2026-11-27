#!/usr/bin/env python3
# >>> deposit header
# Purpose: Alignment quality control on the CAMKK and LKB1 kinase-domain alignments, run and read before any
#          topology is read, so that an inflated alignment can never silently define a classification. Declares
#          the occupancy trigger that decides whether the reference-occupancy mask is applied.
# Methods: Orthology calls
# Inputs:  `<target>_kinase_domain.aln.faa` and `<target>_kinase_domain_manifest.tsv` for CAMKK and LKB1;
#          the CAMKK membership-tree crop `CAMKK/CAMKK_crop.faa`, for a width comparison only
# Outputs: `02_inputs/alignment_qc.tsv` (one row per target: tip and column counts, gap fractions, reference
#          occupancy, reference coverage, the frozen crop width and the QC verdict)
# Tools:   Python 3 (standard library) and the project module `rbx3_common`; versions: see envs/
# <<< deposit header
# The trigger below is the threshold that `rbx3_mask.py` and `rbx3_camkk_calls.py` both cite as fixed before any
# topology was read: at least half the alignment columns must carry a residue in at least one positive reference,
# or the lane stops and the deterministic reference-occupancy mask is applied. No column is ever trimmed because
# of how the tree looks.
## GOVERNED BY: RBX3_EXECUTION_CURRENT.md r2 QC gate. Alignment quality control, run and read BEFORE any
## topology is read, so an inflated alignment can never silently define a constituent classification. The
## deterministic mask rule and its trigger are declared here, before the numbers are seen: if fewer than half
## the columns carry a residue in any positive reference, the alignment is pathologically inflated and the lane
## stops for masking. No trimming is ever chosen by how the tree looks.
import os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rbx3_common as C  # noqa: E402

SELF = os.path.abspath(__file__)
OCCUPANCY_TRIGGER = 0.50     # declared before reading: fraction of columns occupied by any positive reference
## TEN-TARGET REFIT: the width of the CAMKK membership tree's crop, read from it; continuity comparison only
FROZEN_CROP_WIDTH = len(next(iter(C.read_fasta(os.path.join(C.PHYLO, "CAMKK/CAMKK_crop.faa")).values())))


def main():
    C.verify_freezes()
    out, stop = [], []
    for target, sub in (("CAMKK", "04_camkk_subfamilies"), ("LKB1", "10_validation")):
        alnp = os.path.join(C.X3, sub, "%s_kinase_domain.aln.faa" % target)
        if not os.path.exists(alnp):
            print("skip %s: alignment not written yet" % target)
            continue
        aln = C.read_fasta(alnp)
        man = C.rows(os.path.join(C.X3, sub, "%s_kinase_domain_manifest.tsv" % target))
        tipof = lambda r: "%s__%s" % (r["SPECIES"], r["PROTEIN_ID"].replace(".", "_"))
        posrefs = [tipof(r) for r in man if "REFERENCE" in r["ROLE"] and tipof(r) in aln]
        ncol = len(next(iter(aln.values())))
        if any(len(s) != ncol for s in aln.values()):
            sys.exit("REFUSED: %s alignment is ragged" % target)
        ung = [len(s.replace("-", "")) for s in aln.values()]
        gapf = sorted((s.count("-") / float(ncol)) for s in aln.values())
        q = lambda v, f: v[min(len(v) - 1, int(f * len(v)))]
        occ_all = sum(1 for i in range(ncol) if all(aln[t][i] != "-" for t in posrefs))
        occ_any = sum(1 for i in range(ncol) if any(aln[t][i] != "-" for t in posrefs))
        gap_all = ncol - occ_any
        cov = sorted(len(aln[t].replace("-", "")) / float(ncol) for t in posrefs)
        frac = occ_any / float(ncol)
        verdict = "OK" if frac >= OCCUPANCY_TRIGGER else "STOP_PATHOLOGICAL_INFLATION"
        if verdict != "OK":
            stop.append(target)
        out.append([target, len(aln), ncol, int(statistics.median(ung)), min(ung), max(ung),
                    "%.3f" % statistics.median(gapf), "%.3f" % q(gapf, 0.25), "%.3f" % q(gapf, 0.75),
                    occ_all, occ_any, gap_all, "%.3f" % frac, "%.3f" % cov[0],
                    "%.3f" % statistics.median(cov), len(posrefs),
                    FROZEN_CROP_WIDTH if target == "CAMKK" else "-",
                    "%.2f" % (ncol / float(FROZEN_CROP_WIDTH)) if target == "CAMKK" else "-", verdict])
    C.write("02_inputs/alignment_qc.tsv",
            ["TARGET", "N_TIPS", "ALIGNMENT_COLUMNS", "MEDIAN_UNGAPPED_LENGTH", "MIN_UNGAPPED_LENGTH",
             "MAX_UNGAPPED_LENGTH", "MEDIAN_GAP_FRACTION_PER_TIP", "Q1_GAP_FRACTION_PER_TIP",
             "Q3_GAP_FRACTION_PER_TIP", "N_COLUMNS_OCCUPIED_BY_ALL_REFERENCES",
             "N_COLUMNS_OCCUPIED_BY_ANY_POSITIVE_REFERENCE",
             "N_COLUMNS_WITH_REFERENCE_GAP_IN_ALL_POSITIVE_REFERENCES",
             "FRACTION_COLUMNS_OCCUPIED_BY_ANY_POSITIVE_REFERENCE", "REFERENCE_SEQUENCE_COVERAGE_MIN",
             "REFERENCE_SEQUENCE_COVERAGE_MEDIAN", "N_POSITIVE_REFERENCES", "FROZEN_CROP_WIDTH",
             "WIDTH_RATIO_TO_FROZEN_CROP", "QC_VERDICT"], out, script=SELF,
            notes=["Read before any topology. The declared trigger, fixed before these numbers were seen, is",
                   "that at least %.0f%% of alignment columns must carry a residue in at least one positive"
                   % (OCCUPANCY_TRIGGER * 100),
                   "reference; below that the alignment is pathologically inflated, the lane stops, and a",
                   "deterministic reference-occupancy mask is applied, exactly as the RELAX lane's gap",
                   "inflation was handled. The CAMKK membership-tree crop width of %d columns is a" % FROZEN_CROP_WIDTH,
                   "continuity comparison only: the new alignment is not required to equal it, and no column",
                   "is ever trimmed because of how the tree looks."])
    for r in out:
        print("  %-6s tips=%d cols=%d median_ungapped=%d median_gap=%s ref_occupancy=%s verdict=%s"
              % (r[0], r[1], r[2], r[3], r[6], r[12], r[18]))
    C.record_outputs("r2_qc", SELF)
    if stop:
        sys.exit("REFUSED: alignment QC stops these lanes before classification: %s" % stop)


if __name__ == "__main__":
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("REFUSED: not in Slurm")
    main()
