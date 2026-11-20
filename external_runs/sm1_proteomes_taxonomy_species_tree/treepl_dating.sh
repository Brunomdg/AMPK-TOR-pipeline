# >>> deposit header
# Program:  treePL
# Version:  1.0, built from source on 22 August 2026 from commit 1f5c99f088227c7c8dacdf343e5d1ab23405009b
# Evidence: banner, log, source clone
# Methods:  Supplementary Methods, Proteomes, taxonomy and the species tree
# Source:   the dating job log, which records the stages, their return codes and the chosen smoothing
# Outputs:  one dated tree per rooting and calibration set, under treepl/
# <<< deposit header

# treePL publishes no release number beyond the banner "treePL version 1.0", which its build log records.
# The version of record is therefore that banner plus the commit of the clone that was built.
#
# treePL is driven by a configuration file, and the configuration files are not in the records, so no
# flag string is reconstructed here. What the job log does record is the procedure and its outcome:
#
#   prime      returned 0, settings opt=2, optad=2, optcvad=5
#   cv         returned 0; best smoothing 1e-05 on calibration set A, 1e-04 on set B
#   thorough   returned 0, writing A2_setA.dated.tre and A2_setB.dated.tre
#
# Run on root A2 under both calibration sets. On set A it dropped two of the declared calibrations,
# Archaeplastida and Fungi_crown, because the most recent common ancestor subtends the wrong number of
# tips; on set B it dropped Bangio_Florideo and Eumetazoa for the same reason. One branch was floored in
# each. These are treePL's own drops and are not the drops of the chronos fit of record, which are three
# and are listed in ../../docs/DATING_GAP.md.
#
# Neither treePL tree was selected as the chronogram of record.

treePL <configuration file>
