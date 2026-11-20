# >>> deposit header
# Program:  R with ape chronos
# Version:  R 4.5.1 with the home-folder library, not the R 4.5.3 of the cluster steps
# Evidence: log
# Methods:  Supplementary Methods, Proteomes, taxonomy and the species tree
# Source:   the dating job log and the chronos run log; the script of record is lost
# Outputs:  the chronogram of record, published as species_chronogram_209.treefile
# <<< deposit header

The script of record is lost. It is not reconstructed here.

See `../../docs/DATING_GAP.md` for the two checksums, the calibrations that entered and were dropped, and
what the records do fix.

## What the logs record

The chronogram of record is `chronos/A2_setA.correlated_l0.1.tre`, md5
`301836d4f15956d8b14ffbfe420a64a9`, byte-identical to the published `species_chronogram_209.treefile`.

- Rooting A2, calibration set A, correlated-rates model, smoothing parameter 0.1.
- The smoothing value was chosen by leave-one-calibration-out cross-validation on the correlated model.
- The fit reports `conv=FALSE`. Its root age is 2,299.0 Ma.
- Of the eight calibrations declared in set A, **five entered the fit and three were dropped**, each drop
  logged with the tip counts that caused it.
- The sweep behind it: three rootings by two calibration sets, and within each the correlated, relaxed,
  discrete and clock models over smoothing values 0, 0.1, 1, 10, 100 and 1,000. 86 chronograms were
  compared before one was selected.

No command line survives. The settings above are the run's own log, not a transcription of a script.
