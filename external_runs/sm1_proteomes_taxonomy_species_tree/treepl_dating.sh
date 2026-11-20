# >>> deposit header
# Program:  treePL
# Version:  not recorded
# Evidence: not recorded
# Methods:  Supplementary Methods, Proteomes, taxonomy and the species tree
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the second alternative chronogram
# <<< deposit header

# No command line survives. What the records show is the procedure: prime, then
# cross-validation of the smoothing value from 1e6 to 1e-12 in steps of 0.1, then a
# thorough run, over 7,933 sites. treePL is driven by a configuration file, and that
# file is not in the records; the line above is a placeholder for its invocation, not a
# transcription. The version is not recorded either.

treePL <configuration file>
