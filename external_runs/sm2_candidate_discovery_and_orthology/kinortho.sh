# >>> deposit header
# Program:  KinOrtho
# Version:  not recorded
# Evidence: not recorded
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the AMPK-alpha ortholog sets of the orthology arm
# <<< deposit header

# Run for AMPK-alpha only. The version is not recorded.
# The Supplementary Methods records three variants: the published PF00069 query at E 1e-5,
# the query recut to the declared AMPK-alpha controls, and the BLAST default E 10. Only the
# flag string below survives from the job script; which flags carried each variant's E value
# is not recorded, so the three are named here rather than written as three command lines.

kinortho.py -f <full-length queries> -d <kinase-domain queries> -t 32
