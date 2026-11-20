# >>> deposit header
# Program:  KinOrtho
# Version:  not recorded
# Evidence: not recorded
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the AMPK-alpha ortholog sets of the orthology arm
# <<< deposit header

# Run for AMPK-alpha only. The version is not recorded.

kinortho.py -f <full-length queries> -d <kinase-domain queries> -t 32
