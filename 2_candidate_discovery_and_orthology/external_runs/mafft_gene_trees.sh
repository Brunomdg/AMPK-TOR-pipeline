# >>> deposit header
# Program:  MAFFT
# Version:  7.526
# Evidence: env, banner
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  one alignment per target, for the gene trees
# <<< deposit header

# Line 1 is G-INS-i, as van Dam et al. (2011). Line 2 was used for the one tip set above
# the G-INS-i size limit.

mafft --globalpair --maxiterate 1000 --anysymbol
mafft --auto
