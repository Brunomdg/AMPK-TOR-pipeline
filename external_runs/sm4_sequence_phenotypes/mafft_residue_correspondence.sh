# >>> deposit header
# Program:  MAFFT
# Version:  7.526 (orthofinder3 environment)
# Evidence: env
# Methods:  Supplementary Methods, Sequence phenotypes of the retained subunits
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the coccidian residue-correspondence test
# <<< deposit header

# L-INS-i, for the coccidian residue-correspondence test only.

mafft --localpair --maxiterate 1000
