# >>> deposit header
# Program:  DIAMOND
# Version:  2.2.5
# Evidence: env
# Methods:  Supplementary Methods, Proteomes, taxonomy and the species tree
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  marker hits for the species-tree build
# <<< deposit header

diamond blastp --very-sensitive --evalue 1e-5 --max-target-seqs 5
