# >>> deposit header
# Program:  HyPhy (GARD)
# Version:  2.5.101
# Evidence: banner, env
# Methods:  Supplementary Methods, Phylogenetic analyses
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the recombination screen
# <<< deposit header

# Line 1 is the amino-acid screen, line 2 the codon screen. docs/SOFTWARE.md records the
# analysis flags; --alignment is written here as a placeholder for the input argument.
# The GARD log prints only the analysis version; 2.5.101 comes from the job script's
# banner and the environment record.

hyphy gard --type amino-acid --model LG --rv GDD --rate-classes 4 --alignment <alignment>
hyphy gard --type codon --code Universal --rv GDD --rate-classes 4 --alignment <alignment>
