# >>> deposit header
# Program:  rate4site
# Version:  3.0.0, compiled from source
# Evidence: banner
# Methods:  Supplementary Methods, Sequence phenotypes of the retained subunits
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  per-site rates, for the comparison against burial-matched positions
# <<< deposit header

# Line 1 is empirical Bayes, with the assertion build. Line 2 is maximum likelihood with
# 16 rate categories, with the build without assertions. The reference sequence is
# Homo sapiens. docs/SOFTWARE.md records the mode flags; -s and -t are written here as
# placeholders for the alignment and tree arguments.

rate4site -ib -s <alignment> -t <tree>
rate4site -im -k 16 -s <alignment> -t <tree>
