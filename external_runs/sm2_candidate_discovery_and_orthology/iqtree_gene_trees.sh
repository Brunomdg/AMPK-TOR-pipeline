# >>> deposit header
# Program:  IQ-TREE
# Version:  3.1.3
# Evidence: env, banner, log
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  one gene tree per target, and the site-heterogeneous cross-check
# <<< deposit header

# Line 1 is the ModelFinder search. Line 2 is the site-heterogeneous cross-check.

iqtree3 -m MFP -B 1000 --alrt 1000 --seed 488761235
iqtree3 -m LG+C60+F+G -ft <ModelFinder tree> -B 1000 --alrt 1000 --seed 488761235
