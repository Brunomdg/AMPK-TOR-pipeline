# >>> deposit header
# Program:  OrthoFinder
# Version:  3.1.5 (orthofinder-3.1.5-hdfd78af_0)
# Evidence: env
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the orthogroups of the orthology arm, at every inflation value scored
# <<< deposit header

# Line 1 is the delivered run, with DIAMOND as the search program in groups-only mode.
# Line 2 re-clusters the saved similarity graph at the other MCL inflation values of the
# sweep (1.5, 2.0, 3.0) from the same working folder, which
# 02_candidate_discovery/rl_arm1_selection.py scores beside the delivered run.
# Per the Supplementary Methods, membership was read either from the orthogroups at
# inflation 1.2 or 3.0, or from the similarity graph at a data-derived edge weight.

orthofinder -S diamond --only-groups -I 1.2
orthofinder -b <working folder> --only-groups -I <value>
