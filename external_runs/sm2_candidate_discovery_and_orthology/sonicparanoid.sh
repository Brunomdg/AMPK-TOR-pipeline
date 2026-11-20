# >>> deposit header
# Program:  SonicParanoid
# Version:  2.0.9 (sonicparanoid-2.0.9-py311h93dcfea_0)
# Evidence: env (pinned at creation)
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository; the delivered setting is fixed by two deposited scripts
# Outputs:  the ortholog groups of the orthology arm, at every setting scored
# <<< deposit header

# Line 1 is the delivered clustering: the program's own defaults, arch-merging coverage
# 0.75 and MCL inflation 1.5. 02_candidate_discovery/rl_arm1_selection.py names it
# "cov0.75/I1.5 (CURRENT)" (line 544) and 02_candidate_discovery/rl_invariants.py records
# that the jointly chosen settings use it for six of the eight targets (lines 74 to 78).
# Line 2 is the sweep run that sp2_cov90.sbatch submitted; 0.90 is a sweep value, not the
# delivered one. Line 3 varies MCL inflation (2.0 and 3.0).
# docs/SOFTWARE.md gave 0.90 as the setting as run until this was corrected. The Supplementary
# Methods states the same settings as the scripts: arch-merging coverage 0.75 with inflation
# 1.5, or inflation 2.0 or 3.0. Membership was read from those ortholog groups, except for
# AMPK-alpha, where it was read from the pairwise ortholog table at confidence 1.0.

sonicparanoid -m default -ca -op -ot
sonicparanoid -m default -ca -op -ot --min-arch-merging-cov 0.90
sonicparanoid -m default -ca -op -ot -I <value>
