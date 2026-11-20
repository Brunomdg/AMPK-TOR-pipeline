# >>> deposit header
# Program:  RAxML
# Version:  8.2.12 per its own logs; 8.2.13 is the installed conda package
# Evidence: log
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the RAxML cross-check trees
# <<< deposit header

# The version conflict is unresolved and is recorded as such: the logs of the runs report
# 8.2.12, the conda package installed in the busco6 and orthofinder3 environments is
# 8.2.13. The RAxML replication is not described in the manuscript.

raxmlHPC-PTHREADS -f a -x 488761235 -p 488761235 -N autoMRE -m PROTGAMMAIWAG
