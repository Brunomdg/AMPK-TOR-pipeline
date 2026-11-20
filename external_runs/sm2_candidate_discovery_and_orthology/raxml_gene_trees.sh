# >>> deposit header
# Program:  RAxML
# Version:  8.2.13 (conda package raxml-8.2.13-h7b50bb2_3); its own banner prints 8.2.12
# Evidence: env, banner, log
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the RAxML cross-check trees
# <<< deposit header

# The two version numbers describe one binary. The job scripts call
# <envs>/orthofinder3/bin/raxmlHPC-PTHREADS; that environment's package metadata shows the file is
# owned by raxml-8.2.13-h7b50bb2_3 from bioconda, and running it prints "This is RAxML version 8.2.12
# released by Alexandros Stamatakis on May 2018". Upstream did not change the version string for that
# release, so the run logs say 8.2.12 and the package list says 8.2.13 for the same file.
# The RAxML replication is not described in the manuscript.

raxmlHPC-PTHREADS -f a -x 488761235 -p 488761235 -N autoMRE -m PROTGAMMAIWAG
