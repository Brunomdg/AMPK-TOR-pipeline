# >>> deposit header
# Program:  STRPsearch, with Foldseek
# Version:  Foldseek 10.941cd33; STRPsearch commit ed325a7f6b77578ee96753f31dbdcfda4931cca8 (branch master)
# Evidence: env, source clone
# Methods:  Supplementary Methods, Sequence phenotypes of the retained subunits
# Source:   the run's own receipt file, which records the command as executed, argument by argument
# Outputs:  one structured tandem-repeat result directory per structure
# <<< deposit header

# STRPsearch publishes no release number, so the version is the commit of the clone that ran, read from
# that clone's git refs: ed325a7f6b77578ee96753f31dbdcfda4931cca8, from
# https://github.com/BioComputingUP/STRPsearch.git, branch master.
#
# The command below is transcribed from a run receipt written beside the results, which stores the
# argument vector as executed. It is the command as run, not a reconstruction. The python interpreter is
# the tor_heat_tools environment's, which is also where Foldseek comes from.
#
# --no-chainsaw skips the domain-segmentation step. --max-eval 0.1 is the E-value ceiling.
# --min-height auto leaves the dendrogram cut to the program.

python bin/strpsearch.py query-file <structure.pdb> <output dir> --no-chainsaw --max-eval 0.1 --min-height auto --temp-dir <tmp dir>
