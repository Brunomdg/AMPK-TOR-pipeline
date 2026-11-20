# >>> deposit header
# Program:  IQ-TREE
# Version:  3.1.3
# Evidence: env, banner, log
# Methods:  Supplementary Methods, Proteomes, taxonomy and the species tree
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the 210-taxon species tree, the 209-taxon PMSF tree and the concordance factors
# <<< deposit header

# Line 1 is the partitioned 210-taxon tree, line 2 the 209-taxon PMSF tree,
# line 3 the gene- and site-concordance factors.

iqtree3 -p <partitions> -m LG+G4 -B 1000 -seed 20260817
iqtree3 -m LG+C20+F+G4 -ft <210-taxon tree> -B 1000 -seed 20260817
iqtree3 --gcf <20 marker trees> --scf 100
