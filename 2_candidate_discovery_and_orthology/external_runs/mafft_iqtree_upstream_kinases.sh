# >>> deposit header
# Program:  MAFFT, IQ-TREE
# Version:  7.526, 3.1.3
# Evidence: env, banner
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the masked kinase-domain trees that rbx3_tree_calls.R reads
# <<< deposit header

# The external step between 07_phylogenetic_regressions/upstream_kinases/rbx3_mask.py and
# rbx3_tree_calls.py. docs/CODE_MAP.md records it as MAFFT L-INS-i and IQ-TREE.
# L-INS-i is --localpair --maxiterate 1000; the IQ-TREE line repeats the gene-tree
# settings. Neither line is copied from a surviving job script, and no per-run log for
# this step is recorded, so both are the documented method, not a transcription.

mafft --localpair --maxiterate 1000 <kinase domains>
iqtree3 -m MFP -B 1000 --alrt 1000 --seed 488761235
