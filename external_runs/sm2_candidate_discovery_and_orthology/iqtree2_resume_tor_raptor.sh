# >>> deposit header
# Program:  IQ-TREE
# Version:  2.3.6 (cluster module iqtree2/2.3.6)
# Evidence: job script
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the completed TOR and RAPTOR gene trees
# <<< deposit header

# The TOR and RAPTOR searches were interrupted and resumed from their checkpoints under
# the cluster module iqtree2/2.3.6. docs/SOFTWARE.md records the program, the version and
# the fact of the resume, but no flag string, so the line above is the bare invocation
# and not a transcription. The search settings are those of iqtree_gene_trees.sh.

iqtree2 -s <alignment> --redo-tree-check
