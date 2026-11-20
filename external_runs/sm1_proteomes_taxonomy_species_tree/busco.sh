# >>> deposit header
# Program:  BUSCO
# Version:  6.1.0, lineage eukaryota_odb12.2
# Evidence: env
# Methods:  Supplementary Methods, Proteomes, taxonomy and the species tree
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  one BUSCO run per proteome; summarised in the proteome-completeness table
# <<< deposit header

# One proteome per task, over the 210 proteomes.

busco -m proteins -l eukaryota_odb12.2 --offline
