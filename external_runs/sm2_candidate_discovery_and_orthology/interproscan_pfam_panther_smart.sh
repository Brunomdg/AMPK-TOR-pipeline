# >>> deposit header
# Program:  InterProScan
# Version:  5.78-109.0
# Evidence: probe
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the InterPro annotation table of the domain arm
# <<< deposit header

# The installation no longer exists; its outputs remain. The same annotations are read
# again by the TOR HEAT screens (see ../sm4_sequence_phenotypes/).

interproscan.sh -f TSV -appl Pfam,PANTHER,SMART
