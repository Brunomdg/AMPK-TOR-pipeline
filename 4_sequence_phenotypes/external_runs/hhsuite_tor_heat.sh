# >>> deposit header
# Program:  HH-suite (hhsearch, HHrepID)
# Version:  3.3.0
# Evidence: env
# Methods:  Supplementary Methods, Sequence phenotypes of the present subunits
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  remote-homology and repeat evidence for the TOR N-terminal region
# <<< deposit header

# The databases are SCOP40 (scop40_01Mar17) and PDB70 (pdb70_from_mmcif_220313).
# HHrepID ran with its own calibration files; docs/SOFTWARE.md records no flag string
# for it, so line 3 is the bare invocation and not a transcription.

hhsearch -M first -p 0 -E 100 -Z 100 -B 100 -e 100 -d <scop40_01Mar17>
hhsearch -M first -p 0 -E 100 -Z 100 -B 100 -e 100 -d <pdb70_from_mmcif_220313>
hhrepid -i <input> <HHrepID calibration files>
