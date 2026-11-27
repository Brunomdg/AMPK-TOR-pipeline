# >>> deposit header
# Program:  US-align
# Version:  20241201
# Evidence: env
# Methods:  Supplementary Methods, Sequence phenotypes of the present subunits
# Source:   the deposited script that issues the call
# Outputs:  one US-align report per TOR copy and experimental reference
# <<< deposit header

# Issued from 06_tor_heat/structural_screen/run_direct_alignment.py (line 55), which
# builds the argument vector ['USalign', <query>, <reference>, '-mol', 'prot'].
# This step serves the TOR N-terminal HEAT calls; the Supplementary Methods subsection
# on structure prediction covers the AMPK-gamma pocket, not this alignment.

USalign <model crop> <reference structure> -mol prot
