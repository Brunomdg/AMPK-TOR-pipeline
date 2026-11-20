# >>> deposit header
# Program:  MAFFT
# Version:  7.526 (orthofinder3 environment)
# Evidence: env
# Methods:  Supplementary Methods, Sequence phenotypes of the retained subunits
# Source:   the deposited scripts that issue the call
# Outputs:  the pairwise and joint alignments of the AMPK-gamma, AMPK-alpha and AMPK-beta measurements
# <<< deposit header

# Issued from 05_sequence_phenotypes/rbx1_chain/rbx1_score.py (mafft(), line 57), which
# rbx1_gamma.py, rbx1_alpha.py and rbx2_ampkb_cbm.py all call. The binary is the
# orthofinder3 environment's mafft (rbx1_lib.py line 70, rbx2_common.py line 49).
# docs/SOFTWARE.md omits --quiet; the call site is the authority.

mafft --globalpair --maxiterate 1000 --anysymbol --quiet <input.faa>
