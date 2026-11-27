# >>> deposit header
# Program:  metapredict
# Version:  3.0.2
# Evidence: pinned
# Methods:  Supplementary Methods, Sequence phenotypes of the present subunits
# Source:   the deposited script that issues the call
# Outputs:  the per-residue disorder of the AMPK-alpha C-terminal region
# <<< deposit header

`metapredict.predict_disorder(sequence)`

metapredict is called as a Python library, not from a command line.
05_sequence_phenotypes/rbx1_chain/rbx1_alpha.py imports it (line 77), refuses to run
unless the installed version is the pinned 3.0.2 (lines 78 and 79, against
rbx1_lib.py line 73), and calls `meta.predict_disorder(clean)` at line 123 with no
network argument, so the release default applied. In metapredict 3.0.2
`metapredict/parameters.py` sets `DEFAULT_NETWORK = 'V3'`, so the model was V3.
Non-standard residues are substituted before prediction (rbx1_alpha.py line 39); the
predictor never sees one.
