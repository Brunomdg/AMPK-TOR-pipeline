# >>> deposit header
# Program:  pyKVFinder
# Version:  0.9.5
# Evidence: env
# Methods:  Supplementary Methods, Structure prediction and AMPK-gamma pocket geometry
# Source:   the deposited script that issues the call
# Outputs:  cavity volumes at the three AMP sites
# <<< deposit header

`pyKVFinder.run_workflow(structure, box=<site box>, ligand=<site AMP>, ligand_cutoff=5.0, nthreads=1)`

pyKVFinder is called as a Python library, not from a command line.
11_structure_pocket/panel_86_queries/p45b_colabfold_analysis.py line 92 issues the call,
on the model moved into the 2V8Q frame by the site superposition, with the B-factor
column set to 0. The box is fixed per site in the 2V8Q frame: the site residues'
heavy-atom extent plus a padding of 3.5 A (PAD, line 27, the pyKVFinder residue-box
default). The probe and step settings are pyKVFinder's own defaults.
See ../../docs/constants/ampkg_amp_site_structures.md for the structures involved.
