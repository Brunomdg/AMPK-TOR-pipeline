# >>> deposit header
# Program:  ColabFold
# Version:  1.6.2 (alphafold-colabfold 2.3.18, JAX 0.10.2)
# Evidence: env, project log
# Methods:  Supplementary Methods, Structure prediction and AMPK-gamma pocket geometry
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  25 models per query, read by the pocket-geometry scripts
# <<< deposit header

# Every AMPK-gamma copy, predicted alone. The earlier 86-query panel ran at the same
# settings. Five models by five seeds is the 25 models per query that
# 11_structure_pocket/all_copies/analyze_all_ampkg_structures.py reads.
# No template flag and no relaxation flag was passed, so ColabFold's defaults applied:
# no templates and no Amber relaxation. The deposited analysis scripts read only
# *_unrelaxed_rank_*_seed_*.pdb files (p45b_colabfold_analysis.py line 83,
# analyze_all_ampkg_structures.py line 52), which is the artifact that shows no relaxed
# model was produced.

colabfold_batch --msa-only <queries.fasta> <msa dir>
colabfold_batch --model-type alphafold2_ptm --num-models 5 --num-seeds 5 --num-recycle 12 --random-seed 0 <msa dir> <out dir>
