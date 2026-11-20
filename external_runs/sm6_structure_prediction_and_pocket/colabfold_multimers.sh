# >>> deposit header
# Program:  ColabFold
# Version:  1.6.2 (alphafold-colabfold 2.3.18, JAX 0.10.2)
# Evidence: env, project log
# Methods:  Supplementary Methods, Structure prediction and AMPK-gamma pocket geometry
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the three-chain complexes and the partner-context benchmark
# <<< deposit header

# The same settings as the monomer runs, with the multimer model type.

colabfold_batch --model-type alphafold2_multimer_v3 --num-models 5 --num-seeds 5 --num-recycle 12 --random-seed 0 <msa dir> <out dir>
