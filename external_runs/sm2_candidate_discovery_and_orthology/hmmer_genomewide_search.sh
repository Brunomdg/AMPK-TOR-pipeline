# >>> deposit header
# Program:  HMMER
# Version:  3.4
# Evidence: env, probe
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository; the first command is also issued by a deposited script
# Outputs:  the per-sequence and per-domain tables that every candidate-discovery script reads
# <<< deposit header

# Line 1 is the genome-wide search of all 210 proteomes against the profile panel. The
# same invocation is issued from 02_candidate_discovery/rl_threshold_pilot_pertarget.py
# (run_panel), which passes -E 1e5 --domE 1e5 --max and --tblout; that call site is the
# evidence for the score ceilings. Line 2 builds the contaminant-family profiles. Line 3
# is the genome-grade re-examination of absences.

hmmsearch --cpu ${SLURM_CPUS_PER_TASK} -E 1e5 --domE 1e5 --max --tblout <tbl> --domtblout <domtbl> <panel.hmm> <proteome.faa>
hmmbuild <family.hmm> <family alignment>
phmmer <query.faa> <proteome.faa>
