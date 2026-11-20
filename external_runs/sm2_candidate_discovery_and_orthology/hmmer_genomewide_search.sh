# >>> deposit header
# Program:  HMMER
# Version:  3.4 (Aug 2023)
# Evidence: log, env, probe, script
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   the option-settings line that hmmsearch writes into every one of its own output tables
# Outputs:  the per-sequence and per-domain tables that every candidate-discovery script reads
# <<< deposit header

# Line 1 is the genome-wide search of all 210 proteomes against the profile panel, one task per
# proteome. It is not a transcription: hmmsearch writes its version and its complete option settings
# into the footer of every table it produces, and line 1 is that footer with the absolute paths
# replaced by placeholders. The same invocation is issued from
# 02_candidate_discovery/rl_threshold_pilot_pertarget.py (run_panel, line 180), which passes the same
# ceilings and --max. Line 2 builds the contaminant-family profiles. Line 3 is the genome-grade
# re-examination of absences.

hmmsearch --tblout <tbl> --domtblout <domtbl> -E 1e5 --domE 1e5 --max --cpu 8 <panel_v2.hmm> <proteome.faa>
hmmbuild <family.hmm> <family alignment>
phmmer <query.faa> <proteome.faa>
