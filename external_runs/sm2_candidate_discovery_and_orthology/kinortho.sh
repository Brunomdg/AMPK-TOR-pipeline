# >>> deposit header
# Program:  KinOrtho
# Version:  commit 5c14ed1d0b49b6e98d16b56522e0b886d48b11e1, 26 July 2021
# Evidence: source clone
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   the job scripts, which survive on the cluster that ran them
# Outputs:  the AMPK-alpha ortholog sets of the orthology arm, one Results_*.txt set per arm
# <<< deposit header

# Run for AMPK-alpha only. KinOrtho publishes no release number, so the version is the commit of the
# clone that ran.
#
# Three arms, each over the same 210 proteomes with the same full-length query:
#
#   ARM A  the published PF00069 domain query, with no -E passed, so KinOrtho's own default applied.
#          The arms script's header describes this arm as E=1e-5.
#   ARM B  the domain query recut with the AMPK-alpha control HMM instead of PF00069, at -E 1e-5.
#          The recut runs hmmsearch --cpu 8 --domtblout against AMPKa_controls.hmm first.
#   ARM C  the published PF00069 domain query at -E 1e1, the BLAST default of 10.
#
# The three lines below are the invocations as the job scripts issue them.

python3 kinortho.py -i <proteomes/> -f AMPKa_query_full.fasta -d AMPKa_query_domain.fasta -o out/ -t 32
python3 kinortho.py -i <proteomes/> -f AMPKa_query_full.fasta -d <recut domain query> -o out/ -t 32 -E 1e-5
python3 kinortho.py -i <proteomes/> -f AMPKa_query_full.fasta -d AMPKa_query_domain.fasta -o out/ -t 32 -E 1e1
