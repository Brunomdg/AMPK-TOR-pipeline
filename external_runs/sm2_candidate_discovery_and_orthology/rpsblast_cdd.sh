# >>> deposit header
# Program:  BLAST+ (makeprofiledb, rpsblast)
# Version:  2.17.0
# Evidence: env, probe
# Methods:  Supplementary Methods, Candidate discovery and orthology calls
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  the CDD hit table of the domain arm
# <<< deposit header

makeprofiledb -in <every profile of CDD 3.21>
rpsblast -evalue 0.01 -outfmt '6 qseqid sseqid bitscore evalue qstart qend'
