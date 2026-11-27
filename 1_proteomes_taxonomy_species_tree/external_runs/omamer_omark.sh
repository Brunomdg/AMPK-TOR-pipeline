# >>> deposit header
# Program:  OMAmer, OMArk
# Version:  2.1.2, 0.5.0 (database LUCA.h5)
# Evidence: env
# Methods:  Supplementary Methods, Proteomes, taxonomy and the species tree
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  one OMArk report per proteome
# <<< deposit header

# -t was passed where a taxon id was known. The query and output arguments are
# written here as placeholders: docs/SOFTWARE.md records the flags, not the paths.

omamer search --db LUCA.h5 --query <proteome.faa> --out <search>
omark -f <search> -d LUCA.h5 -t <taxon id>
