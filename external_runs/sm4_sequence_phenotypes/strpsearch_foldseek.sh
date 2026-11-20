# >>> deposit header
# Program:  STRPsearch, with Foldseek
# Version:  Foldseek 10.941cd33; STRPsearch not recorded
# Evidence: env
# Methods:  Supplementary Methods, Sequence phenotypes of the retained subunits
# Source:   docs/SOFTWARE.md, which copied the command line from the job script. The job script itself is not in this repository
# Outputs:  structured tandem-repeat calls on each structure
# <<< deposit header

# Structured tandem-repeat search on each structure, maximum E-value 0.1.
# docs/SOFTWARE.md records the E-value ceiling but no full command line, so the line
# above is the documented setting, not a transcription. The STRPsearch version is not
# recorded.

STRPsearch -i <structure> -e 0.1
