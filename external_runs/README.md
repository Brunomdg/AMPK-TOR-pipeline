# External program runs

Every program the study ran that is not itself deposited as a script. Each file names the program, its
version and the evidence for that version, the Supplementary Methods subsection it serves, and the command
line as run. The Supplementary Methods cite these paths instead of restating a flag string.

The directories follow the eight Supplementary Methods subsections:

| Directory | Supplementary Methods subsection |
|---|---|
| `sm1_proteomes_taxonomy_species_tree/` | Proteomes, taxonomy and the species tree |
| `sm2_candidate_discovery_and_orthology/` | Candidate discovery and orthology calls |
| `sm4_sequence_phenotypes/` | Sequence phenotypes of the retained subunits |
| `sm5_phylogenetic_analyses/` | Phylogenetic analyses |
| `sm6_structure_prediction_and_pocket/` | Structure prediction and AMPK-gamma pocket geometry |

Three subsections have no directory here. Species states and composition, and the comparison with earlier
surveys, ran no external program: they are computed by the deposited scripts. The development of the
*Toxoplasma gondii* TOR-knockout strains is wet-lab work.

## What these files are, and what they are not

**The job scripts themselves are not in this repository.** The command lines below were transcribed into
`docs/SOFTWARE.md` from those job scripts while they still existed, and this directory carries them at one
path per program and step. Each file states its own source:

- *transcribed from `docs/SOFTWARE.md`* means the flag string is the job script's, copied at the time;
- *the deposited script that issues the call* means the command is built in code that is in this
  repository, and the file names the script and line, which is the stronger evidence of the two;
- where neither holds, the file says so and gives the documented procedure rather than a flag string. No
  command line is reconstructed from memory.

Two programs ran with no surviving command line at all: `treePL` and the `chronos` dating of record.
`sm1_proteomes_taxonomy_species_tree/chronos_dating.md` records that its script is lost and points to
`docs/DATING_GAP.md`. Nothing is reconstructed in its place.

Placeholders are written in angle brackets, for example `<proteome.faa>`. They mark arguments the records do
not carry, never a flag whose value is known.

## Running these

These files are a record, not a pipeline. They carry no path defaults, no module loads and no scheduler
directives, because the job scripts that held those are gone. To rerun a step, read the command line here,
install the version in `docs/SOFTWARE.md`, and supply your own inputs.
