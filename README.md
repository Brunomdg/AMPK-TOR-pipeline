# Code for "Divergence of TOR and AMPK regulators across 210 eukaryotes predicts a dispensable TOR in *Toxoplasma gondii*"

Version 1.1.0

Aditya Patwardhan<sup>1,#</sup>, Euwen Brennan<sup>1,#</sup>, Bruno Martorelli Di Genova<sup>1,*</sup>

<sup>1</sup> Department of Microbiology and Molecular Genetics, Larner College of Medicine, University of Vermont,
Burlington, Vermont, USA  
<sup>#</sup> Shared first authors  
<sup>*</sup> Corresponding author: bruno.mdigenova@med.uvm.edu

## What is here

The code of the analyses that are new in this study and that produce the data of its main results. The
directories are the subsections of the Supplementary Methods, so a reader who has the paper open knows where
to look.

| Directory | Supplementary Methods subsection | What the code does |
|---|---|---|
| `1_proteomes_taxonomy_species_tree/` | Proteomes, taxonomy, and the species tree | External programs only. Their command lines are here; the dating of record cannot be reproduced (`docs/dating_gap.md`) |
| `2_candidate_discovery_and_orthology/` | Candidate discovery and orthology calls | The search settings chosen jointly against declared controls and decoys; the gene-tree clade criterion anchored on curated controls, with the two-arm rule where the tree abstains; the re-examination of absences and their clade support; the upstream-kinase classification |
| `3_species_states_and_composition/` | Species states and composition | The candidate table and the 210-species state matrix, the canonical TORC core states, each species' upstream-kinase state, and the coding used by the tree analyses |
| `4_sequence_phenotypes/` | Sequence phenotypes of the present subunits | The four measurements on the present AMPK subunits, the coccidian AMPKγ comparison, and the TOR N-terminal HEAT calls |
| `5_phylogenetic_analyses/` | Phylogenetic analyses | Penalised phylogenetic logistic regression with a multistart search, every ordered pair of the ten targets, the AMPKα tail models, and the ancestral-state reconstructions with the origins of absence |
| `6_structure_prediction_and_pocket/` | Structure prediction and AMPKγ pocket geometry | The domain-preserving prediction queries and the AMP-site geometry of the predicted models |

Two subsections have no directory. The comparison with earlier surveys is computed inside
`4_sequence_phenotypes/ampkg_contact_scores.py`, which reproduces the earlier single-species site-1 scores
under the primary procedure and reports every difference with its cause. The development of the *T. gondii*
TOR-knockout strains is wet-lab work.

Each stage directory holds an `external_runs/` folder with the command line of every external program that
stage ran, one file per program and step. `lib/` at the root holds the modules the candidate-discovery and
orthology scripts share; each stage that has its own helpers keeps them in its own `lib/`.

- `docs/code_map.md` maps every table and figure of the article to this code, and lists each method's scripts
  in the order they run.
- `docs/software.md` gives the version of every program the study used, with the class of evidence for that
  version, and the settings it ran with.
- `docs/constants/` holds the values that decide a call, as files that can be read without opening a script:
  the accessions that count as HEAT evidence, the references and acceptance rule of the HEAT structural
  screen, the AMPKβ reference and scoring rule, the structures behind the AMPKγ AMP-contact positions, and
  the three input variants of the site 1 score.
- `docs/dating_gap.md` states the one step that cannot be reproduced from this repository.
- `docs/deposit_gaps.md` states the other limits of the deposit, including the one script that is not the
  revision which produced the published fits.
- `docs/external_runs.md` explains what the `external_runs/` folders are and what they are not.
- `environment/` holds the complete package lists of the software environments and their creation commands.

## What is not here

- The job scripts of the standard programs. What is here is their command lines, in the `external_runs/`
  folders, with versions and evidence in `docs/software.md`. The scheduler directives, the module loads, the
  path defaults, and the inputs and outputs of those runs are not here.
- The drawing of figures and the typesetting of tables.
- The independent recomputation checks.
- Supplementary analyses that apply standard tests to the outputs of these methods.

They are available from the corresponding author on request.

On 22 September 2026, three AMPKγ copies were recovered from genome-grade proteomes, and the
AMPKγ-dependent results were carried forward on the updated inputs. The reruns that did so are not included
(see `docs/code_map.md`, section 3).

## The scripts

The scripts are the versions that produced the results.
- Each begins with a header giving its purpose, the Materials and Methods subsection it serves, its inputs,
  outputs and programs.
- Each gives the name it ran under on a `Ran as:` line, so that cluster logs and output folders can still be
  traced to it.
- Comments were rewritten in plain language.
- Absolute paths were replaced by variables whose defaults are the original paths. Run without those
  variables set, a script reads and writes where it did.
- Every script finds its own libraries relative to the repository root, so a clone imports without
  configuration.

## Running the scripts

The scripts expect the inputs of the data deposit: [data DOI to be added]. To run one, install the programs
of `docs/software.md`, and set the `TA210_*` variables listed in `docs/code_map.md`, section 3, to local
folders with the layout that section describes.

## Licence and citation

The code is licensed under the MIT licence (`LICENSE`). Please cite the article and this repository
(`CITATION.cff`).
