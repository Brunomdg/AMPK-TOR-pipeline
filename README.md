# Code for "Beyond Presence and Absence: Noncanonical TOR and AMPK Divergence Across 210 Eukaryotes Challenges Mammalian Generalizations"

Version 1.0.0

Aditya Patwardhan<sup>1,#</sup>, Euwen Brennan<sup>1,#</sup>, Bruno Martorelli Di Genova<sup>1,*</sup>

<sup>1</sup> Department of Microbiology and Molecular Genetics, Larner College of Medicine, University of Vermont,
Burlington, Vermont, USA  
<sup>#</sup> Shared first authors  
<sup>*</sup> Corresponding author: bruno.mdigenova@med.uvm.edu

## What is here

This repository holds the code of the eleven methods that are new in this study and that produce the data of
its main results.

| # | Method | Folder | Main result |
|---|---|---|---|
| 1 | Search settings chosen jointly against declared controls and decoys | `02_candidate_discovery/` | the candidate set |
| 2 | Gene-tree clade criterion anchored on curated controls, with the two-arm rule where the tree abstains | `03_orthology_calls/` | every orthology call and the 210-species state matrix (Figure 1) |
| 3 | Re-examination of absences and their clade support | `03_orthology_calls/` | Table 1 |
| 4 | TOR N-terminal HEAT calls from signatures and structural alignment | `06_tor_heat/` | Table 1 and Figure 3 |
| 5 | Sequence measurements of the retained AMPK subunits | `05_sequence_phenotypes/` | Figures 1, 2 and 4 |
| 6 | Measurements of the coccidian AMPKγ copies | `05_sequence_phenotypes/coccidian_remodelling/` | Figure 2 |
| 7 | Structure queries that keep every CBS domain, and AMP-site geometry in the predicted models | `11_structure_pocket/` | Figure 2 |
| 8 | Penalised phylogenetic logistic regression with a multistart search | `07_phylogenetic_regressions/` | Figure 5 |
| 9 | Models of the AMPKα C-terminal tail | `07_phylogenetic_regressions/ctd_cross_feature/` | Figure 4 |
| 10 | Upstream-kinase classification on reference-masked kinase-domain trees | `07_phylogenetic_regressions/upstream_kinases/` | Figure 4 |
| 11 | Ancestral-state reconstructions and origins of absence | `08_ancestral_states/` | Figure 6 |

- `docs/CODE_MAP.md` maps every table and figure of the article to this code, and lists each method's
  scripts in the order they run. Section 6 maps the eight subsections of the Supplementary Methods to the
  scripts that implement them.
- `docs/SOFTWARE.md` gives the version of every program the study used and the settings it ran with.
- `external_runs/` holds the command line of every external program, one file per program and step, each
  naming the Supplementary Methods subsection it serves.
- `docs/constants/` holds the values that decide a call, as files that can be read without opening a script:
  the accessions that count as HEAT evidence, the references and acceptance rule of the HEAT structural
  screen, the AMPKβ reference and scoring rule, the structures behind the AMPKγ AMP-contact positions, and
  the three input variants of the site 1 score.
- `docs/DATING_GAP.md` states the one step that cannot be reproduced from this repository.
- `envs/` holds the complete package lists of the software environments.

`shared/pipeline_lib/` and the other helper modules are the code these scripts import.

## The Methods of the article

The Materials and Methods is split in two. The main-text Methods carries the decision rules, the thresholds
that determine a call, the headline software and the wet-lab protocol. The Supplementary Methods carries the
exact invocations, the reference-site definitions, the sensitivity variants and the model-by-model
specifications, under the same eight headings, and **cites paths in this repository** rather than restating a
flag string. `docs/CODE_MAP.md`, sections 6 to 8, is the index from those headings to the code, the external
command lines and the published constants.

## What is not here

- The job scripts of the standard programs: HMMER, OrthoFinder, SonicParanoid, KinOrtho, rpsblast,
  InterProScan, MAFFT, trimAl, IQ-TREE, RAxML, BUSCO, OMArk, ColabFold, US-align, HH-suite, HyPhy and
  rate4site. What is here is their command lines, in `external_runs/`, transcribed from those job scripts
  while they still existed, with versions and evidence in `docs/SOFTWARE.md`. The scheduler directives, the
  module loads, the path defaults, and the inputs and outputs of those runs are not here. Two command lines
  did not survive at all: the `chronos` dating of record, whose script is lost (`docs/DATING_GAP.md`), and
  treePL, which is driven by a configuration file that is not in the records.
- The drawing of figures and the typesetting of tables.
- Checks, including the independent recomputation of every fitted model, and the corrections and reruns made
  after the main runs.
- Supplementary analyses that apply standard tests to the outputs of these methods.

They are available from the corresponding author on request.

On 22 September 2026, three AMPKγ copies were recovered from genome-grade proteomes, and the AMPKγ-dependent
results were carried forward on the updated inputs. The reruns that did so are not included (see
`docs/CODE_MAP.md`, section 4).

## The scripts

The scripts are the versions that produced the results.
- Each begins with a header giving its purpose, the Materials and Methods subsection it serves, its inputs,
  outputs and programs.
- Comments were rewritten in plain language.
- Absolute paths were replaced by variables whose defaults are the original paths. Run without those
  variables set, a script behaves as it did.
- The code itself is unchanged.
- Three files carry no header, because other scripts check their checksums (`docs/CODE_MAP.md`, section 3).
- Scripts keep the names they ran under; `docs/CODE_MAP.md`, section 5, explains the labels in them.

## Running the scripts

The scripts expect the folder layout of the VACC cluster of the University of Vermont, and the inputs of the
data deposit: [data DOI to be added]. To run one elsewhere, set the `TA210_*` variables listed in
`docs/CODE_MAP.md`, section 4, to local folders with the same layout, and install the programs of
`docs/SOFTWARE.md`.

The two scripts of method 9 cannot run as deposited. The version of the model engine they check was not
preserved.

## Licence and citation

The code is licensed under the MIT licence (`LICENSE`). Please cite the article and this repository
(`CITATION.cff`).
