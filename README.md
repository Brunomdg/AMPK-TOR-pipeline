# Code for "Beyond presence and absence: a survey of TOR and AMPK across 210 eukaryotes reveals unusually long coccidian AMPKγ"

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
| 4 | TOR N-terminal HEAT calls from signatures and structural alignment | `06_tor_heat/` | Table 1 |
| 5 | Sequence measurements of the retained AMPK subunits | `05_sequence_phenotypes/` | Figures 2 and 5 |
| 6 | Measurements of the coccidian AMPKγ copies | `05_sequence_phenotypes/coccidian_remodelling/` | Figure 3 |
| 7 | Structure queries that keep every CBS domain, and AMP-site geometry in the predicted models | `11_structure_pocket/` | Figure 3 |
| 8 | Penalised phylogenetic logistic regression with a multistart search | `07_phylogenetic_regressions/` | Figure 4, Table 2 |
| 9 | Models of the AMPKα C-terminal tail | `07_phylogenetic_regressions/ctd_cross_feature/` | Figure 6 |
| 10 | Upstream-kinase classification on reference-masked kinase-domain trees | `07_phylogenetic_regressions/upstream_kinases/` | Figure 5 |
| 11 | Ancestral-state reconstructions and origins of absence | `08_ancestral_states/` | Figure 7 |

- `docs/CODE_MAP.md` maps every table and figure of the article to this code, and lists each method's
  scripts in the order they run.
- `docs/SOFTWARE.md` gives the version of every program the study used and the settings it ran with.
- `envs/` holds the complete package lists of the software environments.

`shared/pipeline_lib/` and the other helper modules are the code these scripts import.

## What is not here

- The runs of standard programs: HMMER, OrthoFinder, SonicParanoid, KinOrtho, rpsblast, InterProScan, MAFFT,
  trimAl, IQ-TREE, RAxML, BUSCO, OMArk, ColabFold, US-align, HH-suite, HyPhy and rate4site. Their versions,
  and the settings as run, are in `docs/SOFTWARE.md`.
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
