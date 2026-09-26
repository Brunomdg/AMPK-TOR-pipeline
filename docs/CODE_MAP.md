# Code map

This repository holds the code of the eleven methods that are new in this study and that produce the data
of the main results. For each method, section 2 lists the scripts in the order they run. Section 1 goes
table by table and figure by figure through the paper and names the method behind each one's data. Paths are
relative to the repository root. Each script begins with a header giving its purpose, inputs, outputs and
programs.

Everything else the study ran is not included here, and the authors can provide it on request:
- the runs of standard programs, whose versions and settings are in `docs/SOFTWARE.md`;
- the drawing of figures and the typesetting of tables;
- checks, including the independent recomputations of the fitted models;
- corrections and reruns;
- supplementary analyses that apply standard tests to the outputs of these methods.

## 1. Tables and figures

Numbering follows the manuscript.

| Item | Data produced by | Notes |
|---|---|---|
| Table 1 | Methods 3 and 4 | The other rows are not included: the residue-correspondence test, the genome-grade re-examination, the joint loss-direction fits and proteome completeness |
| Table 2 | Method 8, `network10_refit.R` | |
| Figure 1 | Method 2 (the species-state matrix) | Supergroup shares are counts of that matrix |
| Figure 2 | Method 5 | The per-phylum summaries of the measurements are not included |
| Figure 3 | Methods 6 and 7 | The plotted values were prepared by hand from Supplementary Table S28 |
| Figure 4 | Method 8, `rbx2_torc_core.py` and `core18_penalized.R` | |
| Figure 5 | Method 10 (repertoire) and method 5, `rbx2_ampka_activation.py` (activation-site residue) | |
| Figure 6 | Method 9, with the measurements of method 5 | |
| Figure 7 | Method 11 | Panel A summarises sensitivity analyses that are not included |
| Supplementary Figure S1 | Method 2, with the species tree (`docs/SOFTWARE.md`) | |
| Supplementary Figures S2, S3 | Built on method 5 | The descriptive analyses of all copies are not included |
| Supplementary Figure S4 | Method 7 | |
| Supplementary Figures S5, S6 | Not included | The earlier 86-query multimer panel and the partner-context benchmark |
| Supplementary Figure S7 | Methods 1 and 2 | The performance measurements are not included |
| Supplementary Figure S8 | Method 8, `network10_refit.R` | |
| Supplementary Figure S9 | Method 8, `core18_penalized.R` | |
| Supplementary Figure S10 | Method 5, `rbx2_ampka_activation.py` | |
| Supplementary Figure S11 | Built on method 10 | The models are not included |
| Supplementary Figure S12 | Method 9 | The post hoc models are not included |
| Supplementary Figure S13 | Method 11 | The co-transition count is not included |
| Supplementary Figures S14, S15 | Not included | Convergence tests and lifestyle models |
| Supplementary Figure S16 | Method 10 | |
| Supplementary Tables S1, S2 | Methods 1 and 2 | |
| Supplementary Tables S3, S4 | Method 1 | The table builders are not included |
| Supplementary Table S5 | Method 2 (the criterion) | Its measured performance is not included |
| Supplementary Tables S6, S7, S9 | Built on method 2 | Not included |
| Supplementary Table S8 | Not included | The LKB1 candidate screen |
| Supplementary Table S10 | Not included | Candidate replacements |
| Supplementary Table S11 | Method 11 | |
| Supplementary Table S12 | Built on method 11 | Not included |
| Supplementary Table S13 | Not included | Lifestyle models |
| Supplementary Table S14 | Not included | The index of the numbers in the text |
| Supplementary Tables S15, S18 | Method 5, `rbx1_gamma.py` | |
| Supplementary Tables S16, S17 | Built on method 5 | The descriptive analyses of all 239 copies are not included |
| Supplementary Table S19 | Method 5, `rbx1_alpha.py` | |
| Supplementary Tables S20 to S22 | Not included | Built on methods 2 and 5 |
| Supplementary Table S23 | Built on method 5 | The summaries are not included |
| Supplementary Tables S24 to S26 | Not included | The block of additional estimates |
| Supplementary Table S27 | Not included | The comparison with earlier surveys |
| Supplementary Table S28 | Method 6 | |
| Supplementary Table S29 | External programs (BUSCO, OMArk) | Their summary is not included |
| Supplementary Tables S30 to S34 | Not included | Built on methods 2, 4 and 5 with standard tests: Pagel's test, the no-regain comparison, rate4site and the within-clade decomposition |
| Supplementary Table S35 | External program only (GARD) | |
| Supplementary Table S36 | Built on method 2 | The composition counts are not included |
| Supplementary Table S37 | Method 8, `network10_refit.R` | |
| Supplementary Table S38 | Method 5, `rbx2_ampkb_cbm.py` | The models are not included |

## 2. The methods, with their scripts in run order

External programs and the settings they ran with are in `docs/SOFTWARE.md`. Inputs that these scripts read
but do not make will be in the data deposit:
- the proteomes;
- the declared control sequences, decoys and negatives;
- the profile panel;
- the raw search outputs;
- the gene trees and the dated species tree;
- the structure predictions.

### Method 1. Search settings chosen jointly against declared controls and decoys

Before these scripts, the external searches run:
- a genome-wide HMMER search of the 210 proteomes;
- rpsblast against CDD;
- InterProScan;
- OrthoFinder, SonicParanoid and KinOrtho.

1. `02_candidate_discovery/rl_threshold_pilot_pertarget.py`: for each target and profile, the pair of sequence
   and domain score cuts that holds the false-positive rate on the target's same-species decoys at or below 5,
   1, 0.5, 0.1 and 0 percent, with the recall of the target's controls at each.
2. `02_candidate_discovery/rl_r59_union_optimum.py`: the joint choice of one setting per method and target,
   so that together the methods recover every declared control with the fewest declared negatives called.
3. `02_candidate_discovery/rl_per_method_report.py`: each method's calls at the chosen settings, the candidate
   list, and the contaminant screen in which each protein goes to the best-scoring family.

Helpers:
- `02_candidate_discovery/rl_arm1_selection.py`;
- `02_candidate_discovery/rl_invariants.py`, which checks the reference inputs against fixed checksums;
- `shared/pipeline_lib/`: `rl_accession.py`, `rl_declaration.py`, `rl_evidence_roots.py`, `rl_join.py`,
  `rl_targets.py` and `rl_vocab.py`, which hold the target lists, the accession joins and the call labels.

### Method 2. Gene-tree clade criterion and the species-state matrix

Before these scripts, the gene trees are built externally: MAFFT, IQ-TREE with ModelFinder, the PMSF
cross-check, and RAxML.

1. `03_orthology_calls/rl_phylo_subtree.py`: the clade criterion of van Dam et al. (2011), anchored on the
   curated controls, applied to each gene tree.
2. `03_orthology_calls/rl_ambiguous_arm_rule.py`: the two-arm call for candidates the tree left ambiguous.
3. `03_orthology_calls/rl_final_ortholog_table.py`: one call per candidate.
4. `03_orthology_calls/rl_core8_tree_rule.py`: resolves the core candidates that remained open.
5. `03_orthology_calls/rl_pub_t1_t2_active10.py`: the candidate table and the species-state matrix.

Helpers:
- `03_orthology_calls/lkb1_camkk/rl_require_camkk_lkb1_joint_freeze.py` and
  `rl_require_clean_camkk_lkb1_master_inputs.py`, which stop a step until the LKB1 and CAMKK calls are fixed;
- `shared/pipeline_lib/rl_release_guard.py`.

### Method 3. Re-examination of absences

1. `03_orthology_calls/rl_absence_research.py`: rescores every protein of a species against the eight core
   targets from the genome-wide search, looking for missed orthologs.
2. `03_orthology_calls/rl_absence_support.R`: classifies each absence as clade-supported or species-specific,
   on five species-tree topologies and two rootings.

### Method 4. TOR N-terminal HEAT calls

Before these scripts, InterProScan signatures, the AlphaFold models and US-align are run externally.

1. `06_tor_heat/signature_screen/classify_tor_heat.py`: HEAT signatures N-terminal to the FAT domain.
2. `06_tor_heat/signature_screen/refine_nterminal_heat.py`: transfers the 32 N-terminal HEAT repeats of human
   mTOR to every copy.
3. `06_tor_heat/structural_screen/run_direct_alignment.py`: aligns each model's N-terminal crop to the
   experimental references and applies the structural acceptance rule.
4. `06_tor_heat/final_calls/resolve_heat_final.py`: the final present or absent call, with its confidence
   tier.

### Method 5. Sequence measurements of the retained subunits

1. `05_sequence_phenotypes/rbx1_chain/rbx1_gamma.py`: scores the AMP-contact positions of every AMPKγ copy
   against the human reference.
2. `05_sequence_phenotypes/rbx1_chain/rbx1_alpha.py`: the AMPKα C-terminal length and its predicted disorder.
3. `05_sequence_phenotypes/rbx2/rbx2_ampkb_cbm.py`: the AMPKβ carbohydrate-binding module score.
4. `05_sequence_phenotypes/rbx2/rbx2_ampka_activation.py`: the residue at the AMPKα activation site, per copy
   and per species.

Helpers: `05_sequence_phenotypes/rbx1_chain/rbx1_lib.py`, `05_sequence_phenotypes/rbx1_chain/rbx1_score.py`,
`05_sequence_phenotypes/rbx2/rbx2_common.py`.

### Method 6. Coccidian AMPKγ copies

1. `05_sequence_phenotypes/coccidian_remodelling/build_degeneration_metrics.py`: collects the per-copy
   measurements of every AMPKγ copy and compares the coccidian copies with the others. The measurements are
   length, CBS-domain evidence and identity away from the AMP contacts, the site 1 contact score and its
   agreement between alignments, and predicted-structure confidence and RMSD.

### Method 7. Structure queries and AMP-site geometry

1. `11_structure_pocket/all_copies/prepare_ampkg_domain_preserving.py`: builds each copy's prediction
   sequence so that every observed CBS domain is kept.
2. ColabFold runs externally.
3. `11_structure_pocket/all_copies/analyze_all_ampkg_structures.py`: measures the three reference AMP sites
   in all 25 models of every prediction.

Helper: `11_structure_pocket/panel_86_queries/p45b_colabfold_analysis.py`, the measurement library.

### Method 8. Penalised phylogenetic logistic regression

1. `05_sequence_phenotypes/rbx2/rbx2_torc_core.py`: the TORC1 and TORC2 core state of every species.
2. `07_phylogenetic_regressions/core_complements/core18_penalized.R`: AMPK subunit absence on core
   disruption, with the multistart search.
3. `07_phylogenetic_regressions/ten_target_pairs/network10_refit.R`: every ordered pair of the ten targets.

Helpers: `07_phylogenetic_regressions/core_complements/core_multistart_helper.R` and
`07_phylogenetic_regressions/ten_target_pairs/network10_multistart_helper.R`, which hold the grid of starting
values and the local refinement.

### Method 9. The AMPKα tail against the retained modules and the upstream kinases

1. `07_phylogenetic_regressions/ctd_cross_feature/ctdx_cross_ampk.R`: tail length and disorder against
   AMPKγ and AMPKβ absence and the retained AMPKβ module score.
2. `07_phylogenetic_regressions/ctd_cross_feature/ctdx_upstream.R`: the same against absence of the
   upstream kinases.

Helper: `07_phylogenetic_regressions/ctd_cross_feature/ctdx_lib.R`.

**These two scripts cannot run as deposited.** `ctdx_lib.R` checks the checksum of the model engine it
reads, and that version of the engine was not preserved. `ctdx_cross_ampk.R` also reads the record of a
reproduction step that is not included.

### Method 10. Upstream-kinase classification

1. `07_phylogenetic_regressions/upstream_kinases/rbx3_mask.py`: masks each kinase-domain alignment to the
   columns the reference sequences occupy.
2. MAFFT L-INS-i and IQ-TREE run externally.
3. `07_phylogenetic_regressions/upstream_kinases/rbx3_tree_calls.R`: assigns each tip by its placement
   among the references.
4. `07_phylogenetic_regressions/upstream_kinases/rbx3_calls.py` for LKB1 and `rbx3_camkk_calls.py` for the
   CAMKK category: the per-copy calls.
5. `07_phylogenetic_regressions/upstream_kinases/rbx3_architecture.py`: each species' upstream-kinase state.

Helper: `07_phylogenetic_regressions/upstream_kinases/rbx3_common.py`.

### Method 11. Ancestral states and origins of absence

1. `08_ancestral_states/rl_rm_q7_standard.py`: the standard coding of each target. Present is 0, not
   assessable 0.5, absent 1.
2. `08_ancestral_states/repair/rl_rm_q0_asr.R`: the state at both ends of every edge, with corHMM (ER, SYM
   and ARD) and parsimony, on every tree, rooting and taxon subset.
3. `08_ancestral_states/rl_rm_q7_origins.py`: the number of origins of absence in each reconstruction block.

Helpers:
- `08_ancestral_states/repair/rl_rm_asr_lib.R`, the repaired reconstruction library;
- `08_ancestral_states/repair/rl_rm_q7_lib.py`;
- `07_phylogenetic_regressions/replacement_and_lifestyle_inputs/rl_rm_prov.R`, which writes the header lines
  of the outputs.

## 3. Files without a header

These three files are exactly as they ran, because other scripts check their checksums.

| File | What it is |
|---|---|
| `07_phylogenetic_regressions/core_complements/core_multistart_helper.R` | The multistart search of the core-complement fits: a grid of starting values, then local refinement, as described in the Methods |
| `07_phylogenetic_regressions/ctd_cross_feature/ctdx_lib.R` | The shared library of the tail models. It checks every input table and the model engine, reads the engine's functions, and fits the models |
| `11_structure_pocket/panel_86_queries/p45b_colabfold_analysis.py` | The measurement library of the structure predictions. It reads each query's models and measures the three AMP sites; the measurement scripts load it after checking its checksum |

## 4. Paths and the data the three recovered AMPKγ copies changed

**Paths.** The scripts address the folders of the VACC cluster. Each script sets its own defaults, between
two `deposit path defaults` comment lines. To run one elsewhere, set these variables to local folders with
the same layout:

| Variable | Default |
|---|---|
| `TA210_SCRATCH` | `/gpfs2/scratch/bdigenov` |
| `TA210_ENVS` | `/gpfs2/scratch/bdigenov/envs` |
| `TA210_HOME` | `/users/b/d/bdigenov` |
| `TA210_HOME_REAL` | `/gpfs1/home/b/d/bdigenov` |
| `TA210_ONEDRIVE_ALT` | the authors' workstation folder |

Some scripts also read the project's own variables, such as `PROJECT_ROOT` and the root of each analysis
chain. Each script's header states the defaults of these.

**The three recovered copies.** On 22 September 2026 three AMPKγ copies were recovered from genome-grade
proteomes: *Eimeria necatrix* XP_013440160.1, *Eimeria tenella* XP_013230916.1 and *Porospora gigantea* B
XP_068374756.1. The AMPKγ-dependent results were then carried forward on the updated inputs. The scripts
here are the ones that made the results. The reruns that carried the three copies into them are not
included, and the data deposit holds the final tables.

## 5. Labels in file and column names

Scripts keep the names they ran under, because logs and output folders refer to them.

| Label | Meaning |
|---|---|
| `rl_`, `rl_rm_` | Scripts of the orthology and species-state pipeline, and of the ancestral-state analyses |
| `rbx1` | The AMPKγ contact-position and AMPKα tail measurements |
| `rbx2` | The AMPKβ module score, the activation-site residue and the TORC core states |
| `rbx3` | The upstream-kinase calls (LKB1 and the CAMKK category) |
| `ctdx` | Models of the AMPKα C-terminal region (tail) against the other subunits and the upstream kinases |
| `core18` | The core-complement models: 6 primary models plus sensitivities, 18 fits |
| `network10` | The models of all ordered pairs of the ten targets |
| `core8`, `active10` | The eight TOR and AMPK targets; the ten targets including LKB1 and the CAMKK category |
| `T1`, `T2` | The candidate table and the species-state matrix |
| `Q0`, `Q7` | The per-edge reconstructions; the origins and co-transitions of absence |
| `r59` | The rule for the joint choice of search settings |
| `arm1`, `arm2` | The orthology arm (OrthoFinder, SonicParanoid, KinOrtho) and the domain arm (HMMER, CDD, InterProScan) |
| `decoy` | A same-species protein known not to be the target, used to count false nominations |
| `p45` | The structure-prediction step |
