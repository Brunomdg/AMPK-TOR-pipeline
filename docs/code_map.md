# Code map

This repository holds the code of the eleven methods that are new in this study and that produce the data
of the main results. For each method, section 2 lists the scripts in the order they run. Section 1 goes
table by table and figure by figure through the paper and names the method behind each one's data. Section 5
maps the eight Supplementary Methods subsections to the scripts that implement them, section 6 names the
command line of every external program the study ran, and section 7 lists the operative constants published
as files. Paths are relative to the repository root. Each script begins with a header giving its purpose,
inputs, outputs and programs.

The Materials and Methods of the article is split in two: a main-text Methods, which carries the decision
rules, the thresholds that determine a call, the headline software and the wet-lab protocol, and a
Supplementary Methods, which carries the exact invocations, the reference-site definitions, the sensitivity
variants and the model-by-model specifications. It cites no path in this repository: it refers to this code
only as the analysis code named under Data Availability. The directories here are named for its subsections
so that the two can be read side by side.

Everything else the study ran is not included here, and the authors can provide it on request:
- the drawing of figures and the typesetting of tables;
- checks, including the independent recomputations of the fitted models;
- corrections and reruns;
- supplementary analyses that apply standard tests to the outputs of these methods.

The runs of standard programs were in that list until the command lines were deposited. They are now in
`external_runs/`, one file per program and step, with their versions and evidence in `docs/software.md`. The
programs themselves still have to be installed to rerun a step.

## 1. Tables and figures

**Numbering follows the manuscript**, which is authoritative. An earlier draft numbered the supplement
differently, and this table was rebuilt against the submitted legends; a number here means what the
manuscript's legend gives it, not what an earlier version of this file called it.

"Not included" means the code that built that object is not deposited, for the reasons listed above. The
object itself is in the supplement; the data behind it are in the data deposit.

### Main tables and figures

| Item | Data produced by | Notes |
|---|---|---|
| Table 1 | Methods 3 and 4 | The other rows are not included: the residue-correspondence test, the genome-grade re-examination, the joint loss-direction fits and proteome completeness |
| Figure 1 | Method 2 (the species-state matrix), with the measurements of method 5 | Supergroup shares are counts of that matrix |
| Figure 2 | Methods 6 and 7, with the contact scores of method 5 | The plotted values were prepared by hand from the coccidian AMPKγ copy measurements (Table S4); the per-copy descriptive analyses are not included |
| Figure 3 | Method 4 | The HEAT-region fits against TORC partner loss are not included |
| Figure 4 | Method 5, `ampka_activation_site.py` (activation-site residue), method 10 (the repertoire) and method 9 (the tail against the retained module) | |
| Figure 5 | Method 8, `torc_core_states.py` and `core_complement_models.R` | |
| Figure 6 | Method 11 | Ranges are across reconstruction blocks and are never summed or averaged |
| Figure 7 | Not included | The *T. gondii* TOR-knockout strains and their assays. This repository holds no wet-lab code |

### Supplementary tables

| Item | What the legend gives it | Data produced by |
|---|---|---|
| Table S1 | prevalence by target and supergroup | Methods 1 and 2 |
| Table S2 | component composition by phylum | Methods 1 and 2 |
| Table S3 | descriptive sequence comparisons and coverage | Built on method 5; the descriptive analyses of all copies are not included |
| Table S4 | coccidian AMPKγ copy measurements | Method 6, `coccidian_ampkg_metrics.py` |
| Table S5 | comparison with the earlier single-species scoring | Method 5, `ampkg_contact_scores.py`, which reproduces the original `gamma_pocket.tsv` site-1 scores under the primary procedure and reports every difference with its cause |
| Table S6 | AMP-contact scores per copy | Method 5, `ampkg_contact_scores.py` (`AMPKg_contact_copy.tsv`) |
| Table S7 | substitution rate against burial-matched positions | rate4site, run externally (`4_sequence_phenotypes/external_runs/rate4site.sh`); the 46 runs of the comparison against matched null sets are not included. The site 1 set only |
| Table S8 | TOR HEAT architecture metrics | Method 4: `tor_heat_final_calls.py` for the per-copy call over all 246 accepted copies, on the repeat transfer of `tor_heat_repeat_transfer.py`. The filter to one tree-usable copy per species is applied in the figure, not here |
| Table S9 | mapped length and identity per human HEAT repeat | Method 4, `tor_heat_repeat_transfer.py`, which transfers the 32 UniProt-annotated HEAT repeats of human mTOR (P42345) to every copy through the 246-copy alignment |
| Table S10 | correlated evolution of the HEAT region and TORC partner loss | Not included; the 17 correlated-evolution fits on the dated tree, built on the HEAT calls of method 4 and the states of method 2 |
| Table S11 | identity at TOR partner-contacting surfaces | Not included; built on method 4. The contact residues come from the experimental complexes PDB 6BCX (RAPTOR) and 6ZWM (RICTOR), not from this repository |
| Table S12 | recombination screen | External program only (GARD; `5_phylogenetic_analyses/external_runs/hyphy_gard.sh`) |
| Table S13 | component-state and sequence-feature models | Methods 8 and 9: `core_complement_models.R` with the TORC core states of `torc_core_states.py` for the component-state rows, and `ampka_tail_models.R` and `ampka_tail_upstream_models.R` for the AMPKγ contact, AMPKα tail and AMPKβ module rows, the site-set and alignment-procedure comparisons, and the cross-feature sensitivities. The mapping sensitivities are the declared predictor swaps (`docs/constants/site1_score_input_variants.tsv`) |
| Table S14 | within- and between-clade decomposition | Not included; the post hoc decomposition of the method 9 fit of AMPKα tail length on the AMPKβ module score |
| Table S15 | AMPKα copy measurements | Method 5, `ampka_tail.py` |
| Table S16 | AMPKβ reference-position states | Method 5, `ampkb_module_score.py`, whose long table is one row per copy, position and alignment over the seven structural positions; the reference panel and the scoring rule are in `docs/constants/ampkb_cbm_reference_and_score.md`. The models fitted on those states are method 8 and method 9 |
| Table S17 | component absence against subunit absence | Method 8, `core_complement_models.R`: the nine primary comparisons of TOR, RAPTOR or RICTOR absence with AMPK subunit absence, the same fits the manuscript also carries in Table S13. The species counts, odds ratios and the four-state table are not included |
| Table S18 | all ordered comparisons | Method 8, `ordered_pair_models.R`: the 90 primary ordered comparisons and the sensitivity refits |
| Table S19 | rate dependence between absences | Not included; Pagel's test of dependent against independent rates for all 45 unordered pairs on each of the 28 pruned trees, on the method 2 states |
| Table S20 | candidate replacement comparisons | Not included; the candidate replacements |
| Table S21 | reconstructed origins | Method 11, `absence_origins.py` |
| Table S22 | reconstructed co-transitions | Built on method 11; the co-transition count is not included |
| Table S23 | lifestyle models | Not included; the lifestyle models |
| Table S24 | matched-endpoint comparison with earlier surveys | Not included; the comparison with earlier surveys |
| Table S25 | gene-tree criterion by target | Method 2 (the criterion); its measured performance is not included |
| Table S26 | basis of present cells per target | Built on method 2; not included |
| Table S27 | present cells whose orthologs are all model sensitive | Built on method 2; not included |
| Table S28 | every candidate and its orthology call | Methods 1 and 2, `species_state_matrix.py` (the candidate table) |
| Table S29 | search methods by target | Method 1, `per_method_calls_and_screen.py` |
| Table S30 | every LKB1 candidate | Not included; the LKB1 candidate screen |
| Table S31 | orthology calls by line of evidence | Method 2, `ortholog_table.py`, with the two-arm rule of `two_arm_rule.py` |
| Table S32 | proteome completeness | External programs (BUSCO, OMArk; `1_proteomes_taxonomy_species_tree/external_runs/`); their summary is not included |
| Table S33 | species states | Method 2, `species_state_matrix.py` (the 210-species state matrix, with each species' clade label) |
| Table S34 | AMP-contact calls per position | Method 5, `ampkg_contact_scores.py` (`AMPKg_contact_long.tsv`) |
| Table S35 | AMP-contact reference positions | The frozen panel that method 5 reads (`AMPKg_contact_reference.tsv`, three site sets of ten positions each, in the data deposit). The positions are a 4.0 Å consensus over six crystal structures of the mammalian γ1 subunit, computed outside this repository; the six, the consensus and the numbering check are in `docs/constants/ampkg_amp_site_structures.md` |
| Table S36 | per-species input and output for HEAT divergence against partner loss | Not included; the input and output of the test on methods 4 and 2 |
| Table S37 | numerical estimates of the first additional block | Not included; the block of additional estimates |
| Table S38 | joint three-response model | Not included; mvMORPH, run on the authors' workstation (`docs/software.md`) |
| Table S39 | model comparison for loss with and without regain | Not included; the no-regain comparison |
| Table S40 | broad-group comparisons | Not included |
| Table S41 | descendant-minus-sister contrasts | Not included |
| Table S42 | species-tree markers, supermatrix and topology checks | External programs only (DIAMOND, MAFFT, trimAl, IQ-TREE; `1_proteomes_taxonomy_species_tree/external_runs/`) |
| Table S43 | divergence-time calibrations and the chronogram ensemble | External programs and the lost `chronos` script (`docs/dating_gap.md`); LSD2 and treePL command lines are deposited |
| Table S44 | convergence statistics for the AMPKα tail phenotypes | Not included; the convergence tests. The packages are recorded in the `port68_R` environment (windex 2.1.0, convevol 2.2.1) |

### Supplementary figures

| Item | What the legend gives it | Data produced by |
|---|---|---|
| Figure S1 | repertoires across 37 clade blocks with the species tree | Method 2, with the species tree (`docs/software.md`) |
| Figure S2 | contact-position evidence across all AMPKγ | Built on method 5; the descriptive analyses of all copies are not included |
| Figure S3 | convergence tests | Not included; the convergence tests |
| Figure S4 | alignment, domain evidence and partner multiplicity across all AMPKγ copies | Built on methods 5 and 7 |
| Figure S5 | predicted geometry of the three AMP sites | Method 7, `measure_amp_sites.py` |
| Figure S6 | geometry and confidence across the earlier multimer panel | Not included; the earlier 86-query multimer panel |
| Figure S7 | partner context and variability of the human multimer predictions | Not included; the partner-context benchmark |
| Figure S8 | AMPKα activation-loop site | Method 5, `ampka_activation_site.py` |
| Figure S9 | upstream kinase repertoire and TOR state | Method 10 |
| Figure S10 | retained-module length fit and model-sensitivity status | Method 9; the post hoc models are not included |
| Figure S11 | associations of absence across ten targets | Method 8, `ordered_pair_models.R` |
| Figure S12 | core-disruption associations after adjustment | Method 8, `core_complement_models.R` |
| Figure S13 | lifestyle and proteome size as predictors | Not included; the lifestyle models |
| Figure S14 | reconstructed origins and co-transitions | Method 11; the co-transition count is not included |
| Figure S15 | how each orthology call was reached | Methods 1 and 2; the performance measurements are not included |
| Figure S16 | CAMKK-category and LKB1 trees | Method 10 |

## 2. The methods, with their scripts in run order

External programs and the settings they ran with are in `external_runs/`, with their versions and the
evidence for them in `docs/software.md`. Inputs that these scripts read but do not make will be in the data
deposit:
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

Their command lines are in `2_candidate_discovery_and_orthology/external_runs/`.

1. `2_candidate_discovery_and_orthology/search_settings/score_cuts_per_target.py`: for each target and profile, the pair of sequence
   and domain score cuts that holds the false-positive rate on the target's same-species decoys at or below 5,
   1, 0.5, 0.1 and 0 percent, with the recall of the target's controls at each.
2. `2_candidate_discovery_and_orthology/search_settings/enumerate_arm_settings.py`: every OrthoFinder and
   SonicParanoid setting held on disk, with the declared controls each recovers and the declared negatives it
   admits, which is the frontier the next step solves over.
3. `2_candidate_discovery_and_orthology/search_settings/joint_setting_choice.py`: the joint choice of one setting per method and target,
   so that together the methods recover every declared control with the fewest declared negatives called.
4. `2_candidate_discovery_and_orthology/search_settings/per_method_calls_and_screen.py`: each method's calls at the chosen settings, the candidate
   list, and the contaminant screen in which each protein goes to the best-scoring family.

Helpers:
- `2_candidate_discovery_and_orthology/search_settings/select_arm_settings.py`, which also defines the
  targets, controls, decoys and negative controls the other scripts of this stage import;
- `lib/`: `accessions.py`, `control_declaration.py`, `evidence_roots.py`, `join.py`, `targets.py` and
  `call_vocabulary.py`, which hold the target lists, the accession joins and the call labels.

### Method 2. Gene-tree clade criterion and the species-state matrix

Before these scripts, the gene trees are built externally: MAFFT, IQ-TREE with ModelFinder, the PMSF
cross-check, and RAxML (`2_candidate_discovery_and_orthology/external_runs/`).

1. `2_candidate_discovery_and_orthology/orthology_calls/gene_tree_clade_criterion.py`: the clade criterion of van Dam et al. (2011), anchored on the
   curated controls, applied to each gene tree.
2. `2_candidate_discovery_and_orthology/orthology_calls/two_arm_rule.py`: the two-arm call for candidates the tree left ambiguous.
3. `2_candidate_discovery_and_orthology/orthology_calls/ortholog_table.py`: one call per candidate.
4. `2_candidate_discovery_and_orthology/orthology_calls/core_target_tree_rule.py`: resolves the core candidates that remained open.
5. `3_species_states_and_composition/species_state_matrix.py`: the candidate table and the species-state matrix.

Helpers:
- the shared modules of `lib/`.

### Method 3. Re-examination of absences

1. `2_candidate_discovery_and_orthology/orthology_calls/absence_reexamination.py`: rescores every protein of a species against the eight core
   targets from the genome-wide search, looking for missed orthologs.
2. `2_candidate_discovery_and_orthology/orthology_calls/absence_clade_support.R`: classifies each absence as clade-supported or species-specific,
   on five species-tree topologies and two rootings.

### Method 4. TOR N-terminal HEAT calls

Before these scripts, InterProScan signatures, the AlphaFold models and US-align are run externally
(`2_candidate_discovery_and_orthology/external_runs/interproscan_pfam_panther_smart.sh` and
`4_sequence_phenotypes/external_runs/usalign_tor_heat.sh`, `hhsuite_tor_heat.sh`,
`strpsearch_foldseek.sh`).

1. `4_sequence_phenotypes/tor_heat_signature_screen.py`: HEAT signatures N-terminal to the FAT domain.
2. `4_sequence_phenotypes/tor_heat_repeat_transfer.py`: transfers the 32 N-terminal HEAT repeats of human
   mTOR to every copy.
3. `4_sequence_phenotypes/tor_heat_structural_screen.py`: aligns each model's N-terminal crop to the
   experimental references and applies the structural acceptance rule. The references, the negative controls
   and both gates of the rule are in `docs/constants/tor_heat_structural_references.md`.
4. `4_sequence_phenotypes/tor_heat_final_calls.py`: the final present or absent call, with its confidence
   tier. The accessions it accepts as HEAT evidence are in
   `docs/constants/tor_heat_evidence_accessions.tsv`.

### Method 5. Sequence measurements of the retained subunits

1. `4_sequence_phenotypes/ampkg_contact_scores.py`: scores the AMP-contact positions of every AMPKγ copy
   against the human reference.
2. `4_sequence_phenotypes/ampka_tail.py`: the AMPKα C-terminal length and its predicted disorder.
3. `4_sequence_phenotypes/ampkb_module_score.py`: the AMPKβ carbohydrate-binding module score.
4. `4_sequence_phenotypes/ampka_activation_site.py`: the residue at the AMPKα activation site, per copy
   and per species.

Helpers: `4_sequence_phenotypes/lib/phenotype_tables_chain1.py`, `4_sequence_phenotypes/lib/sequence_scoring.py`,
`4_sequence_phenotypes/lib/phenotype_tables_chain2.py`.

### Method 6. Coccidian AMPKγ copies

1. `4_sequence_phenotypes/coccidian_ampkg_metrics.py`: collects the per-copy
   measurements of every AMPKγ copy and compares the coccidian copies with the others. The measurements are
   length, CBS-domain evidence and identity away from the AMP contacts, the site 1 contact score and its
   agreement between alignments, and predicted-structure confidence and RMSD.

### Method 7. Structure queries and AMP-site geometry

1. `6_structure_prediction_and_pocket/build_prediction_queries.py`: builds each copy's prediction
   sequence so that every observed CBS domain is kept.
2. ColabFold runs externally (`6_structure_prediction_and_pocket/external_runs/`).
3. `6_structure_prediction_and_pocket/measure_amp_sites.py`: measures the three reference AMP sites
   in all 25 models of every prediction.

Helper: `6_structure_prediction_and_pocket/lib/amp_site_measurement.py`, the measurement library. The
structures it measures against are in `docs/constants/ampkg_amp_site_structures.md`.

### Method 8. Penalised phylogenetic logistic regression

1. `3_species_states_and_composition/torc_core_states.py`: the TORC1 and TORC2 core state of every species.
2. `5_phylogenetic_analyses/core_complement_models.R`: AMPK subunit absence on core
   disruption, with the multistart search.
3. `5_phylogenetic_analyses/ordered_pair_models.R`: every ordered pair of the ten targets.

Helpers: `5_phylogenetic_analyses/lib/multistart_search.R` and
`5_phylogenetic_analyses/lib/network10_multistart_helper.R`, which hold the grid of starting
values and the local refinement.

### Method 9. The AMPKα tail against the retained modules and the upstream kinases

1. `5_phylogenetic_analyses/ampka_tail_models.R`: tail length and disorder against
   AMPKγ and AMPKβ absence and the retained AMPKβ module score.
2. `5_phylogenetic_analyses/ampka_tail_upstream_models.R`: the same against absence of the
   upstream kinases.

Helper: `5_phylogenetic_analyses/lib/tail_model_lib.R`.

The declared sensitivities that swap one input for another are listed in
`docs/constants/site1_score_input_variants.tsv` (AMPKγ) and
`docs/constants/ampkb_cbm_reference_and_score.md` (AMPKβ).

The model engine these two load is `5_phylogenetic_analyses/lib/model_engine.R`, which was absent from the
first deposit. It is here now, but the revision that survived is not the one the models were pinned against,
and `ampka_tail_models.R` also reads the record of a reproduction step that is not included.
`docs/deposit_gaps.md` states both limits.

### Method 10. Upstream-kinase classification

1. `2_candidate_discovery_and_orthology/upstream_kinases/alignment_occupancy_check.py`: the alignment quality
   control, read before any topology, which declares the occupancy trigger that decides whether the mask below
   is applied.
2. `2_candidate_discovery_and_orthology/upstream_kinases/mask_kinase_alignment.py`: masks each kinase-domain alignment to the
   columns the reference sequences occupy.
3. MAFFT L-INS-i and IQ-TREE run externally
   (`2_candidate_discovery_and_orthology/external_runs/mafft_iqtree_upstream_kinases.sh`).
4. `2_candidate_discovery_and_orthology/upstream_kinases/tree_tip_classification.R`: assigns each tip by its placement
   among the references.
5. `2_candidate_discovery_and_orthology/upstream_kinases/lkb1_calls.py` for LKB1 and `camkk_calls.py` for the
   CAMKK category: the per-copy calls.
6. `3_species_states_and_composition/upstream_kinase_architecture.py`: each species' upstream-kinase state.

Helper: `2_candidate_discovery_and_orthology/upstream_kinases/lib/upstream_kinase_tables.py`.

### Method 11. Ancestral states and origins of absence

1. `3_species_states_and_composition/state_coding_for_tree_analyses.py`: the standard coding of each target. Present is 0, not
   assessable 0.5, absent 1.
2. `5_phylogenetic_analyses/ancestral_states_per_edge.R`: the state at both ends of every edge, with corHMM (ER, SYM
   and ARD) and parsimony, on every tree, rooting and taxon subset.
3. `5_phylogenetic_analyses/absence_origins.py`: the number of origins of absence in each reconstruction block.

Helpers:
- `5_phylogenetic_analyses/lib/ancestral_state_lib.R`, the repaired reconstruction library;
- `5_phylogenetic_analyses/lib/reconstruction_paths.py`;
- `5_phylogenetic_analyses/lib/output_provenance.R`, which writes the header lines
  of the outputs.

## 3. Paths and the data the three recovered AMPKγ copies changed

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

## 4. The names these scripts ran under

The files are named for what the paper calls them. The names they ran under are kept on a `Ran as:` line in
each script's header, and collected here, because the cluster logs and the output folders of the original
runs refer to them.

| Ran as | Now |
|---|---|
| `rl_threshold_pilot_pertarget.py` | `2_candidate_discovery_and_orthology/search_settings/score_cuts_per_target.py` |
| `rl_arm1_all_settings.py` | `2_candidate_discovery_and_orthology/search_settings/enumerate_arm_settings.py` |
| `rl_arm1_selection.py` | `2_candidate_discovery_and_orthology/search_settings/select_arm_settings.py` |
| `rl_r59_union_optimum.py` | `2_candidate_discovery_and_orthology/search_settings/joint_setting_choice.py` |
| `rl_per_method_report.py` | `2_candidate_discovery_and_orthology/search_settings/per_method_calls_and_screen.py` |
| `rl_phylo_subtree.py` | `2_candidate_discovery_and_orthology/orthology_calls/gene_tree_clade_criterion.py` |
| `rl_ambiguous_arm_rule.py` | `2_candidate_discovery_and_orthology/orthology_calls/two_arm_rule.py` |
| `rl_final_ortholog_table.py` | `2_candidate_discovery_and_orthology/orthology_calls/ortholog_table.py` |
| `rl_core8_tree_rule.py` | `2_candidate_discovery_and_orthology/orthology_calls/core_target_tree_rule.py` |
| `rl_absence_research.py` | `2_candidate_discovery_and_orthology/orthology_calls/absence_reexamination.py` |
| `rl_absence_support.R` | `2_candidate_discovery_and_orthology/orthology_calls/absence_clade_support.R` |
| `rbx3_aln_qc.py` | `2_candidate_discovery_and_orthology/upstream_kinases/alignment_occupancy_check.py` |
| `rbx3_mask.py` | `2_candidate_discovery_and_orthology/upstream_kinases/mask_kinase_alignment.py` |
| `rbx3_tree_calls.R` | `2_candidate_discovery_and_orthology/upstream_kinases/tree_tip_classification.R` |
| `rbx3_calls.py` | `2_candidate_discovery_and_orthology/upstream_kinases/lkb1_calls.py` |
| `rbx3_camkk_calls.py` | `2_candidate_discovery_and_orthology/upstream_kinases/camkk_calls.py` |
| `rbx3_common.py` | `2_candidate_discovery_and_orthology/upstream_kinases/lib/upstream_kinase_tables.py` |
| `rl_pub_t1_t2_active10.py` | `3_species_states_and_composition/species_state_matrix.py` |
| `rbx2_torc_core.py` | `3_species_states_and_composition/torc_core_states.py` |
| `rbx3_architecture.py` | `3_species_states_and_composition/upstream_kinase_architecture.py` |
| `rl_rm_q7_standard.py` | `3_species_states_and_composition/state_coding_for_tree_analyses.py` |
| `rbx1_gamma.py` | `4_sequence_phenotypes/ampkg_contact_scores.py` |
| `rbx1_alpha.py` | `4_sequence_phenotypes/ampka_tail.py` |
| `rbx2_ampkb_cbm.py` | `4_sequence_phenotypes/ampkb_module_score.py` |
| `rbx2_ampka_activation.py` | `4_sequence_phenotypes/ampka_activation_site.py` |
| `build_degeneration_metrics.py` | `4_sequence_phenotypes/coccidian_ampkg_metrics.py` |
| `classify_tor_heat.py` | `4_sequence_phenotypes/tor_heat_signature_screen.py` |
| `refine_nterminal_heat.py` | `4_sequence_phenotypes/tor_heat_repeat_transfer.py` |
| `run_direct_alignment.py` | `4_sequence_phenotypes/tor_heat_structural_screen.py` |
| `resolve_heat_final.py` | `4_sequence_phenotypes/tor_heat_final_calls.py` |
| `rbx1_score.py` | `4_sequence_phenotypes/lib/sequence_scoring.py` |
| `rbx1_lib.py` | `4_sequence_phenotypes/lib/phenotype_tables_chain1.py` |
| `rbx2_common.py` | `4_sequence_phenotypes/lib/phenotype_tables_chain2.py` |
| `core18_penalized.R` | `5_phylogenetic_analyses/core_complement_models.R` |
| `network10_refit.R` | `5_phylogenetic_analyses/ordered_pair_models.R` |
| `ctdx_cross_ampk.R` | `5_phylogenetic_analyses/ampka_tail_models.R` |
| `ctdx_upstream.R` | `5_phylogenetic_analyses/ampka_tail_upstream_models.R` |
| `rl_rm_q0_asr.R` | `5_phylogenetic_analyses/ancestral_states_per_edge.R` |
| `rl_rm_q7_origins.py` | `5_phylogenetic_analyses/absence_origins.py` |
| `core_multistart_helper.R`, `network10_multistart_helper.R` | `5_phylogenetic_analyses/lib/multistart_search.R` (one file; the two were the same code) |
| `rbx1_models.R` | `5_phylogenetic_analyses/lib/model_engine.R` |
| `ctdx_lib.R` | `5_phylogenetic_analyses/lib/tail_model_lib.R` |
| `rl_rm_asr_lib.R` | `5_phylogenetic_analyses/lib/ancestral_state_lib.R` |
| `rl_rm_q7_lib.py` | `5_phylogenetic_analyses/lib/reconstruction_paths.py` |
| `rl_rm_prov.R` | `5_phylogenetic_analyses/lib/output_provenance.R` |
| `prepare_ampkg_domain_preserving.py` | `6_structure_prediction_and_pocket/build_prediction_queries.py` |
| `analyze_all_ampkg_structures.py` | `6_structure_prediction_and_pocket/measure_amp_sites.py` |
| `p45b_colabfold_analysis.py` | `6_structure_prediction_and_pocket/lib/amp_site_measurement.py` |
| `rl_targets.py`, `rl_accession.py`, `rl_declaration.py`, `rl_vocab.py`, `rl_evidence_roots.py`, `rl_join.py` | `lib/targets.py`, `lib/accessions.py`, `lib/control_declaration.py`, `lib/call_vocabulary.py`, `lib/evidence_roots.py`, `lib/join.py` |

Column names and output folders are unchanged, and some of them carry the old labels: `core18` is the
core-complement models (6 primary fits plus sensitivities, 18 in all), `network10` the models of all ordered
pairs of the ten targets, `r59` the rule for the joint choice of search settings, `arm1` and `arm2` the
orthology arm (OrthoFinder, SonicParanoid, KinOrtho) and the domain arm (HMMER, CDD, InterProScan), `T1` and
`T2` the candidate table and the species-state matrix, `Q0` and `Q7` the per-edge reconstructions and the
origins of absence, `core8` and `active10` the eight TOR and AMPK targets and the ten including LKB1 and the
CAMKK category, and `decoy` a same-species protein known not to be the target.

## 5. The Supplementary Methods, subsection by subsection

The eight subsections of the Supplementary Methods, with the code that implements each. The mapping is by
content. Where a script's own `# Methods:` header names a different subsection, the header is the earlier
internal label and the row below says so.

### Proteomes, taxonomy and the species tree

No deposited scripts: every step is an external program. The command lines are in
`1_proteomes_taxonomy_species_tree/external_runs/`, covering BUSCO, OMAmer and OMArk, DIAMOND, MAFFT,
trimAl, IQ-TREE, LSD2 and treePL. The dating of record cannot be reproduced; `docs/dating_gap.md` states
what is missing.

### Candidate discovery and orthology calls

| Script | What it does |
|---|---|
| `2_candidate_discovery_and_orthology/search_settings/score_cuts_per_target.py` | the per-target score cuts at each false-positive rate |
| `2_candidate_discovery_and_orthology/search_settings/joint_setting_choice.py` | the joint choice of one setting per method and target |
| `2_candidate_discovery_and_orthology/search_settings/per_method_calls_and_screen.py` | each method's calls, the candidate list and the contaminant screen |
| `2_candidate_discovery_and_orthology/search_settings/select_arm_settings.py` | the OrthoFinder and SonicParanoid setting per target |
| `2_candidate_discovery_and_orthology/orthology_calls/gene_tree_clade_criterion.py` | the gene-tree clade criterion |
| `2_candidate_discovery_and_orthology/orthology_calls/two_arm_rule.py` | the two-arm rule where the tree abstains |
| `2_candidate_discovery_and_orthology/orthology_calls/ortholog_table.py` | one call per candidate |
| `2_candidate_discovery_and_orthology/orthology_calls/core_target_tree_rule.py` | the core candidates that remained open |
| `2_candidate_discovery_and_orthology/orthology_calls/absence_reexamination.py` | the genome-grade re-examination of absences |
| `2_candidate_discovery_and_orthology/orthology_calls/absence_clade_support.R` | clade support for each absence |
| `2_candidate_discovery_and_orthology/upstream_kinases/mask_kinase_alignment.py`, `tree_tip_classification.R`, `lkb1_calls.py`, `camkk_calls.py` | the LKB1 and CAMKK-category calls; their headers say "Orthology calls", and they sit under `07_` only because they share that chain's helpers |
| `lib/` | the target lists, accession joins and call labels |

External programs: `2_candidate_discovery_and_orthology/external_runs/`.

### Species states and composition

| Script | What it does |
|---|---|
| `3_species_states_and_composition/species_state_matrix.py` | the candidate table and the 210-species state matrix |
| `3_species_states_and_composition/torc_core_states.py` | the TORC1 and TORC2 core state of every species; its header says "Phylogenetic regressions", the chain it was written for |
| `3_species_states_and_composition/upstream_kinase_architecture.py` | each species' upstream-kinase state; its header says "Phylogenetic regressions" |
| `3_species_states_and_composition/state_coding_for_tree_analyses.py` | the standard coding of each target for the evolutionary analyses |

No external programs.

### Sequence phenotypes of the present subunits

| Script | What it does |
|---|---|
| `4_sequence_phenotypes/ampkg_contact_scores.py` | the AMPKγ AMP-contact scores, per copy and per position |
| `4_sequence_phenotypes/ampka_tail.py` | the AMPKα C-terminal length and disorder |
| `4_sequence_phenotypes/ampkb_module_score.py` | the AMPKβ carbohydrate-binding module score |
| `4_sequence_phenotypes/ampka_activation_site.py` | the AMPKα activation-site residue; its header says "Phylogenetic regressions" |
| `4_sequence_phenotypes/coccidian_ampkg_metrics.py` | the coccidian AMPKγ comparison |
| `4_sequence_phenotypes/tor_heat_signature_screen.py`, `tor_heat_repeat_transfer.py` | the TOR N-terminal HEAT signature screen and the repeat transfer |
| `4_sequence_phenotypes/tor_heat_structural_screen.py` | the structural acceptance rule |
| `4_sequence_phenotypes/tor_heat_final_calls.py` | the final HEAT call and its confidence tier |
| Helpers | `phenotype_tables_chain1.py`, `sequence_scoring.py`, `phenotype_tables_chain2.py` |

External programs: `4_sequence_phenotypes/external_runs/`. Constants:
`docs/constants/tor_heat_evidence_accessions.tsv`, `docs/constants/ampkb_cbm_reference_and_score.md`,
`docs/constants/site1_score_input_variants.tsv`.

### Phylogenetic analyses

| Script | What it does |
|---|---|
| `5_phylogenetic_analyses/core_complement_models.R` | subunit absence on core disruption |
| `5_phylogenetic_analyses/ordered_pair_models.R` | every ordered pair of the ten targets |
| `5_phylogenetic_analyses/ampka_tail_models.R`, `ampka_tail_upstream_models.R` | the AMPKα tail models |
| `5_phylogenetic_analyses/ancestral_states_per_edge.R` | the state at both ends of every edge |
| `5_phylogenetic_analyses/absence_origins.py` | the origins of absence per reconstruction block |
| Helpers | the two multistart helpers, `tail_model_lib.R`, `ancestral_state_lib.R`, `reconstruction_paths.py`, `output_provenance.R` |

External programs: `5_phylogenetic_analyses/external_runs/` (GARD). The R packages that fit the models
are in `docs/software.md`, not here, because the scripts call them directly.

### Structure prediction and AMPKγ pocket geometry

| Script | What it does |
|---|---|
| `6_structure_prediction_and_pocket/build_prediction_queries.py` | the domain-preserving prediction sequences |
| `6_structure_prediction_and_pocket/measure_amp_sites.py` | the three AMP sites in all 25 models per prediction |
| `6_structure_prediction_and_pocket/lib/amp_site_measurement.py` | the measurement library |

External programs: `6_structure_prediction_and_pocket/external_runs/`. Constants:
`docs/constants/ampkg_amp_site_structures.md`.

### Comparison with earlier surveys

| Script | What it does |
|---|---|
| `4_sequence_phenotypes/ampkg_contact_scores.py` | reproduces the earlier single-species site-1 scores under the primary procedure and reports every difference with its cause (`AMPKg_legacy_reproduction.tsv`, `AMPKg_legacy_pairwise_reproduction_delta.tsv`) |

The matched-endpoint comparison with the earlier surveys is not included. No external programs.

### Development of *Toxoplasma gondii* TOR-knockout strains

Wet-lab work. No code, no external programs. The protocol is in the main-text Methods.

## 6. External program runs

`external_runs/` holds the command line of every program the study ran that is not itself deposited as a
script, one file per program and step, each naming the Supplementary Methods subsection it serves.

| Directory | Supplementary Methods subsection |
|---|---|
| `1_proteomes_taxonomy_species_tree/external_runs/` | Proteomes, taxonomy and the species tree |
| `2_candidate_discovery_and_orthology/external_runs/` | Candidate discovery and orthology calls |
| `4_sequence_phenotypes/external_runs/` | Sequence phenotypes of the present subunits |
| `5_phylogenetic_analyses/external_runs/` | Phylogenetic analyses |
| `6_structure_prediction_and_pocket/external_runs/` | Structure prediction and AMPKγ pocket geometry |

The job scripts themselves are not deposited. Each file states whether its command line was transcribed from
`docs/software.md`, which copied it from the job script while it still existed, or is issued by a deposited
script, which is the stronger evidence of the two. Where neither holds, the file says so and gives the
documented procedure rather than a flag string. `docs/external_runs.md` explains the convention.

One step has no command line at all: the `chronos` dating of record, whose script is lost
(`docs/dating_gap.md`).

## 7. Operative constants

The values that decide a call, published as files so that they can be read without opening a script. The
script named in each row remains the authority.

| File | What it holds | From |
|---|---|---|
| `docs/constants/tor_heat_evidence_accessions.tsv` | the Pfam, SMART, PROSITE, SUPERFAMILY, Gene3D and FunFam accessions that count as HEAT evidence, by class, with the tier rule | `4_sequence_phenotypes/tor_heat_final_calls.py` |
| `docs/constants/ampkb_cbm_reference_and_score.md` | the AMPKβ reference, its numbering, the module envelope, the two site sets and the scoring rule | `4_sequence_phenotypes/ampkb_module_score.py` |
| `docs/constants/ampkg_amp_site_structures.md` | the six crystal structures of the 4.0 Å position consensus, the numbering check, and the separate 5.0 Å geometry measurement against `2V8Q` | the Supplementary Methods for the six; `6_structure_prediction_and_pocket/lib/amp_site_measurement.py` for the measurement |
| `docs/constants/tor_heat_structural_references.md` | the experimental references and negative controls of the TOR HEAT structural screen, and both gates of its acceptance rule | `4_sequence_phenotypes/tor_heat_structural_screen.py` |
| `docs/constants/site1_score_input_variants.tsv` | the three input variants of the site-1 score and the primary input they replace | `5_phylogenetic_analyses/ampka_tail_models.R` |

One of these states a limit of the records rather than a value: the AMPKβ residue positions are read at run
time from a frozen panel that is in the data deposit, not here.
