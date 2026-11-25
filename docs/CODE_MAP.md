# Code map

This repository holds the code of the eleven methods that are new in this study and that produce the data
of the main results. For each method, section 2 lists the scripts in the order they run. Section 1 goes
table by table and figure by figure through the paper and names the method behind each one's data. Section 6
maps the eight Supplementary Methods subsections to the scripts that implement them, section 7 names the
command line of every external program the study ran, and section 8 lists the operative constants published
as files. Paths are relative to the repository root. Each script begins with a header giving its purpose,
inputs, outputs and programs.

The Materials and Methods of the article is split in two: a main-text Methods, which carries the decision
rules, the thresholds that determine a call, the headline software and the wet-lab protocol, and a
Supplementary Methods, which carries the exact invocations, the reference-site definitions, the sensitivity
variants and the model-by-model specifications. The Supplementary Methods cites paths in this repository.

Everything else the study ran is not included here, and the authors can provide it on request:
- the drawing of figures and the typesetting of tables;
- checks, including the independent recomputations of the fitted models;
- corrections and reruns;
- supplementary analyses that apply standard tests to the outputs of these methods.

The runs of standard programs were in that list until the command lines were deposited. They are now in
`external_runs/`, one file per program and step, with their versions and evidence in `docs/SOFTWARE.md`. The
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
| Figure 4 | Method 5, `rbx2_ampka_activation.py` (activation-site residue), method 10 (the repertoire) and method 9 (the tail against the retained module) | |
| Figure 5 | Method 8, `rbx2_torc_core.py` and `core18_penalized.R` | |
| Figure 6 | Method 11 | Ranges are across reconstruction blocks and are never summed or averaged |
| Figure 7 | Not included | The *T. gondii* TOR-knockout strains and their assays. This repository holds no wet-lab code |

### Supplementary tables

| Item | What the legend gives it | Data produced by |
|---|---|---|
| Table S1 | prevalence by target and supergroup | Methods 1 and 2 |
| Table S2 | component composition by phylum | Methods 1 and 2 |
| Table S3 | descriptive sequence comparisons and coverage | Built on method 5; the descriptive analyses of all copies are not included |
| Table S4 | coccidian AMPKγ copy measurements | Method 6, `build_degeneration_metrics.py` |
| Table S5 | comparison with the earlier single-species scoring | Method 5, `rbx1_gamma.py`, which reproduces the original `gamma_pocket.tsv` site-1 scores under the primary procedure and reports every difference with its cause |
| Table S6 | AMP-contact scores per copy | Method 5, `rbx1_gamma.py` (`AMPKg_contact_copy.tsv`) |
| Table S7 | substitution rate against burial-matched positions | rate4site, run externally (`external_runs/sm4_sequence_phenotypes/rate4site.sh`); the 46 runs of the comparison against matched null sets are not included. The site 1 set only |
| Table S8 | TOR HEAT architecture metrics | Method 4: `resolve_heat_final.py` for the per-copy call over all 246 accepted copies, on the repeat transfer of `refine_nterminal_heat.py`. The filter to one tree-usable copy per species is applied in the figure, not here |
| Table S9 | mapped length and identity per human HEAT repeat | Method 4, `refine_nterminal_heat.py`, which transfers the 32 UniProt-annotated HEAT repeats of human mTOR (P42345) to every copy through the 246-copy alignment |
| Table S10 | correlated evolution of the HEAT region and TORC partner loss | Not included; the 17 correlated-evolution fits on the dated tree, built on the HEAT calls of method 4 and the states of method 2 |
| Table S11 | identity at TOR partner-contacting surfaces | Not included; built on method 4. The contact residues come from the experimental complexes PDB 6BCX (RAPTOR) and 6ZWM (RICTOR), not from this repository |
| Table S12 | recombination screen | External program only (GARD; `external_runs/sm5_phylogenetic_analyses/hyphy_gard.sh`) |
| Table S13 | component-state and sequence-feature models | Methods 8 and 9: `core18_penalized.R` with the TORC core states of `rbx2_torc_core.py` for the component-state rows, and `ctdx_cross_ampk.R` and `ctdx_upstream.R` for the AMPKγ contact, AMPKα tail and AMPKβ module rows, the site-set and alignment-procedure comparisons, and the cross-feature sensitivities. The mapping sensitivities are the declared predictor swaps (`docs/constants/site1_score_input_variants.tsv`) |
| Table S14 | within- and between-clade decomposition | Not included; the post hoc decomposition of the method 9 fit of AMPKα tail length on the AMPKβ module score |
| Table S15 | AMPKα copy measurements | Method 5, `rbx1_alpha.py` |
| Table S16 | AMPKβ reference-position states | Method 5, `rbx2_ampkb_cbm.py`, whose long table is one row per copy, position and alignment over the seven structural positions; the reference panel and the scoring rule are in `docs/constants/ampkb_cbm_reference_and_score.md`. The models fitted on those states are method 8 and method 9 |
| Table S17 | component absence against subunit absence | Method 8, `core18_penalized.R`: the nine primary comparisons of TOR, RAPTOR or RICTOR absence with AMPK subunit absence, the same fits the manuscript also carries in Table S13. The species counts, odds ratios and the four-state table are not included |
| Table S18 | all ordered comparisons | Method 8, `network10_refit.R`: the 90 primary ordered comparisons and the sensitivity refits |
| Table S19 | rate dependence between absences | Not included; Pagel's test of dependent against independent rates for all 45 unordered pairs on each of the 28 pruned trees, on the method 2 states |
| Table S20 | candidate replacement comparisons | Not included; the candidate replacements |
| Table S21 | reconstructed origins | Method 11, `rl_rm_q7_origins.py` |
| Table S22 | reconstructed co-transitions | Built on method 11; the co-transition count is not included |
| Table S23 | lifestyle models | Not included; the lifestyle models |
| Table S24 | matched-endpoint comparison with earlier surveys | Not included; the comparison with earlier surveys |
| Table S25 | gene-tree criterion by target | Method 2 (the criterion); its measured performance is not included |
| Table S26 | basis of present cells per target | Built on method 2; not included |
| Table S27 | present cells whose orthologs are all model sensitive | Built on method 2; not included |
| Table S28 | every candidate and its orthology call | Methods 1 and 2, `rl_pub_t1_t2_active10.py` (the candidate table) |
| Table S29 | search methods by target | Method 1, `rl_per_method_report.py` |
| Table S30 | every LKB1 candidate | Not included; the LKB1 candidate screen |
| Table S31 | orthology calls by line of evidence | Method 2, `rl_final_ortholog_table.py`, with the two-arm rule of `rl_ambiguous_arm_rule.py` |
| Table S32 | proteome completeness | External programs (BUSCO, OMArk; `external_runs/sm1_proteomes_taxonomy_species_tree/`); their summary is not included |
| Table S33 | species states | Method 2, `rl_pub_t1_t2_active10.py` (the 210-species state matrix, with each species' clade label) |
| Table S34 | AMP-contact calls per position | Method 5, `rbx1_gamma.py` (`AMPKg_contact_long.tsv`) |
| Table S35 | AMP-contact reference positions | The frozen panel that method 5 reads (`AMPKg_contact_reference.tsv`, three site sets of ten positions each, in the data deposit). The positions are a 4.0 Å consensus over six crystal structures of the mammalian γ1 subunit, computed outside this repository; the six, the consensus and the numbering check are in `docs/constants/ampkg_amp_site_structures.md` |
| Table S36 | per-species input and output for HEAT divergence against partner loss | Not included; the input and output of the test on methods 4 and 2 |
| Table S37 | numerical estimates of the first additional block | Not included; the block of additional estimates |
| Table S38 | joint three-response model | Not included; mvMORPH, run on the authors' workstation (`docs/SOFTWARE.md`) |
| Table S39 | model comparison for loss with and without regain | Not included; the no-regain comparison |
| Table S40 | broad-group comparisons | Not included |
| Table S41 | descendant-minus-sister contrasts | Not included |

### Supplementary figures

| Item | What the legend gives it | Data produced by |
|---|---|---|
| Figure S1 | repertoires across 37 clade blocks with the species tree | Method 2, with the species tree (`docs/SOFTWARE.md`) |
| Figure S2 | contact-position evidence across all AMPKγ | Built on method 5; the descriptive analyses of all copies are not included |
| Figure S3 | convergence tests | Not included; the convergence tests |
| Figure S4 | alignment, domain evidence and partner multiplicity across all AMPKγ copies | Built on methods 5 and 7 |
| Figure S5 | predicted geometry of the three AMP sites | Method 7, `analyze_all_ampkg_structures.py` |
| Figure S6 | geometry and confidence across the earlier multimer panel | Not included; the earlier 86-query multimer panel |
| Figure S7 | partner context and variability of the human multimer predictions | Not included; the partner-context benchmark |
| Figure S8 | AMPKα activation-loop site | Method 5, `rbx2_ampka_activation.py` |
| Figure S9 | upstream kinase repertoire and TOR state | Method 10 |
| Figure S10 | retained-module length fit and model-sensitivity status | Method 9; the post hoc models are not included |
| Figure S11 | associations of absence across ten targets | Method 8, `network10_refit.R` |
| Figure S12 | core-disruption associations after adjustment | Method 8, `core18_penalized.R` |
| Figure S13 | lifestyle and proteome size as predictors | Not included; the lifestyle models |
| Figure S14 | reconstructed origins and co-transitions | Method 11; the co-transition count is not included |
| Figure S15 | how each orthology call was reached | Methods 1 and 2; the performance measurements are not included |
| Figure S16 | CAMKK-category and LKB1 trees | Method 10 |

## 2. The methods, with their scripts in run order

External programs and the settings they ran with are in `external_runs/`, with their versions and the
evidence for them in `docs/SOFTWARE.md`. Inputs that these scripts read but do not make will be in the data
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

Their command lines are in `external_runs/sm2_candidate_discovery_and_orthology/`.

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
cross-check, and RAxML (`external_runs/sm2_candidate_discovery_and_orthology/`).

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

Before these scripts, InterProScan signatures, the AlphaFold models and US-align are run externally
(`external_runs/sm2_candidate_discovery_and_orthology/interproscan_pfam_panther_smart.sh` and
`external_runs/sm4_sequence_phenotypes/usalign_tor_heat.sh`, `hhsuite_tor_heat.sh`,
`strpsearch_foldseek.sh`).

1. `06_tor_heat/signature_screen/classify_tor_heat.py`: HEAT signatures N-terminal to the FAT domain.
2. `06_tor_heat/signature_screen/refine_nterminal_heat.py`: transfers the 32 N-terminal HEAT repeats of human
   mTOR to every copy.
3. `06_tor_heat/structural_screen/run_direct_alignment.py`: aligns each model's N-terminal crop to the
   experimental references and applies the structural acceptance rule. The references, the negative controls
   and both gates of the rule are in `docs/constants/tor_heat_structural_references.md`.
4. `06_tor_heat/final_calls/resolve_heat_final.py`: the final present or absent call, with its confidence
   tier. The accessions it accepts as HEAT evidence are in
   `docs/constants/tor_heat_evidence_accessions.tsv`.

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
2. ColabFold runs externally (`external_runs/sm6_structure_prediction_and_pocket/`).
3. `11_structure_pocket/all_copies/analyze_all_ampkg_structures.py`: measures the three reference AMP sites
   in all 25 models of every prediction.

Helper: `11_structure_pocket/panel_86_queries/p45b_colabfold_analysis.py`, the measurement library. The
structures it measures against are in `docs/constants/ampkg_amp_site_structures.md`.

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

The declared sensitivities that swap one input for another are listed in
`docs/constants/site1_score_input_variants.tsv` (AMPKγ) and
`docs/constants/ampkb_cbm_reference_and_score.md` (AMPKβ).

**These two scripts cannot run as deposited.** `ctdx_lib.R` checks the checksum of the model engine it
reads, and that version of the engine was not preserved. `ctdx_cross_ampk.R` also reads the record of a
reproduction step that is not included.

### Method 10. Upstream-kinase classification

1. `07_phylogenetic_regressions/upstream_kinases/rbx3_mask.py`: masks each kinase-domain alignment to the
   columns the reference sequences occupy.
2. MAFFT L-INS-i and IQ-TREE run externally
   (`external_runs/sm2_candidate_discovery_and_orthology/mafft_iqtree_upstream_kinases.sh`).
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

## 6. The Supplementary Methods, subsection by subsection

The eight subsections of the Supplementary Methods, with the code that implements each. The mapping is by
content. Where a script's own `# Methods:` header names a different subsection, the header is the earlier
internal label and the row below says so.

### Proteomes, taxonomy and the species tree

No deposited scripts: every step is an external program. The command lines are in
`external_runs/sm1_proteomes_taxonomy_species_tree/`, covering BUSCO, OMAmer and OMArk, DIAMOND, MAFFT,
trimAl, IQ-TREE, LSD2 and treePL. The dating of record cannot be reproduced; `docs/DATING_GAP.md` states
what is missing.

### Candidate discovery and orthology calls

| Script | What it does |
|---|---|
| `02_candidate_discovery/rl_threshold_pilot_pertarget.py` | the per-target score cuts at each false-positive rate |
| `02_candidate_discovery/rl_r59_union_optimum.py` | the joint choice of one setting per method and target |
| `02_candidate_discovery/rl_per_method_report.py` | each method's calls, the candidate list and the contaminant screen |
| `02_candidate_discovery/rl_arm1_selection.py` | the OrthoFinder and SonicParanoid setting per target |
| `02_candidate_discovery/rl_invariants.py` | the checksum gate on the reference inputs |
| `03_orthology_calls/rl_phylo_subtree.py` | the gene-tree clade criterion |
| `03_orthology_calls/rl_ambiguous_arm_rule.py` | the two-arm rule where the tree abstains |
| `03_orthology_calls/rl_final_ortholog_table.py` | one call per candidate |
| `03_orthology_calls/rl_core8_tree_rule.py` | the core candidates that remained open |
| `03_orthology_calls/rl_absence_research.py` | the genome-grade re-examination of absences |
| `03_orthology_calls/rl_absence_support.R` | clade support for each absence |
| `03_orthology_calls/lkb1_camkk/` (2 scripts) | the freeze gates on the LKB1 and CAMKK calls |
| `07_phylogenetic_regressions/upstream_kinases/rbx3_mask.py`, `rbx3_tree_calls.R`, `rbx3_calls.py`, `rbx3_camkk_calls.py` | the LKB1 and CAMKK-category calls; their headers say "Orthology calls", and they sit under `07_` only because they share that chain's helpers |
| `shared/pipeline_lib/` | the target lists, accession joins and call labels |

External programs: `external_runs/sm2_candidate_discovery_and_orthology/`.

### Species states and composition

| Script | What it does |
|---|---|
| `03_orthology_calls/rl_pub_t1_t2_active10.py` | the candidate table and the 210-species state matrix |
| `05_sequence_phenotypes/rbx2/rbx2_torc_core.py` | the TORC1 and TORC2 core state of every species; its header says "Phylogenetic regressions", the chain it was written for |
| `07_phylogenetic_regressions/upstream_kinases/rbx3_architecture.py` | each species' upstream-kinase state; its header says "Phylogenetic regressions" |
| `08_ancestral_states/rl_rm_q7_standard.py` | the standard coding of each target for the evolutionary analyses |

No external programs.

### Sequence phenotypes of the retained subunits

| Script | What it does |
|---|---|
| `05_sequence_phenotypes/rbx1_chain/rbx1_gamma.py` | the AMPKγ AMP-contact scores, per copy and per position |
| `05_sequence_phenotypes/rbx1_chain/rbx1_alpha.py` | the AMPKα C-terminal length and disorder |
| `05_sequence_phenotypes/rbx2/rbx2_ampkb_cbm.py` | the AMPKβ carbohydrate-binding module score |
| `05_sequence_phenotypes/rbx2/rbx2_ampka_activation.py` | the AMPKα activation-site residue; its header says "Phylogenetic regressions" |
| `05_sequence_phenotypes/coccidian_remodelling/build_degeneration_metrics.py` | the coccidian AMPKγ comparison |
| `06_tor_heat/signature_screen/classify_tor_heat.py`, `refine_nterminal_heat.py` | the TOR N-terminal HEAT signature screen and the repeat transfer |
| `06_tor_heat/structural_screen/run_direct_alignment.py` | the structural acceptance rule |
| `06_tor_heat/final_calls/resolve_heat_final.py` | the final HEAT call and its confidence tier |
| Helpers | `rbx1_lib.py`, `rbx1_score.py`, `rbx2_common.py` |

External programs: `external_runs/sm4_sequence_phenotypes/`. Constants:
`docs/constants/tor_heat_evidence_accessions.tsv`, `docs/constants/ampkb_cbm_reference_and_score.md`,
`docs/constants/site1_score_input_variants.tsv`.

### Phylogenetic analyses

| Script | What it does |
|---|---|
| `07_phylogenetic_regressions/core_complements/core18_penalized.R` | subunit absence on core disruption |
| `07_phylogenetic_regressions/ten_target_pairs/network10_refit.R` | every ordered pair of the ten targets |
| `07_phylogenetic_regressions/ctd_cross_feature/ctdx_cross_ampk.R`, `ctdx_upstream.R` | the AMPKα tail models |
| `08_ancestral_states/repair/rl_rm_q0_asr.R` | the state at both ends of every edge |
| `08_ancestral_states/rl_rm_q7_origins.py` | the origins of absence per reconstruction block |
| Helpers | the two multistart helpers, `ctdx_lib.R`, `rl_rm_asr_lib.R`, `rl_rm_q7_lib.py`, `rl_rm_prov.R` |

External programs: `external_runs/sm5_phylogenetic_analyses/` (GARD). The R packages that fit the models
are in `docs/SOFTWARE.md`, not here, because the scripts call them directly.

### Structure prediction and AMPKγ pocket geometry

| Script | What it does |
|---|---|
| `11_structure_pocket/all_copies/prepare_ampkg_domain_preserving.py` | the domain-preserving prediction sequences |
| `11_structure_pocket/all_copies/analyze_all_ampkg_structures.py` | the three AMP sites in all 25 models per prediction |
| `11_structure_pocket/panel_86_queries/p45b_colabfold_analysis.py` | the measurement library |

External programs: `external_runs/sm6_structure_prediction_and_pocket/`. Constants:
`docs/constants/ampkg_amp_site_structures.md`.

### Comparison with earlier surveys

| Script | What it does |
|---|---|
| `05_sequence_phenotypes/rbx1_chain/rbx1_gamma.py` | reproduces the earlier single-species site-1 scores under the primary procedure and reports every difference with its cause (`AMPKg_legacy_reproduction.tsv`, `AMPKg_legacy_pairwise_reproduction_delta.tsv`) |

The matched-endpoint comparison with the earlier surveys is not included. No external programs.

### Development of *Toxoplasma gondii* TOR-knockout strains

Wet-lab work. No code, no external programs. The protocol is in the main-text Methods.

## 7. External program runs

`external_runs/` holds the command line of every program the study ran that is not itself deposited as a
script, one file per program and step, each naming the Supplementary Methods subsection it serves.

| Directory | Supplementary Methods subsection |
|---|---|
| `external_runs/sm1_proteomes_taxonomy_species_tree/` | Proteomes, taxonomy and the species tree |
| `external_runs/sm2_candidate_discovery_and_orthology/` | Candidate discovery and orthology calls |
| `external_runs/sm4_sequence_phenotypes/` | Sequence phenotypes of the retained subunits |
| `external_runs/sm5_phylogenetic_analyses/` | Phylogenetic analyses |
| `external_runs/sm6_structure_prediction_and_pocket/` | Structure prediction and AMPKγ pocket geometry |

The job scripts themselves are not deposited. Each file states whether its command line was transcribed from
`docs/SOFTWARE.md`, which copied it from the job script while it still existed, or is issued by a deposited
script, which is the stronger evidence of the two. Where neither holds, the file says so and gives the
documented procedure rather than a flag string. `external_runs/README.md` explains the convention.

One step has no command line at all: the `chronos` dating of record, whose script is lost
(`docs/DATING_GAP.md`).

## 8. Operative constants

The values that decide a call, published as files so that they can be read without opening a script. The
script named in each row remains the authority.

| File | What it holds | From |
|---|---|---|
| `docs/constants/tor_heat_evidence_accessions.tsv` | the Pfam, SMART, PROSITE, SUPERFAMILY, Gene3D and FunFam accessions that count as HEAT evidence, by class, with the tier rule | `06_tor_heat/final_calls/resolve_heat_final.py` |
| `docs/constants/ampkb_cbm_reference_and_score.md` | the AMPKβ reference, its numbering, the module envelope, the two site sets and the scoring rule | `05_sequence_phenotypes/rbx2/rbx2_ampkb_cbm.py` |
| `docs/constants/ampkg_amp_site_structures.md` | the six crystal structures of the 4.0 Å position consensus, the numbering check, and the separate 5.0 Å geometry measurement against `2V8Q` | the Supplementary Methods for the six; `11_structure_pocket/panel_86_queries/p45b_colabfold_analysis.py` for the measurement |
| `docs/constants/tor_heat_structural_references.md` | the experimental references and negative controls of the TOR HEAT structural screen, and both gates of its acceptance rule | `06_tor_heat/structural_screen/run_direct_alignment.py` |
| `docs/constants/site1_score_input_variants.tsv` | the three input variants of the site-1 score and the primary input they replace | `07_phylogenetic_regressions/ctd_cross_feature/ctdx_cross_ampk.R` |

One of these states a limit of the records rather than a value: the AMPKβ residue positions are read at run
time from a frozen panel that is in the data deposit, not here.
