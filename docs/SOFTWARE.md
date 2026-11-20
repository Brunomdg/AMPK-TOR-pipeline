# Software and versions

This page lists every program and package the analyses used, with its version and the evidence for that
version. Where a program was run by a job script that is not part of this repository, the settings given
are copied from that job script's command line. Those command lines are also deposited one per program and
step in `external_runs/`, where each file names the Supplementary Methods subsection it serves; this page
remains the record of versions and evidence.

Where a deposited script issues the call itself, the script is the authority and this page follows it.

The analyses ran on the VACC cluster of the University of Vermont (Red Hat Enterprise Linux 9), from conda
environments that the job scripts call by absolute path, for example `${TA210_ENVS}/ms210_phase2/bin/Rscript`.
A few steps ran on the authors' workstation; they are marked below.

The complete package lists are in `envs/`:

| File | Content |
|---|---|
| `envs/PACKAGES.tsv` | Every conda, Python and R package of every environment, read on 25 September 2026 after the last analysis: environment, kind, name, version, build, channel, and the R build line |
| `envs/PACKAGES_20260913.tsv` | The `conda list` and `pip freeze` exports of the environments that existed on 13 September 2026, while the analyses were running |
| `envs/CONDA_HISTORY.txt` | Each environment's own log of its `conda create` and `conda install` commands, with dates |
| `envs/R_SESSIONS.tsv` | A summary of the 494 `sessionInfo()` records that the R steps wrote beside their outputs |

A package installed after a given run also appears in the 25 September lists. The dated history shows when
each package arrived.

Evidence for a version:
- **log**: the run's own log;
- **env**: the environment record in `envs/`;
- **banner**: printed or fixed by the job script;
- **pinned**: the script stops unless it finds this version;
- **sessionInfo**: the R session records;
- **probe**: a one-off software check, not a per-run log;
- **script**: the job script or a deposited script fixes it;
- **source clone**: the commit of the clone that was built or run, for a tool with no release number.

## Environments

| Environment | Main contents |
|---|---|
| `ms210_phase2` | R 4.5.3 with phylolm 2.6.5, ape 5.8-1, phytools 2.5-2, geiger 2.0.11, corHMM 2.8, phangorn 2.12.1; Python 3.12 with pandas, numpy, Biopython 1.88 |
| `busco6` | BUSCO 6.1.0, HMMER 3.4, MAFFT 7.526, BLAST+ 2.17.0, DIAMOND 2.2.5, RAxML 8.2.13, FastTree 2.2.0 |
| `orthofinder3` | OrthoFinder 3.1.5, IQ-TREE 3.1.3, MAFFT 7.526, HMMER 3.4, DIAMOND 2.2.4, MCL 22.282, RAxML 8.2.13 |
| `ms210` | IQ-TREE 3.1.3, trimAl 1.5.1, HyPhy 2.5.101, metapredict 3.0.2 |
| `sonic_sp2` | SonicParanoid 2.0.9, DIAMOND 2.1.9, MCL 14.137, BLAST+ 2.15.0 |
| `omark` | OMArk 0.5.0, OMAmer 2.1.2 |
| `dating` | LSD2 2.2.3 |
| `tor_heat_tools` | US-align 20241201, HH-suite 3.3.0, Foldseek 10.941cd33 |
| `tgampk_plm` | the same package set as `ms210` |
| `colabfold` | ColabFold 1.6.2 (alphafold-colabfold 2.3.18, JAX 0.10.2) |
| `p68_codon` | HyPhy 2.5.101, IQ-TREE 3.1.3, MAFFT 7.526 |
| `p68_seq` | Python 3.11 with Biopython 1.88 and FreeSASA 2.2.1 |
| `p68_kvf` | Python 3.11 with pyKVFinder 0.9.5, Biopython, numpy and pandas |
| `port68_R` | R 4.5.3 with windex 2.1.0, convevol 2.2.1, phylolm 2.6.5, ape 5.8-1 |

Two cluster modules were also loaded: `miniforge/25.11.0-py3.12`, and `iqtree2/2.3.6` for resuming two
interrupted gene-tree searches.

## Programs, by step of the Materials and Methods

The headings below are the steps of the analysis. They are not one to one with the eight subsections of the
Supplementary Methods; `docs/CODE_MAP.md`, section 6, gives that mapping, and `external_runs/` is organised
by it.

### Proteomes and the species tree

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| BUSCO | 6.1.0, lineage eukaryota_odb12.2 | env | `busco -m proteins -l eukaryota_odb12.2 --offline`, one proteome per task |
| OMAmer, OMArk | 2.1.2, 0.5.0 (database LUCA.h5) | env | `omamer search --db LUCA.h5`, then `omark -f <search> -d LUCA.h5`, with `-t <taxon id>` where one was known |
| DIAMOND | 2.2.5 | env | species-tree markers: `diamond blastp --very-sensitive --evalue 1e-5 --max-target-seqs 5` |
| MAFFT | 7.526 | env, banner | species-tree markers: `mafft --auto` |
| trimAl | 1.5.1 | env, banner | `trimal -automated1` |
| IQ-TREE | 3.1.3 | env, banner, log | 210-taxon tree: `-p <partitions> -m LG+G4 -B 1000 -seed 20260817`. 209-taxon PMSF tree: `-m LG+C20+F+G4 -ft <210-taxon tree> -B 1000 -seed 20260817`. Concordance factors: `--gcf <20 marker trees> --scf 100` |
| R (system installation) with ape `chronos` | R 4.5.1 with the home-folder library, **not** the R 4.5.3 of the cluster steps | log | penalised-likelihood chronograms under the correlated, relaxed, discrete and clock models, with the calibrations of the Methods. The script of record is lost and the run of record reports `conv=FALSE`; `docs/DATING_GAP.md` states what is missing |
| LSD2 | 2.2.3 (conda package `lsd` 2.2.3, build `hdfd78af_3`, in the `dating` environment) | env | alternative dating: `lsd2 -i <rooted tree> -d <dates> -s 7933 -z 0`. That environment holds no R, which is consistent with the `chronos` dating having run under the system R |
| treePL | 1.0, built from source on 22 August 2026 from commit `1f5c99f088227c7c8dacdf343e5d1ab23405009b` | banner, log, source clone | alternative dating: prime, then cross-validation of the smoothing value, then a thorough run, on root A2 under both calibration sets. The build log records the banner `treePL version 1.0`; the clone records the commit |

### Candidate discovery

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| OrthoFinder | 3.1.5 | env | `orthofinder -S diamond --only-groups -I 1.2`; the saved similarity graph was re-clustered at further MCL inflation values with `-b <working folder> --only-groups -I <value>` |
| SonicParanoid | 2.0.9 | env (pinned at creation), script | Delivered clustering: `sonicparanoid -m default -ca -op -ot`, at the program's own defaults, arch-merging coverage 0.75 and MCL inflation 1.5. Sweeps: `--min-arch-merging-cov 0.90` (the run submitted by `sp2_cov90.sbatch`) and `-I` at 2.0 and 3.0 |
| KinOrtho | commit `5c14ed1d0b49b6e98d16b56522e0b886d48b11e1`, 26 July 2021; the project publishes no release number | source clone | AMPKα only: `kinortho.py -f <full-length queries> -d <kinase-domain queries> -t 32`, in three variants (the published PF00069 query at E 1e-5, the query recut to the declared AMPKα controls, and the BLAST default E 10) |
| MCL | 22.282 (with OrthoFinder), 14.137 (with SonicParanoid) | env | called by OrthoFinder, SonicParanoid and KinOrtho |
| HMMER | 3.4 (Aug 2023) | log, env, probe, script | genome-wide search of all 210 proteomes against the profile panel: `hmmsearch --tblout <tbl> --domtblout <domtbl> -E 1e5 --domE 1e5 --max --cpu 8 <panel_v2.hmm> <proteome.faa>`. Every per-species output carries its own `# Version: 3.4 (Aug 2023)` header and the complete option-settings line, so the version and the settings are recorded per run, not inferred. The ceilings and `--max` also appear at the call site in `rl_threshold_pilot_pertarget.py` (line 180). `hmmbuild` for the contaminant-family profiles. `phmmer` in the genome-grade re-examination |
| BLAST+ (`rpsblast`, `makeprofiledb`) | 2.17.0 | env, probe | `makeprofiledb` over every profile of CDD 3.21, then `rpsblast -evalue 0.01 -outfmt '6 qseqid sseqid bitscore evalue qstart qend'` |
| InterProScan | 5.78-109.0 | probe, script | `interproscan.sh -i <input.faa> -f TSV -appl Pfam,PANTHER,SMART`. The version is also pinned in the job scripts themselves, which set the install directory to `interproscan-5.78-109.0`. The installation no longer exists; its outputs remain |
| FastTree | 2.2.0 | env | the MAFFT and FastTree tree method of the earlier candidate build |

### Orthology calls

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| MAFFT | 7.526 | env, banner | gene trees: `mafft --globalpair --maxiterate 1000 --anysymbol` (G-INS-i, as van Dam et al. 2011), and `mafft --auto` for a tip set above the G-INS-i size limit |
| IQ-TREE | 3.1.3 | env, banner, log | gene trees: `-m MFP -B 1000 --alrt 1000 --seed 488761235`. Site-heterogeneous cross-check: `-m LG+C60+F+G -ft <ModelFinder tree> -B 1000 --alrt 1000 --seed 488761235` |
| IQ-TREE | 2.3.6 (cluster module) | job script | resumed the interrupted TOR and RAPTOR gene-tree searches from their checkpoints |
| RAxML | 8.2.13, conda package `raxml-8.2.13-h7b50bb2_3` from bioconda. Its own banner prints "RAxML version 8.2.12 released by Alexandros Stamatakis on May 2018", because upstream did not change the version string in that release | env, banner, log | `raxmlHPC-PTHREADS -f a -x 488761235 -p 488761235 -N autoMRE -m PROTGAMMAIWAG`, from the `orthofinder3` environment. The two numbers describe one binary: the package metadata of the environment that owns `raxmlHPC-PTHREADS` says 8.2.13, and that same file prints 8.2.12 |

### Sequence phenotypes and structure

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| MAFFT | 7.526 (`orthofinder3` environment) | env, script | AMPKγ contact positions, AMPKα tail and AMPKβ module: `mafft --globalpair --maxiterate 1000 --anysymbol --quiet` (built in `rbx1_score.py`, line 57). Coccidian residue-correspondence test: `mafft --localpair --maxiterate 1000` |
| metapredict | 3.0.2 | pinned, script | `predict_disorder` called with no model argument, so the release default network applied: **V3**, which `metapredict/parameters.py` sets as `DEFAULT_NETWORK` in 3.0.2. Called from `rbx1_alpha.py`, which refuses to run on any other version |
| ColabFold | 1.6.2 | env, project log | every AMPKγ copy, predicted alone: `colabfold_batch --model-type alphafold2_ptm --num-models 5 --num-seeds 5 --num-recycle 12 --random-seed 0`, with multiple sequence alignments from `colabfold_batch --msa-only`. The earlier 86-query panel: the same settings. The three-chain complexes and the partner-context benchmark: the same settings with `--model-type alphafold2_multimer_v3`. No template flag and no relaxation flag was passed, so the defaults applied: no templates, no Amber relaxation. The analysis scripts read only `*_unrelaxed_rank_*_seed_*.pdb`, which is the artifact for the second of those |
| pyKVFinder | 0.9.5 | env | cavity volumes at the AMP sites, through the measurement library `p45b_colabfold_analysis.py` |
| rate4site | 3.0.0, compiled from source | banner | empirical Bayes `-ib` with the assertion build; maximum likelihood `-im -k 16` with the build without assertions; reference sequence *Homo sapiens* |

### TOR N-terminal HEAT repeats

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| US-align | 20241201 | env | `USalign <model crop> <reference> -mol prot` (in `run_direct_alignment.py`) |
| HH-suite (`hhsearch`, HHrepID) | 3.3.0 | env | `hhsearch -M first -p 0 -E 100 -Z 100 -B 100 -e 100` against SCOP40 (`scop40_01Mar17`); HHrepID with its own calibration files; PDB70 (`pdb70_from_mmcif_220313`) |
| STRPsearch, with Foldseek | Foldseek 10.941cd33; STRPsearch commit `ed325a7f6b77578ee96753f31dbdcfda4931cca8` (BioComputingUP/STRPsearch, branch `master`); the project publishes no release number | env, source clone | `strpsearch.py query-file <structure> <out dir> --no-chainsaw --max-eval 0.1 --min-height auto --temp-dir <tmp>`, recorded per structure in that run's own receipt file |

### Phylogenetic regressions, reconstructions and other tests

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| R | 4.5.3 | sessionInfo | all cluster R steps; every record lists ape 5.8-1 and phylolm 2.6.5 |
| phylolm | 2.6.5 | sessionInfo, env | `phyloglm(method = "logistic_MPLE", btol = 30)` with the multistart search of the Methods; `phylolm(model = "lambda")`, and `model = "BM"` as a sensitivity |
| corHMM | 2.8 | env | `corHMM(rate.cat = 1, node.states = "marginal")` under the ER, SYM and ARD models |
| phangorn | 2.12.1 | env | `ancestral.pars(type = "MPR")`, the parsimony cross-check |
| phytools | 2.5-2 | env | `fitPagel` and `fitMk` in the tests of correlated evolution and of loss direction |
| windex, convevol | 2.1.0, 2.2.1 | env | Wheatsheaf index and Ct1 in the convergence tests |
| HyPhy (GARD) | 2.5.101 | banner, env | amino acids: `gard --type amino-acid --model LG --rv GDD --rate-classes 4`; codons: `gard --type codon --code Universal --rv GDD --rate-classes 4` |
| mvMORPH | 1.2.2 | sessionInfo (workstation) | the joint three-response model of Supplementary Table S38, run on the authors' workstation under R 4.6.1. The manuscript does not describe this model and rests no claim on it. The record is the workstation's own `sessionInfo()`; `envs/R_SESSIONS.tsv` holds the cluster records only, all of them R 4.5.3, and carries no row for this run |
| Python | 3.12 (`ms210_phase2`), 3.11 (`ms210`, `p68_*`) | env | the Python steps; standard library plus pandas, numpy and Biopython 1.88 |

## Online resources queried

- NCBI Taxonomy, through Entrez EFetch, on 19 September 2026.
- UniProt and the AlphaFold Protein Structure Database, for the exact-sequence TOR models.
- The InterPro API, for the curated Pfam and PANTHER profiles.
- The RCSB Protein Data Bank, for the experimental control structures.

## Notes on particular versions

Every program on this page has a stated version resting on at least one record. These are the entries whose
record needs a word of explanation.

- **RAxML, 8.2.13 and 8.2.12 are the same binary.** The `orthofinder3` environment's package metadata shows
  that `raxml-8.2.13-h7b50bb2_3` owns the `raxmlHPC-PTHREADS` that ran, and that file's own banner prints
  "RAxML version 8.2.12 released by Alexandros Stamatakis on May 2018". Upstream did not change the version
  string for that release, so the run logs and the package list disagree by one patch number while describing
  one file. 8.2.13 is the package; 8.2.12 is what the program calls itself.
- **Three tools publish no release number**, so their version is the commit that was built or run: KinOrtho,
  treePL (which also prints `treePL version 1.0`) and STRPsearch. Each commit is recorded in its row above,
  taken from the clone the analysis used.
- **metapredict's network is V3.** The call passes no model argument, and metapredict 3.0.2 sets
  `DEFAULT_NETWORK = 'V3'` in `metapredict/parameters.py`, so the release default applied.
- **HyPhy 2.5.101** rests on the job script's banner and the environment record; the GARD log itself prints
  only the analysis version, which is a property of GARD rather than of HyPhy.
- **InterProScan's installation no longer exists.** Its outputs remain, and the version is pinned in two
  places that do survive: the software check and the install directory named in the job scripts.

## The clade assignment

Clades are read from the NCBI Taxonomy lineage of each species, retrieved through Entrez EFetch on
19 September 2026, which gives 31 clades over the 210 species. Three species that NCBI Taxonomy nests within
Apicomplexa are set aside as Squirmida on published phylogenomic evidence, which takes Apicomplexa from 41
species to 38 and the roster from 31 clades to 32; every analysis that groups by clade is reported both ways.
The per-species assignment is published in full in the supplementary table of species states, with a label
for all 210 species and no blanks, and the working records keep the pre-curation and post-curation
assignments side by side so that the three changes stay visible.

## The one thing the records do not fix

**The dating script of record.** It was edited after the run that produced the chronogram, and the surviving
copy's md5 differs from the one the job log stamps. Everything else about that run is recorded: the rooting,
the calibration set with the five calibrations that entered and the three that were dropped, the model and
smoothing sweep, the convergence flag, the selection rule and the identity of the resulting file.
`docs/DATING_GAP.md` states it in full and the Supplementary Methods cites that file.
