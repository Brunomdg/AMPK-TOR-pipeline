# Software and versions

This page lists every program and package the analyses used, with its version and the evidence for that
version. Where a program was run by a job script that is not part of this repository, the settings given
are copied from that job script's command line.

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
- **not recorded**: no record was found.

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

### Proteomes and the species tree

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| BUSCO | 6.1.0, lineage eukaryota_odb12.2 | env | `busco -m proteins -l eukaryota_odb12.2 --offline`, one proteome per task |
| OMAmer, OMArk | 2.1.2, 0.5.0 (database LUCA.h5) | env | `omamer search --db LUCA.h5`, then `omark -f <search> -d LUCA.h5`, with `-t <taxon id>` where one was known |
| DIAMOND | 2.2.5 | env | species-tree markers: `diamond blastp --very-sensitive --evalue 1e-5 --max-target-seqs 5` |
| MAFFT | 7.526 | env, banner | species-tree markers: `mafft --auto` |
| trimAl | 1.5.1 | env, banner | `trimal -automated1` |
| IQ-TREE | 3.1.3 | env, banner, log | 210-taxon tree: `-p <partitions> -m LG+G4 -B 1000 -seed 20260817`. 209-taxon PMSF tree: `-m LG+C20+F+G4 -ft <210-taxon tree> -B 1000 -seed 20260817`. Concordance factors: `--gcf <20 marker trees> --scf 100` |
| R (system installation) with ape `chronos` | R 4.5.1 with the home-folder library | log | penalised-likelihood chronograms under the correlated, relaxed, discrete and clock models, with the calibrations of the Methods |
| LSD2 | 2.2.3 | env | alternative dating: `lsd2 -i <rooted tree> -d <dates> -s 7933 -z 0` |
| treePL | not recorded | not recorded | alternative dating: prime, then cross-validation of the smoothing value from 1e6 to 1e-12 (step 0.1), then a thorough run, 7,933 sites |

### Candidate discovery

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| OrthoFinder | 3.1.5 | env | `orthofinder -S diamond --only-groups -I 1.2`; the saved similarity graph was re-clustered at further MCL inflation values with `-b <working folder> --only-groups -I <value>` |
| SonicParanoid | 2.0.9 | env (pinned at creation) | `sonicparanoid -m default -ca -op -ot --min-arch-merging-cov 0.90`, with MCL inflation varied by `-I` |
| KinOrtho | not recorded | not recorded | AMPKα only: `kinortho.py -f <full-length queries> -d <kinase-domain queries> -t 32` |
| MCL | 22.282 (with OrthoFinder), 14.137 (with SonicParanoid) | env | called by OrthoFinder, SonicParanoid and KinOrtho |
| HMMER | 3.4 | env, probe | genome-wide search of all 210 proteomes against the profile panel: `hmmsearch -E 1e5 --domE 1e5 --max --tblout --domtblout`. `hmmbuild` for the contaminant-family profiles. `phmmer` in the genome-grade re-examination |
| BLAST+ (`rpsblast`, `makeprofiledb`) | 2.17.0 | env, probe | `makeprofiledb` over every profile of CDD 3.21, then `rpsblast -evalue 0.01 -outfmt '6 qseqid sseqid bitscore evalue qstart qend'` |
| InterProScan | 5.78-109.0 | probe | `interproscan.sh -f TSV -appl Pfam,PANTHER,SMART`. The installation no longer exists; its outputs remain |
| FastTree | 2.2.0 | env | the MAFFT and FastTree tree method of the earlier candidate build |

### Orthology calls

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| MAFFT | 7.526 | env, banner | gene trees: `mafft --globalpair --maxiterate 1000 --anysymbol` (G-INS-i, as van Dam et al. 2011), and `mafft --auto` for a tip set above the G-INS-i size limit |
| IQ-TREE | 3.1.3 | env, banner, log | gene trees: `-m MFP -B 1000 --alrt 1000 --seed 488761235`. Site-heterogeneous cross-check: `-m LG+C60+F+G -ft <ModelFinder tree> -B 1000 --alrt 1000 --seed 488761235` |
| IQ-TREE | 2.3.6 (cluster module) | job script | resumed the interrupted TOR and RAPTOR gene-tree searches from their checkpoints |
| RAxML | 8.2.12 per its own logs (8.2.13 is the installed conda package) | log | `raxmlHPC-PTHREADS -f a -x 488761235 -p 488761235 -N autoMRE -m PROTGAMMAIWAG` |

### Sequence phenotypes and structure

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| MAFFT | 7.526 (`orthofinder3` environment) | env | AMPKγ contact positions, AMPKα tail and AMPKβ module: `mafft --globalpair --maxiterate 1000 --anysymbol` (in `rbx1_score.py`). Coccidian residue-correspondence test: `mafft --localpair --maxiterate 1000` |
| metapredict | 3.0.2 | pinned | default predictor, `predict_disorder`; the model identifier it resolved to was not recorded |
| ColabFold | 1.6.2 | env, project log | every AMPKγ copy, predicted alone: `colabfold_batch --model-type alphafold2_ptm --num-models 5 --num-seeds 5 --num-recycle 12 --random-seed 0`, with multiple sequence alignments from `colabfold_batch --msa-only`. The earlier 86-query panel: the same settings. The three-chain complexes and the partner-context benchmark: the same settings with `--model-type alphafold2_multimer_v3` |
| pyKVFinder | 0.9.5 | env | cavity volumes at the AMP sites, through the measurement library `p45b_colabfold_analysis.py` |
| rate4site | 3.0.0, compiled from source | banner | empirical Bayes `-ib` with the assertion build; maximum likelihood `-im -k 16` with the build without assertions; reference sequence *Homo sapiens* |

### TOR N-terminal HEAT repeats

| Program | Version | Evidence | Settings as run |
|---|---|---|---|
| US-align | 20241201 | env | `USalign <model crop> <reference> -mol prot` (in `run_direct_alignment.py`) |
| HH-suite (`hhsearch`, HHrepID) | 3.3.0 | env | `hhsearch -M first -p 0 -E 100 -Z 100 -B 100 -e 100` against SCOP40 (`scop40_01Mar17`); HHrepID with its own calibration files; PDB70 (`pdb70_from_mmcif_220313`) |
| STRPsearch, with Foldseek | Foldseek 10.941cd33; STRPsearch not recorded | env | structured tandem-repeat search on each structure, maximum E-value 0.1 |

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
| mvMORPH | 1.2.2 | sessionInfo | the joint three-response model of Supplementary Table S26, run on the authors' workstation under R 4.6.1 |
| Python | 3.12 (`ms210_phase2`), 3.11 (`ms210`, `p68_*`) | env | the Python steps; standard library plus pandas, numpy and Biopython 1.88 |

## Online resources queried

- NCBI Taxonomy, through Entrez EFetch, on 19 September 2026.
- UniProt and the AlphaFold Protein Structure Database, for the exact-sequence TOR models.
- The InterPro API, for the curated Pfam and PANTHER profiles.
- The RCSB Protein Data Bank, for the experimental control structures.

## What the records do not show

- The versions of KinOrtho, treePL and STRPsearch, and the metapredict model identifier.
- The criteria that delimit the remaining clades.
- **Versions resting on a single source:**
  - InterProScan, HMMER, rpsblast and MAFFT versions come from one software check, not from each run's log;
    the environment records agree for HMMER, BLAST+ and MAFFT.
  - HyPhy 2.5.101 comes from the job script's banner and the environment; the GARD log prints only the
    analysis version.
- **RAxML:** its logs report 8.2.12, while the installed package is 8.2.13.
- **Species-tree dating:** the script was edited after the run that produced the chronogram. The job log stamps
  its md5 as `2171c20ab00817e67bf598fffdf2d48a`, and the surviving copy differs. The chronos log of the
  chronogram of record reports `conv=FALSE`, and two of the eight calibrations were dropped from that fit
  because they are not clades of this tree.
