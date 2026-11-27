# The species-tree dating: what is recorded, and the one thing that is not

The chronogram of record is identified exactly, with its inputs, its settings and the rule that selected it.
What is missing is the text of the script that produced it. This page states both. Nothing here is
reconstructed. The Supplementary Methods cites this file.

## The chronogram of record

| Item | Value |
|---|---|
| File | `chronos/A2_setA.correlated_l0.1.tre`, published as `species_chronogram_209.treefile` |
| md5 | `301836d4f15956d8b14ffbfe420a64a9` (the two files are byte-identical) |
| Rooting | A2 |
| Calibration set | A, the secondary calibrations |
| Model | correlated rates, smoothing parameter 0.1 |
| How 0.1 was chosen | leave-one-calibration-out cross-validation on the correlated model, root A2, set A |
| Reference for comparison | `chronos/A2_setA.correlated_l1.tre` |
| Chronograms compared | 86 |
| Root age | 2,299.0 Ma |
| Convergence | **`conv=FALSE`** |

The run swept three rooting positions by two calibration sets, and within each the correlated, relaxed,
discrete and clock models over smoothing values 0, 0.1, 1, 10, 100 and 1,000. LSD2 and treePL were run beside
it on the same rootings and sets. The selection among all of them is recorded, not reconstructed.

## The calibrations: eight declared, five entered, three dropped

Calibration set A declares eight calibrations, every one a uniform minimum and maximum interval, frozen before
any dating run. **Five entered the fit of record and three were dropped.**

| Calibration | Entered | Node | Interval (Ma) |
|---|---|---|---|
| LECA_root | yes | the root, 209 tips | 1,958 to 2,386 |
| Cryptophyta | yes | 2 tips | 440 to 1,658 |
| Ochrophyta | yes | 5 tips | 622 to 1,298 |
| Haptophyta | yes | 2 tips | 579 to 1,943 |
| Myzozoa | yes | 60 tips | 696 to 1,520 |
| Archaeplastida | **dropped** | MRCA subtends 31 tips, 29 expected; extra *Cryptomonas paramecium* and *Guillardia theta* | 1,807 to 2,137 |
| Rhodophyta_crown | **dropped** | MRCA subtends 4 tips, 5 expected; missing *Galdieria sulphuraria* | 1,732 to 1,984 |
| Fungi_crown | **dropped** | MRCA subtends 19 tips, 22 expected; missing *Paraphelidium tribonemae*, *Paramicrosporidium saccamoebae* and *Rozella allomycis* | 1,078 to 1,759 |

A calibration whose MRCA tip set is not exactly the expected lineage is dropped and listed; it is never
applied to a deeper node. That rule is stated in the frozen calibration file itself.

**Five, not six.** An earlier version of this page, and of `docs/software.md`, said that six of the eight
entered and two were dropped. That is the count for **treePL**, a different method in the same run, which
dropped Archaeplastida and Fungi_crown only. The `chronos` fit of record dropped three.

## What is missing

The script that produced the chronogram was edited after the run.

| | md5 |
|---|---|
| Stamped in the job log of the run | `2171c20ab00817e67bf598fffdf2d48a` |
| The surviving copy | `6b815479f610ccfe1bca2c126126faf9` |

The two differ, so the file that still exists is not the file that ran. The surviving copy is not deposited
here in its place, because it is not the code that produced the result, and the script of record is not
reconstructed.

This is the only part of the dating that cannot be reproduced from a record. The inputs, the rooting, the
calibration set and its drops, the model, the smoothing sweep, the selection rule and the resulting file are
all recorded.

## Supporting records

| What | Where |
|---|---|
| The stamped md5, with the other scripts of the same job | the job log of the dating run |
| The calibration file, frozen before any run, with every bound, its source and the quoted passage | `calibrations.tsv` beside the chronogram |
| The applied and dropped calibrations, with the tip-count reason for each drop | the chronos run log |
| The full model and smoothing sweep, with log likelihood, root age and convergence flag per fit | the same log |
| The comparison of all 86 chronograms and the selection of the record | the chronogram report |

Those records are on the cluster that ran the analysis, not in this repository. The data deposit holds the
chronogram itself.

## The alternative datings

Both ran on the same rootings and calibration sets, and both survive as command lines in
`1_proteomes_taxonomy_species_tree/external_runs/`.

- **LSD2** 2.2.3: six runs, three rootings by two sets, all returning successfully.
- **treePL** 1.0: prime, then cross-validation of the smoothing value, then a thorough run, on root A2 under
  both sets. On set A it dropped two calibrations, Archaeplastida and Fungi_crown; on set B it dropped
  Bangio_Florideo and Eumetazoa. Its cross-validation chose a smoothing of 1e-05 on set A.

Neither was selected as the chronogram of record: against the reference both fall well short of the
agreement benchmark, with coefficients of determination near 0.1 to 0.3 against a benchmark of 0.94.
