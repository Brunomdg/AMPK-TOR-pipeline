# The species-tree dating: the one reproduction gap

The chronogram of record cannot be reproduced from this repository. This page states exactly what is missing
and what the records do show. Nothing here is reconstructed. The Supplementary Methods cites this file.

## The script of record is lost

The job log of the run that produced the chronogram stamps the script's md5 as:

```text
2171c20ab00817e67bf598fffdf2d48a
```

The script was edited after that run. **The surviving copy's md5 differs from the stamped one**, which is how
the loss was found: the file that still exists is not the file that ran.

The surviving copy is not in this repository, and neither `docs/SOFTWARE.md` nor any other record here carries
its md5, so that value is not stated on this page. The script of record is not reconstructed, and the surviving
copy is not deposited in its place, because it is not the code that produced the result.

## The run of record did not converge

The `chronos` log of the chronogram of record reports:

```text
conv=FALSE
```

## Two of eight calibrations did not enter the fit

Eight calibrations were declared. **Six entered the fit.** Two were dropped because the nodes they constrain
are not clades of this topology, so there was no node for the calibration to apply to.

## What the records do show

| Item | Value | Evidence |
|---|---|---|
| Program | `chronos`, from ape | the job log |
| R | 4.5.1, with the home-folder library, not the R 4.5.3 of the cluster steps | the job log |
| Models fitted | correlated, relaxed, discrete and clock | the job log |
| Calibrations | those stated in the Methods; six of the eight declared entered the fit | the job log |
| Command line | none survives | - |

The two alternative datings do survive as command lines, in
`external_runs/sm1_proteomes_taxonomy_species_tree/`: LSD2 2.2.3 in `lsd2_dating.sh`, and treePL in
`treepl_dating.sh`, whose version is not recorded either.

## What this means for a reader

The dated species tree is an input to the ancestral-state reconstructions and to the phylogenetic regressions,
which are deposited in full and run on the tree as given. The tree itself is in the data deposit. What cannot
be redone from this repository is the dating step that produced it.
