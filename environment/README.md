# The software environments

Fourteen conda environments ran the analyses, on Red Hat Enterprise Linux 9. The files here are the record:

| File | What it is |
|---|---|
| `ms210_phase2.yml` | An installable specification of the environment that runs most of the deposited Python and R, generated from the package list below. The one file to start from |
| `packages_after_analysis.tsv` | Every conda, Python and R package of every environment, read on 25 September 2026 after the last analysis: environment, kind, name, version, build, channel |
| `packages_during_analysis.tsv` | The `conda list` and `pip freeze` exports of the nine environments that existed on 13 September 2026, while the analyses were running. Not a superseded copy of the file above: a different format, and the contemporaneous state |
| `conda_history.txt` | Each environment's own log of its `conda create` and `conda install` commands, with dates |
| `r_sessions.tsv` | A summary of the 494 `sessionInfo()` records the R steps wrote beside their outputs, which tie each R output to the package versions in force when it was written |

## Recreating an environment

Nine of the fourteen creation commands survive in `conda_history.txt` and are repeated here with the
cluster prefix replaced. They install what was asked for, not what was resolved: the versions that actually
ran are in `packages_after_analysis.tsv`, and for an exact rebuild that file is the authority.

### `busco6`

```
conda create -y -p <prefix>/busco6 -c conda-forge -c bioconda busco
```

### `colabfold`

```
conda create -y -p <prefix>/colabfold python=3.11 pip
```

### `dating`

```
conda create -y -p <prefix>/dating -c conda-forge -c bioconda lsd nlopt
```

### `ms210`

```
conda create -y -p <prefix>/ms210 -c conda-forge -c bioconda iqtree>=2.2 trimal pal2nal hyphy>=2.5 python=3.11
```

### `ms210_phase2`

```
conda create -y -p <prefix>/ms210_phase2 -c conda-forge -c bioconda python=3.12 numpy pandas scipy statsmodels scikit-learn networkx matplotlib seaborn biopython dendropy pyarrow umap-learn r-base=4.5 r-ape r-phytools r-geiger r-caper r-castor r-phylolm r-nlme r-mass r-ggplot2 r-data.table r-igraph r-segmented r-circlize r-nloptr r-expm r-phangorn r-gensa r-numderiv r-corpcor r-nnet r-viridis r-svglite bioconductor-ggtree bioconductor-complexheatmap
```

### `omark`

No creation command survives. Its packages and versions are the rows of
`packages_after_analysis.tsv` where the first column is `omark`.

### `orthofinder3`

```
conda create -y -p <prefix>/orthofinder3 -c conda-forge -c bioconda python=3.12 orthofinder>=3.1 mcl
```

### `p68_codon`

No creation command survives. Its packages and versions are the rows of
`packages_after_analysis.tsv` where the first column is `p68_codon`.

### `p68_kvf`

No creation command survives. Its packages and versions are the rows of
`packages_after_analysis.tsv` where the first column is `p68_kvf`.

### `p68_seq`

No creation command survives. Its packages and versions are the rows of
`packages_after_analysis.tsv` where the first column is `p68_seq`.

### `port68_R`

No creation command survives. Its packages and versions are the rows of
`packages_after_analysis.tsv` where the first column is `port68_R`.

### `sonic_sp2`

```
conda create -y -p <prefix>/sonic_sp2 -c conda-forge -c bioconda python=3.11 sonicparanoid=2.0.9 mmseqs2=13.45111 diamond=2.1.9 blast=2.15.0 mcl=14.137 seqkit
```

### `tgampk_plm`

```
conda create -y -p <prefix>/tgampk_plm -c conda-forge -c bioconda iqtree>=2.2 trimal pal2nal hyphy>=2.5 python=3.11
```

### `tor_heat_tools`

```
conda create -y -p <prefix>/tor_heat_tools -c conda-forge -c bioconda python=3.11 hhsuite foldseek usalign biopython numpy scipy pandas matplotlib seaborn typer rich requests pip
```
