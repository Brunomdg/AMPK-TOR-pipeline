# >>> deposit header
# Purpose: For one ordered pair of the ten targets (chosen by `SLURM_ARRAY_TASK_ID`, 0 to 89), fits the phylogenetic
#          logistic regression of response absence on predictor absence under both codings of the matrix, with and
#          without proteome size.
# Methods: Phylogenetic regressions
# Inputs:  `Q0_REPLACEMENT_per_species.tsv` (primary coding) and `Q0_REPLACEMENT_per_species__R37coding.tsv`
#          (sensitivity coding); the dated species tree (209 tips); `proteome_index.tsv` (sequence counts);
#          `MULTI_RECORD_SPECIES.tsv` (proteomes with more than one sequence per gene)
# Outputs: per pair, `MODELS.tsv` (four fits), one `model_NN.rds` per fit, `INPUT_MD5.tsv` (MD5 of every input,
#          including this script and its helper) and `DONE.txt`
# Tools:   R; phylolm 2.6.5 (through `lib/multistart_search.R`); ape (version: see envs/)
# <<< deposit header
## deposit path defaults: begin
.TA210_ENVS <- local({v <- Sys.getenv("TA210_ENVS"); if (nzchar(v)) v else "/gpfs2/scratch/bdigenov/envs"})
.TA210_SCRATCH <- local({v <- Sys.getenv("TA210_SCRATCH"); if (nzchar(v)) v else "/gpfs2/scratch/bdigenov"})
## deposit path defaults: end
## deposit repo root: begin
.REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(f[1]) else "."
  normalizePath(file.path(d, ".."))
})
## deposit repo root: end
## Associations of absence among all ten targets, as 90 directed pairs; this family was added on 18 September
## 2026 and is separate from the nine-test family of component-to-subunit models.
## Predictor and response are both coded absence = 1, presence = 0, so a positive coefficient means that the two
## absences go together.
suppressPackageStartupMessages({library(ape);library(phylolm)})
S <- paste0(.TA210_SCRATCH, '/ms210_phase2/TOR-AMPK_orthologs_10target_R1')
P <- file.path(S,'20_port_68taxa_analyses')
Q <- file.path(S,'14_remodeling/Q0_replacement')
targets <- c('TOR','RAPTOR','RICTOR','LST8','SIN1','AMPKa','AMPKb','AMPKg','LKB1','CAMKK')
pairs <- expand.grid(PREDICTOR=targets,RESPONSE=targets,stringsAsFactors=FALSE)
pairs <- pairs[pairs$PREDICTOR!=pairs$RESPONSE,];rownames(pairs)<-NULL
stopifnot(nrow(pairs)==90,!anyDuplicated(pairs))
id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
stopifnot(!is.na(id),id>=0,id<90)
set.seed(488761235)
D <- file.path(P,'02_derived/network10_loss_20260918',sprintf('pair_%02d',id))
dir.create(D,recursive=TRUE,showWarnings=FALSE)
source(file.path(.REPO_ROOT, '5_phylogenetic_analyses', 'lib', 'multistart_search.R'))
rd <- function(p)read.delim(p,comment.char='#',quote='',stringsAsFactors=FALSE,check.names=FALSE)
treepath <- paste0(.TA210_SCRATCH, '/concord/ms210/P04b_chronogram/species_chronogram_209.treefile')
indexpath <- path.expand('~/ref/proteomes/proteome_index.tsv')
multipath <- file.path(S,'05_method_and_rules/MULTI_RECORD_SPECIES.tsv')
files <- c(TRAITS=file.path(Q,'Q0_REPLACEMENT_per_species.tsv'),R37=file.path(Q,'Q0_REPLACEMENT_per_species__R37coding.tsv'))
inputpaths <- c(files,treepath,indexpath,multipath,file.path(.REPO_ROOT,'5_phylogenetic_analyses/ordered_pair_models.R'),file.path(.REPO_ROOT,'5_phylogenetic_analyses/lib/multistart_search.R'))
hashes <- tools::md5sum(inputpaths)
write.table(data.frame(FILE=names(hashes),MD5=unname(hashes)),file.path(D,'INPUT_MD5.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
tree <- read.tree(treepath)
idx <- read.delim(indexpath,quote='',comment.char='',stringsAsFactors=FALSE)
idx$SPECIES <- sub('\\.faa$','',idx$file)
stopifnot(nrow(idx)==210,!anyDuplicated(idx$SPECIES))
nrec <- setNames(as.numeric(idx$records),idx$SPECIES)
multi <- rd(multipath);multi <- multi$SPECIES[multi$VERDICT!='GENE_COLLAPSED']
pr <- pairs[id+1,];rows <- list()
for(coding in names(files))for(withcov in c(FALSE,TRUE)){
 raw <- rd(files[[coding]])
 stopifnot(nrow(raw)==210,!anyDuplicated(raw$SPECIES),all(paste0(targets,'_absent_p1') %in% names(raw)))
 x <- raw[[paste0(pr$PREDICTOR,'_absent_p1')]];y <- raw[[paste0(pr$RESPONSE,'_absent_p1')]]
 stopifnot(all(x %in% c(0,.5,1)),all(y %in% c(0,.5,1)))
 known <- x!=.5 & y!=.5;ontree <- raw$SPECIES %in% tree$tip.label
 covok <- !(raw$SPECIES %in% multi) & is.finite(nrec[raw$SPECIES]) & nrec[raw$SPECIES]>0
 ok <- known & ontree & if(withcov)covok else TRUE
 d <- data.frame(SPECIES=raw$SPECIES[ok],predictor_absent=as.integer(x[ok]),response_absent=as.integer(y[ok]))
 if(withcov)d$log10_n_proteins <- log10(nrec[d$SPECIES])
 t2 <- keep.tip(tree,d$SPECIES);d <- d[match(t2$tip.label,d$SPECIES),];rownames(d)<-d$SPECIES
 r <- data.frame(PREDICTOR=pr$PREDICTOR,RESPONSE=pr$RESPONSE,CODING=coding,COVARIATE=if(withcov)'log10_n_proteins' else 'none',N=nrow(d),EXCLUDED_UNCERTAIN=sum(!known),EXCLUDED_OFF_TREE=sum(known & !ontree),EXCLUDED_COVARIATE=if(withcov)sum(known & ontree & !covok) else 0,COEF=NA_real_,SE=NA_real_,CI95_LO=NA_real_,CI95_HI=NA_real_,NOMINAL_P=NA_real_,ALPHA=NA_real_,ALPHA_WARN=NA_integer_,NONPHYLO_COEF=NA_real_,NONPHYLO_P=NA_real_,EMPTY_2X2_CELL=any(table(factor(d$predictor_absent,levels=0:1),factor(d$response_absent,levels=0:1))==0),STATUS='NOT_ESTIMABLE',RDS_FILE=NA_character_)
 if(nrow(d)>=25 && length(unique(d$predictor_absent))==2 && length(unique(d$response_absent))==2){
  f <- if(withcov)response_absent~predictor_absent+log10_n_proteins else response_absent~predictor_absent
  fit <- fit_multiple(f,d,t2);cf<-summary(fit)$coefficients['predictor_absent',]
  r$COEF<-cf[1];r$SE<-cf[2];r$CI95_LO<-cf[1]-qnorm(.975)*cf[2];r$CI95_HI<-cf[1]+qnorm(.975)*cf[2];r$NOMINAL_P<-cf[4];r$ALPHA<-fit$alpha;r$ALPHA_WARN<-fit$alphaWarn;r$STATUS<-'OK'
  r$RDS_FILE<-sprintf('model_%02d.rds',fit_counter)
  gl<-tryCatch(suppressWarnings(glm(f,data=d,family=binomial())),error=function(e)NULL)
  if(!is.null(gl)){gc<-summary(gl)$coefficients['predictor_absent',];r$NONPHYLO_COEF<-gc[1];r$NONPHYLO_P<-gc[4]}
 }
 rows[[length(rows)+1]]<-r
}
stopifnot(identical(hashes,tools::md5sum(inputpaths)))
out<-file.path(D,'MODELS.tsv');stopifnot(!file.exists(out))
write.table(do.call(rbind,rows),out,sep='\t',quote=FALSE,row.names=FALSE)
writeLines('COMPLETE: four model specifications; family adjustment deferred until all 90 pairs are validated.',file.path(D,'DONE.txt'))
