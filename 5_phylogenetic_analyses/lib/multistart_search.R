# >>> deposit header
# Purpose: The multistart search of the penalised phylogenetic logistic fits: a grid of starting values, then
#          shrinking local moves ranked by penalised log-likelihood, with the search history and every warning
#          stored beside each fit. This is the search the Data Availability statement names.
# Methods: Phylogenetic analyses
# Inputs:  the model formula, the species data and the pruned tree, passed in by the calling script, which also
#          sets the output folder `D` before sourcing this file
# Outputs: one `model_NN.rds` per fit in `D` (baseline fit, best fit, search history, attempts, warnings) and
#          `sessionInfo.txt`
# Tools:   R; phylolm 2.6.5 (phyloglm, logistic_MPLE, btol = 30; the version is checked here)
# <<< deposit header
# Ran as: core_multistart_helper.R for the core-complement fits and network10_multistart_helper.R for the
# ordered-pair fits. The two were the same code: the only differences were this header and a hardcoded output
# folder on the first line of the core copy. `D` now always comes from the caller.
if (!exists("D")) stop("REFUSED: the calling script must set D, the output folder, before sourcing this file")
dir.create(D,recursive=TRUE,showWarnings=FALSE)
stopifnot(as.character(packageVersion('phylolm'))=='2.6.5')
capture.output(sessionInfo(),file=file.path(D,'sessionInfo.txt'))
fit_counter<-0L
fit_multiple<-function(f,d,t2){
 fit_counter<<-fit_counter+1L
 dest<-file.path(D,sprintf('model_%02d.rds',fit_counter));stopifnot(!file.exists(dest))
 warnings_seen<-character()
 fit<-function(b=NULL,a=NULL){
  args<-list(formula=f,data=d,phy=t2,method='logistic_MPLE',btol=30)
  if(!is.null(b))args$start.beta<-unname(b)
  if(!is.null(a))args$start.alpha<-a
  tryCatch(withCallingHandlers(do.call(phylolm::phyloglm,args),warning=function(w){warnings_seen<<-unique(c(warnings_seen,conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e)NULL)
 }
 baseline<-fit();stopifnot(!is.null(baseline),is.finite(baseline$penlogLik))
 best<-baseline;attempts<-1L
 consider<-function(b,a){g<-fit(b,a);attempts<<-attempts+1L;if(!is.null(g)&&isTRUE(g$convergence==0)&&is.finite(g$penlogLik)&&g$penlogLik>best$penlogLik+1e-8)best<<-g}
 np<-length(coef(best));grid<-expand.grid(rep(list(c(-6,-3,-1,0,1,3,6)),np))
 for(i in seq_len(nrow(grid)))for(a in c(.001,.01,.1,1))consider(as.numeric(grid[i,]),a)
 step<-2;settled<-FALSE;history<-data.frame()
 for(round in 1:60){
  center<-coef(best);a<-best$alpha;before<-best$penlogLik
  offsets<-expand.grid(rep(list(c(-step,-step/2,0,step/2,step)),np))
  for(i in seq_len(nrow(offsets)))for(am in c(.25,.5,1,2,4))consider(center+as.numeric(offsets[i,]),a*am)
  gain<-best$penlogLik-before;history<-rbind(history,data.frame(round=round,step=step,gain=gain))
  if(gain<=1e-4){settled<-TRUE;break};step<-step*.65
 }
 saveRDS(list(formula=f,data=d,tree=t2,baseline=baseline,best=best,settled=settled,history=history,attempts=attempts,warnings=warnings_seen),dest)
 stopifnot(settled,best$convergence==0)
 cat('MODEL',fit_counter,'N',nrow(d),'GAIN',best$penlogLik-baseline$penlogLik,'ALPHA_WARN',best$alphaWarn,'\n');flush.console()
 best
}
