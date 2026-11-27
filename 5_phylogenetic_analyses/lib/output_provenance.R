# >>> deposit header
# Purpose: Shared R helper that returns the header lines written at the top of each output of the
#          candidate-replacement regressions and the ancestral-state reconstructions: script and its MD5
#          checksum, seed, R and package versions, the MD5 checksum of each input, and the time.
# Methods: several steps (shared helper)
# Inputs:  the calling script's path, a named vector or list of its input paths, and the seed
# Outputs: none written; returns the `##` lines
# Tools:   R (tools::md5sum); reports the versions of ape, corHMM, phangorn and phytools it finds; versions: see envs/
# <<< deposit header
# Ran as: rl_rm_prov.R
## Shared header lines for every output under `14_remodeling`.
## Why: an output that cannot name the version of its inputs cannot be traced to them (the trait tables
## changed three times on 9 September 2026). So every output names the script that wrote it and that
## script's MD5 checksum, the MD5 checksum of every input it read, the seed, the software versions and
## the time. Reproducing a number then means re-running one command against named inputs, not
## reconstructing what was on disk at the time.
prov_lines <- function(script, inputs, seed) {
  md5 <- function(p) tryCatch(tools::md5sum(path.expand(p))[[1]],
                              error = function(e) "UNREADABLE")
  L <- c(sprintf("## PRODUCED BY %s  md5=%s", basename(script), md5(script)),
         sprintf("## SEED %s", as.character(seed)),
         sprintf("## R %s; %s", getRversion(),
                 paste(sapply(c("ape", "corHMM", "phangorn", "phytools"), function(p)
                   tryCatch(sprintf("%s %s", p, as.character(packageVersion(p))),
                            error = function(e) sprintf("%s absent", p))),
                   collapse = "; ")))
  for (nm in names(inputs)) L <- c(L, sprintf("## INPUT %s = %s  md5=%s",
                                              nm, inputs[[nm]], md5(inputs[[nm]])))
  c(L, sprintf("## GENERATED %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
}
