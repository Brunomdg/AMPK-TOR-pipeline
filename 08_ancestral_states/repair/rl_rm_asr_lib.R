#!/usr/bin/env Rscript
# >>> deposit header
# Purpose: Repaired library of ancestral-state reconstruction functions (corHMM marginal reconstruction and
#          parsimony) shared by the loss-count script and the per-edge reconstruction script.
# Methods: Ancestral-state reconstructions
# Inputs:  none of its own; the caller passes the tree, the species names and the coded tip states
# Outputs: none of its own (functions only)
# Tools:   R with ape, corHMM, phangorn and phytools attached by the caller; R 4.5.3, ape 5.8.1, corHMM 2.8,
#          phangorn 2.12.1 and phytools 2.5.2 in the ms210_phase2 environment (noted on 9 September 2026)
# <<< deposit header
## Shared functions for the count of independent losses of each component (rl_rm_asr.R) and for the per-edge
## reconstruction used by the replacement analysis and the co-transition analysis (rl_rm_q0_asr.R).
##
## Contents. The five functions that rl_rm_asr.R defined inline (tipstate, is_amb, count_origins, fit_one,
## pars_floor), copied rather than rewritten, with the two corrections marked below, plus four helpers:
## a2_outgroup() (the outgroup taken from the reference tree rooted_A2), root_tree() (the rooting step,
## returning the rooted tree and its ROOT_MODE label as a value rather than as a log line), tree_md5() (the
## checksum of the pruned tree that every per-edge row carries), and edge_states(), which returns what
## fit_one and pars_floor compute but do not return: the hardened state at both ends of every edge.
##
## Why per-edge states. The replacement analysis classifies each edge on which a core component is lost by
## what a candidate replacement did on the same edge and below it, on the same pruned tree. It therefore needs
## the per-edge states behind the loss counts, produced by the same code and not by a re-implementation.
## When the functions were first moved here, the move was checked by requiring that rl_rm_asr__refactored.R
## (rl_rm_asr.R with these definitions replaced by a source() of this file, nothing else changed) write
## ORIGIN_COUNTS.tsv with data rows identical to those of rl_rm_asr.R. The two corrections below were made
## later, when the reconstructions were repaired, so results from this library differ from the original.
##
## Ambiguity. An ambiguous tip (coded 0.5) is written "0&1": a flat prior, kept in the likelihood with equal
## weight on both states and excluded from origin counts (fit_one codes such a tip NA and count_origins uses
## na.rm). Under parsimony (MPR), edge_states treats a flat-prior tip exactly as pars_floor does (after the
## second correction, both return it as NA), so both give the same parsimony origin count for a block.
##
## This file defines functions only. It loads no data, writes nothing, and expects ape, corHMM, phangorn and
## phytools to have been attached by the caller.

## ---- copied from rl_rm_asr.R (version of 2026-09-09, md5 f8f359b24918098ad64c5eb9edd77a17), then corrected ----
## First correction (ambiguity code). The earlier code wrote an ambiguous tip as "1&2". The traits carry the
## observed states "0" and "1", and corHMM numbers observed labels from 1, so "1&2" asked for ambiguity between
## the second state and a third state that does not exist. Because "0" means present here, no ambiguous tip
## (for example one in a proteome where absence cannot be assessed) could be present, and a phantom third state
## was fitted. "0&1" is the intended flat prior over the two observed states and gives the same likelihood as "?".
tipstate <- function(p1) ifelse(p1 == 1, "1", ifelse(p1 == 0, "0", "0&1"))
is_amb   <- function(s) s == "0&1"

count_origins <- function(phy, st) {
  e <- phy$edge
  sum(st[e[, 1]] == 0 & st[e[, 2]] == 1, na.rm = TRUE)
}

fit_one <- function(phy, sp, x, model) {
  fit <- try(corHMM(phy = phy, data = data.frame(sp = sp, t = x,
                                                 stringsAsFactors = FALSE),
                    rate.cat = 1, model = model, node.states = "marginal",
                    nstarts = 0, n.cores = 1), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  nd <- apply(fit$states, 1, function(r) if (r[2] > 0.5) 1L else 0L)
  tp <- x[match(phy$tip.label, sp)]
  st <- c(ifelse(is_amb(tp), NA_integer_, as.integer(tp)), nd)
  list(aicc = fit$AICc, origins = count_origins(phy, st))
}

pars_floor <- function(phy, sp, x) {
  tp <- x[match(phy$tip.label, sp)]
  ch <- ifelse(is_amb(tp), "?", as.character(tp))
  pd <- phyDat(matrix(ch, ncol = 1, dimnames = list(phy$tip.label, NULL)),
               type = "USER", levels = c("0", "1"), ambiguity = "?")
  a <- try(ancestral.pars(phy, pd, type = "MPR"), silent = TRUE)
  if (inherits(a, "try-error")) return(NA_integer_)
  m  <- do.call(rbind, lapply(a, function(z) z[1, ]))
  nd <- apply(m, 1, function(r) if (r[2] > r[1]) 1L else 0L)
  ## Second correction (unobserved tips under parsimony). MPR leaves a "?" tip at 50/50, which fails
  ## r[2] > r[1] and so was hardened to 0, present. A tip whose state was never observed must not be resolved
  ## to a state, so it is set back to NA, as the corHMM path already does and as the Methods describe.
  nd[which(is_amb(tp))] <- NA_integer_
  count_origins(phy, nd)
}
## ---- end of the functions copied from rl_rm_asr.R ----

## The outgroup that rl_rm_asr.R derives from rooted_A2 (its lines from "rc <- ..." to
## "OUTGROUP <- ..."), with the reference tree as the argument: the smaller of the two
## clades below the root.
a2_outgroup <- function(ref) {
  rc  <- ref$edge[ref$edge[, 1] == Ntip(ref) + 1, 2]
  d1  <- if (rc[1] <= Ntip(ref)) ref$tip.label[rc[1]] else extract.clade(ref, rc[1])$tip.label
  d2  <- if (rc[2] <= Ntip(ref)) ref$tip.label[rc[2]] else extract.clade(ref, rc[2])$tip.label
  if (length(d1) <= length(d2)) d1 else d2
}

## The rooting step of rl_rm_asr.R, with its three branches and their ROOT_MODE labels
## unchanged. `mode` is "auto" or "midpoint" (what the job asks for); `tn` names the tree
## in the SKIP line. Returns list(phy, rootmode), or NULL when the tree could not be rooted,
## in which case the same SKIP line as in rl_rm_asr.R is printed.
root_tree <- function(phy, tn, mode, outgroup) {
  og  <- intersect(outgroup, phy$tip.label)
  if (mode == "midpoint") {
    rootmode <- "midpoint_control"
    phy <- tryCatch(midpoint.root(phy), error = function(e) phy)
  } else if (length(og) > 0 && length(og) < Ntip(phy) && is.monophyletic(phy, og)) {
    rootmode <- "outgroup"
    phy <- tryCatch(root(phy, outgroup = og, resolve.root = TRUE),
                    error = function(e) phy)
  } else {
    rootmode <- "midpoint_forced_outgroup_not_monophyletic"
    phy <- tryCatch(midpoint.root(phy), error = function(e) phy)
  }
  if (!is.rooted(phy)) {
    cat(sprintf("  SKIP %s/%s: could not root\n", tn, rootmode))
    return(NULL)
  }
  list(phy = phy, rootmode = rootmode)
}

## md5 checksum of the Newick text of a (pruned) tree: carried by every per-edge row, so that
## a core-component block and a replacement block are compared only when they were
## reconstructed on the same tree. Base R only (tools::md5sum on a temporary file).
tree_md5 <- function(phy) {
  f <- tempfile(fileext = ".nwk")
  write.tree(phy, file = f)
  m <- unname(tools::md5sum(f))
  unlink(f)
  m
}

## Per-edge states. For model ER, SYM or ARD this is fit_one's reconstruction (corHMM
## marginal, internal nodes hardened at 0.5, flat-prior tips NA); for model "MPR" it is
## pars_floor's (phangorn MPR, each node resolved by r[2] > r[1], flat-prior tips NA).
## Returns NULL when the fit fails, exactly where fit_one/pars_floor return NULL/NA.
## Otherwise list(aicc, states, edges): `states` is indexed by ape node number (tips
## 1..Ntip, then internal nodes), so count_origins(phy, states) is the origin count for
## the same block as reported by the loss-count script; `edges` has one row per row of
## phy$edge.
edge_states <- function(phy, sp, x, model) {
  tp <- x[match(phy$tip.label, sp)]
  if (model == "MPR") {
    ch <- ifelse(is_amb(tp), "?", as.character(tp))
    pd <- phyDat(matrix(ch, ncol = 1, dimnames = list(phy$tip.label, NULL)),
                 type = "USER", levels = c("0", "1"), ambiguity = "?")
    a <- try(ancestral.pars(phy, pd, type = "MPR"), silent = TRUE)
    if (inherits(a, "try-error")) return(NULL)
    m  <- do.call(rbind, lapply(a, function(z) z[1, ]))
    st <- apply(m, 1, function(r) if (r[2] > r[1]) 1L else 0L)
    ## Second correction, the same tie as in pars_floor. This tie is what produced spurious regains: an
    ## absent parent above a tip hardened to present reads as a reacquisition of the whole protein.
    st[which(is_amb(tp))] <- NA_integer_
    aicc <- NA_real_
  } else {
    fit <- try(corHMM(phy = phy, data = data.frame(sp = sp, t = x,
                                                   stringsAsFactors = FALSE),
                      rate.cat = 1, model = model, node.states = "marginal",
                      nstarts = 0, n.cores = 1), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    nd <- apply(fit$states, 1, function(r) if (r[2] > 0.5) 1L else 0L)
    st <- c(ifelse(is_amb(tp), NA_integer_, as.integer(tp)), nd)
    aicc <- fit$AICc
  }
  if (length(st) != Ntip(phy) + Nnode(phy))
    stop("REFUSED: ", length(st), " states for ", Ntip(phy) + Nnode(phy),
         " nodes under model ", model)
  e <- phy$edge
  tip <- e[, 2] <= Ntip(phy)
  list(aicc = aicc, states = unname(st),
       edges = data.frame(PARENT_NODE = e[, 1], CHILD_NODE = e[, 2],
                          PARENT_STATE = unname(st[e[, 1]]),
                          CHILD_STATE = unname(st[e[, 2]]),
                          IS_TIP = ifelse(tip, "yes", "no"),
                          TIP_LABEL = ifelse(tip, phy$tip.label[pmin(e[, 2], Ntip(phy))], "-"),
                          stringsAsFactors = FALSE))
}
