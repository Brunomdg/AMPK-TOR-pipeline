#!/usr/bin/env Rscript
# >>> deposit header
# Purpose: Reconstructs, for one trait under one coding, the ancestral state at both ends of every edge in each
#          of 112 blocks (tree, rooting, taxon subset, model); this repaired copy reads the repaired library.
# Methods: Ancestral-state reconstructions
# Inputs:  one per-species coding table (`Q0_REPLACEMENT_per_species.tsv`, its second-coding version, or the file
#          named by Q0_ASR_STANDARD); five species trees (`rooted_A2`, the chronogram `species_chronogram_209`,
#          `p04_pmsf_209_constr`, `guide_209`, `p04_au_209`)
# Outputs: `<trait>__<coding>.edges.tsv` and `<trait>__<coding>.blocks.tsv` in `Q0_replacement/asr`, or in the
#          folder named by Q0_ASR_OUT
# Tools:   R 4.5.3 with ape 5.8.1, corHMM 2.8, phangorn 2.12.1 and phytools 2.5.2 (the ms210_phase2 environment,
#          versions as noted on 9 September 2026 in its job script)
# <<< deposit header
## deposit path defaults: begin
.TA210_ENVS <- local({v <- Sys.getenv("TA210_ENVS"); if (nzchar(v)) v else "/gpfs2/scratch/bdigenov/envs"})
.TA210_SCRATCH <- local({v <- Sys.getenv("TA210_SCRATCH"); if (nzchar(v)) v else "/gpfs2/scratch/bdigenov"})
## deposit path defaults: end
## Per-edge ancestral states for the replacement analysis: the reconstruction behind the test that separates
## decay (a candidate replacement is lost as well) from reassignment (the replacement is kept where the core
## component was lost).
##
## What it does. For one `_p1` trait column of the per-species table under one coding, it reconstructs the
## state at both ends of every edge of every block that rl_rm_asr.R enumerates for the loss counts: five
## topologies, two rooting criteria where the outgroup allows it (seven combinations of tree and rooting),
## four taxon subsets (all, no_Apicomplexa, no_SAR, genome_only), and four models (the three corHMM transition
## models, hardened at 0.5, and MPR parsimony): 7 x 4 x 4 = 112 blocks. It writes
##   14_remodeling/Q0_replacement/asr/<trait>__<coding>.edges.tsv    one row per edge per block
##   14_remodeling/Q0_replacement/asr/<trait>__<coding>.blocks.tsv   one row per block
## The blocks file carries STATUS (OK or FIT_FAILED), the tip counts, AICc, the origin count (the number
## rl_rm_asr.R reports as ORIGINS_MODEL or ORIGINS_PARSIMONY for the same block, which rl_rm_q0_check.py
## compares on AMPKg_lost against gamma_lost) and the md5 checksum of the pruned tree.
##
## Codings. The primary coding is Q0_REPLACEMENT_per_species.tsv (present 0, uncertain 0.5, absent 1). The
## second coding, Q0_REPLACEMENT_per_species__R37coding.tsv, is identical except that every proteome in which
## absence cannot be assessed is coded 0.5 on every target, present calls included. A third, standard coding
## can be enabled (see Q0_ASR_STANDARD below). An ambiguous (0.5) tip is kept in the likelihood as a flat
## prior and excluded from origin counts.
##
## Why per edge and not per count. The loss counts say how many times a component was lost; the replacement
## analysis asks, for each of those losses, what the candidate replacement was doing on the same edge and below
## it. That can be read only from the states at both ends of every edge on the same pruned tree, so the
## functions of rl_rm_asr.R are reused through rl_rm_asr_lib.R (here its repaired copy) and the hardened states
## are written out instead of being summed. The classification (REPLACEMENT_ABSENT_BEFORE, CO_LOSS_SAME_EDGE,
## CO_LOSS_DESCENDANT, REASSIGNMENT_COMPATIBLE) is made by rl_rm_q0_replacement.py --stage pairs, which reads
## the core and the replacement edge files of one coding and refuses unless their pruned trees are identical
## (same PRUNED_TREE_MD5 and edge list).
##
## One trait and coding per process, each with its own lock. The lock of rl_rm_asr.R (`.asr.lock`) is never
## touched; this script locks `.q0asr.<trait>.<coding>.lock` inside the output folder, so array tasks can run
## side by side and beside a rerun of rl_rm_asr.R.
##
## Index mapping (single-argument form, used by the array job rl_rm_q0_asr.sbatch). Let n be the number of
## `_p1` columns in the per-species header, in header order. Index i selects trait (i mod n) + 1, under the
## primary coding for i < n and the second coding for n <= i < 2n. An index at or beyond 2n (for example in a
## partial run with fewer traits) prints NOTHING TO DO and exits 0, so an array of fixed size also serves any
## smaller n (the 20 tasks of rl_rm_q0_asr.sbatch cover up to ten traits).
##
## Sign. A `<T>_absent` trait reconstructs the absence of a replacement, so its state 0 is the replacement
## present and a 0 -> 1 edge on it is a loss of the replacement, the event compatible with decay.
##
## Usage:  Rscript rl_rm_q0_asr.R <index>              one array task, see the index mapping above
##         Rscript rl_rm_q0_asr.R <trait> <coding>     for example `AMPKg_absent STANDARD` (as asr_repair_run.sbatch
##                                                     calls it); coding is `TRAITS` (primary), `R37` (second) or
##                                                     `STANDARD` (standard)
## Env:    PROJECT_ROOT (project folder), CONCORD_MS210 (folder of the species trees, default the concord path),
##         Q0_ASR_STANDARD and Q0_ASR_OUT (see below)

suppressMessages({library(ape); library(corHMM); library(phangorn); library(phytools)})

SEED <- 488761235
set.seed(SEED)

ROOT <- Sys.getenv("PROJECT_ROOT",
  paste0(.TA210_SCRATCH, "/ms210_phase2/TOR-AMPK_orthologs_10target_R1"))
SCR  <- file.path(ROOT, "14_remodeling", "scripts")
SELF <- file.path(ROOT, "19_asr_repair", "scripts", "rl_rm_q0_asr.R")
source(file.path(SCR, "rl_rm_prov.R"))
## The library is read from its repaired copy; rl_rm_prov.R, above, is unchanged and still read from SCR.
source(file.path(ROOT, "19_asr_repair", "scripts", "rl_rm_asr_lib.R"))
CONCORD <- Sys.getenv("CONCORD_MS210", paste0(.TA210_SCRATCH, "/concord/ms210"))
CH  <- file.path(CONCORD, "P04b_chronogram")
P04 <- file.path(CONCORD, "P04_tree")
Q0  <- file.path(ROOT, "14_remodeling", "Q0_replacement")
OUT <- file.path(Q0, "asr")
CODINGS <- c(TRAITS = "Q0_REPLACEMENT_per_species.tsv",
             R37    = "Q0_REPLACEMENT_per_species__R37coding.tsv")
## A third coding, the standard coding (added on 13 September 2026: one coding rule and one `<T>_absent`
## name for every trait), whose per-species file is in Q7_alltargets. It is enabled only by setting
## Q0_ASR_STANDARD to that file's absolute path, and its outputs go to Q0_ASR_OUT when that is set. The index
## mapping still covers only the first two codings, so the standard coding needs the two-argument form.
if (nzchar(Sys.getenv("Q0_ASR_STANDARD"))) CODINGS <- c(CODINGS, STANDARD = Sys.getenv("Q0_ASR_STANDARD"))
if (nzchar(Sys.getenv("Q0_ASR_OUT"))) OUT <- Sys.getenv("Q0_ASR_OUT")
coding_path <- function(coding) if (startsWith(CODINGS[[coding]], "/")) CODINGS[[coding]] else file.path(Q0, CODINGS[[coding]])
MODELS <- c("ER", "SYM", "ARD")
TREES  <- c(rooted_A2   = file.path(CH, "rooted_A2.treefile"),
            chronogram  = file.path(CH, "species_chronogram_209.treefile"),
            pmsf_constr = file.path(P04, "p04_pmsf_209_constr.treefile"),
            guide_209   = file.path(P04, "guide_209.treefile"),
            au_209      = file.path(P04, "p04_au_209.treefile"))
## trees that are also run midpoint-rooted, so that rooting varies independently of topology
MIDPOINT_ALSO <- c("rooted_A2", "pmsf_constr")

read_species <- function(coding) {
  p <- coding_path(coding)
  if (!file.exists(p))
    stop("REFUSED: ", p, " is missing; run rl_rm_q0_replacement.py --stage traits first")
  d <- read.delim(p, comment.char = "#", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) != 210) stop("REFUSED: ", p, " has ", nrow(d), " rows, not 210")
  d
}
trait_columns <- function(d) sub("_p1$", "", grep("_p1$", names(d), value = TRUE))
banner_of <- function(p, needle) {
  L <- readLines(p)
  L[startsWith(L, "##") & grepl(needle, L, fixed = TRUE)]
}

## ---- choose the trait and the coding ----
args <- commandArgs(trailingOnly = TRUE)
d0 <- read_species("TRAITS")
traits <- trait_columns(d0)
n <- length(traits)
if (n < 1) stop("REFUSED: no _p1 column in the per_species header")
if (length(args) == 1 && grepl("^[0-9]+$", args[1])) {
  i <- as.integer(args[1])
  if (i >= 2 * n) {
    cat(sprintf("NOTHING TO DO: index %d is at or beyond 2 x %d traits in the per_species header\n", i, n))
    quit(status = 0)
  }
  TRAIT  <- traits[(i %% n) + 1]
  CODING <- names(CODINGS)[(i %/% n) + 1]
} else if (length(args) == 2) {
  TRAIT  <- args[1]
  CODING <- args[2]
} else {
  stop("REFUSED: usage is `rl_rm_q0_asr.R <index>` or `rl_rm_q0_asr.R <trait> <coding>`")
}
if (!CODING %in% names(CODINGS))
  stop("REFUSED: coding ", CODING, " is not one of ", paste(names(CODINGS), collapse = ", "))
if (CODING == "STANDARD") traits <- trait_columns(read_species("STANDARD"))
if (!TRAIT %in% traits)
  stop("REFUSED: trait ", TRAIT, " is not a _p1 column of the per_species header: ",
       paste(traits, collapse = ", "))
## The CAMKK and LKB1 absence traits are reconstructed only after rl_require_camkk_lkb1_joint_freeze.py
## confirms that the calls for both kinases are final (required since 16 September 2026); the other traits
## are unaffected.
if (TRAIT %in% c("CAMKK_absent", "LKB1_absent") && system2(paste0(.TA210_ENVS, "/ms210_phase2/bin/python3"), paste0(.TA210_SCRATCH, "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/scripts/rl_require_camkk_lkb1_joint_freeze.py")) != 0)
  stop("REFUSED: CAMKK/LKB1 joint freeze is not complete; ", TRAIT, " is not reconstructed")

SPECIESF <- coding_path(CODING)
d <- read_species(CODING)
col <- paste0(TRAIT, "_p1")
if (!col %in% names(d)) stop("REFUSED: ", SPECIESF, " lacks column ", col)
if (!all(d[[col]] %in% c(0, 0.5, 1)))
  stop("REFUSED: ", col, " carries a value outside {0, 0.5, 1}")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOCK <- file.path(OUT, sprintf(".q0asr.%s.%s.lock", TRAIT, CODING))
if (file.exists(LOCK))
  stop("REFUSED: ", LOCK, " exists, another rl_rm_q0_asr.R is running this trait and coding (pid ",
       tryCatch(readLines(LOCK)[1], error = function(e) "unknown"),
       "). If that pid is not alive the lock is stale (a killed job): remove it and re-run.")
writeLines(as.character(Sys.getpid()), LOCK)
## The lock is removed by the `finally` of the tryCatch that encloses all the work below.
## on.exit() does not do this at the top level of an Rscript (in ten test runs on 9 September
## 2026 each run left its lock behind), so it is not used here.
tryCatch({

cat(sprintf("Q0 ASR trait=%s coding=%s  (%d species, %d traits in the header)\n",
            TRAIT, CODING, nrow(d), n))

ref <- read.tree(TREES[["rooted_A2"]])
OUTGROUP <- a2_outgroup(ref)

jobs <- list()
for (tn in names(TREES)) {
  jobs[[length(jobs) + 1]] <- list(tn = tn, mode = "auto")
  if (tn %in% MIDPOINT_ALSO)
    jobs[[length(jobs) + 1]] <- list(tn = tn, mode = "midpoint")
}

edge_rows <- list()
block_rows <- list()
for (jb in jobs) {
  tn <- jb$tn
  phy <- tryCatch(read.tree(TREES[[tn]]), error = function(e) NULL)
  if (is.null(phy)) { cat(sprintf("  SKIP %s: tree unreadable\n", tn)); next }
  phy <- keep.tip(phy, intersect(phy$tip.label, d$SPECIES))
  rt <- root_tree(phy, tn, jb$mode, OUTGROUP)
  if (is.null(rt)) next
  phy <- rt$phy
  rootmode <- rt$rootmode
  cat(sprintf("  %-12s %s\n", tn, rootmode))
  dd <- d[match(phy$tip.label, d$SPECIES), ]
  x <- tipstate(dd[[col]])
  subsets <- list(all            = rep(TRUE, nrow(dd)),
                  no_Apicomplexa = dd$CLADE != "Apicomplexa",
                  no_SAR         = dd$SUPERGROUP != "SAR",
                  genome_only    = dd$DATA_GRADE == "G")
  for (sn in names(subsets)) {
    k <- subsets[[sn]]
    p2 <- keep.tip(phy, phy$tip.label[k]); x2 <- x[k]; s2 <- dd$SPECIES[k]
    tm <- tree_md5(p2)
    for (md in c(MODELS, "MPR")) {
      es <- edge_states(p2, s2, x2, md)
      ok <- !is.null(es)
      block_rows[[length(block_rows) + 1]] <- data.frame(
        TREE = tn, ROOT_MODE = rootmode, TRAIT = TRAIT, CODING = CODING, SUBSET = sn,
        MODEL = md, STATUS = if (ok) "OK" else "FIT_FAILED",
        N_TIP = Ntip(p2), N_STATE1 = sum(x2 == "1"), N_STATE0 = sum(x2 == "0"),
        N_FLAT_PRIOR = sum(is_amb(x2)),
        AICC = if (ok && !is.na(es$aicc)) round(es$aicc, 2) else NA_real_,
        ORIGINS = if (ok) count_origins(p2, es$states) else NA_integer_,
        PRUNED_TREE_MD5 = tm, stringsAsFactors = FALSE)
      if (!ok) { cat(sprintf("    FIT_FAILED %s %s %s\n", tn, sn, md)); next }
      edge_rows[[length(edge_rows) + 1]] <- data.frame(
        TREE = tn, ROOT_MODE = rootmode, SUBSET = sn, MODEL = md, es$edges,
        PRUNED_TREE_MD5 = tm, stringsAsFactors = FALSE)
    }
  }
}

if (length(block_rows) == 0) stop("REFUSED: no block was enumerated (no tree could be read or rooted)")
blocks <- do.call(rbind, block_rows)

missing <- banner_of(SPECIESF, "TARGETS_MISSING:")
prov <- prov_lines(SELF, c(per_species = SPECIESF,
                           rooted_A2 = TREES[["rooted_A2"]], chronogram = TREES[["chronogram"]],
                           pmsf_constr = TREES[["pmsf_constr"]], guide_209 = TREES[["guide_209"]],
                           au_209 = TREES[["au_209"]]), SEED)
sign_line <- paste0(
  "## SIGN INVERSION (DECISIONS I7): a `<T>_absent` trait is the ABSENCE of a replacement; its 0 state is the ",
  "replacement present, the reassignment-compatible pattern against a core loss, and its 0 -> 1 edge the ",
  "decay-compatible one.")

write_out <- function(f, banner, tab) {
  ftmp <- paste0(f, ".tmp")
  con <- file(ftmp, "w")
  writeLines(banner, con)
  close(con)
  write.table(tab, ftmp, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "NA")
  if (!file.rename(ftmp, f)) stop("REFUSED: could not rename ", ftmp, " over ", f)
}

fb <- file.path(OUT, sprintf("%s__%s.blocks.tsv", TRAIT, CODING))
write_out(fb, c(
  sprintf("## Q0 ASR BLOCKS for trait %s under the %s coding: one row per (TREE, ROOT_MODE, SUBSET, MODEL) block rl_rm_asr.R enumerates for Q2, on the same pruned trees.", TRAIT, CODING),
  "## STATUS OK means the reconstruction produced per-edge states (in the matching .edges.tsv); FIT_FAILED means corHMM or MPR failed and the block has no edges and is counted, never filled.",
  "## ORIGINS counts parent-0 to child-1 edges on the hardened states: for ER/SYM/ARD the number Q2 reports as ORIGINS_MODEL, for MPR its ORIGINS_PARSIMONY, for the same block of the same trait.",
  "## N_FLAT_PRIOR tips are coded 1&2 per R5: kept in the likelihood with equal weight, excluded from origin counts.",
  sign_line, missing, prov), blocks)
cat(sprintf("wrote %s  (%d blocks, %d OK)\n", fb, nrow(blocks), sum(blocks$STATUS == "OK")))

if (length(edge_rows) == 0)
  stop("REFUSED: every block FIT_FAILED for ", TRAIT, "/", CODING, "; no edges to write (blocks file kept as the record)")
edges <- do.call(rbind, edge_rows)
fe <- file.path(OUT, sprintf("%s__%s.edges.tsv", TRAIT, CODING))
write_out(fe, c(
  sprintf("## Q0 PER-EDGE ANCESTRAL STATES for trait %s under the %s coding: the state at both ends of every edge of every (TREE x ROOT_MODE x SUBSET x MODEL) block rl_rm_asr.R enumerates for Q2, on the same pruned trees.", TRAIT, CODING),
  "## PARENT_STATE/CHILD_STATE: corHMM marginal reconstruction hardened at 0.5 for ER/SYM/ARD (a flat-prior tip is NA: kept in the likelihood, excluded from every count, R5); phangorn MPR for MODEL=MPR (a flat-prior tip resolved as pars_floor resolves it). A 0 -> 1 edge is a LOSS of the component on a `_lost` trait and a loss of the REPLACEMENT on an `_absent` trait.",
  "## Node numbers are ape's for the pruned tree of the block (tips 1..N_TIP, root N_TIP+1); PRUNED_TREE_MD5 is the md5 of that tree's Newick text, and rl_rm_q0_replacement.py refuses to compare a core block with a replacement block whose md5 or edge list differs.",
  "## ROOT_MODE separates rooting from topology as in Q2: rooted_A2 and pmsf_constr run both outgroup-rooted and midpoint-rooted; guide_209 and au_209 are midpoint-rooted of necessity.",
  sign_line, missing, prov), edges)
cat(sprintf("wrote %s  (%d edge rows)\n", fe, nrow(edges)))

}, finally = unlink(LOCK))
